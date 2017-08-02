from prospect.models import model_setup
from prospect.io import read_results
import os, sys, time
import prosp_dutils_orig as prosp_dutils  # original file without bunches of print statements
import numpy as np
from copy import copy
from astropy import constants
import argparse
import matplotlib.pyplot as plt
import pickle
import prospect.io.read_results as bread
import matplotlib.pyplot as plt
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log
np.errstate(invalid='ignore')


# FROM sfh_output.py
def maxprob_model(sample_results, sps):
    ### grab maximum probability, plus the thetas that gave it
    maxprob = sample_results['flatprob'].max()
    maxtheta = sample_results['flatchain'][sample_results['flatprob'].argmax()]

    ### ensure that maxprob stored is the same as calculated now
    current_maxprob = prosp_dutils.test_likelihood(sps, sample_results['model'], sample_results['obs'], maxtheta,
                                                   sample_results['run_params']['param_file'])

    print('Best-fit lnprob currently: {0}'.format(current_maxprob))
    print('Best-fit lnprob during sampling: {0}'.format(maxprob))

    return maxtheta, maxprob


def sample_flatchain(chain, lnprob, parnames, ir_priors=True, include_maxlnprob=True, nsamp=2000):
    '''
    CURRENTLY UNDER DEVELOPMENT
    goal: sample the flatchain in a smart way
    '''

    ##### use randomized, flattened chain for posterior draws
    # don't allow draws which are outside the priors
    good = np.isfinite(lnprob) == True

    ### cut in IR priors
    if ir_priors:
        gamma_idx = parnames == 'duste_gamma'
        umin_idx = parnames == 'duste_umin'
        qpah_idx = parnames == 'duste_qpah'
        gamma_prior = 0.15
        umin_prior = 15
        qpah_prior = 7
        in_priors = (chain[:, gamma_idx] < gamma_prior) & (chain[:, umin_idx] < umin_prior) & (
        chain[:, qpah_idx] < qpah_prior) & good[:, None]
        if in_priors.sum() < 2 * nsamp:
            print('Insufficient number of samples within the IR priors! Not applying IR priors.')
        else:
            good = in_priors

    sample_idx = np.random.choice(np.where(good)[0], nsamp)

    ### include maxlnprob?
    if include_maxlnprob:
        sample_idx[0] = lnprob.argmax()

    return sample_idx


def set_sfh_time_vector(sample_results, ncalc):
    # if parameterized, calculate linearly in 100 steps from t=0 to t=tage
    # if nonparameterized, calculate at bin edges.
    if 'tage' in sample_results['model'].theta_labels():
        nt = 100
        idx = np.array(sample_results['model'].theta_labels()) == 'tage'
        maxtime = np.max(sample_results['flatchain'][:ncalc, idx])
        t = np.linspace(0, maxtime, num=nt)
    elif 'agebins' in sample_results['model'].params:
        in_years = 10 ** sample_results['model'].params['agebins'] / 1e9
        t = np.concatenate((np.ravel(in_years) * 0.9999, np.ravel(in_years) * 1.001))
        t.sort()  # samples on either side of the bin edges
        t = t[1:-1]  # remove points older than oldest bin, younger than youngest bin
        t = np.clip(t, 1e-3, np.inf)  # nothing younger than 1 Myr!
    else:
        sys.exit('ERROR: not sure how to set up the time array here!')
    return t


def calc_extra_quantities(sample_results, ncalc=2000, **kwargs):  # ncalc=3000
    ''''
    CALCULATED QUANTITIES
    model nebular emission line strength
    model star formation history parameters (ssfr,sfr,half-mass time)
    '''

    # different options for what to calculate
    # speedup is measured in runtime, where runtime = ncalc * model_call
    opts = {
        'restframe_optical_photometry': False,  # currently deprecated! but framework exists in
                                                # restframe_optical_properties
        'ir_priors': False,  # no cost
        'measure_spectral_features': False,  # cost = 2 runtimes
        'mags_nodust': False  # cost = 1 runtime
    }
    if kwargs:
        for key in kwargs.keys():
            opts[key] = kwargs[key]

    parnames = np.array(sample_results['model'].theta_labels())

    ##### describe number of components in Prospector model [legacy]
    sample_results['ncomp'] = np.sum(['mass' in x for x in sample_results['model'].theta_labels()])

    ##### array indexes over which to sample the flatchain
    sample_idx = sample_flatchain(sample_results['flatchain'], sample_results['flatprob'], parnames,
                                  ir_priors=opts['ir_priors'], include_maxlnprob=True, nsamp=ncalc)

    ##### initialize output arrays for SFH + emission line posterior draws
    half_time, sfr_10, sfr_100, ssfr_100, stellar_mass, ssfr_10, ssfr_full = [np.zeros(shape=(ncalc)) for i in range(7)]

    ##### set up time vector for full SFHs
    t = set_sfh_time_vector(sample_results, ncalc)  # returns array of len=18
    intsfr = np.zeros(shape=(t.shape[0], ncalc))

    ##### initialize sps, calculate maxprob
    # also confirm probability calculations are consistent with fit
    sps = model_setup.load_sps(**sample_results['run_params'])
    maxthetas, maxprob = maxprob_model(sample_results, sps)

    ##### set up model flux vectors
    mags = np.zeros(shape=(len(sample_results['obs']['filters']), ncalc))
    try:
        wavelengths = sps.wavelengths
    except AttributeError:
        wavelengths = sps.csp.wavelengths
    spec = np.zeros(shape=(wavelengths.shape[0], ncalc))

    ##### modify nebular status to ensure emission line production
    # don't cache, and turn on
    if sample_results['model'].params['add_neb_emission'] == 2:
        sample_results['model'].params['add_neb_emission'] = np.array([True])
    sample_results['model'].params['nebemlineinspec'] = np.array([True])

    loop = 0
    ######## posterior sampling #########
    for jj, idx in enumerate(sample_idx):

        ##### model call, to set parameters
        thetas = copy(sample_results['flatchain'][idx])
        spec[:, jj], mags[:, jj], sm = sample_results['model'].mean_model(thetas, sample_results['obs'], sps=sps)

        ##### extract sfh parameters
        # pass stellar mass to avoid extra model call
        sfh_params = prosp_dutils.find_sfh_params(sample_results['model'], thetas, sample_results['obs'], sps, sm=sm)

        ##### calculate SFH
        intsfr[:, jj] = prosp_dutils.return_full_sfh(t, sfh_params)

        ##### solve for half-mass assembly time
        half_time[jj] = prosp_dutils.halfmass_assembly_time(sfh_params)

        ##### calculate time-averaged SFR
        sfr_10[jj] = prosp_dutils.calculate_sfr(sfh_params, 0.01, minsfr=-np.inf, maxsfr=np.inf)  # avg over 10 Myr
        sfr_100[jj] = prosp_dutils.calculate_sfr(sfh_params, 0.1, minsfr=-np.inf, maxsfr=np.inf)  # avg over 100 Myr

        ##### calculate mass, sSFR
        stellar_mass[jj] = sfh_params['mass']
        ssfr_10[jj] = sfr_10[jj] / stellar_mass[jj]
        ssfr_100[jj] = sfr_100[jj] / stellar_mass[jj]
        ssfr_full = intsfr[:, jj] / stellar_mass[jj]

        loop += 1
        print('loop', loop)

    #### QUANTILE OUTPUTS #
    extra_output = {}

    ##### CALCULATE Q16,Q50,Q84 FOR EXTRA PARAMETERS
    extra_flatchain = np.dstack(
        (half_time, sfr_10, sfr_100, ssfr_10, ssfr_100, stellar_mass))[0]

    #### EXTRA PARAMETER OUTPUTS
    extras = {'flatchain': extra_flatchain,
              'parnames': np.array(
                  ['half_time', 'sfr_10', 'sfr_100', 'ssfr_10', 'ssfr_100', 'stellar_mass']),
              'sfh': intsfr,
              't_sfh': t,
              'ssfr': ssfr_full}
    extra_output['extras'] = extras

    #### BEST-FITS
    bfit = {'maxprob_params': maxthetas,
            'maxprob': maxprob,
            'sfh': intsfr[:, 0],
            'half_time': half_time[0],
            'sfr_10': sfr_10[0],
            'sfr_100': sfr_100[0],
            'spec': spec[:, 0],
            'mags': mags[:, 0]}
    extra_output['bfit'] = bfit

    ##### spectral features
    return extra_output


def update_all(runname, **kwargs):
    '''
    change some parameters, need to update the post-processing?
    run this!
    '''
    filebase, parm_basename, ancilname = prosp_dutils.generate_basenames(runname)
    for param in parm_basename:
        post_processing(param, **kwargs)


def post_processing(out_file, filename, **kwargs):
    '''
    Driver. Loads output, runs post-processing routine.
    '''

    sample_results, powell_results, model = read_results.results_from(out_file)

    ### create flatchain, run post-processing
    sample_results['flatchain'] = prosp_dutils.chop_chain(sample_results['chain'], **sample_results['run_params'])
    sample_results['flatprob'] = prosp_dutils.chop_chain(sample_results['lnprobability'],
                                                         **sample_results['run_params'])
    extra_output = calc_extra_quantities(sample_results, **kwargs)

    with open(filename, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(extra_output, newfile, pickle.HIGHEST_PROTOCOL)


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False


# FROM quickgrab.py
def printer(out_file, corner=False, trace=False):
    objname = ''
    count = 0
    field = ''
    for i in kwargs['outname']:
        if i == '_':
            count += 1
        elif count == 0:
            objname += i
        elif count == 1:
            field += i
        elif count == 2:
            break
    print(field)

    res, pr, mod = bread.results_from(out_file)

    # PRINT TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
    tracefig, prob = bread.param_evol(res)  # print tracefig, store probability
    if trace:
        plt.show()

    # FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
    # print(prob)
    print('max', prob.max())
    row = prob.argmax() / len(prob[0])
    col = prob.argmax() - row * len(prob[0])
    walker, iteration = row, col
    print(walker, iteration)

    # PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
    if corner:
        cornerfig = bread.subtriangle(res, start=400, thin=5, show_titles=True)
        plt.show()
    # For FAST: truths = [mass, age, tau, dust2] (for 1824: [9.78, 0.25, -1., 0.00])
    return objname, field, res, mod, walker, iteration


def sed(objname, field, res, mod, walker, iteration, param_file, base, print_it=False):
    # PRINT MODEL SED FOR GALAXY
    # We need the correct sps object to generate models
    sargv = sys.argv
    argdict = {'param_file': param_file}
    clargs = model_setup.parse_args(sargv, argdict=argdict)
    run_params = model_setup.get_run_params(argv=sargv, **clargs)
    sps = model_setup.load_sps(**run_params)


    # GET MODELED SPECTRA AND PHOTOMETRY
    # These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
    spec, phot, mfrac = mod.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)

    # PLOT SPECTRUM
    wave = [f.wave_effective for f in res['obs']['filters']]
    wave = np.asarray(wave)

    # CHANGING OBSERVED TO REST FRAME WAVELENGTH
    if field == 'cdfs':
        datname = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
    elif field == 'cosmos':
        datname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'  # main catalog
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog
    elif field == 'uds':
        datname = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.zout'

    with open(datname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(datname, comments='#', delimiter=' ', dtype=dtype)

    with open(zname, 'r') as fz:
        hdr_z = fz.readline().split()
    dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
    zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)

    idx = dat['id'] == objname  # array filled: False when dat['id'] != objname, True when dat['id'] == objname
    zred = zout['z_spec'][idx][0]  # z = z_spec
    if zred == -99:
        zred = zout['z_peak'][idx][0]  # if z_spec does not exist, z = z_phot
    print('redshift', zred)

    wave_rest = []  # REST FRAME WAVELENGTH
    for i in range(len(wave)):
        wave_rest.append(wave[i]/(1 + zred))  # 1 + z = l_obs / l_emit --> l_emit = l_obs / (1 + z)

    # OUTPUT SED results to files
    write_res = objname + '_' + field + '_' + base + '_res_out.pkl'  # results
    with open(write_res, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(res, newfile, pickle.HIGHEST_PROTOCOL)  # res includes res['obs']['maggies'] and ...['maggies_unc']
    write_sed = objname + '_' + field + '_' + base + '_sed_out.pkl'  # model sed
    with open(write_sed, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(phot, newfile, pickle.HIGHEST_PROTOCOL)
    write_restwave = objname + '_' + field + '_' + base + '_restwave_out.pkl'  # rest frame wavelengths
    with open(write_restwave, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(wave_rest, newfile, pickle.HIGHEST_PROTOCOL)
    write_spec = objname + '_' + field + '_' + base + '_spec_out.pkl'  # spectrum
    with open(write_spec, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(spec, newfile, pickle.HIGHEST_PROTOCOL)
    write_sps = objname + '_' + field + '_' + base + '_spswave_out.pkl'  # wavelengths that go with spectrum
    with open(write_sps, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(sps.wavelengths, newfile, pickle.HIGHEST_PROTOCOL)

    # OUTPUT CHI_SQ results to files
    chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
    write_chisq = objname + '_' + field + '_' + base + '_chisq_out.pkl'
    with open(write_chisq, 'wb') as newfile:  # 'wb' because binary format
        pickle.dump(chi_sq, newfile, pickle.HIGHEST_PROTOCOL)
    write_justchi = objname + '_' + field + '_' + base + '_justchi_out.pkl'
    with open(write_justchi, 'wb') as newfile:
        pickle.dump((res['obs']['maggies'] - phot) / res['obs']['maggies_unc'], newfile, pickle.HIGHEST_PROTOCOL)

    if print_it:
        plt.plot(sps.wavelengths, spec)
        plt.xlabel('Wavelength [angstroms]')
        plt.title(str(objname) + ' spec')
        plt.show()

        # HOW CONVERGED IS THE CODE?? LET'S FIND OUT!
        parnames = np.array(res['model'].theta_labels())
        fig, kl_ax = plt.subplots(1, 1, figsize=(7, 7))
        for i in xrange(parnames.shape[0]):
            kl_ax.plot(res['kl_iteration'], np.log10(res['kl_divergence'][:, i]), 'o', label=parnames[i], lw=1.5,
                       linestyle='-', alpha=0.6)

        kl_ax.set_ylabel('log(KL divergence)')
        kl_ax.set_xlabel('iteration')
        # kl_ax.set_xlim(0, nsteps*1.1)

        kl_div_lim = res['run_params'].get('convergence_kl_threshold', 0.018)
        kl_ax.axhline(np.log10(kl_div_lim), linestyle='--', color='red', lw=2, zorder=1)

        kl_ax.legend(prop={'size': 5}, ncol=2, numpoints=1, markerscale=0.7)
        plt.title(str(objname) + ' kl')
        plt.show()


if __name__ == "__main__":

    ### don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--parfile')
    parser.add_argument('--outname')

    files = {'outname': '', 'parfile': ''}

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]
        files[key] = kwargs[key]

    outname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + kwargs['outname']
    obj = ''
    field = ''
    base1 = ''
    count = 0
    for i in kwargs['outname']:
        if i == '_':
            count += 1
        elif count == 0:
            obj += i
        elif count == 1:
            field += i
        elif count == 2:
            base1 += i
        elif count == 3:
            break

    base = base1
    # write = obj + '_' + field + '_' + base + '_sfh_out.pkl'
    write = obj + '_' + field + '_' + base + '_extra_out.pkl'

    out_file = files['outname']
    param_file = files['parfile']
    print(param_file, out_file)

    objname, field, res, mod, walker, iteration = printer(out_file)

    sed(objname, field, res, mod, walker, iteration, param_file, base)

    print(kwargs)
    post_processing(out_file=outname, filename=write, **kwargs)

'''
RUNNING WITH:
python output_all.py --outname=1824_cosmos_decoupled_n1200_1495729287_mcmc.h5 --parfile=eelg_multirun_params.py

# Takes ~8 mins to run
'''
