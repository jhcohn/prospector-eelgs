from prospect.models import model_setup
from prospect.io import read_results
import os, sys, time
import prosp_dutils_orig as prosp_dutils  # original file without bunches of print statements
# import prosp_dutils
import numpy as np
from copy import copy
from astropy import constants
import argparse
import matplotlib.pyplot as plt
import pickle


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


def calc_extra_quantities(sample_results, ncalc=3000, **kwargs):
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
    half_time, sfr_10, sfr_100, ssfr_100, stellar_mass, ssfr_10 = [np.zeros(shape=(ncalc)) for i in range(6)]

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
              't_sfh': t}
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

    '''
    # NEW PLOTTING
    fig = plt.figure()
    sfh_ax = fig.add_axes([0.15, 0.15, 0.6, 0.6], zorder=32)
    import boop
    boop.add_sfh_plot([extra_output], fig, main_color=['black'], ax_inset=sfh_ax, text_size=3, lw=3)  # lw=5
    plt.show()
    # NEW PLOTTING
    plt.plot(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'])
    plt.ylabel(r'M$_\odot$ yr$^{-1}$')
    plt.xlabel('Gyr')
    plt.show()
    '''

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False


if __name__ == "__main__":

    ### don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    # parser.add_argument('parfile', type=str)
    parser.add_argument('--outname')
    # parser.add_argument('--measure_spectral_features', type=str2bool)
    # parser.add_argument('--mags_nodust', type=str2bool)
    # parser.add_argument('--ir_priors', type=str2bool)
    parser.add_argument('--ncalc', type=int)

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    outname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + kwargs['outname']
    obj = ''
    count = 0
    field = ''
    for i in kwargs['outname']:
        if i == '_':
            count += 1
        if count == 0:
            obj += i
        elif count == 1:
            field += i
        elif count == 2:
            break

    write = obj + field + '_sfh_out2.pkl'

    print(kwargs)
    post_processing(out_file=outname, filename=write, **kwargs)

'''
RUNNING WITH:

python sfh_output.py --outname=6459_multirun_commentedcontinuum_1498501917_mcmc.h5 --ncalc=2000
--measure_spectral_features=False --mags_nodust=False parfile=eelg_multirun_params.py

python sfh_output.py --outname=6459_cosmos_multirun_commentedcontinuum_1498501917_mcmc.h5 --ncalc=2000
'''
