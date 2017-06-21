from prospect.models import model_setup
from prospect.io import read_results
# import os, prosp_dutils, hickle, sys, time  # I probably don't need hickle
import os, sys, time  # I probably don't need hickle
import prosp_dutils_orig as prosp_dutils  # import prosp_dutils
import numpy as np
from copy import copy
from astropy import constants
import argparse
import matplotlib.pyplot as plt


def maxprob_model(sample_results, sps):  # GOTEM
    ### grab maximum probability, plus the thetas that gave it
    maxprob = sample_results['flatprob'].max()
    maxtheta = sample_results['flatchain'][sample_results['flatprob'].argmax()]

    ### ensure that maxprob stored is the same as calculated now
    current_maxprob = prosp_dutils.test_likelihood(sps, sample_results['model'], sample_results['obs'], maxtheta,
                                                   sample_results['run_params']['param_file'])

    print('Best-fit lnprob currently: {0}'.format(current_maxprob))
    print('Best-fit lnprob during sampling: {0}'.format(maxprob))

    return maxtheta, maxprob


def sample_flatchain(chain, lnprob, parnames, ir_priors=True, include_maxlnprob=True, nsamp=2000):  # GOTEM?
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


def set_sfh_time_vector(sample_results, ncalc):  # GOTEM
    # if parameterized, calculate linearly in 100 steps from t=0 to t=tage
    # if nonparameterized, calculate at bin edges.
    if 'tage' in sample_results['model'].theta_labels():
        nt = 100
        idx = np.array(sample_results['model'].theta_labels()) == 'tage'
        maxtime = np.max(sample_results['flatchain'][:ncalc, idx])
        t = np.linspace(0, maxtime, num=nt)
    elif 'agebins' in sample_results['model'].params:
        # print(sample_results['model'].params['agebins'], 'bins')  # PRINTER (bins[0] = (0, 8.))
        # print(len(sample_results['model'].params['agebins']), 'lenbins')  # PRINTER (lenbins=5)
        in_years = 10 ** sample_results['model'].params['agebins'] / 1e9
        # print(in_years, 'in_years')  # PRINTER  (inyears[0] = [[1e-9, 1e-1],[1e-1, 2.14515035e-1],... )
        # print(len(in_years), 'lenyears')  # PRINTER (lenyears = 5)
        t = np.concatenate((np.ravel(in_years) * 0.9999, np.ravel(in_years) * 1.001))
        t.sort()
        # print(t, 'concat, ravel, sort')  # PRINTER (1e-9 * 0.9999, 1e-9*1.001, 1e-1*0.9999, 1e-1*0.9999, 1e-1*1.001,)
        # print(len(t), 'len_t0')  # PRINTER (concat, ravel, sort takes each 2 things in all 5 bins, multiplies by
        # 0.9999 and 1.001, then appends them to one single list and sorts them in increasing order --> (len=20); for
        # the purpose of sampling on either side of the bin edges
        t = t[1:-1]  # remove points older than oldest bin, younger than youngest bin
        # print(t, 'oldest/youngest removed')  # PRINTER (removes first, last entries)
        # print(len(t), 'len_t1')  # PRINTER (len=18)
        t = np.clip(t, 1e-3, np.inf)  # nothing younger than 1 Myr!
        # print(t, 'clipped')  # PRINTER (smallest entry was 1e-9, now is 1e-3)
        # print(len(t), 'len_t2')  # PRINTER (len=18)
    else:
        sys.exit('ERROR: not sure how to set up the time array here!')
    return t


def calc_extra_quantities(sample_results, ncalc=3000, **kwargs):  # GOTEM??
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
        'measure_spectral_features': True,  # cost = 2 runtimes
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
    '''
    ### DOES THIS DO ANYTHING? # LOL NOPE
    sargv = sys.argv
    argdict = {'param_file': 'eelg_emission_params.py'}
    clargs = model_setup.parse_args(sargv, argdict=argdict)
    run_params = model_setup.get_run_params(argv=sargv, **clargs)
    sps = model_setup.load_sps(**run_params)
    ### DOES THIS DO ANYTHING?
    '''

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
        # t1 = time.time()

        ##### model call, to set parameters
        thetas = copy(sample_results['flatchain'][idx])
        spec[:, jj], mags[:, jj], sm = sample_results['model'].mean_model(thetas, sample_results['obs'], sps=sps)

        ##### extract sfh parameters
        # pass stellar mass to avoid extra model call
        sfh_params = prosp_dutils.find_sfh_params(sample_results['model'], thetas, sample_results['obs'], sps, sm=sm)
        # find sfh params: returns out: a dictionary containing 'sfr_fraction_full', 'mass_fraction', 'mformed', and
        #   'mass'; where 'mass_fraction' = 'sfr_fraction_full' * time_per_bin / 'sfr_fraction_full'.sum(); and where
        #   'mass_formed' = 'mass' / stellar_mass

        # print(t, 't')  # PRINTER (what t always is)
        # print(prosp_dutils.return_full_sfh(t, sfh_params), 'intsfr[:, jj]')  # PRINTER (array of all 1 value)

        ##### calculate SFH
        print(0.5)  # PRINTER
        print(prosp_dutils.return_full_sfh(t, sfh_params), 'what this do')
        intsfr[:, jj] = prosp_dutils.return_full_sfh(t, sfh_params)
        print(1)  # PRINTER
        # return_full_sfh: returns for i in len(tcalc): {where tcalc = t - np.max(10 ** sfh_params['agebins']) / 1e9
        #   tcalc = tcalc[tcalc < 0] * -1} intsfr[i] = calculate_sfr(sfh_params, deltat, tcalc=tcalc[i], **kwargs)
        # calculate_sfr: returns sfr = integrate_sfh(tcalc - timescale, tcalc, sfh_params) * sfh_params['mformed'].sum
        #   / (timescale * 1e9); and clips sfr to account for minsfr and maxsfr; sfr = units (?) * (mass OR none) / yr
        # integrate_sfh (for npSFH): returns tot_mformed = sum(weights / time_per_bin) * sfh_params['mass_fraction'];
        #   linearizes bins: to_linear_bins = 10**sfh_params['agebins'] / 1e9
        #   time_per_bin = to_linear_bins[:, 1] - to_linear_bins[:, 0]; time_bins = max(to_linear_bins) - to_linear_bins
        #   clips times outside SFH bins
        #   weights: initialized to same shape as sfh_params['mass_fraction']
        #   for bins inside t1,t2 boundaries: weights[i] = time_per_bin[i]
        #   for edge cases: weights[i] = t2 - time_bins[i, 1] or weights[i] = time_bins[i, 0] - t1
        #   calculates tot_mformed --> tot_mformed units: (time) / (time) * sfh_params['mass_fraction']
        #   --> calculate_sfr units: sfh_params['mass_fraction'] * sfh_params['mformed'] / time = none * mass / time?
        # print(len(intsfr[:, jj]), 'len')  # len 18 (column len is 18)
        # print(len(intsfr[jj, :]), 'lenny')  # len 2000 (row len is 2000)
        # --> there are 2000 columns
        '''
        print(intsfr[0, 0], '0, 0')
        print(intsfr[1, 0], '1, 0')
        print(intsfr[9, 0], '9, 0')
        # PRINTER the three above are ALWAYS the same 0.2365572
        print(intsfr[15, 0], '15, 0')
        print(intsfr[17, 0], '18, 0')
        # print(intsfr[:, 0], 'hello')  # same as 'hi-there'
        # print(intsfr[0, 1], '0, 1')
        # print(intsfr[9, 1], '9, 1')
        # PRINTER the two above are always the same
        '''

        ##### solve for half-mass assembly time
        '''  # PRINTER TEMPORARY COMMENT
        half_time[jj] = prosp_dutils.halfmass_assembly_time(sfh_params)
        '''

        ##### calculate time-averaged SFR
        sfr_10[jj] = prosp_dutils.calculate_sfr(sfh_params, 0.01, minsfr=-np.inf, maxsfr=np.inf)  # avg over 10 Myr
        print(2)
        sfr_100[jj] = prosp_dutils.calculate_sfr(sfh_params, 0.1, minsfr=-np.inf, maxsfr=np.inf)  # avg over 100 Myr
        print(3)

        '''
        # print(t, 't')  # PRINTER
        print(intsfr[:, jj], 'sfh')  # PRINTER
        print(sfr_100[jj], '100')  # PRINTER
        '''

        ##### calculate mass, sSFR
        stellar_mass[jj] = sfh_params['mass']
        ssfr_10[jj] = sfr_10[jj] / stellar_mass[jj]
        ssfr_100[jj] = sfr_100[jj] / stellar_mass[jj]

        loop += 1
        print('loop', loop)

    print(intsfr[:, 0], 'hi-there')
    print(t, 'tea')
    '''
    (array([ 0.23655724,  0.23655724,  0.23655724,  0.23655724,  0.23655724,
        0.23655724,  0.23655724,  0.23655724,  0.23655724,  0.23655724,
        0.23655724,  0.23655724,  0.23655724,  0.23655724,  0.23655724,
        0.23655724,  0.23655724,  0.        ]), 'hi-there')
    (array([  1.00000000e-03,   9.99900000e-02,   9.99900000e-02,
         1.00100000e-01,   1.00100000e-01,   2.14493583e-01,
         2.14493583e-01,   2.14729550e-01,   2.14729550e-01,
         4.60120984e-01,   4.60120984e-01,   4.60627168e-01,
         4.60627168e-01,   9.87028689e-01,   9.87028689e-01,
         9.88114529e-01,   9.88114529e-01,   2.11732493e+00]), 'tea')
    '''


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


# def post_processing(param_name, outname=None, **kwargs):
def post_processing(param_name, **kwargs):
    '''
    Driver. Loads output, runs post-processing routine.
    '''
    print(outname)

    sample_results, powell_results, model = read_results.results_from(outname)

    ### create flatchain, run post-processing
    sample_results['flatchain'] = prosp_dutils.chop_chain(sample_results['chain'], **sample_results['run_params'])
    sample_results['flatprob'] = prosp_dutils.chop_chain(sample_results['lnprobability'],
                                                         **sample_results['run_params'])
    extra_output = calc_extra_quantities(sample_results, **kwargs)

    print(len(extra_output['extras']['sfh'][0]), '2000?')  # YAY

    '''
    print('len flatchain', len(extra_output['extras']['flatchain']))  # 2000 :D
    print('index 0', extra_output['extras']['flatchain'][0])  # should return set of params (half time, sfr10, etc)?
        # Returns array with len 6: (half_time [units?], sfr10, sfr100, ssfr10, ssfr100, stellar_mass)?
    print('thing', extra_output['extras']['flatchain'])  # Returns 2000 arrays each with len(6)
    '''

    print(extra_output['bfit']['maxprob_params'])  # [10.34, 0.33, 0.59, 0.0023, 0.03, 0.0095, 0.69, -2, -0.2, -1]
    # emission_decoupled_recent: [10.32, 0.32, 0.61, 0.0246, 0.00155, 0.015, 0.7, -2, -0.2, -1]
    # 10 thetas are: logmass, sfr_frac1, sfr_frac2, sfr_frac3, sfr_frac4, sfr_frac5, dust2, logzsol, gas_logz, gas_logu
    # print('max', len(extra_output['bfit']['maxprob_params']))  # 10
    print(extra_output['bfit']['half_time'])  # -1.0216794961; emission_decoupled_recent: -1.03132084372
    print(extra_output['bfit']['sfr_10'])  # 66.4716596583; emission_decoupled_recent: 66.2740632613
    print(extra_output['bfit']['sfr_100'])  # 114.443681589; emission_decoupled_recent: 121.448368493

    # plt.plot(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'])
    # plt.plot(extra_output['extras']['t_sfh'], extra_output['extras']['sfh'])

    plt.plot(extra_output['extras']['t_sfh'], extra_output['bfit']['sfh'])
    plt.ylabel(r'M$_\odot$ yr$^{-1}$')
    plt.xlabel('Gyr')
    plt.show()


def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False


if __name__ == "__main__":

    ### don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('parfile', type=str)
    parser.add_argument('--outname')
    parser.add_argument('--measure_spectral_features', type=str2bool)
    parser.add_argument('--mags_nodust', type=str2bool)
    parser.add_argument('--ir_priors', type=str2bool)
    parser.add_argument('--ncalc', type=int)

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    #
    # out_file = "1824_decoupled_n1200_1495729287_mcmc.h5"
    # out_file = "1824_emission_decoupled_recent_1497037005_mcmc.h5"
    # out_file = "1824_emission_newagelims_1497289273_mcmc.h5"
    outname = kwargs['outname']
    '''
    res, pr, mod = read_results.results_from(outname)
    sample_results = res
    '''
    #

    print(kwargs)
    post_processing(outname, **kwargs)

'''
RUNNING WITH:

python only_sfh.py --outname=/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/1824_emission_newagelims_1497289273_mcmc.h5
    --ncalc=2000 --measure_spectral_features=True --mags_nodust=False --ir_priors=False parfile=eelg_emission_params.py

'''
