import numpy as np
from scipy.optimize import brentq
from prosp_dutils import integrate_sfh, chop_chain, find_sfh_params
# from prospector_io import load_prospector_data
import prospect.io.read_results as bread
from prospect.models import model_setup
import os
import copy
# import sys


def sfh_half_time(x, sfh_params, c):
    """
    wrapper for use with halfmass assembly time
    """
    # check for nonparametric
    sf_start = sfh_params['sf_start']
    if sfh_params['sf_start'].shape[0] == 0:
        sf_start = 0.0
    return integrate_sfh(sf_start, x, sfh_params) - c


def halfmass_assembly_time(sfh_params, c=0.5):
    try:
        half_time = brentq(sfh_half_time, 0, 14, args=(sfh_params, c), rtol=1.48e-08, maxiter=1000)
    except ValueError:
        # big problem
        warnings.warn("You've passed SFH parameters that don't allow t_half to be calculated. Check for bugs.",
                      UserWarning)
        half_time = np.nan

    # define age of galaxy
    tgal = sfh_params['tage']
    if tgal.shape[0] == 0:
        tgal = np.max(10 ** sfh_params['agebins'] / 1e9)

    return tgal - half_time


# DOES SETUP
def sample_posterior(outname=None, shortname=None, mass_folder=None):
    # I/O
    # paramfile = model_setup.import_module_from_file(param_name)
    # outname = paramfile.run_params['outfile']
    # sample_results, powell_results, model, eout = load_prospector_data(outname, hdf5=True, load_extra_output=True)
    sample_results, powell_results, model = bread.results_from(outname)

    # create useful quantities
    sample_results['flatchain'] = chop_chain(sample_results['chain'], **sample_results['run_params'])
    sample_results['flatprob'] = chop_chain(sample_results['lnprobability'], **sample_results['run_params'])

    sps = model_setup.load_sps(**sample_results['run_params'])
    # sps = paramfile.load_sps(**sample_results['run_params'])
    # obs = paramfile.load_obs(**sample_results['run_params'])

    # sample from posterior
    nsamp = 3000
    good = np.isfinite(sample_results['flatprob']) == True
    sample_idx = np.random.choice(np.where(good)[0], nsamp)

    # define outputs
    mfrac = np.linspace(0, 0.95, 20)
    mfrac_out = np.zeros(shape=(nsamp, mfrac.shape[0]))
    for jj, idx in enumerate(sample_idx):
        print(jj)
        ##### model call, to set parameters
        thetas = copy.copy(sample_results['flatchain'][idx])
        spec, mags, sm = sample_results['model'].mean_model(thetas, sample_results['obs'], sps=sps)

        ##### extract sfh parameters
        sfh_params = find_sfh_params(sample_results['model'], thetas,
                                     sample_results['obs'], sps, sm=sm)

        for ii, m in enumerate(mfrac):
            mfrac_out[jj, ii] = halfmass_assembly_time(sfh_params, c=m)

    # fixing negatives
    mfrac_out = np.clip(mfrac_out, 0.0, np.inf)
    # write out
    out = np.percentile(mfrac_out, [50, 84, 16], axis=0)
    with open('out/' + mass_folder + shortname + 'mass.txt', 'w') as f:
        f.write('# mass_fraction median_time err_up err_down\n')
        for ii in range(out.shape[1]):
            f.write("{:.2f}".format(mfrac[ii]) + ' ' + "{:.3f}".format(out[0, ii]) + ' ' + "{:.3f}".format(
                out[1, ii]) + ' ' + "{:.3f}".format(out[2, ii]) + ' ')
            f.write('\n')

'''
base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/out_noelg/' # out_fixedmet/
out = base + '7065_uds_noelg_1506704017_mcmc.h5'  # '16067_cosmos_fixedmet_1506462907_mcmc.h5'
print(out)
sample_posterior(outname=out, shortname='7065_uds_noelg_', mass_folder='mass_noelg/')
'''

folders = [['out_fixedmet/', 'mass_fixedmet/'], ['out_noelg/', 'mass_noelg/'],
           ['out_eth/', 'mass_eth/'], ['out_nth/', 'mass_nth/']]

for set in folders:

    outs = []
    simples = []
    base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + set[0]
    for filename in os.listdir(base):
        if filename.endswith("_mcmc.h5"):
            outs.append(os.path.join(base, filename))
            simple = ''
            count = 0

            for i in filename:
                if i == '_':
                    count += 1
                if count <= 2:
                    simple += i
            simples.append(simple)
        else:
            pass

    print(simples)  # prints ID_field_paramfile

    for i in range(len(outs)):
        print(set[0], 'galaxy', i)
        sample_posterior(outname=outs[i], shortname=simples[i], mass_folder=set[1])

# '''
