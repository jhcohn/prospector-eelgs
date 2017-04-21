#!/usr/local/bin/python

import time, sys, os
import numpy as np
np.errstate(invalid='ignore')

from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log


# --------------
# Read command line arguments
# --------------
sargv = sys.argv
argdict = {'param_file': ''}
clargs = model_setup.parse_args(sargv, argdict=argdict)
run_params = model_setup.get_run_params(argv=sargv, **clargs)

# --------------
# Globals
# --------------
# SPS Model instance as global
sps = model_setup.load_sps(**run_params)
# GP instances as global
spec_noise, phot_noise = model_setup.load_gp(**run_params)
# Model as global
global_model = model_setup.load_model(**run_params)
# Obs as global
global_obs = model_setup.load_obs(**run_params)

# -------------

# ADDED LINES FOR FAST TEST
import matplotlib.pyplot as plt
mu, phot, x = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)  # don't do chi^2 from initial guess
errs = np.log10(global_obs['maggies'])-np.log10(global_obs['maggies']-global_obs['maggies_unc'])

print(global_model.initial_theta, 'model')  # check initial theta
print(global_model.theta_labels(), 'labels')  # check theta labels

y = global_obs['maggies']
yerr = global_obs['maggies_unc']
# print(len(y), 'y')
# print(len(yerr), 'yerr')
# print(len(phot), 'phot')
# print(len(global_obs['wave_effective']), 'wavel')

mean = ((global_obs['maggies']-phot)/global_obs['maggies']).mean()
ratio = []
for i in range(len(global_obs['maggies'])):
    ratio.append((global_obs['maggies'] / phot))
# print(ratio, 'ratio')
# print(mean, 'pop')
plt.subplot(111, xscale="log", yscale="log")
plt.errorbar(global_obs['wave_effective'], y, yerr=yerr, marker='o', linestyle='', color='b')
plt.plot(global_obs['wave_effective'], phot, 'o', color='r')  # label='Model at {},{}'.format(walker, iteration), color='r')
plt.show()

# print(global_obs['wave_effective'], 'pop2')
# too_blue = np.log10(global_obs['wave_effective']) < 3.5
# print(too_blue, 'pop3')

chi_sq = ((global_obs['maggies'] - phot) / global_obs['maggies_unc']) ** 2
plt.plot(global_obs['wave_effective'], chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()

'''
plt.plot(global_obs['wavelength'], chi_sq, 'o', color='b')
# plt.plot(np.log10(global_obs['wave_effective']), chi_sq, 'o', color='b')
plt.show()
'''

'''
# wave effective:
array([  4335.6664887 ,   4810.70051748,   7614.53955113,   4259.01530918,
         4846.20888187,   5059.78255171,   5258.12330727,   6229.33890863,
         7072.07696982,   7358.07099729,   6223.2124034 ,   3787.81454827,
         5424.0997122 ,   6480.44749396,   8787.97811387,   9077.42082014,
        17023.22462635,  15509.19183943,  10538.94027717,  11430.97032625,
        12862.76187733,  21414.83623539,   1190.87503643,   2098.91357473,
        12425.7548425 ,  13830.70557458,  15324.46877726,   5813.0166945 ,
         7972.96049318,   1256.76676717,   1639.73547792,   2139.5815672 ,
         1027.19474025,  35378.51012199,  44780.55099406,  56962.1965764 ,
        77980.20467193])
'''

'''
plt.errorbar(np.log10(global_obs['wave_effective']), np.log10(global_obs['maggies']),yerr=errs,fmt='o',color='red')
plt.errorbar(np.log10(global_obs['wave_effective']), np.log10(phot),fmt='o',color='blue')

plt.show()
'''

# -------------

# -----------------
# LnP function as global
# ------------------

def lnprobfn(theta, model=None, obs=None, verbose=run_params['verbose']):
    """Given a parameter vector and optionally a dictionary of observational
    data and a model object, return the ln of the posterior. This requires that
    an sps object (and if using spectra and gaussian processes, a GP object) be
    instantiated.

    :param theta:
        Input parameter vector, ndarray of shape (ndim,)

    :param model:
        bsfh.sedmodel model object, with attributes including ``params``, a
        dictionary of model parameters.  It must also have ``prior_product()``,
        and ``mean_model()`` methods defined.

    :param obs:
        A dictionary of observational data.  The keys should be
          *``wavelength``
          *``spectrum``
          *``unc``
          *``maggies``
          *``maggies_unc``
          *``filters``
          * and optional spectroscopic ``mask`` and ``phot_mask``.

    :returns lnp:
        Ln posterior probability.
    """
    if model is None:
        model = global_model
    if obs is None:
        obs = global_obs

    lnp_prior = model.prior_product(theta)
    if np.isfinite(lnp_prior):
        # Generate mean model
        t1 = time.time()
        try:
            mu, phot, x = model.mean_model(theta, obs, sps=sps)
        except(ValueError):
            return -np.infty
        d1 = time.time() - t1

        # Noise modeling
        if spec_noise is not None:
            spec_noise.update(**model.params)
        if phot_noise is not None:
            phot_noise.update(**model.params)
        vectors = {'spec': mu, 'unc': obs['unc'],
                   'sed': model._spec, 'cal': model._speccal,
                   'phot': phot, 'maggies_unc': obs['maggies_unc']}

        # Calculate likelihoods
        t2 = time.time()
        lnp_spec = lnlike_spec(mu, obs=obs, spec_noise=spec_noise, **vectors)
        lnp_phot = lnlike_phot(phot, obs=obs, phot_noise=phot_noise, **vectors)
        d2 = time.time() - t2
        if verbose:
            write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

        return lnp_prior + lnp_phot + lnp_spec
    else:
        return -np.infty


def chisqfn(theta, model, obs):
    """Negative of lnprobfn for minimization, and also handles passing in
    keyword arguments which can only be postional arguments when using scipy
    minimize.
    """
    return -lnprobfn(theta, model=model, obs=obs)

# -----------------
# MPI pool.  This must be done *after* lnprob and
# chi2 are defined since slaves will only see up to
# sys.exit()
# ------------------
try:
    from emcee.utils import MPIPool
    pool = MPIPool(debug=False, loadbalance=True)
    if not pool.is_master():
        # Wait for instructions from the master process.
        pool.wait()
        sys.exit(0)
except(ImportError, ValueError):
    pool = None
    print('Not using MPI')


def halt(message):
    """Exit, closing pool safely.
    """
    print(message)
    try:
        pool.close()
    except:
        pass
    sys.exit(0)

# --------------
# Master branch
# --------------

if __name__ == "__main__":

    # --------------
    # Setup
    # --------------
    rp = run_params
    rp['sys.argv'] = sys.argv
    # Use the globals
    model = global_model
    obsdat = global_obs
    chi2args = [None, None]
    postkwargs = {}

    # make zeros into tiny numbers
    initial_theta = model.rectify_theta(model.initial_theta)
    if rp.get('debug', False):
        halt('stopping for debug')

    # Try to set up an HDF5 file and write basic info to it
    outroot = "{0}_{1}".format(rp['outfile'], int(time.time()))
    try:
        hfilename = outroot + '_mcmc.h5'
        hfile = h5py.File(hfilename, "a")
        print("Writing to file {}".format(hfilename))
        write_results.write_h5_header(hfile, run_params, model)
        write_results.write_obs_to_h5(hfile, obs)
    except:
        hfile = None
        
    # -----------------------------------------
    # Initial guesses using powell minimization
    # -----------------------------------------
    if bool(rp.get('do_powell', True)):
        if rp['verbose']:
            print('minimizing chi-square...')
        ts = time.time()
        powell_opt = {'ftol': rp['ftol'], 'xtol': 1e-6, 'maxfev': rp['maxfev']}
        powell_guesses, pinit = fitting.pminimize(chisqfn, initial_theta,
                                                  args=chi2args, model=model,
                                                  method='powell', opts=powell_opt,
                                                  pool=pool, nthreads=rp.get('nthreads', 1))
        best = np.argmin([p.fun for p in powell_guesses])
        initial_center = fitting.reinitialize(powell_guesses[best].x, model,
                                              edge_trunc=rp.get('edge_trunc', 0.1))
        initial_prob = -1 * powell_guesses[best]['fun']
        pdur = time.time() - ts
        if rp['verbose']:
            print('done Powell in {0}s'.format(pdur))
            print('best Powell guess:{0}'.format(initial_center))
    else:
        powell_guesses = None
        pdur = 0.0
        initial_center = initial_theta.copy()
        initial_prob = None

    # -------
    # Sample
    # -------
    if rp['verbose']:
        print('emcee sampling...')
    tstart = time.time()
    out = fitting.run_emcee_sampler(lnprobfn, initial_center, model,
                                    postkwargs=postkwargs, initial_prob=initial_prob,
                                    pool=pool, hdf5=hfile, **rp)
    esampler, burn_p0, burn_prob0 = out
    edur = time.time() - tstart
    if rp['verbose']:
        print('done emcee in {0}s'.format(edur))

    # -------------------------
    # Output pickles (and HDF5)
    # -------------------------
    write_results.write_pickles(rp, model, obsdat, esampler, powell_guesses,
                                outroot=outroot, toptimize=pdur, tsample=edur,
                                sampling_initial_center=initial_center,
                                post_burnin_center=burn_p0,
                                post_burnin_prob=burn_prob0)
    if hfile is None:
        hfile = hfilename

    write_results.write_hdf5(hfile, rp, model, obsdat, esampler, powell_guesses,
                             toptimize=pdur, tsample=edur,
                             sampling_initial_center=initial_center,
                             post_burnin_center=burn_p0,
                             post_burnin_prob=burn_prob0)

    halt('Finished')
