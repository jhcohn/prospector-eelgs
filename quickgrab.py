import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np
# import prosp_dutils  # for SFH --> SFR(t)
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log
import sys
np.errstate(invalid='ignore')

out_file = "1824_logzsol_0_1496245637_mcmc.h5"
param_file = 'eelg_emission_params.py'

objname = ''
for i in out_file:
    if i == '_':
        break
    objname += i  # object id is named at the start of each out_file

res, pr, mod = bread.results_from(out_file)
# print(res['chain'], 'chain')

# ''' #

# PRINT TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
tracefig, prob = bread.param_evol(res)  # print tracefig, store probability
plt.show()

# FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
# print(prob)
print('max', prob.max())
row = prob.argmax() / len(prob[0])
col = prob.argmax() - row * len(prob[0])
walker, iteration = row, col
print(walker, iteration)

# ''' #

# PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
cornerfig = bread.subtriangle(res, start=400, thin=5, show_titles=True)
plt.show()
# For FAST: truths = [mass, age, tau, dust2] (for 1824: [9.78, 0.25, -1., 0.00])

# ''' #

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

mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
print(mean, 'mean')  # print normalized mean difference between model and observations

# PLOT SPECTRUM
wave = [f.wave_effective for f in res['obs']['filters']]
wave = np.asarray(wave)
print('len', len(sps.wavelengths), len(spec))

plt.plot(sps.wavelengths, spec)
plt.xlabel('Wavelength [angstroms]')
plt.show()


# '''
# HOW CONVERGED IS THE CODE?? LET'S FIND OUT!
parnames = np.array(res['model'].theta_labels())
fig, kl_ax = plt.subplots(1, 1, figsize=(7, 7))
for i in xrange(parnames.shape[0]):
    kl_ax.plot(res['kl_iteration'], np.log10(res['kl_divergence'][:, i]),
               'o', label=parnames[i], lw=1.5, linestyle='-', alpha=0.6)

kl_ax.set_ylabel('log(KL divergence)')
kl_ax.set_xlabel('iteration')
# kl_ax.set_xlim(0, nsteps*1.1)

kl_div_lim = res['run_params'].get('convergence_kl_threshold', 0.018)
kl_ax.axhline(np.log10(kl_div_lim), linestyle='--', color='red', lw=2, zorder=1)

kl_ax.legend(prop={'size': 5}, ncol=2, numpoints=1, markerscale=0.7)
plt.show()
# '''

# CHANGING OBSERVED TO REST FRAME WAVELENGTH
datname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'  # main catalog
zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog

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
print(zred)

wave_rest = []  # REST FRAME WAVELENGTH
for i in range(len(wave)):
    wave_rest.append(wave[i]/(1 + zred))  # 1 + z = l_obs / l_emit --> l_emit = l_obs / (1 + z)


# PLOT MODEL SED BEST FIT, INPUT PHOT
yerr = res['obs']['maggies_unc']

plt.subplot(111, xscale="log", yscale="log")
plt.errorbar(wave_rest, res['obs']['maggies'], yerr=yerr, marker='o', linestyle='', color='b',
             label='Observed photometry')
plt.plot(wave_rest, phot, 'o', label='Model at {},{}'.format(walker, iteration), color='r')
plt.legend(loc="best", fontsize=20)
plt.plot(sps.wavelengths, spec, color='b', alpha=0.5)
plt.xlabel('Rest frame wavelength [angstroms]')
plt.ylabel('Maggies')
plt.show()


# PLOT CHI_SQ BESTFIT
chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
plt.plot(wave_rest, chi_sq, 'o', color='b')
plt.xlabel('Rest frame wavelength [angstroms]')
plt.ylabel(r'$\chi^2$')
plt.show()
# ''' #
