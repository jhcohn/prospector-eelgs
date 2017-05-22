import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np

filename = "1824_emission1_nompirun_1495123382_mcmc.h5"

objname = str(1824)  # objname = str(1824)  # = filename[0:6]; if objname[-1] == '_': objname = filename[0:5]

res, pr, mod = bread.results_from(filename)
print(res['chain'], 'chain')

# ''' #

# PRINT TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
tracefig, prob = bread.param_evol(res)  # print tracefig, store probability
plt.show()

# FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
# print(prob)
print('max', prob.max())
row = prob.argmax() / len(prob[0])
col = prob.argmax() - row * len(prob[0])
print('rc', row, col)
walk, iter = row, col
# print(walk, iter)

# ''' #

# PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
cornerfig = bread.subtriangle(res, start=400, thin=5, show_titles=True)
plt.show()
# For FAST: truths = [mass, age, tau, dust2] (for 1824: [9.78, 0.25, -1., 0.00])

# ''' #

# PRINT MODEL SED FOR GALAXY
# We need the correct sps object to generate models
### SPS ###
import sys
np.errstate(invalid='ignore')

from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

sargv = sys.argv
argdict = {'param_file': 'eelg_emission_params.py'}
clargs = model_setup.parse_args(sargv, argdict=argdict)
run_params = model_setup.get_run_params(argv=sargv, **clargs)
sps = model_setup.load_sps(**run_params)

# ''' #

# GET MODELED SPECTRA AND PHOTOMETRY
# These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
spec, phot, mfrac = mod.mean_model(res['chain'][walk, iter, :], obs=res['obs'], sps=sps)

# PLOT MODEL SED
wave = [f.wave_effective for f in res['obs']['filters']]

wave = np.asarray(wave)

print('len', len(sps.wavelengths), len(spec))

# PRINT SPECTRUM
print(sps.wavelengths)
plt.plot(sps.wavelengths, spec)
plt.xlabel('Wavelength [angstroms]')
plt.show()

mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
print(mean, 'mean')

yerr = res['obs']['maggies_unc']

# CHANGING OBSERVED TO REST FRAME WAVELENGTH
datname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'  # main catalog
zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'

with open(datname, 'r') as f:
    hdr = f.readline().split()
dtype = np.dtype([(hdr[1],'S20')] + [(n, np.float) for n in hdr[2:]])
dat = np.loadtxt(datname, comments = '#', delimiter=' ', dtype = dtype)

with open(zname, 'r') as fz:
    hdr_z = fz.readline().split()
dtype_z = np.dtype([(hdr_z[1],'S20')] + [(n, np.float) for n in hdr_z[2:]])
zout = np.loadtxt(zname, comments = '#', delimiter=' ', dtype = dtype_z)

idx = dat['id'] == objname
zred = zout['z_spec'][idx][0]  # GET REDSHIFT (z_spec, or if no z_spec exists, z_phot)
if zred == -99:
    zred = zout['z_peak'][idx][0]
print(zred)
wave_rest = []  # REST FRAME WAVELENGTH
for i in range(len(wave)):
    wave_rest.append(wave[i]/(1 + zred))

# PLOT MODEL SED BEST FIT, INPUT PHOT
plt.subplot(111, xscale="log", yscale="log")
plt.errorbar(wave_rest, res['obs']['maggies'], yerr=yerr, marker='o', linestyle='', color='b', label='Observed photometry')
plt.plot(wave_rest, phot, 'o', label='Model at {},{}'.format(walk, iter), color='r')
plt.legend(loc="best", fontsize=20)
plt.plot(sps.wavelengths, spec, color='b', alpha=0.5)
# plt.xlabel('Observed wavelength [angstroms]')
plt.xlabel('Rest frame wavelength [angstroms]')
plt.ylabel('Maggies')
plt.show()

# PLOT CHI_SQ BESTFIT
chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
plt.plot(wave_rest, chi_sq, 'o', color='b')
# plt.plot(np.log10(global_obs['wave_effective']), chi_sq, 'o', color='b')
plt.xlabel('Rest frame wavelength [angstroms]')
plt.ylabel(r'$\chi^2$')
plt.show()
# ''' #
