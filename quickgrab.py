import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np

# filename = "mpirun_fixedconvergence_1824_1494274211_mcmc.h5"
filename = "mpirun_cvgtest_1824_1494290102_mcmc.h5"  # 101, 102, 103

res, pr, mod = bread.results_from(filename)

print(res['chain'], 'chain')

# ''' #

# PRINT TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
tracefig, prob = bread.param_evol(res)  # print tracefig, store probability
plt.show()

# FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
# print(prob)
print('max', prob.max())
for i in range(len(prob)):
    for j in range(len(prob[i])):
        if prob[i][j] == prob.max():
            print('max', prob[i][j], i, j)  # i, j = walker, iteration
            walk, iter = i, j  # this will re-write if no single global max exists, and return just the last max
        # if prob[i][j] == prob.min():
        #     print('min', prob[i][j], i, j)
print(len(prob))

# ''' #

# PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
cornerfig = bread.subtriangle(res, start=400, thin=5, show_titles=True)  # add truths back in
plt.show()
# For FAST: truths = [mass, age, tau, dust2] (for 1824: [9.78, 0.25, -1., 0.00])

# ''' #

# PRINT MODEL SED FOR GALAXY
# We need the correct sps object to generate models
from prospect.sources import FastStepBasis  # NEW FastStepBasis, changed from CSPBasis
sps = FastStepBasis(**res['run_params'])  # NEW FastSepBasis

walker, iteration = walk, iter  # BUCKET switch to using: walk, iter variables
# Get the modeled spectra and photometry.
# These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
spec, phot, mfrac = mod.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)
# Plot the model SED
import matplotlib.pyplot as pl
wave = [f.wave_effective for f in res['obs']['filters']]

wave = np.asarray(wave)

mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
print(mean, 'mean')

yerr = res['obs']['maggies_unc']

# PLOT MODEL SED BEST FIT, INPUT PHOT
pl.subplot(111, xscale="log", yscale="log")
pl.errorbar(wave, res['obs']['maggies'], yerr=yerr, marker='o', linestyle='', color='b', label='Observed photometry')
pl.plot(wave, phot, 'o', label='Model at {},{}'.format(walker, iteration), color='r')
pl.legend(loc='best', fontsize=20)
pl.xlabel('Wavelength')
pl.ylabel('Maggies')
pl.show()

# PLOT CHI_SQ BESTFIT
chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
plt.plot(wave, chi_sq, 'o', color='b')
# plt.plot(np.log10(global_obs['wave_effective']), chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()
# ''' #
