import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np

filename = "quicktest_1493325578_mcmc.h5"
'''
Testing 'chain' keyword error
print(filename.split('.')[-1] == 'h5', 'true')
if filename.split('.')[-1] == 'h5':
    print('yep')
    res = bread.read_hdf5(filename)

print(res)
'''

res, pr, mod = bread.results_from(filename)

id = 1824

# ''' #

# FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
tracefig, prob = bread.param_evol(res)
print(prob)
print('max', prob.max())
for i in range(len(prob)):
    for j in range(len(prob[i])):
        if prob[i][j] == prob.max():
            print('max', prob[i][j], i, j)
        if prob[i][j] == prob.min():
            print('min', prob[i][j], i, j)  # i, j = walker, iteration
print(len(prob))

# ''' #

# PRINTS TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
# THEN PRINTS CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
# tracefig = bread.param_evol(res)
plt.show()
# cornerfig = bread.subtriangle(res, start=0, thin=5)
if id == 1824:
    cornerfig = bread.subtriangle(res, start=400, thin=5, show_titles=True)  # add truths back in

# For FAST: truths = mass, age, tau, dust2 (for 1824: [9.78, 0.25, -1., 0.00])
plt.show()

# ''' #

# PRINTS MODEL SED FOR OBJECT
# We need the correct sps object to generate models
from prospect.sources import CSPBasis  # NEW FastStepBasis
sps = CSPBasis(**res['run_params'])  # NEW FastSepBasis

walker, iteration = 96, 493
# Get the modeled spectra and photometry.
# These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
spec, phot, mfrac = mod.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)
# Plot the model SED
import matplotlib.pyplot as pl
wave = [f.wave_effective for f in res['obs']['filters']]

wave = np.asarray(wave)

mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
print(mean)

y = res['obs']['maggies']
yerr = res['obs']['maggies_unc']

pl.subplot(111, xscale="log", yscale="log")
pl.errorbar(wave, y, yerr=yerr, marker='o', linestyle='', color='b')
pl.plot(wave, phot, 'o', label='Model at {},{}'.format(walker, iteration), color='r')
pl.show()

pl.loglog(wave, res['obs']['maggies'], 'o', label='Observations')  # '-o'
pl.loglog(wave, phot, 'o', label='Model at {},{}'.format(walker, iteration))  # '-o'
pl.loglog(wave, res['obs']['maggies_unc'], 'o', label='Uncertainties')
pl.ylabel("Maggies")
pl.show()


chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
plt.plot(wave, chi_sq, 'o', color='b')
# plt.plot(np.log10(global_obs['wave_effective']), chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()
# '''