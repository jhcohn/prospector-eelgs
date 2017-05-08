import numpy as np
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

# -------------
# ADDED LINES FOR FAST TEST
import matplotlib.pyplot as plt
mu, phot, x = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)
errs = np.log10(global_obs['maggies'])-np.log10(global_obs['maggies']-global_obs['maggies_unc'])

# print(global_model.initial_theta, 'model')  # check initial theta
# print(global_model.theta_labels(), 'labels')  # check theta labels

y = global_obs['maggies']
yerr = global_obs['maggies_unc']
# print(len(y), 'y')
# print(len(yerr), 'yerr')
# print(len(phot), 'phot')
# print(len(global_obs['wave_effective']), 'wavel')

# CALCULATE MEAN DIFFERENCE BETWEEN INPUT PHOTOMETRY AND INITIAL GUESS MODEL
mean = ((global_obs['maggies']-phot)/global_obs['maggies']).mean()
ratio = []
for i in range(len(global_obs['maggies'])):
    ratio.append((global_obs['maggies'] / phot))
print(ratio, 'ratio')
print(mean, 'mean')

# PRINT INPUT PHOTOMETRY SED AND INITIAL GUESS MODEL SED
plt.subplot(111, xscale="log", yscale="log")
plt.errorbar(global_obs['wave_effective'], y, yerr=yerr, marker='o', linestyle='', color='b')
plt.plot(global_obs['wave_effective'], phot, 'o', color='r')
# label='Model at {},{}'.format(walker, iteration), color='r')
plt.show()

# print(global_obs['wave_effective'], 'pop2')
# too_blue = np.log10(global_obs['wave_effective']) < 3.5
# print(too_blue, 'pop3')

# PRINT INITIAL GUESS CHI_SQ
chi_sq = ((global_obs['maggies'] - phot) / global_obs['maggies_unc']) ** 2
plt.plot(global_obs['wave_effective'], chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()
