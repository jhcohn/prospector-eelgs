import numpy as np
import sys
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

sargv = sys.argv
argdict = {'param_file': 'eelg_addmass_params.py'}
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

mu, phot, x = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)
errs = np.log10(global_obs['maggies'])-np.log10(global_obs['maggies']-global_obs['maggies_unc'])

# print(global_model.initial_theta, 'model')  # check initial theta
# print(global_model.theta_labels(), 'labels')  # check theta labels

y = global_obs['maggies']
ujy = y * 3631 * 10 ** 6
yerr = 3631 * 10**6 * global_obs['maggies_unc']

phot *= 3631 * 10 ** 6  # uJy

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
plt.errorbar(global_obs['wave_effective'], ujy, yerr=yerr, marker='o', linestyle='', color='b')
plt.plot(global_obs['wave_effective'], phot, 'o', color='r')
# label='Model at {},{}'.format(walker, iteration), color='r')
plt.show()

# print(global_obs['wave_effective'], 'pop2')
# too_blue = np.log10(global_obs['wave_effective']) < 3.5
# print(too_blue, 'pop3')

'''
# PRINT INITIAL GUESS CHI_SQ
chi_sq = ((global_obs['maggies'] - phot) / global_obs['maggies_unc']) ** 2
plt.plot(global_obs['wave_effective'], chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()
'''

def get_names(field):
    photname = None
    zname = None

    if field == 'cosmos':
        photname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'

    elif field == 'cdfs':
        photname = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'

    elif field == 'uds':
        photname = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.awk.zout'

    return photname, zname


photname, zname = get_names('cosmos')  # for 12105

# OPEN FILE, LOAD DATA
with open(photname, 'r') as f:
    hdr = f.readline().split()
dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
dat = np.loadtxt(photname, comments='#', delimiter=' ', dtype=dtype)

# EXTRACT FILTERS, FLUXES, ERRORS FOR OBJECT
obj_idx = (dat['id'] == '12105')  # for 12105
cos_filternames = ['B', 'G', 'I', 'IA427', 'IA484', 'IA505', 'IA527', 'IA624', 'IA709', 'IA738', 'R', 'U', 'V', 'Rp',
                   'Z', 'Zp', 'Hl', 'Hs', 'J1', 'J2', 'J3', 'Ks', 'NB118', 'NB209', 'F125W', 'F140W', 'F160W', 'F606W',
                   'F814W', 'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'UVISTA_Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']

flux = np.squeeze([dat[obj_idx]['f_' + f] for f in cos_filternames])
print(len(flux), len(phot))
tot = []
diff = []
for j in range(len(flux)):
    tot.append((float(flux[j]) * (10 ** - 0.44)) + phot[j])
    diff.append(((float(flux[j]) * 10 ** -0.44) - phot[j]) / (float(flux[j]) * 10 ** -0.44))

print(diff)
newphot = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/newphot_21442'
with open(newphot, 'w+') as new:
    new.write('# id ')
    for i in range(len(cos_filternames)):
        new.write('f_' + cos_filternames[i] + ' ')
    new.write('\n')
    new.write('12105 ')
    for k in range(len(tot)):
        new.write(str(tot[k]) + ' ')
    new.write('\n')

phot2 = []
with open(newphot, 'r') as newp:
    for line in newp:
        if line[0] != '#':
            cols = line.split()
            for l in range(len(cols)):
                if l == 0:
                    pass
                else:
                    phot2.append(cols[l])
print(np.squeeze(phot2))
