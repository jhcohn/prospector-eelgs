import numpy as np
import sys
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

obj = 11063  # 15124  # 20366  # 12533  # 11462  # 12552  # 12105  # 21442
field = 'cosmos'  # 'cdfs'  # 'cosmos'  # 'cdfs'
sfh_style = '15_05_17'  # =0.15 mfrac in most recent bin, 0.05 mfrac in second most recent bin, logzsol of -1.7

sargv = sys.argv
argdict = {'param_file': 'eelg_sfhtest_params.py'}
clargs = model_setup.parse_args(sargv, argdict=argdict)
run_params = model_setup.get_run_params(argv=sargv, **clargs)
print('runpars')
print(run_params)

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

print(global_model.initial_theta)
print(global_obs)

# -------------
# ADDED LINES FOR SFH TEST
import matplotlib.pyplot as plt

mu, phot, x = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)
errs = np.log10(global_obs['maggies'])-np.log10(global_obs['maggies']-global_obs['maggies_unc'])
# WHY IS INITIAL GUESS MODEL NOT CHANGING?!

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
# mean = ((global_obs['maggies']-phot)/global_obs['maggies']).mean()
# ratio = []
# for i in range(len(global_obs['maggies'])):
#     ratio.append((global_obs['maggies'] / phot))
# print(ratio, 'ratio')
# print(mean, 'mean')

# PRINT INPUT PHOTOMETRY SED AND INITIAL GUESS MODEL SED
plt.subplot(111, xscale="log", yscale="log")
plt.errorbar(global_obs['wave_effective'], ujy, yerr=yerr, marker='o', linestyle='', color='b')
plt.plot(global_obs['wave_effective'], phot, 'o', color='r')
# label='Model at {},{}'.format(walker, iteration), color='r')
plt.show()


def get_names(field):
    photname = None
    zname = None
    filters = None

    if field == 'cosmos':
        photname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'
        filters = ['B', 'G', 'I', 'IA427', 'IA484', 'IA505', 'IA527', 'IA624', 'IA709', 'IA738', 'R', 'U', 'V', 'Rp',
                   'Z', 'Zp', 'Hl', 'Hs', 'J1', 'J2', 'J3', 'Ks', 'NB118', 'NB209', 'F125W', 'F140W', 'F160W', 'F606W',
                   'F814W', 'UVISTA_J', 'UVISTA_H', 'UVISTA_Ks', 'UVISTA_Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']
    elif field == 'cdfs':
        photname = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
        filters = ['B', 'I', 'R', 'U', 'V', 'Z', 'Hs', 'Hl', 'J1', 'J2', 'J3', 'Ks', 'KsHI', 'NB118', 'NB209', 'F098M',
                   'F105W', 'F125W', 'F140W', 'F160W', 'F814W', 'IA484', 'IA527', 'IA574', 'IA598', 'IA624', 'IA651',
                   'IA679', 'IA738', 'IA767', 'IA797', 'IA856', 'WFI_V', 'WFI_Rc', 'WFI_U38', 'tenisK', 'IRAC_36',
                   'IRAC_45', 'IRAC_58', 'IRAC_80']
    elif field == 'uds':
        photname = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.awk.zout'
        filters = ['u', 'B', 'V', 'R', 'i', 'z', 'J1', 'J2', 'J3', 'Hs', 'Hl', 'Ks', 'J', 'H', 'K', 'KsHI', 'F125W',
                   'F140W', 'F160W', 'F606W', 'F814W', 'Y', 'IRAC_36', 'IRAC_45', 'IRAC_58', 'IRAC_80']

    return photname, zname, filters


photname, zname, filters = get_names(field)  # FIELD

# OPEN FILE, LOAD DATA
with open(photname, 'r') as f:
    hdr = f.readline().split()
dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
dat = np.loadtxt(photname, comments='#', delimiter=' ', dtype=dtype)

# EXTRACT FILTERS, FLUXES, ERRORS FOR OBJECT
obj_idx = (dat['id'] == str(obj))  # OBJID

# flux = np.squeeze([dat[obj_idx]['f_' + f] for f in filters])
# print(len(flux), len(phot))

# ASSUME UNCERTAINTIES SIMILAR TO TYPICAL UNCERTAINTIES
unc = np.squeeze([dat[obj_idx]['e_' + f] for f in filternames])
# tot = []
# diff = []
for j in range(len(phot)):
    tot.append(phot[j])
    # diff.append(((float(flux[j]) * 10 ** -0.44) - phot[j]) / (float(flux[j]) * 10 ** -0.44))

# print(diff)
newphot = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_' + sfh_style  # str(obj)
with open(newphot, 'w+') as new:
    new.write('# id ')
    for i in range(len(filters)):
        new.write('f_' + filters[i] + ' ')
    new.write('\n')
    new.write(sfh_style + ' ')  # (str(obj) + ' ')
    for k in range(len(tot)):
        new.write(str(tot[k]) + ' ')
    new.write('\n')

print('tot', tot)
print('generated phot', phot)  # should be identical to tot
# print('cat flux', flux)

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
print(np.squeeze(phot2))  # should be identical to phot and tot?
