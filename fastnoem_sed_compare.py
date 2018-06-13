import os
import make_all_plots as map
import pickle
import numpy as np
import print_sfh
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
import stellar_ages as sa
import get_mass_dust as gmd

'''
home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
ff = 'fast_fits/'
cdfs_04 = ff + 'cdfs_Z004_EL/'
cdfs_02 = ff + 'cdfs_Z02_EL/'
uds_04 = ff + 'uds_Z004_EL/'
uds_02 = ff + 'uds_Z02_EL/'
cosmos_04 = ff + 'cosmos_Z004_EL/'
cosmos_02 = ff + 'cosmos_Z02_EL/'
cdfs_pre = 'cdfs.v1.6.9_'
cosmos_pre = 'cosmos.v1.3.6_'
uds_pre = 'uds.v1.5.8_'

eelgs = 1  # True=1
# folders = ['pkl_efifty/', 'pkl_nvary/']  # ['pkl_efifty', 'pkl_nvary']  # ['pkl_evar', 'pkl_nvary']
# folders = ['evar2_pkls/', 'pkl_ecorr/']
# folders = ['pkl_ecorr/', 'pkl_ncorr/']
# folder = 'pkl_efifty/'
folder = 'pkl_efico/'
# out_fold = 'out_efifty/'
out_fold = 'out_efico/'
# file = 'fifty'  # 'fifty'  # 'vary'  # 'newmask'
file = 'fico'
# base = ['fifty', 'vary']  # fifty, vary
base = ['fico', 'fico']
fs = 20  # 30
fs_text = 30
textx = 1100
texty = 150  # 500

objs = []
fields = []
count = 0

eelg_objs, eelg_fields, lbg_objs, lbg_fields = sa.get_gal_lists(base, objlists=True)
full_set = [[eelg_objs, eelg_fields], [lbg_objs, lbg_fields]]
'''

'''
for place in [cdfs_04, cosmos_04, uds_04]:
    for file in os.listdir(place):
        for i in range(len(eelg_objs)):

for i in range(len(eelg_objs)):
    print(eelg_objs[i])
'''
obj = 20366
field = 'cdfs'
reg = 'fico'
fast = 'fastnoem'
pkls = ['pkl_efico/', 'pkl_efastnoem/']
compares = [reg, fast]
colors = ['purple', 'r']
labels = [r'Prospector', r'FAST']

fig = plt.figure(1, figsize=(12, 12))
gs = gridspec.GridSpec(1, 1)  # , height_ratios=[3, 1])
#    gs = gridspec.GridSpec(1, 2)
# prepares fig to hold two sets side-by-side of: two axes, one atop other, size ratio 3:1
ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis
# ax3 = plt.subplot(gs[1], sharex=ax1)  # ax3 = smaller lower axis
ax1.set_yscale("log")
ax1.set_xscale("log")
fs = 20  # 30
fs_text = 30
fs_ticks = 25
textx = 1100
texty = 60  # 500

for comp in range(len(compares)):
    folder = pkls[comp]

    pre = folder + str(obj) + '_' + field + '_' + compares[comp] + '_'
    print(pre)
    base = '_out.pkl'
    extra = pre + 'extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res = pre + 'res' + base
    sed = pre + 'sed' + base
    restwave = pre + 'restwave' + base
    spec = pre + 'spec' + base
    spswave = pre + 'spswave' + base
    chisq = pre + 'chisq' + base
    justchi = pre + 'justchi' + base

    files = [extra, res, sed, restwave, spec, spswave, chisq, justchi]
    print(files)
    print(os.path.exists(files[1]))
    # map.all_plots(files, obj, field, file, loc='upper left', sfh=False, curves=False, sep_uvj=False)

    if os.path.exists(files[1]):
        print('exist pkl')

        # PLOT SEPARATE UVJ
        with open(files[1], 'rb') as res:
            results = pickle.load(res)

        with open(files[2], 'rb') as sed:
            sed = pickle.load(sed)

        with open(files[3], 'rb') as restwave:
            wave_rest = pickle.load(restwave)

        with open(files[4], 'rb') as spec:
            spec = pickle.load(spec)

        with open(files[5], 'rb') as spswave:
            sps_wave = pickle.load(spswave)

        with open(files[7], 'rb') as justchi:
            chi = pickle.load(justchi)

        # MASK
        # create mask
        mask = results['obs']['phot_mask']  # mask out -99.0 values
        wave_rest = np.asarray(wave_rest)
        # no longer need this once I use new output that creates wave_rest as an array
        # ly_mask = (1180 < wave_rest) & (wave_rest < 1260)  # mask out ly-alpha values
        # mask[ly_mask] = False  # combine the masks!

        # apply mask
        phot = results['obs']['maggies'][mask]
        unc = results['obs']['maggies_unc'][mask]
        sed = sed[mask]
        wave_rest = wave_rest[mask]
        chi = chi[mask]

        # CONVERTING TO microjansky
        factor = 3631 * 10 ** 6
        res_jan = []
        err_jan = []
        sed_jan = []
        for i in range(len(wave_rest)):
            res_jan.append(phot[i] * factor)
            err_jan.append(unc[i] * factor)
            sed_jan.append(sed[i] * factor)
        spec_jan = []
        for i in range(len(spec)):
            spec_jan.append(spec[i] * factor)

        ax1.plot(sps_wave, spec_jan, lw=2, color=colors[comp], alpha=0.5, label=labels[comp] + r' Model Spectrum')  # plot spectrum
        # ax1.plot(wave_rest, sed_jan, 'D', color='k', markerfacecolor='None', markersize=10,
        #          markeredgewidth=1.25, markeredgecolor='k', label=labels[comp] + r' Model Photometry')  # plot best fit model
        if comp == 1:
            ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='k',
                         label=r'Observed Photometry')  # plot observations

ax1.set_ylabel(r'Flux [$\mu$Jy]', fontsize=fs_text)  # 30

# ax1.text(textx, texty, field.upper() + '-' + str(obj), fontsize=20)  # 700, 600
ax1.legend(numpoints=1, loc='upper left', prop={'size': 20})  # , line2) ... , r'$\chi$']
xmin = 10 ** 3
xmax = 2.5 * 10 ** 4
ymin = 10**-2  # 10 ** -1  # 6*10**-2
ymax = 50  # 10**4
ax1.set_xlim(xmin, xmax)  # 700, xmax
ax1.set_ylim(ymin, ymax)
ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
ax1.set_xticks([10 ** 3, 2 * 10 ** 3, 5 * 10 ** 3, 10 ** 4, 2 * 10 ** 4])  # technically works
ax1.set_xticklabels([r'$10^3$', r'$2\times10^3$', r'$5 \times 10^3$', r'$10^4$', r'$2\times10^4$'], size=fs_ticks)
ax1.set_yticks([10 ** -2, 10**-1, 10**0, 10 ** 1])  # technically works  , 10 ** 2
ax1.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$'], size=fs_ticks)  # , r'$10^2$'
# xlabel!
fig.text(0.5, 0.03, r'Wavelength (Rest) [$\rm \AA$]', ha='center', va='bottom', fontsize=fs_text)  # 30
plt.show()
