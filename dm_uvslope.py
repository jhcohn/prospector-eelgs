import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log
import sys
import argparse
import pickle
import os
import glob
import copy
from astropy import constants
import eelg_fifty_params as pfile
import get_mass_dust as gmd
import stellar_ages as sa
from matplotlib import rc
from scipy import stats
import fast_compare as fc
from astropy.table import Table
import eelg_fifty_params as efp
from scipy.optimize import curve_fit
np.errstate(invalid='ignore')
# NOTE: check if OIII / OII high for SFGs because OII essentially just noisy around 0?


def calc_uv_slope(fileset, plot=False):
    # PLOT SEPARATE UVJ
    with open(fileset[1], 'rb') as res:
        results = pickle.load(res)

    with open(fileset[2], 'rb') as sed:
        sed = pickle.load(sed)

    with open(fileset[3], 'rb') as restwave:
        wave_rest = pickle.load(restwave)

    with open(fileset[4], 'rb') as spec:
        spec = pickle.load(spec)

    with open(fileset[5], 'rb') as spswave:
        sps_wave = pickle.load(spswave)

    with open(fileset[7], 'rb') as justchi:
        chi = pickle.load(justchi)

    # MASK
    mask = results['obs']['phot_mask']  # mask out -99.0 values
    wave_rest = np.asarray(wave_rest)
    # no longer need this once I use new output that creates wave_rest as an array
    # ly_mask = (1180 < wave_rest) & (wave_rest < 1260)  # mask out ly-alpha values
    # mask[ly_mask] = False  # combine the masks!

    # apply mask
    phot = results['obs']['maggies']# [mask]
    sed = sed[mask]
    wave_rest = wave_rest[mask]

    c_ang = 3. * 10 ** 18  # speed of light, angstroms/s
    slope_waves = []
    slope_flux = []
    for i in range(len(wave_rest)):
        if 1500 < wave_rest[i] < 2600:
            slope_waves.append(wave_rest[i])
            # slope_flux.append(sed[i] * 3631 * 10**6.44)  # maggies to F_nu
            slope_flux.append(sed[i] * 3631 * 10 ** 6.44 * c_ang / (wave_rest[i])**2)  # F_nu to F_lam!

    sorted_waves = sorted(slope_waves)
    sorted_fluxes = []
    idx = 0
    while idx < len(sorted_waves):
        # print('hm', idx, len(sorted_waves))
        for j in range(len(slope_waves)):
            if slope_waves[j] == sorted_waves[idx]:
                sorted_fluxes.append(slope_flux[j])
                idx += 1
                if idx >= len(sorted_waves):
                    break
    # print(slope_waves, sorted_waves)
    # print(slope_flux, sorted_fluxes)

    if len(slope_waves) <= 2:
        return [99., 99., 99.]
    else:
        popt, pcov = curve_fit(powerlaw, sorted_waves, sorted_fluxes, maxfev=1500)
        # print(popt, 'look!')
        if plot:
            plt.plot(sorted_waves, sorted_fluxes, 'bo', label='Photometry')
            plt.plot(sorted_waves, powerlaw(sorted_waves, *popt), 'k')
            plt.show()

        return popt


def powerlaw(x, m, a, b):
    return a + b*x ** m


def dmass(obj_e, field_e, order, three_masses, folder='out_efico/', f_ind=0, v_met=False, quiet=True):
    field_dict = {}
    for i in range(len(obj_e)):
        field_dict[obj_e[i]] = field_e[i]

    base_lab = r'M$_{\rm P}$ / M$_{\rm F}$ '
    three_labels = [base_lab + r'(Z$_{\rm F}$ = Z$_{\odot}$, with emission lines)',
                    base_lab + r'(Z$_{\rm F}$ = Z$_{\odot}/5$, with emission lines)',
                    base_lab + r'(FAST with ZFOURGE catalog parameters)']

    use_this = three_masses[f_ind]
    use_lbl = three_labels[f_ind]

    dictionary = fc.get_fast(use_this)
    mass_dict = fc.compare_gmd(dictionary, folder, quiet=quiet)
    # print(mass_dict)

    fast, prosp, xratio, xfield, mratio, mets = [], [], [], [], [], []
    for key in mass_dict:  # mass_diff[key][0] is the FAST mass, mass_diff[key][1] is the Prospector mass
        # print(mass_dict[key][0], mass_dict[key][1])
        fast.append(float(mass_dict[key][0]))  # FAST
        prosp.append(float(mass_dict[key][1]))  # Prospector
        mratio.append((10 ** float(mass_dict[key][1])) / (10 ** float(mass_dict[key][0])))
        xratio.append(key)
        for f_key in field_dict:
            if int(key) == int(f_key):
                xfield.append(field_dict[f_key])

    from collections import OrderedDict
    new_mass_order = []
    for ln in range(len(order)):
        for id in range(len(xratio)):
            if xratio[id] == str(order[ln]):
                new_mass_order.append(mratio[id])
    print(new_mass_order)
    new_mass_order = list(OrderedDict.fromkeys(new_mass_order))
    print(new_mass_order)

    return new_mass_order, use_lbl, mass_dict, len(new_mass_order), len(xratio)


def density_estimation(m1, m2, xs=[8.5, 11.5], ys=[0., 700.]):
    X, Y = np.mgrid[xs[0]:xs[1]:100j, ys[0]:ys[1]:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def plotter(x, y, color, norm=7., xs=[8.5, 11.5], ys=[0., 700.], scat=True):
    X, Y, Z = density_estimation(x, y, xs=xs, ys=ys)
    # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
    levels = np.arange(np.amax(Z) / norm, np.amax(Z), np.amax(Z) / norm) + (np.amax(Z) / norm)
    if color == 'b':
        cmap = 'Blues'
    else:
        cmap = 'Purples'
    plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]
    # plt.contour(X, Y, Z, levels=levels, cmap=cmap, lw=3)#, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]


if __name__ == "__main__":

    path = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    pkls = path + 'pkl_efico/'
    spkls = path + 'pkl_nfico/'
    e_out = path + 'out/out_efico/'
    s_out = path + 'out/out_nfico/'

    # same for EELGs, SFGs
    key = 'e'
    ratio = [False, True]
    idx = 5  # [NII]6585, [OII]3726, 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    print('getting lists')
    ee_order, ee_fd, sf_order, sf_fd = sa.get_gal_lists(base=['fico', 'fico'], objlists=True, normal=True)

    # different for EELGs, SFGs
    colors = ['purple', 'b']  # eelg, sfg

    # f_ind=1 --> FAST with Z = Z_sol / 5, with emission lines
    home = '/home/jonathan/mz_files/'
    three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat',
                    home + 'Comp_10_zm_ZFOURGE.dat']
    three_smasses = [home + 'Comp_00_zm_EL_Z002.dat', home + 'Comp_00_zm_EL_Z004.dat',
                     home + 'Comp_00_zm_ZFOURGE.dat']

    print('mass ratios')
    mass_ratio, mass_label, mdict, nmo, xr = dmass(obj_e=ee_order, field_e=ee_fd, order=ee_order, three_masses=three_masses,
                                          folder='out_efico/', f_ind=1)

    print('smass ratios')
    smass_ratio, smass_label, smdict, snmo, sxr = dmass(obj_e=sf_order, field_e=sf_fd, order=sf_order, three_masses=three_smasses,
                                             folder='out_nfico/', f_ind=1)
    labs = [r'UV slope [$\beta$]', mass_label]

    slopes = []
    for ee in range(len(ee_order)):
        obj = ee_order[ee]
        field = ee_fd[ee]
        pre = 'pkl_efico/' + str(obj) + '_' + field + '_fico'
        base = '_out.pkl'
        extra1 = pre + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
        res1 = pre + '_res' + base
        sed1 = pre + '_sed' + base
        restwave1 = pre + '_restwave' + base
        spec1 = pre + '_spec' + base
        spswave1 = pre + '_spswave' + base
        chisq1 = pre + '_chisq' + base
        justchi1 = pre + '_justchi' + base
        fileset = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
        # print(path + pre)
        if os.path.exists(path + extra1):
            # print('me!')
            popt = calc_uv_slope(fileset)
            print(popt[0])
            slopes.append(popt[0])

    sslopes = []
    sgals = []
    num = 0
    for sf in range(len(sf_order)):
        obj = sf_order[sf]
        field = sf_fd[sf]
        pre = 'pkl_nfico/' + str(obj) + '_' + field + '_fico'
        base = '_out.pkl'
        extra1 = pre + '_extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
        res1 = pre + '_res' + base
        sed1 = pre + '_sed' + base
        restwave1 = pre + '_restwave' + base
        spec1 = pre + '_spec' + base
        spswave1 = pre + '_spswave' + base
        chisq1 = pre + '_chisq' + base
        justchi1 = pre + '_justchi' + base
        fileset = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
        if os.path.exists(path + extra1):
            num += 1
            sgals.append(obj)
            # print('sme!')
            popt = calc_uv_slope(fileset)
            print(popt[0])
            sslopes.append(popt[0])

    sgs = []
    print(len(smdict))
    for key in smdict:
        match = False
        for sg in range(len(sgals)):
            if int(key) == int(sgals[sg]):
                sgs.append(sg)
                match = True
        if not match:
            print(key, 'no match')
    print(sorted(sgs))

    new_sslopes = []
    new_smr = []
    for slp in range(len(sslopes)):
        if -3. <= sslopes[slp] <= 0.5:
            new_sslopes.append(sslopes[slp])
            new_smr.append(smass_ratio[slp])
    smass_ratio = new_smr
    sslopes = new_sslopes

    print(len(sslopes), len(smass_ratio))
    print(len(slopes), len(mass_ratio))
    print(num)
    print(snmo, sxr)
    print(np.percentile(slopes, [16., 50., 84.]))
    print(np.percentile(sslopes, [16., 50., 84.]))
    print(np.percentile(mass_ratio, [16., 50., 84.]))
    print(np.percentile(smass_ratio, [16., 50., 84.]))
    # PLOT PARAMS
    contours = True

    if contours:
        plotter(sslopes, smass_ratio, 'b', xs=[-3., -1.0], ys=[1., 10.])
        plotter(slopes, mass_ratio, 'purple', xs=[-3, -1.5], ys=[1., 14.])
        print(min(sslopes), max(sslopes), max(slopes), min(slopes))
    else:
        plt.scatter(sslopes, smass_ratio, color='b')
        plt.scatter(slopes, mass_ratio, color='purple')
    fs_text = 30
    fs = 20
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}
    xmin, xmax = -2.5, -1.  # 1.2 * max([max(slope), max(sslope)])  # 25.
    ymin, ymax = 0., 16.
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    plt.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.xlabel(labs[0], fontsize=fs_text, **font)  # 20
    plt.ylabel(labs[1], fontsize=fs_text, **font)  # 20
    plt.show()


'''
(array([ 9.57749403,  9.76314425,  9.87482784]), 'e-mass-matched-mass')
(array([  9.77374964,   9.92011118,  10.03033199]), 's-mass-matched-mass')
(array([ 0.11200133,  0.29464107,  0.38801465]), 'e-mass-matched-dust')
(array([ 0.3373157 ,  0.57326892,  0.78298057]), 's-mass-matched-dust')
(array([-0.34760772, -0.30852319, -0.2841036 ]), 'e-mass-matched-met')
(array([-1.32048439, -0.68648973, -0.15814167]), 's-mass-matched-met')
(array([ 1.57868734,  2.86966937,  3.06445666]), 'e-mass-matched-ssfr')
(array([ 0.50228701,  1.20524879,  2.77549461]), 's-mass-matched-ssfr')
'''