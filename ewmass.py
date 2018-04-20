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
np.errstate(invalid='ignore')
# NOTE: check if OIII / OII high for SFGs because OII essentially just noisy around 0?


def list_ssfr(file):
    '''
    :param file: .txt file listing galaxy id followed by the galaxy's ssfr
    :return:
    '''
    galaxies = []
    counter = 0
    with open(file, 'r') as ssfr:
        for line in ssfr:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                galaxies.append(float(cols[1]) * 10 ** 9)  # ssfrs are stored in 2nd column in file; convert to Gyr^-1
                counter += 1
    return galaxies


def list_ewflux(file, line_idx, ratio=False):
    '''
    lines in files are listed in order: '[NII]6585', '[OII]3726', r'H$\alpha$6563', '[OII]3729', '[OIII]4960',
    '[OIII]5007', r'H$\beta$4861',  r'H$\delta$4102'

    :param file: .txt file listing galaxy id followed by the galaxy's emission line widths OR line fluxes
    :param line_idx: idx 0 - 7 for the desired line EW or flux, OR [idx1, idx2] of lines for desired flux ratio
    :return:
    '''
    galaxies = []
    counter = 0
    ids = []
    with open(file, 'r') as ewflux:
        for line in ewflux:
            if line[0] != '#':
                cols = line.split()
                # print(cols[0])
                ids.append(cols[0])
                if not ratio:
                    galaxies.append(float(cols[line_idx + 1]))  # 0th line_idx is ID
                else:
                    print(line_idx, line_idx[0])
                    print(cols)
                    galaxies.append(float(cols[line_idx[0] + 1]) / float(cols[line_idx[1] + 1]))
                counter += 1

    return galaxies, ids


def preplot(order=None, keys='e', file=None, line_idx=5, line_idx2=[5, 3], ratio=[False, False], outfold=None,
            fs_text=30, fs=20, color='purple', font={'fontname': 'Times'}, types=False):
    '''

    :param keys:
    :param files:
    :param line_idx: single idx for desired line flux or EW, OR [idx1, idx2] for desired line ratios
    :param line_idx2: (if no ssfr) single idx for desired line flux or EW, OR [idx1, idx2] for desired line ratios
    :param ratio:
    :param fs_text:
    :param fs:
    :param font:
    :return:
    '''

    lines = ['[NII]6585', '[OII]3726', r'H$\alpha$6563', '[OII]3729', '[OIII]4960', '[OIII]5007', r'H$\beta$4861',
             r'H$\delta$4102']
    # names = '[NII]6585', '[OII]3726', 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    if keys != 's':
        ew, ew_order = list_ewflux(file, line_idx=line_idx, ratio=ratio[0])  # line1
        # xy[1] = list_ewflux(files[1], line_idx=line_idx2, ratio=ratio[1])  # line2
        labels = [r'Stellar Mass [$\log_{10}(\textrm{M} / \textrm{M}_\odot$)]',
                  lines[line_idx] + r' Equivalent Width [\AA]']
        # labels[1] = lines[line_idx2[0]] + '/' + lines[line_idx2[1]]
    else:  # want ssfr
        ew = list_ssfr(file)  # ssfr
        labels = [r'Stellar Mass [log_{10}(M / M$_\odot$)]', 'SSFR [Gyr$^{-1}$]']

    out = outfold
    gals1 = []
    for file in os.listdir(out):
        if file.endswith(".h5"):
            gals1.append(file)
    gals = []
    relevent_idcs = []
    for ln in range(len(order)):
        for id in range(len(gals1)):
            if gals1[id].startswith(str(order[ln])):
                if order[ln] == 11462 or order[ln] == 21076 or order[ln] == 21442:
                    relevent_idcs.append(ln)
                gals.append(gals1[id])
    # print('1', gals)  # print('2', order)  # print('3', ew_order)  # yay is good!

    mass = []
    for i in range(len(gals)):
        mass.append(gmd.printer(out + gals[i], percs=True)[0])  # median mass

    if types:
        return ew, mass, labels, relevent_idcs
    else:
        return ew, mass, labels


def density_estimation(m1, m2, xs=[8.5, 11.5], ys=[0., 700.]):
    X, Y = np.mgrid[xs[0]:xs[1]:100j, ys[0]:ys[1]:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def plotter(x, y, color, errs=True, fs_text=30, fs=20, norm=7., xs=[8.5, 11.5], ys=[0., 700.], scat=True):
    # font={'fontname': 'Times'},
    if scat:
        plt.scatter(x=x, y=y, color=color, s=fs * 2)

    if errs:
        '''
        errx = np.percentile(x, [16., 50., 84.])
        erry = np.percentile(y, [16., 50., 84.])
        plt.errorbar(x=[errx[1]], y=[erry[1]], xerr=[[errx[1] - errx[0]], [errx[2] - errx[1]]],
                     yerr=[[erry[1] - erry[0]], [erry[2] - erry[1]]], color=color, fmt='*')
        '''
        X, Y, Z = density_estimation(x, y, xs=xs, ys=ys)
        # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
        levels = np.arange(np.amax(Z) / norm, np.amax(Z), np.amax(Z) / norm) + (np.amax(Z) / norm)
        if color == 'b':
            cmap = 'Blues'
        else:
            cmap = 'Purples'
        plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]
        # plt.contour(X, Y, Z, levels=levels, cmap=cmap, lw=3)#, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]


def plotter_temp(x, y, color, types=None, errs=True, fs_text=30, fs=20, norm=7., xs=[8.5, 11.5], ys=[0., 700.], scat=True):
    # font={'fontname': 'Times'},
    if scat:
        if types is not None:
            plt.scatter(x=x, y=y, color=color, marker=types, s=fs * 8)
        else:
            plt.scatter(x=x, y=y, color=color, s=fs * 2)

    if errs:
        '''
        errx = np.percentile(x, [16., 50., 84.])
        erry = np.percentile(y, [16., 50., 84.])
        plt.errorbar(x=[errx[1]], y=[erry[1]], xerr=[[errx[1] - errx[0]], [errx[2] - errx[1]]],
                     yerr=[[erry[1] - erry[0]], [erry[2] - erry[1]]], color=color, fmt='*')
        '''
        X, Y, Z = density_estimation(x, y, xs=xs, ys=ys)
        # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
        levels = np.arange(np.amax(Z) / norm, np.amax(Z), np.amax(Z) / norm) + (np.amax(Z) / norm)
        if color == 'b':
            cmap = 'Blues'
        else:
            cmap = 'Purples'
        plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]
        # plt.contour(X, Y, Z, levels=levels, cmap=cmap, lw=3)#, alpha=0.5)  # [0.1, 0.2, 0.5, 1., 25.]


def mass_match_plot(x1, y1, color, order, lims=[9.5, 10.], errs=True, outfold=None, ssfr_file=None, fs_text=30, fs=20,
                    font={'fontname': 'Times'}):
    x = []
    y = []
    do_these = []
    for i in range(len(x1)):
        if lims[0] <= x1[i] <= lims[1]:  # 9.55 <= x1[i] <= 9.9:
            do_these.append(order[i])
            x.append(x1[i])
            y.append(y1[i])

    plt.scatter(x=x, y=y, color=color, s=fs * 2)

    if errs:
        '''
        errx = np.percentile(x, [16., 50., 84.])
        erry = np.percentile(y, [16., 50., 84.])
        plt.errorbar(x=[errx[1]], y=[erry[1]], xerr=[[errx[1] - errx[0]], [errx[2] - errx[1]]],
                     yerr=[[erry[1] - erry[0]], [erry[2] - erry[1]]], color=color, fmt='*')
        '''
        X, Y, Z = density_estimation(x, y)
        # levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
        levels = np.arange(np.amax(Z) / 7., np.amax(Z), np.amax(Z) / 7.) + (np.amax(Z) / 7.)
        if color == 'b':
            cmap = 'Blues'
        else:
            cmap = 'Purples'
        plt.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]

    out = outfold
    gals1 = []
    for file in os.listdir(out):
        if file.endswith(".h5"):
            gals1.append(file)
    gals = []
    for ln in range(len(do_these)):
        for id in range(len(gals1)):
            if gals1[id].startswith(str(do_these[ln])):
                gals.append(gals1[id])
    mass = []
    dust = []
    met = []
    for i in range(len(gals)):
        get = gmd.printer(out + gals[i], percs=True)  # median mass, dust, gasmet, metallicity
        mass.append(get[0])
        dust.append(get[1])
        met.append(get[2])

    ssfrs = []
    counter = 0
    with open(ssfr_file, 'r') as ssfr:
        for line in ssfr:
            print(line)
            if line[0] != '#':
                cols = line.split()
                print(cols)
                for dt in range(len(do_these)):
                    if cols[0].startswith(str(do_these[dt])):
                        ssfrs.append(float(cols[1]) * 10 ** 9)  # ssfrs are stored in 2nd column in file; convert to Gyr^-1
                        counter += 1

    return mass, dust, met, ssfrs, do_these


def return_fracs(gals, order, pklfolder='pkl_efico/', outfold='out_efico/', base='fico', ls=False):
    import ssfr_now as sfnow
    eelgs, lbgs1 = sa.get_gal_lists(base=[base, base], normal=True)
    git = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    out = git + 'out/' + outfold
    f_folder = git + pklfolder
    gals1 = gals
    if ls:
        gals1 = []
        for i in range(len(lbgs1)):
            file = f_folder + lbgs1[i] + '_extra_out.pkl'
            if os.path.exists(file):
                for ga in gals:
                    print(str(ga))
                    if file.startswith(f_folder + str(ga)):
                        gals1.append(lbgs1[i])
    # print(lbgs1)
    print(len(gals1), 'look')

    '''
    ordered = []
    for go in range(len(order)):
        for ga in range(len(gals1)):
            if gals1[ga].startswith(order[go]):
                ordered.append(gals1[ga])
    '''
    ordered = gals1

    time = 0.5 * 10 ** 8
    sfh = []
    highfracs = []
    # print('masses', masses)
    frac = np.zeros(shape=(len(ordered), 10**3))
    for ord in range(len(ordered)):
        # print(ordered[ord])
        # file = f_folder + str(glxy) + '*_extra_out.pkl'
        for pkls in os.listdir(f_folder):
            if pkls.startswith(str(ordered[ord])) and pkls.endswith('_extra_out.pkl'):
                file = f_folder + pkls
                sfh.append(sfnow.get_sfh(file)[0])  # 1st bin
    # for i in range(len(ordered)):
        # print(ord)
        file_i = None
        for outs in os.listdir(out):
            if outs.startswith(str(ordered[ord])) and outs.endswith('.h5'):
                file_i = out + outs
        get_e = gmd.printer(file_i, percs=False, quiet=True)
        for k in range(10 ** 3):
            frac[ord, k] = sfh[ord][np.random.randint(len(sfh[ord]))] * time / (10 ** get_e[np.random.randint(len(get_e))])

    # print(frac)
    all_true = []
    for fr in range(len(frac)):
        for k in range(10**3):
            all_true.append(frac[fr, k])
    # print(all_fracs)

    return np.percentile(all_true, [16., 50., 84.])


def dmass(obj_e, field_e, order, three_masses, folder='out_efico/', f_ind=0, v_met=False):
    field_dict = {}
    for i in range(len(obj_e)):
        field_dict[obj_e[i]] = field_e[i]

    base_lab = r'M$_{\rm P}$ / M$_{\rm F}$ '
    three_labels = [base_lab + r'(Z$_{\rm F}$ = Z$_{\odot}$, with emission lines)',
                    base_lab + r'(Z$_{\rm F}$ = Z$_{\odot}/5$, with emission lines)',
                    base_lab + r'(FAST with ZFOURGE catalog parameters)']
    colors = ['r', 'b', 'purple']
    shapes = ['o', 's', 'v']

    use_this = three_masses[f_ind]
    use_lbl = three_labels[f_ind]

    dictionary = fc.get_fast(use_this)
    if v_met:
        mass_dict, met_dict = fc.compare_gmd(dictionary, folder, include_met=v_met)
    else:
        mass_dict = fc.compare_gmd(dictionary, folder)
    # print(mass_dict)

    fast, prosp, xratio, xfield, mratio, mets = [], [], [], [], [], []
    for key in mass_dict:  # mass_diff[key][0] is the FAST mass, mass_diff[key][1] is the Prospector mass
        # print(mass_dict[key][0], mass_dict[key][1])
        fast.append(float(mass_dict[key][0]))  # FAST
        prosp.append(float(mass_dict[key][1]))  # Prospector
        mratio.append((10 ** float(mass_dict[key][1])) / (10 ** float(mass_dict[key][0])))
        xratio.append(key)
        if v_met:
            mets.append(met_dict[key])
        for f_key in field_dict:
            if int(key) == int(f_key):
                xfield.append(field_dict[f_key])

    new_met_order = []
    new_mass_order = []
    for ln in range(len(order)):
        for id in range(len(xratio)):
            if xratio[id].startswith(str(order[ln])):
                new_mass_order.append(mratio[id])
                if v_met:
                    new_met_order.append(mets[id])
    # PLOT!
    print(len(xfield), len(mratio))
    print(xfield)
    '''
    for l in range(len(xfield)):
        if xfield[l] == 'cosmos':
            xfield[l] = 'cos'
    fs_text = 30  # 30
    fs = 20
    xlabs = []
    for j in range(len(xratio)):
        xlabs.append(str(xfield[j]).upper() + xratio[j])
    # plt.axhline(y=np.median(mratio), color=colors[i])
    # xspace = np.linspace(0, 19, len(ratio))
    # plt.scatter(xspace, mratio, marker=shapes[i], color=colors[i], label=labels[i], s=fs)
    # plt.xticks(xspace, xlabs, rotation=60)  # 'vertical')
    print(xlabs)

    plt.ylabel(r'Stellar mass ratio (Prospector / FAST)', fontsize=fs_text)
    # plt.xlabel(r'Galaxy IDs', fontsize=fs_text)
    plt.axhline(y=1., color='k', linestyle='--')
    plt.ylim(0, 17)  # 37)
    tick_f = 15
    # plt.yscale('log')
    plt.legend(numpoints=1, loc='upper left', prop={'size': 20})

    plt.tick_params('x', length=3, width=1, which='both', labelsize=tick_f)
    plt.tick_params('y', length=3, width=0.5, which='both', labelsize=tick_f)

    plt.show()
    '''
    if v_met:
        return new_mass_order, use_lbl, new_met_order
    else:
        return new_mass_order, use_lbl


def agewrite(gal_order, fields, name1='eelg', base='fico', pfolder='pkl_efico/'):
    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + pfolder
    agelist = None

    with open(name1 + '_ages.txt', 'w+') as write_ages:
        write_ages.write('# ID Age[Gyr]\n')

        for gl in range(len(gal_order)):
            agelist = []
            file = pkls + str(gal_order[gl]) + '_' + fields[gl] + '_' + base + '_extra_out.pkl'
            if os.path.exists(file):
                import eelg_fifty_params as params
                model = params.load_model(objname=str(gal_order[gl]), field=fields[gl])
                agebins = model.params['agebins']
                # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
                temp = sa.randraw_sfr_perc(file)
                temp_draw = sa.smooth(temp[0])
                this_age = sa.stellar_age(temp_draw, agebins) / 1e9
                agelist.append(this_age)
                print(str(gal_order[gl]) + ' ' + str(this_age) + '\n')
                write_ages.write(str(gal_order[gl]) + ' ' + str(this_age) + '\n')
            else:
                print('missing ' + file)  # print galaxy if pkls don't exist for it

    return agelist


def ageread(gal_order, agefile='eelg_ages.txt'):
    agefile = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + agefile

    agelist = []
    gals = []
    with open(agefile, 'r') as read_ages:
        for line in read_ages:
            if line[0] != '#':
                cols = line.split()
                gals.append(cols[0])
                agelist.append(cols[1])

    new_age_order = []
    for go in range(len(gal_order)):
        for ga in range(len(gals)):
            if int(gals[ga]) == int(gal_order[go]):
                new_age_order.append(float(agelist[ga]))

    return new_age_order


def d_age(gal_order, field_order, agefile='eelg_ages.txt'):
    new_age_order = ageread(gal_order, agefile)

    home = '/home/jonathan/'
    field_files = [home + 'cdfs/cdfs.v1.6.9_EL_Z004.awk.fout', home + 'cosmos/cosmos.v1.3.6_EL_Z004.awk.fout',
                   home + 'uds/uds.v1.5.8_EL_Z004.awk.fout']

    fast_ages = []
    for fld in range(len(field_order)):
        if field_order[fld] == 'cdfs':
            field_file = field_files[0]
        elif field_order[fld] == 'cosmos':
            field_file = field_files[1]
        else:
            field_file = field_files[2]

        with open(field_file, 'r') as ffile:
            for line in ffile:
                cols = line.split()
                if cols[0] == str(gal_order[fld]):
                    fast_ages.append(10**float(cols[4])/ 10**9)

    print(len(fast_ages), len(field_order), len(new_age_order))
    age_ratio = []
    for na in range(len(new_age_order)):
        age_ratio.append(new_age_order[na] / fast_ages[na])

    return age_ratio


if __name__ == "__main__":

    do_mass = 0  # do_mass --> dMass_v_ew, else ew_v_mass\
    mass_match = 0
    dm_met = 0
    dm_age = 0
    dm_dage = 1

    path = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    e_out = path + 'out/out_efico/'
    s_out = path + 'out/out_nfico/'
    e_fs = {'s': path + 'eelg_ssfrs.txt', 'e': path + 'ews_eelgs.txt', 'f': path + 'fluxes_eelgs.txt'}
    s_fs = {'s': path + 'sfg_ssfrs.txt', 'e': path + 'ews_sfgs.txt', 'f': path + 'fluxes_sfgs.txt'}

    # same for EELGs, SFGs
    key = 'e'
    ratio = [False, True]
    idx = 5  # [NII]6585, [OII]3726, 'Halpha6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'Hbeta4861', 'Hdelta4102'

    ee_order, ee_fd, sf_order, sf_fd = sa.get_gal_lists(base=['fico', 'fico'], objlists=True, normal=True)

    # different for EELGs, SFGs
    colors = ['purple', 'b']  # eelg, sfg
    eelg_file = e_fs[key]
    sfg_file = s_fs[key]

    # eew, emass, labs = preplot(order=ee_order, keys=key, file=eelg_file, line_idx=idx, ratio=ratio, outfold=e_out)
    eew, emass, labs, r_idcs = preplot(order=ee_order, keys=key, file=eelg_file, line_idx=idx, ratio=ratio,
                                       outfold=e_out, types=True)
    print('hi')
    sew, smass, labs = preplot(order=sf_order, keys=key, file=sfg_file, line_idx=idx, ratio=ratio, outfold=s_out)
    print(len(eew), len(emass), len(sew), len(smass))

    sfg_x = []
    sfg_y = []
    for thing in range(len(smass)):
        if sew[thing] <= 10**-16 or smass[thing] <= 10**-16:
            # print('me')
            pass
        else:
            sfg_x.append(smass[thing])
            sfg_y.append(sew[thing])
    smass = sfg_x
    sew = sfg_y
    print(len(eew), len(emass), len(sew), len(smass))

    if do_mass:
        # f_ind=1 --> FAST with Z = Z_sol / 5, with emission lines
        home = '/home/jonathan/mz_files/'
        three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat',
                        home + 'Comp_10_zm_ZFOURGE.dat']
        three_smasses = [home + 'Comp_00_zm_EL_Z002.dat', home + 'Comp_00_zm_EL_Z004.dat',
                         home + 'Comp_00_zm_ZFOURGE.dat']

        mass_ratio, mass_label = dmass(obj_e=ee_order, field_e=ee_fd, order=ee_order, three_masses=three_masses,
                                       folder='out_efico/', f_ind=1)

        smass_ratio, smass_label = dmass(obj_e=sf_order, field_e=sf_fd, order=sf_order, three_masses=three_smasses,
                                         folder='out_nfico/', f_ind=1)
        plotter(eew, mass_ratio, color=colors[0], errs=False)
        plotter(sew, smass_ratio, color=colors[1], errs=False)
        xmin, xmax = 0., 1.2 * max([max(eew), max(sew)])  # 25.
        ymin, ymax = 0., 16.
        labs = [labs[1], mass_label]

    elif dm_met:
        # f_ind=1 --> FAST with Z = Z_sol / 5, with emission lines
        home = '/home/jonathan/mz_files/'
        three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat',
                        home + 'Comp_10_zm_ZFOURGE.dat']
        three_smasses = [home + 'Comp_00_zm_EL_Z002.dat', home + 'Comp_00_zm_EL_Z004.dat',
                         home + 'Comp_00_zm_ZFOURGE.dat']

        mass_ratio, mass_label, mets = dmass(obj_e=ee_order, field_e=ee_fd, order=ee_order, three_masses=three_masses,
                                             folder='out_efico/', f_ind=1, v_met=True)

        smass_ratio, smass_label, smets = dmass(obj_e=sf_order, field_e=sf_fd, order=sf_order,
                                                three_masses=three_smasses, folder='out_nfico/', f_ind=1, v_met=True)

        print('mass ratio', np.percentile(mass_ratio, [16., 50., 84.]))
        print('smass ratio', np.percentile(smass_ratio, [16., 50., 84.]))
        print(max(mets), min(mets), max(smets), min(smets), 'look!')
        plotter(smets, smass_ratio, color=colors[1], errs=True, xs=[-2, 0.5], ys=[0., 14.], norm=10.)
        plotter(mets, mass_ratio, color=colors[0], errs=True, xs=[-2., 0.5], ys=[0., 14.], norm=10.)
        # allx = np.concatenate((mets, smets), axis=0)
        # ally = np.concatenate((mass_ratio, smass_ratio), axis=0)
        # plt.plot(np.unique(mets), np.poly1d(np.polyfit(mets, mass_ratio, 2))(np.unique(mets)), color='purple',
        #          linestyle='--')
        # plt.plot(np.unique(smets), np.poly1d(np.polyfit(allx, ally, 1))(np.unique(smets)), color='blue',
        #          linestyle='--')

        xmin, xmax = -2., 0.5  #  1.2 * max([max(eew), max(sew)])  # 25.
        ymin, ymax = 0., 14.
        labs = ['$\log_{10}$(Z / Z$_{\odot}$)', mass_label]
        # plt.xscale('log')

    elif dm_age:
        # agelist = agewrite(ee_order, ee_fd, name1='eelg', base='fico', pfolder='pkl_efico/')
        # sagelist = agewrite(sf_order, sf_fd, name1='sfg', base='fico', pfolder='pkl_nfico/')
        agelist = ageread(ee_order, agefile='eelg_ages.txt')
        sagelist = ageread(sf_order, agefile='sfg_ages.txt')

        home = '/home/jonathan/mz_files/'
        three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat',
                        home + 'Comp_10_zm_ZFOURGE.dat']
        three_smasses = [home + 'Comp_00_zm_EL_Z002.dat', home + 'Comp_00_zm_EL_Z004.dat',
                         home + 'Comp_00_zm_ZFOURGE.dat']

        mass_ratio, mass_label = dmass(obj_e=ee_order, field_e=ee_fd, order=ee_order, three_masses=three_masses,
                                       folder='out_efico/', f_ind=1, v_met=False)

        smass_ratio, smass_label = dmass(obj_e=sf_order, field_e=sf_fd, order=sf_order,
                                         three_masses=three_smasses, folder='out_nfico/', f_ind=1, v_met=False)

        print(len(sagelist), len(smass_ratio))  # 128, 128
        print(len(agelist), len(mass_ratio))  # 19, 19
        # print(max(agelist), max(sagelist), min(agelist), min(sagelist))
        fastfits = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'\
                   + 'fast_fits_properties/'

        print('smass', np.percentile(smass_ratio, [16., 50., 84.]))
        print('mass', np.percentile(mass_ratio, [16., 50., 84.]))

        plotter(sagelist, smass_ratio, color=colors[1], errs=True, xs=[0., 3.], ys=[0., 14.], norm=10.)
        plotter(agelist, mass_ratio, color=colors[0], errs=True, xs=[0., 3.], ys=[0., 14.], norm=10.)
        # allx = np.concatenate((mets, smets), axis=0)
        # ally = np.concatenate((mass_ratio, smass_ratio), axis=0)
        # plt.plot(np.unique(mets), np.poly1d(np.polyfit(mets, mass_ratio, 2))(np.unique(mets)), color='purple',
        #          linestyle='--')
        # plt.plot(np.unique(smets), np.poly1d(np.polyfit(allx, ally, 1))(np.unique(smets)), color='blue',
        #          linestyle='--')

        xmin, xmax = 0., 3.0  #  1.2 * max([max(eew), max(sew)])  # 25.
        ymin, ymax = 0., 14.
        labs = ['Stellar age [Gyr]', mass_label]

    elif dm_dage:
        # agelist = agewrite(ee_order, ee_fd, name1='eelg', base='fico', pfolder='pkl_efico/')
        # sagelist = agewrite(sf_order, sf_fd, name1='sfg', base='fico', pfolder='pkl_nfico/')
        ageratio = d_age(ee_order, ee_fd, agefile='eelg_ages.txt')
        sageratio = d_age(sf_order, sf_fd, agefile='sfg_ages.txt')

        home = '/home/jonathan/mz_files/'
        three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat',
                        home + 'Comp_10_zm_ZFOURGE.dat']
        three_smasses = [home + 'Comp_00_zm_EL_Z002.dat', home + 'Comp_00_zm_EL_Z004.dat',
                         home + 'Comp_00_zm_ZFOURGE.dat']

        mass_ratio, mass_label = dmass(obj_e=ee_order, field_e=ee_fd, order=ee_order, three_masses=three_masses,
                                       folder='out_efico/', f_ind=1, v_met=False)

        smass_ratio, smass_label = dmass(obj_e=sf_order, field_e=sf_fd, order=sf_order,
                                         three_masses=three_smasses, folder='out_nfico/', f_ind=1, v_met=False)

        print(len(sageratio), len(smass_ratio))  # 128, 128
        print(len(ageratio), len(mass_ratio))  # 19, 19
        # print(max(agelist), max(sagelist), min(agelist), min(sagelist))
        fastfits = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'\
                   + 'fast_fits_properties/'

        print('smass', np.percentile(smass_ratio, [16., 50., 84.]))
        print('mass', np.percentile(mass_ratio, [16., 50., 84.]))
        xs_age = [0., 200.]
        ys_age = [0., 14.]
        xmin, xmax = 0., 1.2 * max([max(ageratio), max(sageratio)])  # 25.
        ymin, ymax = 0., 14.

        logify = False
        if logify:
            sageratio = [np.log10(x) for x in sageratio]
            ageratio = [np.log10(x) for x in ageratio]
            smass_ratio = [np.log10(x) for x in smass_ratio]
            mass_ratio = [np.log10(x) for x in mass_ratio]
            xs_age = [-1., 2.5]
            ys_age = [-1., 1.5]
            xmin, xmax = 0.5, 2.5  # 1.2 * max([max(ageratio), max(sageratio)])  # 25.
            ymin, ymax = 0.2, 1.2  # 1.2 * max([max(mass_ratio), max(smass_ratio)])

        plotter(sageratio, smass_ratio, color=colors[1], errs=True, xs=xs_age, ys=ys_age, norm=7., scat=False)
        plotter(ageratio, mass_ratio, color=colors[0], errs=True, xs=xs_age, ys=ys_age, norm=7., scat=False)
        print(np.percentile(ageratio, [16., 50., 84.]))
        print(np.percentile(mass_ratio, [16., 50., 84.]))

        labs = ['Prospector age / Fast age', mass_label]

    elif mass_match:
        xmin, xmax = 9.4, 10.1  # 1.2 * max([max(emass), max(smass)])  # 50.  # 10 ** 3  # 700.
        ymin, ymax = 0., 1.2 * max([max(eew), max(sew)])  # 25.
        ethese = mass_match_plot(emass, eew, color=colors[0], order=ee_order, outfold=e_out, ssfr_file=e_fs['s'])
        sthese = mass_match_plot(smass, sew, color=colors[1], order=sf_order, outfold=s_out, ssfr_file=s_fs['s'])

        do_these = ethese[4]
        sdo_these = sthese[4]

        efrac = return_fracs(do_these, ee_order, pklfolder='pkl_efico/', outfold='out_efico/', ls=False)
        sfrac = return_fracs(sdo_these, sf_order, pklfolder='pkl_nfico/', outfold='out_nfico/', ls=True)
        print('efrac', efrac)
        print('sfrac', sfrac)

        print(len(ethese))
        print(len(sthese))
        mass_percs = np.percentile(ethese[0], [16., 50, 84.])
        print(mass_percs, 'e-mass-matched-mass')
        smass_percs = np.percentile(sthese[0], [16., 50, 84.])
        print(smass_percs, 's-mass-matched-mass')
        dust_percs = np.percentile(ethese[1], [16., 50, 84.])
        print(dust_percs, 'e-mass-matched-dust')
        sdust_percs = np.percentile(sthese[1], [16., 50, 84.])
        print(sdust_percs, 's-mass-matched-dust')
        met_percs = np.percentile(ethese[2], [16., 50, 84.])
        print(met_percs, 'e-mass-matched-met')
        smet_percs = np.percentile(sthese[2], [16., 50, 84.])
        print(smet_percs, 's-mass-matched-met')
        ssfr_percs = np.percentile(ethese[3], [16., 50, 84.])
        print(ssfr_percs, 'e-mass-matched-ssfr')
        sssfr_percs = np.percentile(sthese[3], [16., 50, 84.])
        print(sssfr_percs, 's-mass-matched-ssfr')
    else:
        xmin, xmax = 8.5, 11.5  # 1.2 * max([max(emass), max(smass)])  # 50.  # 10 ** 3  # 700.
        ymin, ymax = 0., 1.2 * max([max(eew), max(sew)])  # 25.
        # plotter(emass, eew, color=colors[0])
        # plotter(smass, sew, color=colors[1])

        fmts = 'o'
        x1 = []
        y1 = []
        for r_idc in r_idcs:
            x1.append(emass[r_idc])
            y1.append(eew[r_idc])
            fmt1 = '*'

        plotter(emass, eew, color=colors[0])
        plotter(smass, sew, color=colors[1])
        plotter_temp(x1, y1, color=colors[0], types=fmt1, errs=False)

        sfrac = return_fracs(sf_order, sf_order, pklfolder='pkl_nfico/', outfold='out_nfico/', ls=True)
        efrac = return_fracs(ee_order, ee_order, pklfolder='pkl_efico/', outfold='out_efico/', ls=False)
        print('efrac', efrac)
        print('sfrac', sfrac)

    # PLOT PARAMS
    fs_text = 30
    fs = 20
    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)
    font = {'fontname': 'Times'}
    log = False
    if log:
        plt.yscale('log')
        plt.xscale('log')
        xmin, xmax = 2., 2*10**2  # 10**3
        ymin, ymax = 1., 20  # 3*10
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