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
np.errstate(invalid='ignore')

pc = 3.085677581467192e18  # in cm
dfactor_10pc = 4 * np.pi * (10 * pc) ** 2
to_ergs = 3631e-23  # maggies to Jansky --> 3631, Jasnky to erg/s/cm^2/Hz -->10**-23


def smooth_spectrum(lam, spec, sigma, minlam=0.0, maxlam=1e50):

    '''
    ripped from Charlie Conroy's smoothspec.f90
    the 'fast way'
    integration is truncated at +/-4*sigma
    '''
    c_kms = 2.99e5
    int_trunc=4
    spec_out = copy.copy(spec)

    for ii in xrange(len(lam)):
        if lam[ii] < minlam or lam[ii] > maxlam:
            spec_out[ii] = spec[ii]
            continue

        dellam = lam[ii]*(int_trunc*sigma/c_kms+1)-lam[ii]
        integrate_lam = (lam > lam[ii]-dellam) & (lam < lam[ii]+dellam)

        if np.sum(integrate_lam) <= 1:
            spec_out[ii] = spec[ii]
        else:
            vel = (lam[ii]/lam[integrate_lam]-1)*c_kms
            func = 1/np.sqrt(2*np.pi)/sigma * np.exp(-vel**2/2./sigma**2)  # s/km
            dx = np.abs(np.diff(vel)) # we want absolute value
            func = func / np.trapz(func, dx=dx)  # (km/s)^-1 / integral([km/s]^-1 wrt km/s)
            spec_out[ii] = np.trapz(func * spec[integrate_lam], dx=dx)  # integral(above * spec_flam, wrt km/s)

    return spec_out


def ems(param_file, out_file, objname='21442', field='cdfs', enames=None):
    res, pr, model = bread.results_from(out_file)

    # get lnprob, based on bread.param_evol()
    chain = res['chain'][..., 0:, :]
    lnprob = res['lnprobability'][..., 0:]
    # deal with single chain (i.e. nested sampling) results
    if len(chain.shape) == 2:
        lnprob = lnprob[None, ...]
    # tracefig, prob = bread.param_evol(res)  # store probability
    # plt.show()
    print('max', lnprob.max())
    row = lnprob.argmax() / len(lnprob[0])
    col = lnprob.argmax() - row * len(lnprob[0])
    walker, iteration = row, col
    print(walker, iteration)

    # Get emission lines!
    # We need the correct sps object to generate models
    sargv = sys.argv
    argdict = {'param_file': param_file}
    clargs = model_setup.parse_args(sargv, argdict=argdict)
    run_params = model_setup.get_run_params(argv=sargv, **clargs)
    sps = model_setup.load_sps(**run_params)
    print('sps')
    # spec, mags, sm = model.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)  # spec [maggies/Hz]
    print('mean model')
    w = sps.wavelengths

    ### save redshift, lumdist
    z = model.params.get('zred', np.array(0.0))
    lumdist = model.params.get('lumdist', np.array(0.0))
    nebinspec = model.params.get('nebemlineinspec', True)
    model.params['zred'] = np.array(0.0)
    if lumdist:
        model.params['lumdist'] = np.array(1e-5)
    if nebinspec == False:
        model.params['nebemlineinspec'] = True

    ### if we want restframe optical photometry, generate fake obs file
    ### else generate NO obs file (don't do extra filter convolutions if not necessary)
    obs = {'filters': [], 'wavelength': None}
    ### calculate SED. comes out as maggies per Hz, @ 10pc
    spec, mags, sm = model.mean_model(res['chain'][walker, iteration, :], obs=obs, sps=sps)  # maggies/Hz at 10pc
    w = sps.wavelengths

    ### reset model
    model.params['zred'] = z
    if lumdist:
        model.params['lumdist'] = lumdist
    if nebinspec == False:
        model.params['nebemlineinspec'] = False

    spec *= dfactor_10pc / constants.L_sun.cgs.value * to_ergs  # spec * cm^2 * (s/erg) * erg/maggies = s*cm^2 / Hz
    # note erg = [g cm^2 / s^2]
    # according to measure_restframe_properties(), above line converts to Lsun / Hz
    to_flam = 3e18 / w ** 2  # for f_nu in erg s^-1 Hz^-1 cm^-2: to_flam [(Ang / s^-1) * (1 / Ang^2)]
    spec_flam = spec * to_flam  # erg cm^-2 s^-1 ang^-1

    smooth_spec = smooth_spectrum(w, spec_flam, 250.0, minlam=3e3, maxlam=7e3)

    ### load fsps emission line list
    loc = os.getenv('SPS_HOME') + '/data/emlines_info.dat'
    dat = np.loadtxt(loc, delimiter=',', dtype={'names': ('lam', 'name'), 'formats': ('f16', 'S40')})
    print('env')

    print(type(enames))
    ### define emission lines
    # legacy code compatible
    if type(enames) == bool:
        lines = np.array(['Hdelta', 'Hbeta', '[OIII]1', '[OIII]2', 'Halpha', '[NII]'])
        fsps_name = np.array(['H delta 4102', 'H beta 4861', '[OIII]4960', '[OIII]5007', 'H alpha 6563', '[NII]6585'])
    else:
        lines = enames
        fsps_name = enames

    print(lines, enames)  # None, None
    # lines = np.array(['Hdelta', 'Hbeta', '[OIII]1', '[OIII]2', 'Halpha', '[NII]'])
    # fsps_name = np.array(['H delta 4102', 'H beta 4861', '[OIII]4960', '[OIII]5007', 'H alpha 6563', '[NII]6585'])
    lines = np.array(['[OII]1', '[OII]2', 'Hdelta', 'Hbeta', '[OIII]1', '[OIII]2', 'Halpha', '[NII]'])
    fsps_name = np.array(['[OII]3726', '[OII]3729', 'H delta 4102', 'H beta 4861', '[OIII]4960', '[OIII]5007',
                          'H alpha 6563', '[NII]6585'])
    ##### measure emission line flux + EQW
    out = {}
    # print(1, sps.emline_wavelengths)  # big long array, 900 to 6*1e6
    for jj in xrange(len(lines)):

        # if we don't do nebular emission, zero this out
        if not hasattr(sps, 'get_nebline_luminosity'):
            print('hi')
            out[lines[jj]] = {'flux': 0.0, 'eqw': 0.0}
            continue

        ### calculate luminosity (in Lsun)
        idx = fsps_name[jj] == dat['name']
        # L_so = 3.84*10**33 erg/s
        print(sps.params['mass'].sum())
        eflux = float(sps.get_nebline_luminosity[idx]*sps.params['mass'].sum())  # L_sun
        elam = float(sps.emline_wavelengths[idx])  # Angstroms!
        print(sps.get_nebline_luminosity[idx])
        print(eflux, 'e luminosity')  # typically 10**40 to 10**42
        print(elam, 'elam')  # Angstroms

        # simple continuum estimation
        tidx = np.abs(sps.wavelengths-elam) < 100  # True for sps.wavelengths within 100 Angstroms of elam wavelengths
        eqw = eflux / np.median(smooth_spec[tidx])  # (erg/s) / (erg cm^-2 s^-1 Ang^-1) = Ang * cm**2

        eflux *= 3.84*10**33  # erg/s (Lum, not flux)

        out[lines[jj]] = {'flux': eflux, 'eqw': eqw}

    return out, fsps_name


if __name__ == "__main__":
    # no args
    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    fs_text = 30
    fs_ticks = 25
    m_size = 10

    param_file = 'eelg_fifty_params.py'

    sets = [['Comp_10.dat', 'out_efico/', 'purple'], ['lbg_ids1', 'out_nfico/', 'blue']]
    # dataset = ['Comp_10.dat', 'lbg_ids1']
    # foldset = ['out_efico/', 'out_nfico/']
    for set in sets:
        data = set[0]
        fold = set[1]
        color = set[2]  # for plotting

        with open(data, 'r') as comp:
            objs = []
            fs = []
            for line in comp:
                if line[0] == '#':
                    pass
                else:
                    cols = line.split()
                    if int(cols[0]) - 200000 > 0:
                        objs.append(str(int(cols[0]) - 200000))
                        fs.append('uds')
                    elif int(cols[0]) - 100000 > 0:
                        objs.append(str(int(cols[0]) - 100000))
                        fs.append('cosmos')
                    else:
                        objs.append(str(int(cols[0])))
                        fs.append('cdfs')
        # 12105, 11462, 12533, 12552, 12903, 14808, 15124, 17189, 17342, 18561, 18742, 21076, 21442, 22768, 11063,
        # 17423, 8787, 15462
        lines = np.zeros(shape=(len(objs), 8))  # [OII]1, [OII]2, Hdelta, Hbeta, [OIII]1, [OIII]2, Halpha, [NII] = 8
        width = np.zeros(shape=(len(objs), 8))
        for i in range(len(objs)):
            for infile in glob.glob(os.path.join('/home/jonathan/.conda/envs/snowflakes/lib/python2.7/' +
                                                 'site-packages/prospector/git/out/' + fold, objs[i] + '*.h5')):
                out_file = infile
                print(param_file, out_file)
                true_field = fs[i]
                true_obj = objs[i]

                out, line_names = ems(param_file, out_file, objname=true_obj, field=true_field)
                print(out)

            fluxes = []
            eqws = []
            for idx in out:  # for each line in out
                fluxes.append(out[idx]['flux'])  # that line's flux
                print(fluxes[-1], out[idx], out[idx]['flux'])
                eqws.append(out[idx]['eqw'])  # that line's eqw
                print(eqws[-1])
            lines[i, :] = fluxes
            print('fluxes', fluxes)
            width[i, :] = eqws
            # print('lines', lines)
        print('lines', lines)
        print('width', width)
        ratio = []
        for i in range(len(objs)):
            ratio.append((lines[i][4] + lines[i][5]) / (lines[i][0] + lines[i][1]))  # OIII doublet / OII doublet

        print(ratio, 'ratio')
        print(np.percentile(ratio, [16., 50., 84.]))

        names = ['[NII]6585', '[OII]3726', 'H alpha 6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'H beta 4861',
                 'H delta 4102']
        # names = ['[NII]', '[OII]1', 'Halpha', '[OII]2', '[OIII]1', '[OIII]2', 'Hbeta', 'Hdelta']
        for j in range(len(width)):  # for each galaxy
            Ha = width[j, 2]
            OIII = width[j, 5]
            ax1.plot(Ha, OIII, color=color, marker='o', markersize=m_size)
        percs_Ha = np.percentile(width[:, 2], [16., 50., 84.])
        percs_OIII = np.percentile(width[:, 5], [16., 50., 84.])
        h_err = [[percs_Ha[1] - percs_Ha[0]], [percs_Ha[2] - percs_Ha[1]]]
        o_err = [[percs_OIII[1] - percs_OIII[0]], [percs_OIII[2] - percs_OIII[1]]]
        print(h_err)
        print(o_err)
        ax1.errorbar(x=percs_Ha[1], y=percs_OIII[1], yerr=o_err, xerr=h_err, color=color, markerfacecolor='None',
                     markeredgecolor=color, markeredgewidth=1.25, linewidth=2, markersize=m_size)

    ax1.set_xlabel(r'H$\alpha$ Equivalent Width [$\rm \AA$]', fontsize=fs_text)
    ax1.set_ylabel(r'[OIII] Equivalent Width [$\rm \AA$]', fontsize=fs_text)
    ax1.set_xlim(0., 10**3)
    ax1.set_xticks([0, 200, 400, 600, 800, 1000])  # technically works
    ax1.set_xticklabels([r'$0$', r'$200$', r'$400$', r'$600$', r'$800$', r'$1000$'], size=fs_ticks)
    ax1.set_yticks([0, 200, 400, 600, 800, 1000])  # technically works
    ax1.set_yticklabels([r'$0$', r'$200$', r'$400$', r'$600$', r'$800$', r'$1000$'], size=fs_ticks)
    ax1.set_ylim(0., 10**3)
    plt.show()


'''
RUNNING WITH:
python plot_ews.py
'''