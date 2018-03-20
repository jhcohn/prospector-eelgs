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

'''
After performing a model call, the SPS object has the emission line information associated with that theta vector stored
as the following:

sps.ssp.emline_luminosity = Luminosity in Lsun
sps.ssp.emline_wavelengths = wavelengths in Angstroms

I'm 95% sure these are rest-frame luminosities (a good check: change the redshift, then do the mean_model() call again,
these values should NOT change).

There's an FSPS line list which you can load and use to identify these lines (Halpha, etc) with the following code:

import os
loc = os.getenv('SPS_HOME')+'/data/emlines_info.dat'
dat = np.loadtxt(loc, delimiter=',', dtype = {'names':('lam','name'),'formats':('f16','S40')})


def measure_emlines(smooth_spec,sps,enames=None):
    """ emission line fluxes are part of SPS output now. this is
    largely present to measure the continuum for EQW calculations
    """

    ### load fsps emission line list
    loc = os.getenv('SPS_HOME')+'/data/emlines_info.dat'
    dat = np.loadtxt(loc, delimiter=',',
                     dtype = {'names':('lam','name'),'formats':('f16','S40')})

    ### define emission lines
    # legacy code compatible
    if type(enames) == bool:
        lines = np.array(['Hdelta','Hbeta','[OIII]1','[OIII]2','Halpha','[NII]'])
        fsps_name = np.array(['H delta 4102','H beta 4861','[OIII]4960','[OIII]5007','H alpha 6563','[NII]6585'])
    else:
        lines = enames
        fsps_name = enames

    ##### measure emission line flux + EQW
    out = {}
    for jj in xrange(len(lines)):

        # if we don't do nebular emission, zero this out
        if not hasattr(sps, 'get_nebline_luminosity'):
            out[lines[jj]] = {'flux':0.0,'eqw':0.0}
            continue

        ### calculate luminosity (in Lsun)
        idx = fsps_name[jj] == dat['name']
        eflux = float(sps.get_nebline_luminosity[idx]*sps.params['mass'].sum())
        elam = float(sps.emline_wavelengths[idx])

        # simple continuum estimation
        tidx = np.abs(sps.wavelengths-elam) < 100
        eqw = eflux/np.median(smooth_spec[tidx])

        out[lines[jj]] = {'flux':eflux,'eqw':eqw}

    return out
'''
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
    # don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--parfile')
    parser.add_argument('--outname')

    eels = 0  # 1 = EELGs, 0 = SFGs

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    files = {'outname': '', 'parfile': ''}
    for key in kwargs.keys():
        files[key] = kwargs[key]

    if files['outname'] == 'all':
        if eels:
            data = 'Comp_10.dat'  # 'lbg_ids1'  # 'Comp_10.dat'
            fold = 'out_efico/'  # 'out_nfico/'  # 'out_efico/'
        else:
            data = 'lbg_ids1'  # 'lbg_ids1'  # 'Comp_10.dat'
            fold = 'out_nfico/'  # 'out_nfico/'  # 'out_efico/'
        with open(data, 'r') as comp:
            e_objs = []
            e_fs = []
            for line in comp:
                if line[0] == '#':
                    pass
                else:
                    cols = line.split()
                    if int(cols[0]) - 200000 > 0:
                        e_objs.append(str(int(cols[0]) - 200000))
                        e_fs.append('uds')
                    elif int(cols[0]) - 100000 > 0:
                        e_objs.append(str(int(cols[0]) - 100000))
                        e_fs.append('cosmos')
                    else:
                        e_objs.append(str(int(cols[0])))
                        e_fs.append('cdfs')
        # 12105, 11462, 12533, 12552, 12903, 14808, 15124, 17189, 17342, 18561, 18742, 21076, 21442, 22768, 11063,
        # 17423, 8787, 15462
        lines = np.zeros(shape=(len(e_objs), 8))  # [OII]1, [OII]2, Hdelta, Hbeta, [OIII]1, [OIII]2, Halpha, [NII] = 8
        width = np.zeros(shape=(len(e_objs), 8))
        for i in range(len(e_objs)):
            for infile in glob.glob(os.path.join('/home/jonathan/.conda/envs/snowflakes/lib/python2.7/' +
                                                 'site-packages/prospector/git/out/' + fold, e_objs[i] + '*.h5')):
                out_file = infile
                param_file = files['parfile']
                print(param_file, out_file)
                true_field = e_fs[i]
                true_obj = e_objs[i]

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
        for i in range(len(e_objs)):
            ratio.append((lines[i][4] + lines[i][5]) / (lines[i][0] + lines[i][1]))  # OIII doublet / OII doublet

        print(ratio, 'ratio')
        print(np.percentile(ratio, [16., 50., 84.]))

        fs_text = 15

        fig = plt.figure()
        ax1 = plt.subplot(3, 3, 1)
        ax2 = plt.subplot(3, 3, 2)
        ax3 = plt.subplot(3, 3, 3)
        ax4 = plt.subplot(3, 3, 4)
        ax5 = plt.subplot(3, 3, 5)
        ax6 = plt.subplot(3, 3, 6)
        ax7 = plt.subplot(3, 3, 7)
        ax8 = plt.subplot(3, 3, 8)

        axes = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]

        names = ['[NII]6585', '[OII]3726', 'H alpha 6563', '[OII]3729', '[OIII]4960', '[OIII]5007', 'H beta 4861',
                 'H delta 4102']
        # names = ['[NII]', '[OII]1', 'Halpha', '[OII]2', '[OIII]1', '[OIII]2', 'Hbeta', 'Hdelta']
        for i in range(len(names)):  # for each line
            percs = np.percentile(width[:, i], [16., 50., 84.])
            axes[i].hist(width[:, i], bins=10, histtype='step', weights=[1. / len(width[:, i])] * len(width[:, i]),
                         normed=False, color='k', lw=2, label=names[i])  # plot hist: this line for all galaxies
            axes[i].axvline(x=percs[1], color='k', linestyle='--', lw=2, label='Median')
            axes[i].legend(numpoints=1, loc='upper right', prop={'size': fs_text})
            axes[i].set_xlabel(r'Equivalent Width [$\rm \AA$]', fontsize=fs_text)
            if i == 3:
                axes[i].set_ylabel(r'Fraction', fontsize=fs_text)
            axes[i].set_xlim(0., 10**3)
            axes[i].set_ylim(0., 0.70)
        plt.show()
        plt.close()

    else:
        out_file = files['outname']
        param_file = files['parfile']
        print(param_file, out_file)

        out = ems(param_file, out_file)
        print(out)

'''
RUNNING WITH:
python emlines.py --outname=17342_cdfs_evar2_1518635733_mcmc.h5 --parfile=eelg_varymet_params.py

{'[NII]': {'flux': 29985579.87637427, 'eqw': 4.3604602484966949e+19}, 'Halpha': {'flux': 893603473.0211837, 'eqw':
1.2987595240584614e+21}, '[OIII]1': {'flux': 359786255.8684004, 'eqw': 2.6482665204910591e+20},
'[OIII]2': {'flux': 1085450530.7686644, 'eqw': 8.1265190031614712e+20},
'Hbeta': {'flux': 288758510.9395125, 'eqw': 2.0257050676411587e+20}, 'Hdelta': {'flux': 73832632.02298185,
'eqw': 3.2602449229181383e+19}}

'''