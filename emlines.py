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
to_ergs = 3631e-23


def smooth_spectrum(lam,spec,sigma,
                    minlam=0.0,maxlam=1e50):

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
            func = 1/np.sqrt(2*np.pi)/sigma * np.exp(-vel**2/2./sigma**2)
            dx=np.abs(np.diff(vel)) # we want absolute value
            func = func / np.trapz(func,dx=dx)
            spec_out[ii] = np.trapz(func*spec[integrate_lam],dx=dx)

    return spec_out


def ems(param_file, out_file, enames=None):
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
    spec, mags, sm = model.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)
    print('mean model')
    w = sps.wavelengths
    spec *= dfactor_10pc / constants.L_sun.cgs.value * to_ergs
    to_flam = 3e18 / w ** 2
    spec_flam = spec * to_flam
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
    # print(sps.get_nebline_luminosity)  # 1e-1 to 1e-10, mostly 1e-2 to 1e-6 (big long array)
    # print(1, sps.emline_wavelengths)  # big long array, 900 to 6*1e6
    for jj in xrange(len(lines)):

        # if we don't do nebular emission, zero this out
        if not hasattr(sps, 'get_nebline_luminosity'):
            print('hi')
            out[lines[jj]] = {'flux': 0.0, 'eqw': 0.0}
            continue

        ### calculate luminosity (in Lsun)
        idx = fsps_name[jj] == dat['name']
        # print(sps.params['mass'], 'mass')  # e.g. BIG23104018432.0, 10554496000.0 (~1e10)
        print(sps.params['mass'].sum())
        eflux = float(sps.get_nebline_luminosity[idx]*sps.params['mass'].sum() * 3.848*10**33)  # erg/s per L_sol
        elam = float(sps.emline_wavelengths[idx])
        print(sps.get_nebline_luminosity[idx])  # 1e-1 to 1e-10, typically ~1e-2 to 1e-6
        print(eflux, 'eflux')
        print(elam, 'elam')

        # simple continuum estimation
        tidx = np.abs(sps.wavelengths-elam) < 100
        eqw = eflux / np.median(smooth_spec[tidx])

        out[lines[jj]] = {'flux': eflux, 'eqw': eqw}

    return out


if __name__ == "__main__":
    # don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--parfile')
    parser.add_argument('--outname')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    files = {'outname': '', 'parfile': ''}
    for key in kwargs.keys():
        files[key] = kwargs[key]

    if files['outname'] == 'all':
        data = 'lbg_ids1'  # 'lbg_ids1'  # 'Comp_10.dat'
        fold = 'out_nvary/'  # 'out_nvary/'  # 'out_evar/'
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

                out = ems(param_file, out_file)
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