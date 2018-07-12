import numpy as np
import sys
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log


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

# newphot = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_all'
# photbase = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_'
'''
photbase = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/tphot_'  # 'eephot_'
'''
photbase = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/eephot2_'
photc = photbase + 'cosmos'
photu = photbase + 'uds'
photf = photbase + 'cdfs'

photname, zname, cfilters = get_names('cosmos')  # FIELD
uphotname, uzname, ufilters = get_names('uds')  # FIELD
fphotname, fzname, ffilters = get_names('cdfs')  # FIELD
with open(photc, 'w+') as new:
    new.write('# id ')
    for i in range(len(cfilters)):
        new.write('f_' + cfilters[i] + ' ')
    new.write('\n')
with open(photu, 'w+') as new:
    new.write('# id ')
    for i in range(len(ufilters)):
        new.write('f_' + ufilters[i] + ' ')
    new.write('\n')
with open(photf, 'w+') as new:
    new.write('# id ')
    for i in range(len(ffilters)):
        new.write('f_' + ffilters[i] + ' ')
    new.write('\n')
# GENERATE PHOTOMETRY FROM ALL THE PARAM FILES I JUST MADE WITH GENERATE_PARAM.PY

'''
for x in range(100):  # 200
'''
for x in range(100):
    # parfile = 'newtest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    # 'eetest/sfhtest_'
    '''
    parfile = 'ctest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    '''
    parfile = 'eetest2/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    pset = None
    field = ''
    with open(parfile, 'r') as pfile:
        for line in pfile:
            if line.startswith('# Codename: '):
                pset = line[12:]
            elif line.startswith('# Identity: '):
                counterl = 0
                for l in line:
                    if l == ' ' or l == '_':
                        counterl += 1
                    elif counterl == 2:
                        field += l
    print('pset', pset)
    print('field', field)
    print(parfile, x)

    if field == 'cosmos':
        usephot = photc
    elif field == 'cdfs':
        usephot = photf
    elif field == 'uds':
        usephot = photu
    else:
        print('WARNING', 'field = ' + field)

    sargv = sys.argv
    argdict = {'param_file': parfile}
    clargs = model_setup.parse_args(sargv, argdict=argdict)
    run_params = model_setup.get_run_params(argv=sargv, **clargs)
    print('runpars')
    # --------------
    # Globals
    # --------------
    # SPS Model instance as global
    sps = model_setup.load_sps(**run_params)
    # GP instances as global
    # spec_noise, phot_noise = model_setup.load_gp(**run_params)
    # Model as global
    global_model = model_setup.load_model(**run_params)
    # Obs as global
    global_obs = model_setup.load_obs(**run_params)

    mu, phot, other = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)

    phot *= 3631 * 10 ** 6.44  # maggies to catalog fluxes, so I can use the cat fluxes --> maggies code in param files

    # print(diff)
    with open(usephot, 'a') as new:
        new.write(str(x) + ' ')
        for k in range(len(phot)):
            new.write(str(phot[k]) + ' ')
        new.write('\n')
        print('look!')

print('done!')
