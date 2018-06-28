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


obj = 11063  #12533  # 15124  # 20366  # 12533  # 11462  # 12552  # 12105  # 21442
field = 'cosmos'  # 'cdfs'# 'cosmos'  # 'cdfs'  # 'cosmos'  # 'cdfs'
'''
sfh_style = str(obj) + '_set8'
set9 = 'set9_94_04_17_68_07'
set8 = 'set8_99_03_006_30_05'  # cdfs 12533
set7 = 'set7_94_03_17_66_05'
set6 = 'set6_94_03_17_15_05'
set5 = 'set5_97_04_20_20_05'
set4 = 'set4_95_03_17_60_05'
set3 = 'set3_95_00_10_40_05'
set2 = 'set2_10_03_17_30_10'  # log(M)=10, dust=0.3, log(Z_sol)=-1.7, sfr_frac_msotrecent=0.30, sfr_fracsecond=0.10
set1 = 'set1_9_03_17_90_05'
# =logmass of 9., dust2 of 0.3, logzsol of -1.7, 0.90 mfrac in most recent bin, 0.05 mfrac in second most recent bin

pset = set8
'''

masses = np.linspace(9., 9.8, 5)  # [9. 9.1 ... 9.8]  (...9); with 5 instead: 9., 9.2, 9.4, 9.6, 9.8
dusts = np.linspace(0.15, 0.55, 3)  # [0.1 0.2 ... 0.6]; with 5 instead (0.15:0.55): 0.15, 0.25, ..., 0.55
# 3: 0.15, 0.35, 0.55
mets = np.linspace(-2., -1.5, 3)  # [-2. -1.75 -1.5]; (with5: [-2 -1.75 .. -1.])
sfr1 = np.linspace(0.05, 0.45, 5)  # [0.1 0.2 ... 0.6]  # 5-45%
# sfr2 = np.linspace(0.05, 0.2, 4)  # [0.05 0.1 0.15 0.2]  # NOTE: keep sfr2 at 0.05
newphot = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_all' #\
#           + str(int(masses[m] * 10)) + '_' + str(int(dusts[d] * 100)) + '_' + str(int(mets[z] * -10)) \
#           + '_' + str(int(sfr1[f] * 100))

photname, zname, filters = get_names(field)  # FIELD

with open(newphot, 'w+') as new:
    new.write('# id ')
    for i in range(len(filters)):
        new.write('f_' + filters[i] + ' ')
    new.write('\n')

# GENERATE PHOTOMETRY FROM ALL THE PARAM FILES I JUST MADE WITH GENERATE_PARAM.PY
for m in range(len(masses)):
    print(m, 'm!')
    for d in range(len(dusts)):
        for z in range(len(mets)):
            for f in range(len(sfr1)):
                pset = str(int(masses[m] * 10)) + '_' + str(int(dusts[d] * 100)) + '_' + str(int(mets[z] * -10)) \
                       + '_' + str(int(sfr1[f] * 100))
                print(pset)
                parfile = 'testpars/sfhtest_' + pset + '_params.py'

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

                mu, phot, x = global_model.mean_model(global_model.initial_theta, global_obs, sps=sps)

                phot *= 3631 * 10 ** 6  # maggies to uJy

                # print(diff)
                with open(newphot, 'a') as new:
                    new.write(pset + ' ')  # new.write(sfh_style + ' ')  # (str(obj) + ' ')
                    for k in range(len(phot)):
                        new.write(str(phot[k]) + ' ')
                    new.write('\n')

print('done!')
# print('generated phot', phot)
# print(str(obj) + ' flux', flux)
# print('cat flux', flux)

'''
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
'''
