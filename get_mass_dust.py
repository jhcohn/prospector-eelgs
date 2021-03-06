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
np.errstate(invalid='ignore')


def md(sample_results, start=0, thin=1, percs=True, sfhtest=False, masstest=False, fast=False, five=False, quiet=False,
       draw1=False, fixmet=False):
    """Make a triangle plot of the (thinned, latter) samples of the posterior
    parameter space.  Optionally make the plot only for a supplied subset of
    the parameters.

    :param start:
        The iteration number to start with when drawing samples to plot.

    :param thin:
        The thinning of each chain to perform when drawing samples to plot.
    """
    # pull out the parameter names and flatten the thinned chains
    try:
        parnames = np.array(sample_results['theta_labels'])
    except(KeyError):
        parnames = np.array(sample_results['model'].theta_labels())
    flatchain = sample_results['chain'][:, start::thin, :]
    flatchain = flatchain.reshape(flatchain.shape[0] * flatchain.shape[1],
                                  flatchain.shape[2])

    # print(parnames)
    # print(np.percentile(flatchain[:, 0], [16, 50, 84]), len(flatchain[0]))  len(flatchain[0]) = 8 = len(parnames)
    # print(parnames) = [u'logmass' u'sfr_fraction_1' u'sfr_fraction_2' u'sfr_fraction_3' u'sfr_fraction_4'
    # u'sfr_fraction_5' u'dust2' u'logzsol' u'gas_logz']
    # IF METALLICITY: print(parnames) = [... u'dust2' u'logzsol' u'gaslogz']
    if percs:
        mass = np.percentile(flatchain[:, 0], [16, 50, 84])[1]
        dust = np.percentile(flatchain[:, 6], [16, 50, 84])[1]
        metal = None
        gasmet = None
        if parnames[7] == u'logzsol':
            metal = np.percentile(flatchain[:, 7], [16, 50, 84])[1]
            gasmet = np.percentile(flatchain[:, -1], [16, 50, 84])[1]
        ret = [mass, dust, metal, gasmet]
    elif fixmet:
        mass = np.random.choice(flatchain[:, 0], size=(10**3))
        dust = np.random.choice(flatchain[:, 6], size=(10**3))
        ret = np.asarray([mass, dust])
    elif sfhtest and draw1:
        mass = np.random.choice(flatchain[:, 0], size=(10**3))
        ssfr1 = np.random.choice(flatchain[:, 1], size=(10**3))
        ssfr2 = np.random.choice(flatchain[:, 2], size=(10**3))
        dust = np.random.choice(flatchain[:, 6], size=(10**3))
        gasmet = np.random.choice(flatchain[:, -1], size=(10**3))
        metal = np.random.choice(flatchain[:, 7], size=(10**3))
        ret = np.asarray([mass, dust, metal, gasmet, ssfr1, ssfr2])
    elif masstest:
        mass = np.percentile(flatchain[:, 0], [16, 50, 84])[1]
        dust = np.percentile(flatchain[:, 6], [16, 50, 84])[1]
        gasmet = np.percentile(flatchain[:, -1], [16, 50, 84])[1]
        metal = np.percentile(flatchain[:, 7], [16, 50, 84])[1]
        sfh = []
        ret = [mass]
        for i in [1, 2, 3, 4, 5]:
            sfh.append(np.percentile(flatchain[:, i], [16, 50, 84])[1])
            ret.append(sfh[i - 1])
        ret.append(dust)
        ret.append(metal)
        ret.append(gasmet)
    elif draw1 and fast:
        mass = np.random.choice(flatchain[:, 0], size=(10**3))
        tage = np.random.choice(flatchain[:, 1], size=(10**3))
        logtau = np.random.choice(flatchain[:, 2], size=(10**3))
        dust = np.random.choice(flatchain[:, 3], size=(10**3))
        ret = np.asarray([mass, tage, logtau, dust])
    elif draw1:
        mass = np.random.choice(flatchain[:, 0], size=(10**3))
        dust = np.random.choice(flatchain[:, 6], size=(10**3))
        gasmet = np.random.choice(flatchain[:, -1], size=(10**3))
        metal = np.random.choice(flatchain[:, 7], size=(10**3))
        ret = np.asarray([mass, dust, metal, gasmet])
    elif fast:
        mass = np.percentile(flatchain[:, 0], [16, 50, 84])[1]
        tage = np.percentile(flatchain[:, 1], [16, 50, 84])[1]
        logtau = np.percentile(flatchain[:, 2], [16, 50, 84])[1]
        dust = np.percentile(flatchain[:, 3], [16, 50, 84])[1]
        ret = np.asarray([mass, tage, logtau, dust])
    elif five:
        mass = np.percentile(flatchain[:, 0], [16, 50, 84])[1]
        dust = np.percentile(flatchain[:, 6], [16, 50, 84])[1]
        ret = np.asarray([mass, dust])
    else:
        mass = flatchain[:, 0]
        if not quiet:
            print(len(mass), 'mass')
        dust = flatchain[:, 6]
        gasmet = flatchain[:, -1]
        metal = None
        if parnames[7] == u'logzsol':
            metal = flatchain[:, 7]
        ret = mass
    return ret


def printer(out_file, percs=True, sfhtest=False, masstest=False, fast=False, five=False, quiet=False, draw1=False,
            fixmet=False):
    if not quiet:
        print(out_file)
        print(masstest)
    res, pr, mod = bread.results_from(out_file)
    # ''' #

    # PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
    return md(res, start=-1000, thin=5, percs=percs, sfhtest=sfhtest, masstest=masstest, fast=fast, five=five,
              quiet=quiet, draw1=draw1, fixmet=fixmet)  # -650
    # set start by when kl converges!  # returns mass, dust, stellar Z, gas Z


if __name__ == "__main__":

    corr = 0
    fico = 0
    news = 0
    tenmet = 0
    fast = 0
    masstest = 0
    fivemet = 0
    tt = 0
    sfhtest = 0
    tenmet = 1

    massmatch = 0

    if corr:
        folders = ['out_ecorr/', 'out_ncorr/']
        pars = ['eelg_varymet_params.py', 'eelg_varymet_params.py']
        base = ['corr', 'corr']
    elif fico:
        folders = ['out_efico/', 'out_nfico/']
        pars = ['eelg_fifty_params.py', 'eelg_fifty_params.py']
        base = ['fico', 'fico']
    elif tenmet:
        folders = ['out_etenmet/', 'out_nfico/']
        pars = ['eelg_tenmet_params.py', 'eelg_fifty_params.py']
        base = ['tenmet', 'fico']
    elif fivemet:
        folders = ['out_efivemet/', 'out_efico/']
        pars = ['eelg_fivemet_params.py', 'eelg_fifty_params.py']
        base = ['fivemet', 'fico']
    elif tenmet:
        folders = ['out_etenmet/', 'out_nfico/']
        pars = ['eelg_tenmet_params.py', 'eelg_fifty_params.py']
        base = ['tenmet', 'fico']
    elif fast:
        folders = ['out_efast/', 'out_nfico/']
        pars = ['eelg_fast_params.py', 'eelg_fifty_params.py']
        base = ['fast', 'fico']
    elif news:
        folders = ['out_efico/', 'out_nnewsfg/']
        pars = ['eelg_fifty_params.py', 'eelg_fifty_params.py']
        base = ['fico', 'newsfg']
    elif masstest:
        folders = ['out_masstest/', 'out_efico/']
        pars = ['eelg_masstest_params.py', 'eelg_fifty_params.py']
        base = ['masstest', 'fico']
        obj = '12105'  # '12552'  # '12105'  # '21442'
    elif sfhtest:
        folders = ['out_simsfh/', 'out_sfhtest/']
        pars = ['eelg_simsfh_params.py', 'eelg_test0_params.py']
        base = ['9_03_17_15_05', 'sfhtest']
        obj = '11063'

    eelgs1 = []
    oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[0]
    for file in os.listdir(oute):
        if masstest or sfhtest:
            if file.endswith(".h5") and file.startswith(obj):
                eelgs1.append(file)
        else:
            if file.endswith(".h5"):
                eelgs1.append(file)

    # LIST GALAXIES IN EELG COMPOSITE SED
    eelg_list = open('eelg_specz_ids1', 'r')
    comp = []
    for line in eelg_list:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            if int(cols[0]) - 200000 > 0:
                comp.append(str(int(cols[0]) - 200000) + '_uds_' + base[0])  # base[1] = noelg (or nother)
            elif int(cols[0]) - 100000 > 0:
                comp.append(str(int(cols[0]) - 100000) + '_cosmos_' + base[0])  # base[1] = noelg (or nother)
            else:
                comp.append(str(int(cols[0])) + '_cdfs_' + base[0])  # base[1] = noelg (or nother)
    eelg_list.close()

    # MAKE SURE ALL THE GALAXIES IN THE COMPOSITE SED HAVE BEEN RUN THROUGH PROSPECTOR
    eelgs = []
    count = 0
    for file in eelgs1:
        for i in range(len(comp)):
            if file.startswith(comp[i]):
                count += 1
                print(file)
                eelgs.append(file)
    print(count)

    # LIST GALAXIES IN SFG COMPOSITE SED
    lbgs = []
    outl = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[1]
    countl = 0
    if masstest or sfhtest:
        for file in os.listdir(outl):
            if file.endswith(".h5") and file.startswith(obj):
                lbgs.append(file)
                countl += 1
    else:
        for file in os.listdir(outl):
            if file.endswith(".h5"):
                lbgs.append(file)
                countl += 1
    print(countl)

    # IF COMPARING MASS-MATCHED SUBSAMPLES
    if massmatch:
        mass1, dust1, met1, ssfrs1 = [], [], [], []
        counter = 0
        counter2 = 0
        lims = [9.7, 9.9]  # set mass range on which to compare (ranges are ~equivalent at 9.7 < log(M) < 9.9)
        for i in range(len(eelgs)):
            # print(eelgs[i], 'look!')
            get = printer(oute + eelgs[i], percs=True)  # median mass, dust, gasmet, metallicity

            # IF mass is within limits (i.e.if part of mass-matched subsample)
            if lims[0] <= get[0] <= lims[1]:  # 9.55 <= x1[i] <= 9.9:
                counter += 1
                mass1.append(get[0])
                dust1.append(get[1])
                met1.append(get[2])

                with open('eelg_ssfrs.txt', 'r') as ssfr:
                    for line in ssfr:
                        if line[0] != '#':
                            cols = line.split()
                            # print(cols[0], eelgs[i])
                            if eelgs[i].startswith(cols[0]):
                                # print(cols[0], 'hey!')
                                ssfrs1.append(float(cols[1]) * 10 ** 9)  # ssfrs are stored in 2nd column in file
                            # convert to Gyr^-1
                            # counter += 1
        mass2, dust2, met2, ssfrs2 = [], [], [], []
        for i in range(len(lbgs)):
            get = printer(outl + lbgs[i], percs=True)  # median mass, dust, gasmet, metallicity

            if lims[0] <= get[0] <= lims[1]:  # 9.55 <= x1[i] <= 9.9:
                counter2 += 1
                mass2.append(get[0])
                dust2.append(get[1])
                met2.append(get[2])

                with open('sfg_ssfrs.txt', 'r') as ssfr:
                    for line in ssfr:
                        if line[0] != '#':
                            cols = line.split()
                            if lbgs[i].startswith(cols[0]):
                                ssfrs2.append(float(cols[1]) * 10 ** 9)  # ssfrs are stored in 2nd column in file
                            # convert to Gyr^-1
                            # counter += 1

        print(counter)
        print('mass1', np.percentile(mass1, [16., 50., 84.]))
        print('dust1', np.percentile(dust1, [16., 50., 84.]))
        print('met1', np.percentile(met1, [16., 50., 84.]))
        print('ssfrs1', np.percentile(ssfrs1, [16., 50., 84.]))
        print(counter2)
        print('mass2', np.percentile(mass2, [16., 50., 84.]))
        print('dust2', np.percentile(dust2, [16., 50., 84.]))
        print('met2', np.percentile(met2, [16., 50., 84.]))
        print('ssfrs2', np.percentile(ssfrs2, [16., 50., 84.]))

    else:  # ELSE IF NOT MASS-MATCH
        if fivemet:  # if using fixed met: Z_sol / 5
            get_e = np.zeros(shape=(2, len(eelgs)))  # 2 rows (dust, mass), each row as long as eelgs
        else:
            get_e = np.zeros(shape=(4, len(eelgs)))  # 4 rows (dust, mass, gaslogz, logzsol), each row as long as eelgs
        onedraw = np.zeros(shape=(4, len(eelgs), 10**3))  # *10**3))
        for i in range(len(eelgs)):
            if os.path.exists(oute + eelgs[i]):
                if fast:
                    print('hmm')
                    get_e[:, i] = printer(oute + eelgs[i], percs=False, fast=fast)
                    print('hmm')
                elif fivemet:
                    get_e[:, i] = printer(oute + eelgs[i], percs=False, fast=fast, five=fivemet)
                else:
                    get_e[:, i] = printer(oute + eelgs[i])
                    print(eelgs[i], get_e[0, i])
                    # for k in range(10**3):
                    onedraw[:, i, :] = printer(oute + eelgs[i], percs=False, fast=fast, five=fivemet, draw1=True)

        masse = np.percentile(get_e[0], [16., 50., 84.])
        duste = np.percentile(get_e[1], [16., 50., 84.])
        print(masse)
        print(duste)
        mete = np.percentile(get_e[2], [16., 50., 84.])
        gasmete = np.percentile(get_e[3], [16., 50., 84.])

        # '''
        get_l = np.zeros(shape=(4, len(lbgs)))  # 4 rows (dust, mass, gaslogz, logzsol), each row as long as lbgs
        onedraw_l = np.zeros(shape=(4, len(lbgs), 10**3))  # * 10**3))
        for i in range(len(lbgs)):
            if os.path.exists(outl + lbgs[i]):
                get_l[:, i] = printer(outl + lbgs[i])
                # for k in range(10**3):
                onedraw_l[:, i, :] = printer(outl + lbgs[i], percs=False, draw1=True)

        means = np.zeros(shape=(4, 10**3))
        for m in range(4):
            for k in range(10**3):
                means[m, k] = np.mean(onedraw[m, :, k])
        means_l = np.zeros(shape=(4, 10**3))
        for m in range(4):
            for k in range(10**3):
                means_l[m, k] = np.mean(onedraw_l[m, :, k])

        mass = np.percentile(get_l[0], [16., 50., 84.])
        dust = np.percentile(get_l[1], [16., 50., 84.])
        met = np.percentile(get_l[2], [16., 50., 84.])
        gasmet = np.percentile(get_l[3], [16., 50., 84.])
        # '''

        print('masse', masse)
        print('duste', duste)
        print('mete', mete)
        print('gasmete', gasmete)
        # '''
        print('mass', mass)
        print('dust', np.percentile(get_l[1], [16., 50., 84.]))
        print('met', np.percentile(get_l[2], [16., 50., 84.]))
        print('gasmet', np.percentile(get_l[3], [16., 50., 84.]))
        # '''

        print('errors on means:')
        params = ['mass', 'dust', 'met', 'gasmet']
        for eom in range(4):
            print(params[eom] + 'e', np.percentile(means[eom, :], [16., 50., 84.]))
        for eoml in range(4):
            print(params[eoml] + 'l', np.percentile(means_l[eoml, :], [16., 50., 84.]))
        #print('masse', np.percentile(onedraw[0,:,:], [16., 50., 84.]))
        #print('duste', np.percentile(onedraw[1,:,:], [16., 50., 84.]))
        #print('mete', np.percentile(onedraw[2,:,:], [16., 50., 84.]))
        #print('gasmete', np.percentile(onedraw[3,:,:], [16., 50., 84.]))
        #e_on_mean = np.zeros(shape=(4, len(eelgs), 3))
        #print('massl', np.percentile(onedraw_l[0,:,:], [16., 50., 84.]))
        #print('dustl', np.percentile(onedraw_l[1,:,:], [16., 50., 84.]))
        #print('metl', np.percentile(onedraw_l[2,:,:], [16., 50., 84.]))
        #print('gasmetl', np.percentile(onedraw_l[3,:,:], [16., 50., 84.]))

'''
RUNNING WITH:
python get_mass_dust.py

For dust2 values: convert to A_V by multiplying by 1.86

# TENMET
('masse', array([ 8.88606308,  9.2459774 ,  9.6234944 ]))
('duste', array([ 0.11443192,  0.28719273,  0.39732726]))
('mete', array([-0.96755197, -0.90360114, -0.80587448]))
('gasmete', array([-0.39200348, -0.33855876, -0.30372495]))
('mass', array([  9.90069426,  10.14381647,  10.44132971]))
('dust', array([ 0.33243034,  0.55183363,  0.76055876]))
('met', array([-1.86135451, -1.71223199, -1.26867414]))
('gasmet', array([-1.2800012 , -0.62593657, -0.00709531]))




# MASSTEST
# 12105 --> less mass?!
# 21442 --> slightly less mass, much more dust
# 12552 --> 1e8.5 more mass, much more dust, higher metallicity




# FIXEDMET:
('masse', array([  9.7307283 ,   9.98464298,  10.25774393]))
('duste', array([ 0.16936105,  0.30668478,  0.4184224 ])
('mass', array([  9.98807951,  10.25705242,  10.46080221]))
('dust', array([ 0.27746924,  0.491308  ,  0.58165465]))

# FOR VARY:
('masse', array([  9.91871967,  10.15146065,  10.3976491 ]))
('mass', array([ 10.24227583,  10.51341248,  10.74591787]))
('duste', array([ 0.11079067,  0.26369962,  0.43304285]))
('dust', array([ 0.24633007,  0.46240632,  0.62991244]))
('mete', array([-1.82372492, -1.66893744, -1.09977767]))
('met', array([-1.8266431 , -1.63226652, -1.35153838]))
('gasmete', array([-0.68674988, -0.38936704, -0.20108394])) 0.41 -0.20 +0.22 [0.206, 0.408, 0.63]
('gasmet', array([-1.31453426, -0.80460706, -0.19746851])) [0.004847 0.16 635] 0.16 -0.11 +0.48

# VARY ONE COMPOSITE!
# C14:
('masse', array([ 10.06390381,  10.25686312,  10.45410967]))
('duste', array([ 0.1176697 ,  0.23385095,  0.43429242]))
('mete', array([-1.80541593, -1.58386713, -1.02036411]))
('gasmete', array([-1.05025488, -0.59004486, -0.00989608]))
# C10:
('masse', array([  9.56806931,   9.9428606 ,  10.26535469]))
('duste', array([ 0.1366495 ,  0.29065425,  0.4745783 ]))
('mete', array([-1.85909033, -1.72481412, -1.57049279]))
('gasmete', array([-0.39605677, -0.34745189, -0.29891892]))
# C04:
('masse', array([  9.91665936,  10.11143303,  10.37669609]))
('duste', array([ 0.1015627 ,  0.27033682,  0.42009685]))
('mete', array([-1.80391939, -1.69474024, -1.36637911]))
('gasmete', array([-0.52040582, -0.37772299, -0.22860636]))

# THIRTY MYR BIN VARY (C_10):
('masse', array([  9.434018  ,   9.83425117,  10.21535215]))
('duste', array([ 0.11624586,  0.30704709,  0.45386506]))
('mete', array([-1.80094258, -1.70408946, -1.46885813]))
('gasmete', array([-0.41292822, -0.36352271, -0.30002578]))
('mass', array([ 10.3132814 ,  10.56044149,  10.85358103]))
('dust', array([ 0.39572945,  0.57010597,  0.75727346]))
('met', array([-1.81587191, -1.64381295, -1.13834254]))
('gasmet', array([-1.22951266, -0.56768566, -0.04120839]))

# NEWMASK:
('masse', array([  9.92217859,  10.17154694,  10.37101849]))
('mass', array([ 10.26151201,  10.4980092 ,  10.73596525]))
('duste', array([ 0.11645889,  0.24099925,  0.40087688]))
('dust', array([ 0.25051906,  0.4485528 ,  0.62015834]))
('mete', array([-1.77119631, -1.63574648, -1.1661375 ]))
('met', array([-1.81943177, -1.64793831, -1.40059553]))
('gasmete', array([-1.08805659, -0.53650928, -0.24590697]))
('gasmet', array([-1.17834177, -0.78646141, -0.26640241]))

# SHORT:
('masse', array([  9.79244015,  10.07043171,  10.32201672]))
('mass', array([ 10.16148235,  10.39132977,  10.62170816]))
('duste', array([ 0.11510793,  0.24030764,  0.37787486]))
('dust', array([ 0.31472595,  0.46629035,  0.63637274]))
('mete', array([-1.77362186, -1.45970088, -0.45076948]))
('met', array([-1.75056641, -1.37157118, -0.88023348]))
('gasmete', array([-1.05382175, -0.49581321, -0.14460007]))
('gasmet', array([-1.197699  , -0.74171874, -0.18949151]))

# VARY C_10, NEW SFG COMP
('masse', array([  9.56806931,   9.9428606 ,  10.26535469]))
('duste', array([ 0.1366495 ,  0.29065425,  0.4745783 ]))
('mete', array([-1.85909033, -1.72481412, -1.57049279]))
('gasmete', array([-0.39605677, -0.34745189, -0.29891892]))
('mass', array([ 10.32643452,  10.55351305,  10.83929956]))
('dust', array([ 0.33653897,  0.55233464,  0.74976853]))
('met', array([-1.82982044, -1.64881372, -1.2083779 ]))
('gasmet', array([-1.26446749, -0.59641001, -0.02850388]))
'''

'''
FICO errors on means:
('masse', array([ 9.37587193,  9.39738364,  9.42036556]))
('duste', array([ 0.30516629,  0.31939399,  0.33271809]))
('mete', array([-1.61562099, -1.55379575, -1.49063667]))
('gasmete', array([-0.37862442, -0.33522604, -0.29446069]))
('massl', array([ 10.10920348,  10.11635971,  10.12342247]))
('dustl', array([ 0.54547819,  0.5549387 ,  0.56425811]))
('metl', array([-1.56044648, -1.53799786, -1.51428261]))
('gasmetl', array([-0.72025215, -0.67690548, -0.63743182]))
'''

'''
GP log[o/h] = -3.3 (8.7-12)
sun log[o/h] = -3.09 (8.91-12)
--> 10^-3.3 / 10^-3.09 =0.616
say 8.35 = 0.25solar --> 8.35 - 12 = -3.65 & 10^-3.65 / 10^-x = 0.25 --> 10^-x = 4*10^-3.65 --> x = 3.05 - sol = 8.95
10^-3.3 / 10^-3.05
'''
