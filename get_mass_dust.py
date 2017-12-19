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


def printer(out_file):
    print(out_file)

    res, pr, mod = bread.results_from(out_file)
    # ''' #

    # PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
    return bread.md(res, start=650, thin=5)  # set start by when kl converges!


if __name__ == "__main__":

    vary = False
    others = False
    short = True
    if vary:
        folders = ['out_evar/', 'out_nvar/']
        pars = ['eelg_varymet_params.py', 'eelg_varymet_params.py']
    elif others:
        folders = ['out_eth/', 'out_nth/']
        pars = ['eelg_thirty_params.py', 'noelg_thirty_params.py']
    elif short:
        folders = ['out_eshort/', 'out_nshort/']
        pars = ['eelg_short_params.py', 'eelg_short_params.py']
    else:
        folders = ['out_fixedmet/', 'out_noelg/']
        pars = ['eelg_fixedmet_params.py', 'noelg_multirun_params.py']

    eelgs = []
    oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[0]
    for file in os.listdir(oute):
        if file.endswith(".h5"):
            eelgs.append(file)

    lbgs = []
    outl = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folders[1]
    for file in os.listdir(outl):
        if file.endswith(".h5"):
            lbgs.append(file)

    get_e = np.zeros(shape=(4, len(eelgs)))  # 4 rows (dust, mass, gaslogz, logzsol), each row as long as eelgs
    for i in range(len(eelgs)):
        get_e[:, i] = printer(oute + eelgs[i])

    get_l = np.zeros(shape=(4, len(lbgs)))  # 4 rows (dust, mass, gaslogz, logzsol), each row as long as lbgs
    for i in range(len(lbgs)):
        get_l[:, i] = printer(outl + lbgs[i])

    print('masse', np.percentile(get_e[0], [16., 50., 84.]))
    print('mass', np.percentile(get_l[0], [16., 50., 84.]))
    print('duste', np.percentile(get_e[1], [16., 50., 84.]))
    print('dust', np.percentile(get_l[1], [16., 50., 84.]))
    print('mete', np.percentile(get_e[2], [16., 50., 84.]))
    print('met', np.percentile(get_l[2], [16., 50., 84.]))
    print('gasmete', np.percentile(get_e[3], [16., 50., 84.]))
    print('gasmet', np.percentile(get_l[3], [16., 50., 84.]))

'''
RUNNING WITH:
python get_mass_dust.py

For dust2 values: convert to A_V by multiplying by 1.86

# NOT VARY:
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
'''

'''
GP log[o/h] = -3.3 (8.7-12)
sun log[o/h] = -3.09 (8.91-12)
--> 10^-3.3 / 10^-3.09 =0.616
say 8.35 = 0.25solar --> 8.35 - 12 = -3.65 & 10^-3.65 / 10^-x = 0.25 --> 10^-x = 4*10^-3.65 --> x = 3.05 - sol = 8.95
10^-3.3 / 10^-3.05
'''
