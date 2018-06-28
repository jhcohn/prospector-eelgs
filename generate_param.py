import numpy as np
import sys
from prospect.models import model_setup
from prospect.io import write_results
from prospect import fitting
from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

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
copyfrom = 'eelg_sfhtest_params.py'

for m in range(len(masses)):
    print(m, ' /5!')
    for d in range(len(dusts)):
        print(d, ' /3!')
        for z in range(len(mets)):
            for f in range(len(sfr1)):
                newpar = 'testpars/sfhtest_' + str(int(masses[m]*10)) + '_' + str(int(dusts[d]*100)) + '_' +\
                         str(int(mets[z]*-10)) + '_' + str(int(sfr1[f]*100)) + '_params.py'
                # e.g. newpar = testpars/sfhtest_90_15_20_5_params.py
                orig = open(copyfrom, 'r')
                writenew = open(newpar, 'w+')
                for line in orig.readlines():  # for each line in original param file
                    towrite = line  # copy the line over, unless I need to change it:
                    if line.startswith("    model_params[n.index('logmass')]['init']"):
                        towrite = "    model_params[n.index('logmass')]['init'] = " + str(masses[m]) + '\n'
                    elif line.startswith("    model_params[n.index('dust2')]['init']"):
                        towrite = "    model_params[n.index('dust2')]['init'] = " + str(dusts[d]) + '\n'
                    elif line.startswith("    model_params[n.index('logzsol')]['init']"):
                        towrite = "    model_params[n.index('logzsol')]['init'] = " + str(mets[z]) + '\n'
                    elif line.startswith("    eelg_frac1 = "):
                        towrite = "    eelg_frac1 = " + str(sfr1[f]) + '\n'
                    writenew.write(towrite)  # + '\n')
                orig.close()
                writenew.close()
print('done!')