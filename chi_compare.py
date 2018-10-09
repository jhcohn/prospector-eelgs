import os
import make_all_plots as map
import pickle
import numpy as np
import print_sfh
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
import stellar_ages as sa

home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
e_only = 1  # True=1
# folders = ['pkl_efifty', 'pkl_nvary']
# file = ['fifty', 'vary']  # 'fix'  # 'fifty'  # 'vary'  # 'newmask'
# folders = ['pkl_efico', 'pkl_esol']
folders = ['pkl_efico', 'pkl_etenmet']
# folders = ['pkl_efifty', 'pkl_evar']  # 'pkl_efix']
# file = ['fico', 'sol']  # ['fifty', 'vary']  # 'fix'  # 'fifty'  # 'vary'  # 'newmask'
file = ['fico', 'tenmet']
ylabs = ['50 Myr bin', 'solar']  # '1/10 met minimum boundary']
# ylabs = ['50 Myr bin', '100 Myr bin']

objs = []
fields = []
count = 0

eelg_objs, eelg_fields, lbg_objs, lbg_fields = sa.get_gal_lists(base=file, objlists=True)
if e_only:
    eelg_objs2, eelg_fields2, lbg_objs, lbg_fields = sa.get_gal_lists(base=file, objlists=True)
    objsets = [[eelg_objs, eelg_fields], [eelg_objs2, eelg_fields2]]
    print(objsets)
else:
    objsets = [[eelg_objs, eelg_fields], [lbg_objs, lbg_fields]]
print(objsets[1])

set_count = 0
waves_cos = []
waves_cdf = []
waves_uds = []
chis_cdf = []
chis_cos = []
chis_uds = []
masks_cos = []
masks_cdf = []
masks_uds = []

lwaves_cos = []
lwaves_cdf = []
lwaves_uds = []
lchis_cdf = []
lchis_cos = []
lchis_uds = []
lmasks_cdf = []
lmasks_cos = []
lmasks_uds = []

for objset in objsets:  # for eelgs, then for lbgs
    objs = objset[0]  # eelg_objs, then lbg_objs
    fields = objset[1]  # eelg_fields, then lbg_fields
    folder = folders[set_count]  # 0 for EELG folder, 1 for LBG folder
    for i in range(len(objs)):  # for each galaxy in eelg_objs, then for each galaxy in lbg_objs
        obj = objs[i]
        field = fields[i]
        pre = folder + '/' + str(obj) + '_' + field + '_' + file[set_count] + '_'
        print(pre)
        base = '_out.pkl'
        extra = pre + 'extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
        res = pre + 'res' + base
        sed = pre + 'sed' + base
        restwave = pre + 'restwave' + base
        spec = pre + 'spec' + base
        spswave = pre + 'spswave' + base
        chisq = pre + 'chisq' + base
        justchi = pre + 'justchi' + base

        files = [extra, res, sed, restwave, spec, spswave, chisq, justchi]
        # print(files)
        # map.all_plots(files, obj, field, file, loc='upper left', sfh=False, curves=False, sep_uvj=False)

        if os.path.exists(files[1]) and os.path.exists(files[3]) and os.path.exists(files[6]):
            print('hi')
            # PLOT SEPARATE UVJ
            with open(files[1], 'rb') as res:
                results = pickle.load(res)

            with open(files[3], 'rb') as restwave:
                wave_rest = pickle.load(restwave)

            with open(files[6], 'rb') as chisq:
                chisq = pickle.load(chisq)

            # MASK
            # create mask
            mask = results['obs']['phot_mask']  # mask out -99.0 values
            wave_rest = np.asarray(wave_rest)

            # apply mask
            # phot = results['obs']['maggies'][mask]
            # unc = results['obs']['maggies_unc'][mask]
            for i in range(len(mask)):
                if not mask[i]:
                    chisq[i] = -99
            # wave_rest = wave_rest[mask]
            # chisq = chisq[mask]
            if set_count == 0:
                if field == 'cdfs':
                    waves_cdf.append(wave_rest)
                    chis_cdf.append(chisq)
                    masks_cdf.append(mask)
                elif field == 'cosmos':
                    waves_cos.append(wave_rest)
                    chis_cos.append(chisq)
                    masks_cos.append(mask)
                elif field == 'uds':
                    waves_uds.append(wave_rest)
                    chis_uds.append(chisq)
                    masks_uds.append(mask)
                else:
                    print('!!!!!broken')

            elif set_count == 1:
                if field == 'cdfs':
                    lwaves_cdf.append(wave_rest)
                    lchis_cdf.append(chisq)
                    lmasks_cdf.append(mask)
                elif field == 'cosmos':
                    lwaves_cos.append(wave_rest)
                    lchis_cos.append(chisq)
                    lmasks_cos.append(mask)
                elif field == 'uds':
                    lwaves_uds.append(wave_rest)
                    lchis_uds.append(chisq)
                    lmasks_uds.append(mask)
                else:
                    print('broken!!!!!!')
    set_count += 1

print(len(waves_cdf), len(waves_cos), len(waves_uds))
chis_cdf = np.asarray(chis_cdf)
chis_cos = np.asarray(chis_cos)
chis_cuds = np.asarray(chis_uds)
all_chis_cdf = np.zeros(shape=(len(waves_cdf[0])))
all_chis_cos = np.zeros(shape=(len(waves_cos[0])))
all_chis_uds = np.zeros(shape=(len(waves_uds[0])))
num_cdf = np.zeros(len(all_chis_cdf))
num_cos = np.zeros(len(all_chis_cos))
num_uds = np.zeros(len(all_chis_uds))
chisq_per_nphot = []
rands = []
for i in range(len(chis_cdf)):  # for each of the 19 galaxies
    newsum = 0
    photsum = 0
    for k in range(len(chis_cdf[i])):
        rands.append(np.random.choice(chis_cdf[i]))
    for j in range(len(chis_cdf[i])):  # for each data point
        if chis_cdf[i][j] >= 0:
            all_chis_cdf[j] += chis_cdf[i][j]
            num_cdf[j] += 1
        if 0 <= chis_cdf[i][j] < 300:
            newsum += chis_cdf[i][j]
            photsum += 1
    chisq_per_nphot.append(newsum / photsum)
for i in range(len(chis_cos)):  # for each of the 19 galaxies
    newsum = 0
    photsum = 0
    for l in range(len(chis_cos[i])):
        rands.append(np.random.choice(chis_cos[i]))
    for j in range(len(chis_cos[i])):  # for each data point
        if chis_cos[i][j] >= 0:
            all_chis_cos[j] += chis_cos[i][j]
            num_cos[j] += 1
        if 0 <= chis_cdf[i][j] < 300:
            newsum += chis_cdf[i][j]
            photsum += 1
    chisq_per_nphot.append(newsum / photsum)
for i in range(len(chis_uds)):  # for each of the 19 galaxies
    newsum = 0
    photsum = 0
    for m in range(len(chis_uds[i])):
        rands.append(np.random.choice(chis_uds[i]))
    for j in range(len(chis_uds[i])):  # for each data point
        if chis_uds[i][j] >= 0:
            all_chis_uds[j] += chis_uds[i][j]
            num_uds[j] += 1
        if 0 <= chis_cdf[i][j] < 300:
            newsum += chis_cdf[i][j]
            photsum += 1
    chisq_per_nphot.append(newsum / photsum)

for i in range(len(all_chis_cdf)):
    all_chis_cdf[i] /= num_cdf[i]
for i in range(len(all_chis_cos)):
    all_chis_cos[i] /= num_cos[i]
for i in range(len(all_chis_uds)):
    all_chis_uds[i] /= num_uds[i]

lchis_cdf = np.asarray(lchis_cdf)
lchis_cos = np.asarray(lchis_cos)
lchis_cuds = np.asarray(lchis_uds)
lall_chis_cdf = np.zeros(shape=(len(lwaves_cdf[0])))
lall_chis_cos = np.zeros(shape=(len(lwaves_cos[0])))
lall_chis_uds = np.zeros(shape=(len(lwaves_uds[0])))
lnum_cdf = np.zeros(len(all_chis_cdf))
lnum_cos = np.zeros(len(all_chis_cos))
lnum_uds = np.zeros(len(all_chis_uds))
lchisq_per_nphot = []
lrands = []
for i in range(len(lchis_cdf)):  # for each of the second set of galaxies
    lnewsum = 0
    lphotsum = 0
    for lk in range(len(lchis_cdf[i])):
        lrands.append(np.random.choice(lchis_cdf[i]))
    for j in range(len(lchis_cdf[i])):
        if lchis_cdf[i][j] >= 0:
            lall_chis_cdf[j] += lchis_cdf[i][j]
            lnum_cdf[j] += 1
        if 0 <= lchis_cdf[i][j] < 300:
            lnewsum += lchis_cdf[i][j]
            lphotsum += 1
    lchisq_per_nphot.append(lnewsum / lphotsum)
for i in range(len(lchis_cos)):  # for each of the second set of galaxies
    lnewsum = 0
    lphotsum = 0
    for ll in range(len(lchis_cos[i])):
        lrands.append(np.random.choice(lchis_cos[i]))
    for j in range(len(lchis_cos[i])):
        if lchis_cos[i][j] >= 0:
            lall_chis_cos[j] += lchis_cos[i][j]
            lnum_cos[j] += 1
        if 0 <= lchis_cdf[i][j] < 300:
            lnewsum += lchis_cdf[i][j]
            lphotsum += 1
    lchisq_per_nphot.append(lnewsum / lphotsum)
for i in range(len(lchis_uds)):  # for each of the second set of galaxies
    lnewsum = 0
    lphotsum = 0
    for lm in range(len(lchis_uds[i])):
        lrands.append(np.random.choice(lchis_uds[i]))
    for j in range(len(lchis_uds[i])):
        if lchis_uds[i][j] >= 0:
            lall_chis_uds[j] += lchis_uds[i][j]
            lnum_uds[j] += 1
        if 0 <= lchis_cdf[i][j] < 300:
            lnewsum += lchis_cdf[i][j]
            lphotsum += 1
    lchisq_per_nphot.append(lnewsum / lphotsum)

print(np.std(chisq_per_nphot))
print(np.std(lchisq_per_nphot))
ran = np.zeros(shape=(10**3, 10**3))
lran = np.zeros(shape=(10**3, 10**3))
means = []
lmeans = []
for j in range(10**3):
    for k in range(10**3):
        ran[j, k] = np.random.choice(chisq_per_nphot)
        lran[j, k] = np.random.choice(lchisq_per_nphot)
    means.append(np.mean(ran[j]))
    lmeans.append(np.mean(lran[j]))
print('errors on means')
print(np.percentile(means, [16., 50., 84.]))
print(np.percentile(lmeans, [16., 50., 84.]))

print(np.percentile(rands, [16., 50., 84.]))
print(np.percentile(lrands, [16., 50., 84.]))
print(np.percentile(chisq_per_nphot, [16., 50., 84.]), 'orig')
print(np.mean(chisq_per_nphot))
print(np.percentile(lchisq_per_nphot, [16., 50., 84.]), 'tenmet')
print(np.mean(lchisq_per_nphot))

delt = np.asarray(lchisq_per_nphot) - np.asarray(chisq_per_nphot)
delt_percs = np.percentile(delt, [16., 50., 84.])
siz = 10**3  # 10**4
delt_ran = np.zeros(shape=(siz, siz))
delt_means = []
for j in range(siz):
    for k in range(siz):
        delt_ran[j, k] = np.random.choice(delt)
    delt_means.append(np.mean(delt_ran[j]))
print(delt)
print(delt_percs[1], "-", delt_percs[1] - delt_percs[0], "+", delt_percs[2] - delt_percs[1])
delt_emeans = np.percentile(delt_means, [16., 50., 84.])
print('errors on mean', delt_emeans[1], '-', delt_emeans[1] - delt_emeans[0], '+', delt_emeans[2] - delt_emeans[1])

for i in range(len(lall_chis_cdf)):
    lall_chis_cdf[i] /= lnum_cdf[i]
for i in range(len(lall_chis_cos)):
    lall_chis_cos[i] /= lnum_cos[i]
for i in range(len(lall_chis_uds)):
    lall_chis_uds[i] /= lnum_uds[i]

plt.scatter(waves_cdf[0], all_chis_cdf, s=30, color='b', marker='D', label=r'varying metallicity (CDFS)')
plt.scatter(waves_cos[0], all_chis_cos, s=30, color='b', marker='o', label=r'varying metallicity (COSMOS)')
plt.scatter(waves_uds[0], all_chis_uds, s=30, color='b', marker='*', label=r'varying metallicity (UDS)')

plt.scatter(lwaves_cdf[0], lall_chis_cdf, s=30, color='r', marker='D', label=r'1/10 metallicity (CDFS)')
plt.scatter(lwaves_cos[0], lall_chis_cos, s=30, color='r', marker='o', label=r'1/10 metallicity (COSMOS)')
plt.scatter(lwaves_uds[0], lall_chis_uds, s=30, color='r', marker='*', label=r'1/10 metallicity (UDS)')

if e_only:
    plt.ylim(0, 40)
    plt.xlim(700, 20000)
else:
    plt.ylim(0, 10**3)
    plt.xlim(0, 25000)

plt.ylabel(r'Mean $\chi^2$', fontsize=30)
plt.xlabel(r'Wavelength (Rest) [$\rm \AA$]', fontsize=30)
plt.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
plt.xscale('log')
plt.legend(numpoints=1, loc='upper left', prop={'size': 20})  # , line2) ... , r'$\chi$']
plt.tick_params('x', length=3, width=1, which='both', labelsize=20)
plt.tick_params('y', length=3, width=0.5, which='both', labelsize=20)
plt.show()

cdf = []
cos = []
uds = []
for i in range(len(all_chis_cdf)):
    print(all_chis_cdf[i], lall_chis_cdf[i], 'me!')
    if np.isnan(all_chis_cdf[i]) or np.isnan(lall_chis_cdf[i]):
        cdf.append(-5)
    else:
        cdf.append(all_chis_cdf[i]/lall_chis_cdf[i])
for i in range(len(all_chis_cos)):
    if np.isnan(all_chis_cos[i]) or np.isnan(lall_chis_cos[i]):
        cos.append(-5)
    else:
        cos.append(all_chis_cos[i]/lall_chis_cos[i])
for i in range(len(all_chis_uds)):
    if np.isnan(all_chis_uds[i]) or np.isnan(lall_chis_uds[i]):
        uds.append(-5)
    else:
        uds.append(all_chis_uds[i]/lall_chis_uds[i])
plt.scatter(waves_cdf[0], cdf, s=40, color='purple', marker='D', label=r'CDFS')
plt.scatter(waves_cos[0], cos, s=40, color='purple', marker='o', label=r'COSMOS')
plt.scatter(waves_uds[0], uds, s=40, color='purple', marker='*', label=r'UDS')

cdf1 = [x for x in cdf if x >= 0]
cos1 = [x for x in cos if x >= 0]
uds1 = [x for x in uds if x >= 0]
total = cdf1 + cos1 + uds1
print(len(total), len(cos1))
print(np.percentile(cdf1, [16., 50., 84.]))
print(np.percentile(cos1, [16., 50., 84.]))
print(np.percentile(uds1, [16., 50., 84.]))
print(np.percentile(total, [16., 50., 84.]))  # ~0.911 for 50_myr_bin / 100_
# myr_bin

plt.ylim(10**-2, 10**2)
plt.xlim(700, 21000)
plt.ylabel(ylabs[0] + r' $\chi^2$ / ' + ylabs[1] + r' $\chi^2$', fontsize=30)  # r'Free $\chi^2$ / Fixed $\chi^2$'
plt.xlabel(r'Wavelength (Rest) [$\rm \AA$]', fontsize=30)
plt.axvspan(4800, 5050, color='k', alpha=0.175)  # 0.2
plt.axhline(y=1, color='k', linestyle='--')
plt.xscale('log')
plt.yscale('log')
plt.legend(numpoints=1, loc='upper right', prop={'size': 20})  # , line2) ... , r'$\chi$']
plt.tick_params('x', length=3, width=1, which='both', labelsize=20)
plt.tick_params('y', length=3, width=0.5, which='both', labelsize=20)
plt.show()

# BINS!
#for i in range(len(waves_cdf[0])):




'''
percs = []
for i in range(len(chis)):
    percs.append(np.percentile(chis[i], [16., 50., 84.]))
print(percs)
print(len(percs))

medperc = np.percentile(percs[:1], [16., 50., 84.])
print(medperc)
'''

'''
# efix
[array([ 0.1393269 ,  0.66732137,  2.82596265]), array([ 0.08514731,  0.49976773,  2.35719781]),
array([ 0.0070303 ,  0.34480457,  3.12724055]), array([ 0.09096948,  1.14187229,  3.12762379]),
array([ 0.1654041 ,  0.81367099,  3.86377822]), array([ 0.03623754,  0.79877345,  4.48364795]),
array([  0.22521677,   4.19483908,  10.35979175]), array([ 0.0544622 ,  0.61057814,  3.50712903]),
array([ 0.04231078,  1.87354846,  5.89014532]), array([  0.21247055,   1.57546266,  11.42379847]),
array([ 0.14261526,  1.13629065,  4.97304168]), array([ 0.01445016,  0.49196858,  4.15568844]),
array([ 0.10912794,  1.13757351,  5.20695445]), array([  0.14121575,   1.94560239,  10.77372169]),
array([  0.11575091,   1.2801173 ,  14.23217029]), array([ 0.03329902,  0.50140011,  9.15408339]),
array([ 0.0234051 ,  0.93095627,  3.44954613]), array([ 0.03573371,  0.58266962,  3.88255853]),
array([ 0.10879415,  1.15576025,  4.16048265])]
19
[ 0.30828513  0.66732137  2.13519744]

# efifty
[array([ 0.04202765,  0.45736303,  1.8723331 ]), array([ 0.00762068,  0.21519007,  1.36598681]),
array([ 0.00687951,  0.23571147,  2.26337515]), array([ 0.0659939 ,  0.84250053,  2.89549189]),
array([ 0.05172658,  0.36957506,  1.6232407 ]), array([ 0.0392396 ,  0.40155736,  2.32734246]),
array([ 0.50113971,  4.1163247 ,  9.52775831]), array([ 0.05899883,  0.54884343,  2.17061003]),
array([ 0.07926935,  0.7476057 ,  3.98363519]), array([ 0.06478105,  1.04571208,  2.739091  ]),
array([ 0.10047842,  0.89106899,  2.87402035]), array([ 0.0432529 ,  0.43763784,  2.14528272]),
array([ 0.03678089,  0.64314915,  4.30212667]), array([ 0.1815536 ,  1.59492966,  6.27511173]),
array([ 0.11764305,  0.74189518,  7.27722902]), array([ 0.06587803,  0.77792715,  5.34056368]),
array([ 0.04309122,  0.42629335,  3.0754503 ]), array([ 0.00761615,  0.43407241,  3.65139897]),
array([ 0.11534648,  0.53949223,  2.51586129])]
19
[ 0.17493497  0.45736303  1.41954268]
'''
# fig = plt.figure()
# fig.text(0.5, 0.03, r'Wavelength (Rest) [$\rm \AA$]', ha='center', va='bottom', fontsize=fs_text)  # 30
# plt.show()

