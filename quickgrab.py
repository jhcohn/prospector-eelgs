import prospect.io.read_results as bread
import matplotlib.pyplot as plt
import numpy as np

# res, pr, mod = bread.results_from("cosmos1824_test1_1484768675_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_test2_1484784141_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_phot_3.0_x2niter_test_1485104536_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_test_with_err_3.0_and_x2niter_1485202527_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_err_3.0_mass8to11_1485212425_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_err_3.0_mass8to11_maxtage10_minZ-2_1485276700_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_err_3.0_mass7to11_minZ-2_1485285100_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_test25_1485314083_mcmc")
# res, pr, mod = bread.results_from("demo_mpirun_4_1485372547_mcmc")
# res, pr, mod = bread.results_from("demo_mpirun_2_hdf5_commented_1485798339_mcmc")
# res, pr, mod = bread.results_from("demo_mpirun_4_hdf5_commented_1485803218_mcmc")
# demo_mpirun_4_hdf5_commented_1485803223_mcmc
# demo_mpirun_4_hdf5_commented_1485803224_mcmc
# res, pr, mod = bread.results_from("demo_mpirun_4_x2niter_hdf5_commented_1485809892_mcmc")
# res, pr, mod = bread.results_from("demo_mpirun_4_23dayslater_1486078742_mcmc")
# demo_mpirun_4_23dayslater_1486078743_mcmc
# res, pr, mod = bread.results_from("is_this_thing_on_1486398386_mcmc")
# res, pr, mod = bread.results_from("cosmos1824_23dayslater_fulltest_1486402078_mcmc")
# cosmos1824_23dayslater_fulltest_1486402081_mcm
# cosmos1824_23dayslater_fulltest_1486402082_mcmc
# cosmos1824_23dayslater_fulltest_1486402084_mcmc
# res, pr, mod = bread.results_from("cosmos1824_23dayslater_fulltest_1486405531_mcmc")  # niter=2048
# res, pr, mod = bread.results_from("cosmos1824_23_niter4096_1486410548_mcmc")  # 546, 547, 548
# res, pr, mod = bread.results_from("cosmos1824_23_lumdist_niter4096_1486421418_mcmc")  # 418, 419
# res, pr, mod = bread.results_from("cosmos1824_23_lumdist_niter4096_newranges_1486427756_mcmc")  # 756, 757, 759
# res, pr, mod = bread.results_from("cosmos1824_23_lumdist_niter4096_newzred_1486485286_mcmc")  # 286, 287, 291, 294
# res, pr, mod = bread.results_from("cosmos1824_23_lumdist_niter10000_newzred_1486492983_mcmc")  #
# res, pr, mod = bread.results_from("1824_niter1000_betterparamfile_1486513658_mcmc")  # 658, 663, 664, 666
# res, pr, mod = bread.results_from("cosmos1824_fullfilterset_1487964018_mcmc")  # 4011, 4017, 4018 FIRST OF THE NEW STUFF
# res, pr, mod = bread.results_from("cosmos1824_fullfilterset_niter1000_1487972522_mcmc")  # 521, 522, 523 FIRST OF THE NEW STUFF
# res, pr, mod = bread.results_from("cosmos_1824_fullfilterset_niter2000_1488215032_mcmc")  # 032, 033 Second OF THE NEW STUFF
# res, pr, mod = bread.results_from("cosmos_1824_themiyafilts_1488771358_mcmc")  # 358, 359 Trying with Themiya's filters
# res, pr, mod = bread.results_from("cosmos_1824_adapted_1488998012_mcmc")  # with nebular emission
# res, pr, mod = bread.results_from("cosmos_1824_fastparamrun_1489033721_mcmc")  # fast_mimic
# res, pr, mod = bread.results_from("cosmos_1824_fastparamrun_newphot_1489259293_mcmc")  # fast_mimic now with good photometry hopefully?
# res, pr, mod = bread.results_from("cosmos_1824_fast_mimic_3_1489530383_mcmc")  # fast_mimic now with almost good filters
# res, pr, mod = bread.results_from("cosmos_1824_fast_mimic_fixed_filts_1489786470_mcmc")  # 470, 473, 474, 475 fast_mimic now with good filters
# res, pr, mod = bread.results_from("cosmos_1824_fastparamrun_newphot_1489259293_mcmc")  # 93, 96, 97, 98 fast_mimic now with good filters (not as good though)
# res, pr, mod = bread.results_from("cosmos_1824_fast_mimic_new_NB_1490463705_mcmc")  # 05, 08, 10 fast_mimic now with better?
# res, pr, mod = bread.results_from("cosmos_2329_fastmimic_normalgal_1491261229_mcmc")  # 29, 30, 32, 33 normal gal 2329
# res, pr, mod = bread.results_from("cosmos_2329_fast_mimic_fixed_1491441386_mcmc")  # 29, 30, 32, 33 normal gal 2329
# res, pr, mod = bread.results_from("joel_fast_mimic_massive_gals_1491671811_mcmc")  # 11, 12, 13, 15 joel galaxy
# res, pr, mod = bread.results_from("2329_fast_mimic_try_1492095912_mcmc")  # 2329 edited flux
# res, pr, mod = bread.results_from("2329_fast_mimic_logtau_1492099753_mcmc")  # 53, 56 2329 edited flux (doesn't work for 1824)
# res, pr, mod = bread.results_from("1824_fast_mimic_logtau_1492108987_mcmc")  # 87, 88 1824 logtau
# res, pr, mod = bread.results_from("2329_fast_mimic_logtau_1492099753_mcmc")  # 0997: 53, 56; 1833: 53, 56; 2329 logtau
# res, pr, mod = bread.results_from("1824_fast_mimic_logtau_forreal_1492209450_mcmc")  # 50, 52 1824 logtau, not tau
res, pr, mod = bread.results_from("2329_logtau_fast_1492441685_mcmc")  # 85, 86 2329 logtau, not tau

id = 1824  # 2329

# res, pr, mod = bread.results_from("cosmos_1824_normed_fullfilts_1488418724_mcmc")  # NORMED (new stuff) 724, 725

''' #
# USE THIS TO FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
tracefig, prob = bread.param_evol(res)
print(prob)
print('max', prob.max())
for i in range(len(prob)):
    for j in range(len(prob[i])):
        if prob[i][j] == prob.max():
            print('max', prob[i][j], i, j)
        if prob[i][j] == prob.min():
            print('min', prob[i][j], i, j)  # i, j = walker, iteration
            # cosmos_1824_fast_mimic_fixed_filts_1489786470_mcmc
            # ('min', 597.91032354824006, 372, 399)
            # ('min', 597.91032354824006, 372, 400)
            # ('max', 613.05532603869699, 569, 495)
            # ('max', 613.05532603869699, 569, 496)
            # ('max', 613.05532603869699, 569, 497)
            # ('max', 613.05532603869699, 569, 498)
            # ('max', 613.05532603869699, 569, 499)
            # ('max', 613.05532603869699, 569, 500)
            # ('max', 613.05532603869699, 569, 501)
            # ('max', 613.05532603869699, 569, 502)
            # ('max', 613.05532603869699, 569, 503)
            # ('max', 613.05532603869699, 569, 504)
            # ('max', 613.05532603869699, 569, 505)
            # ('max', 613.05532603869699, 569, 506)
            # ('max', 613.05532603869699, 569, 507)
            # ('max', 613.05532603869699, 569, 508)
            # ('max', 613.05532603869699, 569, 509)
            # ('max', 613.05532603869699, 569, 510)
            # ('max', 613.05532603869699, 569, 511)
print(len(prob))

'''
# 2329_logtau_fast_1492441685_mcmc
# ('max', 753.44505433937184)
# ('min', 742.04121396763639, 73, 173)
# ('max', 753.44505433937184, 96, 493)



# 1824_fast_mimic_logtau_forreal_1492209450_mcmc
# ('max', 567.96651332519787)
# ('min', 556.26781947574068, 47, 342)
# ('max', 567.96651332519787, 110, 168)


# 2329_fast_mimic_logtau_1492099753_mcmc
#('max', 787.71386581812976)
#('min', 770.78390703984974, 74, 339)
#('min', 770.78390703984974, 74, 340)
#('min', 770.78390703984974, 74, 341)
#('min', 770.78390703984974, 74, 342)
#('min', 770.78390703984974, 74, 343)
#('max', 787.71386581812976, 93, 451)
#('max', 787.71386581812976, 93, 452)
#('max', 787.71386581812976, 93, 453)
#('max', 787.71386581812976, 93, 454)
#('max', 787.71386581812976, 93, 455)
#('max', 787.71386581812976, 93, 456)
#('max', 787.71386581812976, 93, 457)
#('max', 787.71386581812976, 93, 458)
#('max', 787.71386581812976, 93, 459)
#('max', 787.71386581812976, 93, 460)
#('max', 787.71386581812976, 93, 461)
#('max', 787.71386581812976, 93, 462)
#('max', 787.71386581812976, 93, 463)
#('max', 787.71386581812976, 93, 464)
#('max', 787.71386581812976, 93, 465)
#('max', 787.71386581812976, 93, 466)
#('max', 787.71386581812976, 93, 467)
#('max', 787.71386581812976, 93, 468)
#('max', 787.71386581812976, 93, 469)
#('max', 787.71386581812976, 93, 470)
#('max', 787.71386581812976, 93, 471)
#('max', 787.71386581812976, 93, 472)
#('max', 787.71386581812976, 93, 473)
#('max', 787.71386581812976, 93, 474)
#('max', 787.71386581812976, 93, 475)
#('max', 787.71386581812976, 93, 476)
#('max', 787.71386581812976, 93, 477)
#('max', 787.71386581812976, 93, 478)
#('max', 787.71386581812976, 93, 479)
#('max', 787.71386581812976, 93, 480)
#('max', 787.71386581812976, 93, 481)
#('max', 787.71386581812976, 93, 482)
#('max', 787.71386581812976, 93, 483)
#('max', 787.71386581812976, 93, 484)
#('max', 787.71386581812976, 93, 485)
#('max', 787.71386581812976, 93, 486)
#('max', 787.71386581812976, 93, 487)
#('max', 787.71386581812976, 93, 488)
#('max', 787.71386581812976, 93, 489)
#('max', 787.71386581812976, 93, 490)
#('max', 787.71386581812976, 93, 491)
#('max', 787.71386581812976, 93, 492)
#('max', 787.71386581812976, 93, 493)
#('max', 787.71386581812976, 93, 494)
#('max', 787.71386581812976, 93, 495)
#('max', 787.71386581812976, 93, 496)
#('max', 787.71386581812976, 93, 497)
#('max', 787.71386581812976, 93, 498)
#('max', 787.71386581812976, 93, 499)

# 1824_fast_mimic_logtau_1492108987_mcmc
#('max', 567.03876383788281)
#('min', 555.98206153986007, 54, 314)
#('min', 555.98206153986007, 54, 315)
#('min', 555.98206153986007, 54, 316)
#('max', 567.03876383788281, 95, 486)
#('max', 567.03876383788281, 95, 487)
#('max', 567.03876383788281, 95, 488)
#('max', 567.03876383788281, 95, 489)
#('max', 567.03876383788281, 95, 490)
#('max', 567.03876383788281, 95, 491)
#('max', 567.03876383788281, 95, 492)
#('max', 567.03876383788281, 95, 493)
#('max', 567.03876383788281, 95, 494)
#('max', 567.03876383788281, 95, 495)
#('max', 567.03876383788281, 95, 496)
#('max', 567.03876383788281, 95, 497)


# 2329_fast_mimic_try_1492095912_mcmc
#('max', 788.9067752985834)
#('max', 788.9067752985834, 11, 487)
#('max', 788.9067752985834, 11, 488)
#('max', 788.9067752985834, 11, 489)
#('max', 788.9067752985834, 11, 490)
#('max', 788.9067752985834, 11, 491)
#('max', 788.9067752985834, 11, 492)
#('max', 788.9067752985834, 11, 493)
#('max', 788.9067752985834, 11, 494)
#('max', 788.9067752985834, 11, 495)
#('max', 788.9067752985834, 11, 496)
#('max', 788.9067752985834, 11, 497)
#('max', 788.9067752985834, 11, 498)
#('max', 788.9067752985834, 11, 499)
#('min', 746.13350683518991, 130, 0)

# joel_fast_mimic_massive_gals_1491671811_mcmc
#('max', 939.65419929356631)
#('min', 925.88206218900143, 31, 463)
#('min', 925.88206218900143, 31, 464)
#('min', 925.88206218900143, 31, 465)
#('max', 939.65419929356631, 117, 378)
#('max', 939.65419929356631, 117, 379)
#('max', 939.65419929356631, 117, 380)
#('max', 939.65419929356631, 117, 381)



# cosmos_2329_fast_mimic_fixed_1491441386_mcmc
#('max', 750.13134681487418)
#('min', 736.40733909546373, 41, 304)
#('min', 736.40733909546373, 41, 305)
#('max', 750.13134681487418, 50, 397)
#('max', 750.13134681487418, 50, 398)
#('max', 750.13134681487418, 50, 399)
#('max', 750.13134681487418, 50, 400)
#('max', 750.13134681487418, 50, 401)
#('max', 750.13134681487418, 50, 402)
#('max', 750.13134681487418, 50, 403)
#('max', 750.13134681487418, 50, 404)
#('max', 750.13134681487418, 50, 405)
#('max', 750.13134681487418, 50, 406)
#('max', 750.13134681487418, 50, 407)
#('max', 750.13134681487418, 50, 408)
#('max', 750.13134681487418, 50, 409)
#('max', 750.13134681487418, 50, 410)
#('max', 750.13134681487418, 50, 411)


# cosmos_1824_fast_mimic_new_NB_1490463705_mcmc
# ('max', 620.45157572586163, 72, 114)
# ('min', 605.83272113396561, 200, 217)

# 2329! cosmos_2329_fastmimic_normalgal_1491261229_mcmc
# ('max', 782.37213178936054, 182, 359)
# ('max', 782.37213178936054, 182, 360)
# ('max', 782.37213178936054, 182, 361)
# ('max', 782.37213178936054, 182, 362)
# ('max', 782.37213178936054, 182, 363)
# ('max', 782.37213178936054, 182, 364)
# ('max', 782.37213178936054, 182, 365)
# ('max', 782.37213178936054, 182, 366)
# ('max', 782.37213178936054, 182, 367)
# ('max', 782.37213178936054, 182, 368)
# ('max', 782.37213178936054, 182, 369)
# ('max', 782.37213178936054, 182, 370)
# ('max', 782.37213178936054, 182, 371)
# ('max', 782.37213178936054, 182, 372)
# ('min', 766.48433321433322, 341, 417)
# ('min', 766.48433321433322, 341, 418)
# ('min', 766.48433321433322, 341, 419)


# '''
# ''' #
# PRINTS TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
# THEN PRINTS CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
tracefig = bread.param_evol(res)
plt.show()
# cornerfig = bread.subtriangle(res, start=0, thin=5)
if id == 1824:
    # cornerfig = bread.subtriangle(res, start=400, thin=5, truths=[9.78, 0.25, -1., 0.00, 0.1])  # 1824 (tau=0.1, logtau=-1.)
    cornerfig = bread.subtriangle(res, start=400, thin=5, truths=[9.78, 0.25, -1., 0.00], show_titles=True)  # 1824 logtau not tau
elif id == 2329:
    # cornerfig = bread.subtriangle(res, start=400, thin=5, truths=[11.23, 0.25, -1.2, 1.70 / 1.086, 10 ** 7.8 / 1e9])  # 2329 logtau
    cornerfig = bread.subtriangle(res, start=400, thin=5, truths=[11.23, 0.25, -1.2, 1.70 / 1.086], show_titles=True)  # 2329 logtau not tau

# cornerfig = bread.subtriangle(res, start=400, thin=5, truths=[11.61, 10 ** 9.4 / 1e9, 10 ** 8.60 / 1e9, 0.9 / 1.086])  # joel
# truths = mass, age, tau, dust2 (recall A_v / 1.086 = dust2)
# logtau --> truths = mass, age, logtau, A_v, tau
plt.show()

'''
# mpirun -n 4 python prospector.py --param_file=1824_fast_params.py --objid=0 --outfile=cosmos_1824_fast_mimic_fixed_filts --niter=512
'''


# ''' #
# PRINTS MODEL SED FOR OBJECT
# We need the correct sps object to generate models
from prospect.sources import CSPBasis
sps = CSPBasis(**res['run_params'])

# Choose the walker and iteration number
# walker, iteration = 72, 114  # 372, 399  # 49, -3  # 0, -1
# walker, iteration = 41, 304  # 182, 370  # galaxy cosmos 2329
# walker, iteration = 117, 380  # joel galaxy
# walker, iteration = 11, 487
# walker, iteration = 95, 487  # 1824 logtau
# walker, iteration = 93, 455  # 2329 logtau
# walker, iteration = 110, 168  # 1824 logtau
walker, iteration = 96, 493  # 2329 logtau only
# Get the modeled spectra and photometry.
# These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
spec, phot, mfrac = mod.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)
# Plot the model SED
import matplotlib.pyplot as pl
wave = [f.wave_effective for f in res['obs']['filters']]

wave = np.asarray(wave)

mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
print(mean)


print(wave)
print(type(wave))
print(res['obs']['maggies'])
print(type(res['obs']['maggies']))

y = res['obs']['maggies']
yerr = res['obs']['maggies_unc']

pl.subplot(111, xscale="log", yscale="log")
pl.errorbar(wave, y, yerr=yerr, marker='o', linestyle='', color='b')
pl.plot(wave, phot, 'o', label='Model at {},{}'.format(walker, iteration), color='r')
pl.show()

pl.loglog(wave, res['obs']['maggies'], 'o', label='Observations')  # '-o'
pl.loglog(wave, phot, 'o', label='Model at {},{}'.format(walker, iteration))  # '-o'
pl.loglog(wave, res['obs']['maggies_unc'], 'o', label='Uncertainties')
pl.ylabel("Maggies")
pl.show()


chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
plt.plot(wave, chi_sq, 'o', color='b')
# plt.plot(np.log10(global_obs['wave_effective']), chi_sq, 'o', color='b')
plt.xlabel('Wave Effective')
plt.ylabel(r'$\chi^2$')
plt.show()
# '''


'''
WARNING:root:Deprecated keyword argument 'extents'. Use 'range' instead.
Traceback (most recent call last):
  File "quickgrab.py", line 13, in <module>
    cornerfig = bread.subtriangle(res, start=0, thin=5)
  File "/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospect/io/read_results.py", line 325, in subtriangle
    quantiles=[0.16, 0.5, 0.84], range=trim_outliers, **kwargs)
  File "/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/corner/corner.py", line 173, in corner
    "{0}".format, np.arange(len(m))[m]))))
ValueError: It looks like the parameter(s) in column(s) 0 have no dynamic range. Please provide a `range` argument.


list1 = [0,1,2,3]
print(list1[0:2])
'''

'''
PROSPECTOR PROJECT
run prospector on z~3 EELGs and make a letter
'''