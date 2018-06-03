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
np.errstate(invalid='ignore')


def printer(out_file, cvg=-1000, masstest=False, obj_true=None):
    objname = ''
    count = 0
    field = ''
    for i in kwargs['outname']:
        if i == '_':
            count += 1
        elif count == 0:
            objname += i
        elif count == 1:
            field += i
        elif count == 2:
            break
    print(field)

    res, pr, mod = bread.results_from(out_file)
    # print(res['chain'], 'chain')

    # ''' #

    # PRINT TRACE SHOWING HOW ITERATIONS CONVERGE FOR EACH PARAMETER
    bread.param_evol(res)  # print tracefig
    plt.show()
    # plt.savefig('delete' + '_tracefig2.png', bbox_inches='tight')

    # FIND WALKER, ITERATION THAT GIVE MAX (AND MIN) PROBABILITY
    # print(prob)
    prob = res['lnprobability'][..., cvg:]
    print('max', prob.max())
    row = prob.argmax() / len(prob[0])
    col = prob.argmax() - row * len(prob[0])
    walker, iteration = row, col
    print(walker, iteration)

    # ''' #

    truths = None
    print(objname)
    if masstest:
        print('hi?')
        import get_mass_dust as gmd
        oute = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/out_efico/'
        for file in os.listdir(oute):
            if file.endswith(".h5") and file.startswith(obj_true):
                print(file, 'file')
                params = gmd.printer(oute + file, percs=False, masstest=True)
                mass = params[0]
                dust = params[6]
                met = params[7]
                gasmet = params[8]
                sfh = []
                for i in [1, 2, 3, 4, 5]:
                    sfh.append(params[i])
                truths = [np.log10(10**mass + 6.5*10**8), sfh[0], sfh[1], sfh[2], sfh[3], sfh[4], dust, met, gasmet]
    # PRINT CORNERFIG CONTOURS/HISTOGRAMS FOR EACH PARAMETER
    # truths now corrected for masses of efico output + 6.5e8
    # truths = [9.466 + 6.5e8, 0.22, 0.12, 0.22, 0.12, 0.1, 0.0737, -1.79, -0.33]  # 12552
    # truths = [9.83 + 6.5e8, 0.43, 0.08, 0.04, 0.12, 0.1, 0.32, -1.57, -0.35]  # 12105
    # truths = [9.51 + 6.5e8, 0.80, 0.02, 0.01, 0.04, 0.04, 0.54, -1.94, -0.31]  # 21442
    print(truths)
    bread.subtriangle(res, start=cvg, thin=5, truths=truths, show_titles=True)  # set start by when kl converges!
    plt.show()
    # plt.savefig('delete' + '_ctracefig2.png', bbox_inches='tight')
    # For FAST: truths = [mass, age, tau, dust2] (for 1824: [9.78, 0.25, -1., 0.00])
    return objname, field, res, mod, walker, iteration
# ''' #


def sed(objname, field, res, mod, walker, iteration, param_file, **kwargs):
    # PRINT MODEL SED FOR GALAXY
    # We need the correct sps object to generate models
    sargv = sys.argv
    argdict = {'param_file': param_file}
    clargs = model_setup.parse_args(sargv, argdict=argdict)
    run_params = model_setup.get_run_params(argv=sargv, **clargs)
    sps = model_setup.load_sps(**run_params)

    # GET MODELED SPECTRA AND PHOTOMETRY
    # These have the same shape as the obs['spectrum'] and obs['maggies'] arrays.
    spec, phot, mfrac = mod.mean_model(res['chain'][walker, iteration, :], obs=res['obs'], sps=sps)

    mean = ((res['obs']['maggies']-phot)/res['obs']['maggies']).mean()
    print(mean, 'mean')  # print normalized mean difference between model and observations

    # PLOT SPECTRUM
    wave = [f.wave_effective for f in res['obs']['filters']]
    wave = np.asarray(wave)
    print('len', len(sps.wavelengths), len(spec))

    ''' #
    plt.plot(sps.wavelengths, spec)
    plt.xlabel('Wavelength [angstroms]')
    plt.title(str(objname) + ' spec')
    plt.show()

    # ''' #
    ''' #
    # HOW CONVERGED IS THE CODE?? LET'S FIND OUT!
    parnames = np.array(res['model'].theta_labels())
    fig, kl_ax = plt.subplots(1, 1, figsize=(7, 7))
    for i in xrange(parnames.shape[0]):
        kl_ax.plot(res['kl_iteration'], np.log10(res['kl_divergence'][:, i]),
                   'o', label=parnames[i], lw=1.5, linestyle='-', alpha=0.6)

    kl_ax.set_ylabel('log(KL divergence)')
    kl_ax.set_xlabel('iteration')
    # kl_ax.set_xlim(0, nsteps*1.1)

    kl_div_lim = res['run_params'].get('convergence_kl_threshold', 0.018)
    kl_ax.axhline(np.log10(kl_div_lim), linestyle='--', color='red', lw=2, zorder=1)

    kl_ax.legend(prop={'size': 5}, ncol=2, numpoints=1, markerscale=0.7)
    plt.title(str(objname) + ' kl')
    plt.show()
    # ''' #

    # field = kwargs['field']
    # CHANGING OBSERVED TO REST FRAME WAVELENGTH
    if field == 'cdfs':
        datname = '/home/jonathan/cdfs/cdfs.v1.6.11.cat'
        zname = '/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout'
    elif field == 'cosmos':
        datname = '/home/jonathan/cosmos/cosmos.v1.3.8.cat'  # main catalog
        zname = '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout'  # redshift catalog
    elif field == 'uds':
        datname = '/home/jonathan/uds/uds.v1.5.10.cat'
        zname = '/home/jonathan/uds/uds.v1.5.8.zout'

    with open(datname, 'r') as f:
        hdr = f.readline().split()
    dtype = np.dtype([(hdr[1], 'S20')] + [(n, np.float) for n in hdr[2:]])
    dat = np.loadtxt(datname, comments='#', delimiter=' ', dtype=dtype)

    with open(zname, 'r') as fz:
        hdr_z = fz.readline().split()
    dtype_z = np.dtype([(hdr_z[1], 'S20')] + [(n, np.float) for n in hdr_z[2:]])
    zout = np.loadtxt(zname, comments='#', delimiter=' ', dtype=dtype_z)

    idx = dat['id'] == objname  # array filled: False when dat['id'] != objname, True when dat['id'] == objname
    zred = zout['z_spec'][idx][0]  # z = z_spec
    if zred == -99:
        zred = zout['z_peak'][idx][0]  # if z_spec does not exist, z = z_phot
    print(zred)

    wave_rest = []  # REST FRAME WAVELENGTH
    for i in range(len(wave)):
        wave_rest.append(wave[i]/(1 + zred))  # 1 + z = l_obs / l_emit --> l_emit = l_obs / (1 + z)

    # PLOT MODEL SED BEST FIT, INPUT PHOT
    yerr = res['obs']['maggies_unc']
    plt.subplot(111, xscale="log", yscale="log")
    plt.errorbar(wave_rest, res['obs']['maggies'], yerr=yerr, marker='o', linestyle='', color='purple',
                 label='Observed photometry')
    plt.plot(wave_rest, phot, 'D', label='Model', color='b', markerfacecolor='None', markersize=10,
             markeredgewidth=1.25, markeredgecolor='k')  # label='Model at {},{}'.format(walker, iteration)
    plt.legend(loc="best", fontsize=20)
    plt.title(str(field) + '-' + str(objname) + ' SED')
    plt.plot(sps.wavelengths, spec, color='k', alpha=0.5)
    plt.xlabel(r'Wavelength (Rest) [$\AA$]')
    plt.ylabel('Maggies')
    plt.xlim(10**3, 2.5*10**4)
    plt.ylim(10**-5, 4*10**3)
    plt.show()

    # ''' #
    # PLOT CHI_SQ BESTFIT
    chi_sq = ((res['obs']['maggies'] - phot) / res['obs']['maggies_unc']) ** 2
    plt.plot(wave_rest, chi_sq, 'o', color='b')
    plt.title(str(objname) + r' $\chi^2$')
    plt.xlabel('Rest frame wavelength [angstroms]')
    plt.ylabel(r'$\chi^2$')
    plt.show()
    # 4 10-30, 1>70 ;; 3 10-30, 1>40 ;; 4 10-30 ;; 2 10-30, 1>50 ;; 2 10-30

    # 5 10-30, 1>50 ;; 2 10-30, 1>30
    # ''' #


if __name__ == "__main__":
    # don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--parfile')
    parser.add_argument('--outname')

    out_folder = 'out_efast/'  # 'out_etenmet/'  # 'out_masstest/'  # 'out_efico/'

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    files = {'outname': '', 'parfile': ''}
    for key in kwargs.keys():
        files[key] = kwargs[key]

    if files['outname'] == 'all':
        with open('Comp_10.dat', 'r') as comp:
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
        # cvg = [1600, 1900, 1400, 1050, 900, 1400, 800, 1150, 875, 1250, 1300, 1250, 1150, 1300, 1200,
        #        2150, 1600, 1050]
        cvg = []
        # 12105, 11462, 12533, 12552, 12903, 14808, 15124, 17189, 17342, 18561, 18742, 21076, 21442, 22768, 11063,
        # 17423, 8787, 15462
        for i in range(len(e_objs)):
            for infile in glob.glob(os.path.join('/home/jonathan/.conda/envs/snowflakes/lib/python2.7/' +
                                                 'site-packages/prospector/git/out/' + out_folder, e_objs[i] + '*.h5')):
                out_file = infile
                param_file = files['parfile']
                print(param_file, out_file)
                true_field = e_fs[i]
                true_obj = e_objs[i]

                objname, field, res, mod, walker, iteration = printer(out_file, masstest=False, obj_true=true_obj)
                # , cvg=cvg[i])
                sed(true_obj, true_field, res, mod, walker, iteration, param_file, **kwargs)

    else:
        out_file = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'+files['outname']
        param_file = files['parfile']
        print(param_file, out_file)

        objname, field, res, mod, walker, iteration = printer(out_file)

        sed(objname, field, res, mod, walker, iteration, param_file, **kwargs)

'''
RUNNING WITH:
python quickgrab.py --outname=10246_cdfs_multirun_1498677216_mcmc.h5 --parfile=eelg_multirun_params.py
'''