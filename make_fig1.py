import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import print_sfh
import uvj
import widths  # WIDTHS
from matplotlib import gridspec
from mpl_toolkits.axes_grid.inset_locator import inset_axes


def all_plots(fileset, objname, field, loc='upper left'):
    # FIGURES
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 2, height_ratios=[3, 1])
#    gs = gridspec.GridSpec(1, 2)
    # prepares fig to hold two sets side-by-side of: two axes, one atop other, size ratio 3:1

    ax1 = plt.subplot(gs[0])  # ax1 = bigger upper axis (left)
    # ax1.set_title(field[0] + '-' + objname[0])
    # plt.axvline(x=5270, color='k')  # proving no offset in two axes
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    fs = 20  # 30
    textx = 1100
    texty = 300  # 500

    for j in range(len(fileset)):
        print(j)
        # PLOT SEPARATE UVJ
        with open(fileset[j][1], 'rb') as res:
            results = pickle.load(res)

        with open(fileset[j][2], 'rb') as sed:
            sed = pickle.load(sed)

        with open(fileset[j][3], 'rb') as restwave:
            wave_rest = pickle.load(restwave)

        with open(fileset[j][4], 'rb') as spec:
            spec = pickle.load(spec)

        with open(fileset[j][5], 'rb') as spswave:
            sps_wave = pickle.load(spswave)

        with open(fileset[j][7], 'rb') as justchi:
            chi = pickle.load(justchi)

        # MASK
        # create mask
        mask = results['obs']['phot_mask']  # mask out -99.0 values
        wave_rest = np.asarray(wave_rest)
        # no longer need this once I use new output that creates wave_rest as an array
        # ly_mask = (1180 < wave_rest) & (wave_rest < 1260)  # mask out ly-alpha values
        # mask[ly_mask] = False  # combine the masks!

        # apply mask
        phot = results['obs']['maggies'][mask]
        unc = results['obs']['maggies_unc'][mask]
        sed = sed[mask]
        wave_rest = wave_rest[mask]
        chi = chi[mask]

        # CONVERTING TO microjansky
        factor = 3631 * 10 ** 6
        res_jan = []
        err_jan = []
        sed_jan = []
        for i in range(len(wave_rest)):
            res_jan.append(phot[i] * factor)
            err_jan.append(unc[i] * factor)
            sed_jan.append(sed[i] * factor)
        spec_jan = []
        for i in range(len(spec)):
            spec_jan.append(spec[i] * 3631 * 10 ** 6)

        if j == 0:
            ax1.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observed Photometry')  # plot observations
            ax1.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
            ax1.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            ax1.set_ylabel(r'$\mu$Jy', fontsize=fs)  # 30

            # ax1.text(1000, 50, 'EELG', fontsize=20)
            ax1.text(textx, texty, 'COSMOS-1824 (EELG)', fontsize=20)  # 700, 600
            # ax1.set_xlabel(r'Rest frame wavelength')
            ax3 = plt.subplot(gs[2], sharex=ax1)  # ax2 = smaller lower axis
            ax3.plot(wave_rest, chi, 'o', color='k')  # plot chi
            plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
            plt.setp(ax1.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
            yticks = ax3.yaxis.get_major_ticks()  # show ytick labels on lower axis
            yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
            ax3.set_ylabel(r'$\chi$', fontsize=fs)
            # ax3.set_xlabel(r'Rest frame wavelength', fontsize=fs)  # 30
            ax1.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            # plt.subplots_adjust(hspace=.0)

            loc_uvj = 1
            inset_axes(ax1, width="32%", height="28%", loc=loc_uvj)
            # NOTE: height = 0.875 width (3.5 units vs 4 units) 20%/0.875=0.22857 --> if height 20%, set width to 23%
            # create inset axis: width (%), height (inches), location
            # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
            # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
            uvj.uvj_plot(objname[j], field[j], title=False, labels=False, lims=True, size=20, show=False)
            print('uvj1')
            # add uvj plot to inset axis
        elif j == 1:
            print(j)
            ax2 = plt.subplot(gs[1], sharey=ax1)  # ax2 = bigger upper axis (right)
            # plt.axvline(x=5270, color='k')  # proving no offset in two axes
            ax2.set_yscale("log")
            ax2.set_xscale("log")
            # ax2.set_ylim(ymin, ymax)
            # ax2.set_title(field[1] + '-' + objname[1])
            ax2.errorbar(wave_rest, res_jan, yerr=err_jan, marker='o', linestyle='', color='r',
                         label=r'Observed Photometry')  # plot observations
            ax2.plot(wave_rest, sed_jan, 'o', color='b', label=r'Model Photometry')  # plot best fit model
            ax2.plot(sps_wave, spec_jan, color='b', alpha=0.5, label=r'Model Spectrum')  # plot spectrum
            # ax2.set_ylabel(r'$\mu$Jy')

            # ax2.text(1000, 50, 'LBG', fontsize=20)
            ax2.text(textx, texty, 'COSMOS-4708 (LBG)', fontsize=20)  # 1025, 600
            # ax2.set_xlabel(r'Rest frame wavelength')
            ax4 = plt.subplot(gs[3], sharex=ax2, sharey=ax3)
            ax4.plot(wave_rest, chi, 'o', color='k')  # plot chi
            plt.axhline(y=0, color='k')  # plot horizontal line at y=0 on chi plot
            plt.setp(ax2.get_xticklabels(), visible=False)  # hide xtick labels on upper axis
            plt.setp(ax2.get_yticklabels(), visible=False)  # hide ytick labels on upper axis
            plt.setp(ax4.get_yticklabels(), visible=False)  # hide ytick labels on lower axis
            plt.setp(ax4, yticks=[-8, -4, 0, 4, 8], yticklabels=['-8', '-4', '0', '4', '8'])
            # ax4.yaxis.get_major_ticks(numticks=5)  # show ytick labels on lower axis
            # yticks[-1].label1.set_visible(False)  # hide uppermost ytick label on lower axis to prevent overlap
            # ax4.set_ylabel(r'$\chi$')
            # ax4.set_xlabel(r'Rest frame wavelength', fontsize=fs)  # 30
            ax2.legend(numpoints=1, loc=loc, prop={'size': 20})  # , line2) ... , r'$\chi$']
            # plt.subplots_adjust(hspace=.0)

            loc_uvj = 1  # 7
            inset_axes(ax2, width="32%", height="28%", loc=loc_uvj)
            # create inset axis: width (%), height (inches), location
            # loc=1 (upper right), loc=2 (upper left) --> loc=3 (lower left), loc=4 (lower right); loc=7 (center right)
            # https://stackoverflow.com/questions/10824156/matplotlib-legend-location-numbers
            uvj.uvj_plot(objname[j], field[j], title=False, labels=False, lims=True, size=20, show=False)
            print('uvj2')
    print('show')
    plt.subplots_adjust(wspace=.0, hspace=.0)

    xmin = 10**3
    xmax = 2.5*10**4
    ymin = 10**-1  # 6*10**-2
    ymax = 4*10**3  # 10**4
    ax1.set_xlim(xmin, xmax)  # 700, xmax
    ax1.set_ylim(ymin, ymax)
    ax2.set_xlim(xmin, xmax)  # 10**3, xmax

    # xlabel!
    fig.text(0.5, 0.04, r'Rest frame wavelength', ha='center', va='bottom', fontsize=fs)  # 30

    # TICK PARAMS
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    ax2.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax2.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    ax3.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax3.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    ax4.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax4.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj1')
    parser.add_argument('--obj2')
    parser.add_argument('--field1')
    parser.add_argument('--field2')
    parser.add_argument('--base1')
    parser.add_argument('--base2')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    obj1 = kwargs['obj1']
    obj2 = kwargs['obj2']

    field1 = kwargs['field1']
    field2 = kwargs['field2']
    pre1 = 'pkls/' + obj1 + '_' + field1 + '_' + kwargs['base1']
    pre2 = 'nmpkls/' + obj2 + '_' + field2 + '_' + kwargs['base2']

    base1 = '_out.pkl'
    extra1 = pre1 + '_extra' + base1  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res1 = pre1 + '_res' + base1
    sed1 = pre1 + '_sed' + base1
    restwave1 = pre1 + '_restwave' + base1
    spec1 = pre1 + '_spec' + base1
    spswave1 = pre1 + '_spswave' + base1
    chisq1 = pre1 + '_chisq' + base1
    justchi1 = pre1 + '_justchi' + base1

    base2 = '_out.pkl'
    extra2 = pre2 + '_extra' + base2  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res2 = pre2 + '_res' + base2
    sed2 = pre2 + '_sed' + base2
    restwave2 = pre2 + '_restwave' + base2
    spec2 = pre2 + '_spec' + base2
    spswave2 = pre2 + '_spswave' + base2
    chisq2 = pre2 + '_chisq' + base2
    justchi2 = pre2 + '_justchi' + base2

    files1 = [extra1, res1, sed1, restwave1, spec1, spswave1, chisq1, justchi1]
    files2 = [extra2, res2, sed2, restwave2, spec2, spswave2, chisq2, justchi2]
    fileset = [files1, files2]
    objs = [obj1, obj2]
    fields = [field1, field2]

    all_plots(fileset, objs, fields, loc='upper left')

'''
Currently running with:
python make_fig1.py --obj1=1824 --field1=cosmos --base1=fixedmet --obj2=4708 --field2=cosmos --base2=noelg
'''
