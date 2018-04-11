import matplotlib.pyplot as plt
import numpy as np
import sys
import argparse
import os
import stellar_ages as sa
import get_mass_dust as gmd
np.errstate(invalid='ignore')


def get_fast(file, objlist=None):
    '''
    Get fast properties from given fast file

    :param file: output file with properties from FAST
    :return: output, a dictionary with redshift, mass stored for each object
    '''
    if objlist is None:
        objs = []
        masses = []
        zs = []
        # sfrs = []
        with open(file) as tm:
            for line in tm:
                if line[0] != '#':
                    cols = line.split()
                    objs.append(cols[0])
                    zs.append(cols[1])
                    masses.append(cols[2])
                    # if file[-5] == '2':
                        # sfrs.append(cols[3])

        # RENAME AND MAKE IT A DICTIONARY
        for i in range(len(objs)):
            if int(objs[i]) > 200000:
                objs[i] = str(int(objs[i]) - 200000)
            elif int(objs[i]) > 100000:
                objs[i] = str(int(objs[i]) - 100000)
        # make it a dict:
        output = {}
        for i in range(len(objs)):
            output[objs[i]] = [zs[i], masses[i]]  # , sfrs[i]]
    else:
        output = {}
        for thing in objlist:
            with open(file) as tm:
                for line in tm:
                    if line[0] != '#':
                        cols = line.split()
                        if cols[0] == thing:
                            output[thing] = [cols[1], cols[2]]
    return output


def compare_gmd(dict, folder, include_met=False):
    '''
    Compare mass calculated from Prospector with mass calculated from FAST

    :param dict: dictionary, with mass calculated from FAST (output of get_fast() function)
    :param folder: output folder in which Prospector output of relevant galaxies are stored
    :return: mass_diff, which is a dictionary with the FAST mass, Prospector mass stored for each object
    '''
    exists = []
    out = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + folder
    for thing in os.listdir(out):
        if thing.endswith(".h5"):
            exists.append(thing)

    ### NEW COMP
    mass_diff = {}
    met_dict = {}
    dust_dict = {}
    count = 0
    for thing in exists:
        for key in dict:
            if thing.startswith(key):
                get = gmd.printer(out + thing)
                mass_diff[key] = [dict[key][1], get[0]]
                dust_dict = get[1]
                met_dict[key] = get[2]
    print(count)

    if include_met:
        return mass_diff, met_dict  # , dust_dict
    else:
        return mass_diff  # mass_diff[key][0] is the FAST mass, mass_diff[key][1] is the Prospector mass


if __name__ == "__main__":
    fico = 1
    delta = 1

    if fico:
        folders = ['out_efico/', 'out_nfico/']
        pars = ['eelg_fifty_params.py', 'eelg_fifty_params.py']
        base = ['fico', 'fico']
    '''
    elif vary:
        folders = ['out_evar/', 'out_nvary/']
        pars = ['eelg_varymet_params.py', 'eelg_varymet_params.py']
        base = 'vary'
    elif fifty:
        folders = ['out_efifty/', 'out_efifty/']  # out_nfifty
        pars = ['eelg_fifty_params.py', 'eelg_fifty_params.py']
        base = 'fifty'
    elif fix:
        folders = ['out_efix/', 'out_efifty/']  # out_nfifty
        pars = ['eelg_fixedmet_params.py', 'eelg_fifty_params.py']
        base = 'fix'
    elif newu:
        folders = ['out_enewu/', 'out_enewu/']  # out_nfifty
        pars = ['eelg_newu_params.py', 'eelg_newu_params.py']
        base = 'newu'
    elif others:
        folders = ['out_ethvary/', 'out_nthvary/']
        pars = ['eelg_thvary_params.py', 'eelg_thvary_params.py']
        base = 'thvary'
    elif short:
        folders = ['out_eshort/', 'out_nshort/']
        pars = ['eelg_short_params.py', 'eelg_short_params.py']
        base = 'short'
    else:
        folders = ['out_fixedmet/', 'out_noelg/']
        pars = ['eelg_fixedmet_params.py', 'noelg_multirun_params.py']
    '''

    obj_e, field_e, obj_l, field_l = sa.get_gal_lists(base, objlists=True)
    field_dict = {}
    for i in range(len(obj_e)):
        field_dict[obj_e[i]] = field_e[i]

    home = '/home/jonathan/'
    three_masses = [home + 'Comp_10_zm_EL_Z002.dat', home + 'Comp_10_zm_EL_Z004.dat', home + 'Comp_10_zm_ZFOURGE.dat']
    labels = [r'FAST Z = Z$_{\odot}$, with emission lines', r'FAST Z = Z$_{\odot}/5$, with emission lines',
              r'FAST (ZFOURGE catalog parameters)']
    colors = ['r', 'b', 'purple']
    shapes = ['o', 's', 'v']

    for i in range(len(three_masses)):
        dictionary = get_fast(three_masses[i])
        mass_dict = compare_gmd(dictionary, folders[0])
        # print(mass_dict)

        x, y, xratio, xfield, ratio = [], [], [], [], []
        for key in mass_dict:  # mass_diff[key][0] is the FAST mass, mass_diff[key][1] is the Prospector mass
            # print(mass_dict[key][0], mass_dict[key][1])
            x.append(float(mass_dict[key][0]))  # FAST
            y.append(float(mass_dict[key][1]))  # Prospector
            ratio.append((10**float(mass_dict[key][1]))/(10**float(mass_dict[key][0])))
            xratio.append(key)
            for f_key in field_dict:
                if int(key) == int(f_key):
                    xfield.append(field_dict[f_key])

        # PLOT!
        print(len(xfield), len(xratio))
        print(xfield)
        print(xratio, 'hi')
        for l in range(len(xfield)):
            if xfield[l] == 'cosmos':
                xfield[l] = 'cos'
        fs_text = 20  # 30
        fs = 20
        # delta = False  # make this a delta-mass plot
        if delta:
            xlabs = []
            for j in range(len(xratio)):
                xlabs.append(str(xfield[j]).upper() + xratio[j])
            median = False
            line = True
            if median:
                xlabs.append('Median')
                avg = np.median(ratio)
                ratio.append(avg)
                xspace = np.linspace(0, 20, len(ratio))
            if line:
                plt.axhline(y=np.median(ratio), color=colors[i])
            xspace = np.linspace(0, 19, len(ratio))
            plt.scatter(xspace, ratio, marker=shapes[i], color=colors[i], label=labels[i], s=fs)
            plt.xticks(xspace, xlabs, rotation=60)  # 'vertical')
            print(xlabs)
        else:
            m, b = np.polyfit(x, y, 1)
            liney = [m*j + b for j in x]
            plt.scatter(x, y, marker=shapes[i], color=colors[i], label=labels[i], s=fs)
            plt.plot(x, liney, linestyle='-', color=colors[i])
            plt.xlabel(r'M$_{\odot}$, FAST, ' + labels[i], fontsize=fs_text)
    if delta:
        plt.ylabel(r'Stellar mass ratio (Prospector / FAST)', fontsize=fs_text)
        # plt.xlabel(r'Galaxy IDs', fontsize=fs_text)
        plt.axhline(y=1., color='k', linestyle='--')
        plt.ylim(0, 17)  # 37)
        tick_f = 15
        # plt.yscale('log')
        plt.legend(numpoints=1, loc='upper left', prop={'size': 20})

    else:
        plt.xlabel(r'M$_{\odot}$ (FAST)', fontsize=fs_text)
        plt.ylabel(r'M$_{\odot}$ (Prospector)', fontsize=fs_text)
        plt.plot([7., 12], [7., 12], color='k', linestyle='--')
        plt.xlim(7.8, 11.2)  # plt.xlim(7.8, 10.2)
        plt.ylim(7.8, 11.2)  # plt.ylim(9.0, 10.7)
        tick_f = 20
        plt.gca().set_aspect('equal', adjustable='box')
        plt.legend(numpoints=1, loc='upper left', prop={'size': 20})
    plt.tick_params('x', length=3, width=1, which='both', labelsize=tick_f)
    plt.tick_params('y', length=3, width=0.5, which='both', labelsize=tick_f)

    plt.show()
