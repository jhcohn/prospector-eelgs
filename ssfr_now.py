import os
import numpy as np
import pickle
import one_stack as one
import get_mass_dust as gmd


def get_ssfr(file, mass):
    with open(file, 'rb') as exout:
        extra_output = pickle.load(exout)

    return extra_output['bfit']['sfr_100'] / mass

if __name__ == "__main__":

    git = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    boot = 0
    vary = 1
    mask = 0
    others = 0
    short = 0
    if vary:
        base = ['vary', 'vary']
        out_folds = ['out_evar/', 'out_nvar/']
        folders = ['pkl_evar/', 'pkl_nvar/']
        mass = [10.15, 10.51]
        import eelg_varymet_params as e_params
        import eelg_varymet_params as n_params
    elif mask:
        base = ['newmask', 'newmask']
        folders = ['pkl_emask/', 'pkl_nmask/']
        mass = [10.15, 10.51]
        import eelg_newmask_params as e_params
        import eelg_newmask_params as n_params
    elif others:
        base = ['thirty', 'nth']  # base = ['otherbins', 'nother']  # use for otherbins
        folders = ['etpkls/', 'ntpkls/']  # ['opkls/', 'nopkls/']
        mass = [10., 10.3]  # PLACEHOLDER
        import eelg_thirty_params as e_params
        import noelg_thirty_params as n_params
    elif short:
        base = ['short', 'short']
        folders = ['pkl_eshort/', 'pkl_nshort/']
        mass = [10.07, 10.39]
        import eelg_short_params as e_params
        import eelg_short_params as n_params
    else:
        base = ['fixedmet', 'noelg']  # use for fixedmet
        folders = ['pkls/', 'nmpkls/']
        mass = [9.98, 10.26]
        import eelg_fixedmet_params as e_params
        import noelg_multirun_params as n_params

    eelg_list = open('eelg_specz_ids', 'r')
    pkls = git + folders[0]
    eelgs = []
    e_objs = []
    e_fields = []
    for line in eelg_list:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            e_objs.append(cols[1])
            e_fields.append(cols[0])
            eelgs.append(cols[1] + '_' + cols[0] + '_' + base[0])  # base[0] = fixedmet (or otherbins)
    eelg_list.close()

    lbg_list = open('lbg_ids', 'r')
    flist = {}
    lbgs = []
    l_objs = []
    l_fields = []
    l_pkls = git + folders[1]
    for line in lbg_list:
        if int(line) - 200000 > 0:
            flist[str(int(line) - 200000)] = 'uds'
            lbgs.append(str(int(line) - 200000) + '_uds_' + base[1])  # base[1] = noelg (or nother)
            l_objs.append(int(line) - 200000)
            l_fields.append('uds')
        elif int(line) - 100000 > 0:
            flist[str(int(line) - 100000)] = 'cosmos'
            lbgs.append(str(int(line) - 100000) + '_cosmos_' + base[1])
            l_objs.append(int(line) - 100000)
            l_fields.append('cosmos')
        else:
            flist[str(int(line))] = 'cdfs'
            lbgs.append(str(int(line)) + '_cdfs_' + base[1])
            l_objs.append(int(line))
            l_fields.append('cdfs')
    lbg_list.close()

    ssfr = []
    for glxy in eelgs:
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            ssfr.append(get_ssfr(file, mass[0]))
        else:
            print(file)

    ssfr_l = []
    for glxy in lbgs:
        file = l_pkls + glxy + '_extra_out.pkl'

        if os.path.exists(file):
            ssfr_l.append(get_ssfr(file, mass[1]))
        else:
            print(file)

    print(np.percentile(ssfr, [16, 50, 84]))
    print(np.percentile(ssfr_l, [16, 50, 84]))

    newpath = '/home/jonathan/'
    comps = ['Comp_04.dat', 'Comp_10.dat', 'Comp_14.dat']  # 'Comp_12.dat',
    for file in comps:
        print(file)
        objs = []
        fields = []
        ids = open(newpath + file, 'r')  # open(base + 'specz_ids', 'r')
        for line in ids:
            cols = line.split()
            yep = 1
            if line[0] == '#':  # ignore commented lines
                yep = 0  # nope
            if yep:
                if int(cols[0]) - 200000 > 0:
                    objs.append(int(cols[0]) - 200000)
                    fields.append('uds')
                elif int(cols[0]) - 100000 > 0:
                    objs.append(int(cols[0]) - 100000)
                    fields.append('cosmos')
                else:
                    objs.append(int(cols[0]))
                    fields.append('cdfs')
        # print(objs)
        # print(fields)
        ids.close()

        pkls = git + folders[0]

        gals = []
        for i in range(len(objs)):  # for each obj in the given composite
            gals.append(str(objs[i]) + '_' + fields[i] + '_' + base[0])
        print(len(gals))

        tun = True
        if tun:  # flat ssfr prior = 1 / t_univ, based on perc of t_univ values for each population
            pri = [0.41772065 * 1e-9, 0.50135904 * 1e-9, 0.55399038 * 1e-9]  # USE
            pri_l = [0.40128419 * 1e-9, 0.44860297 * 1e-9, 0.56183993 * 1e-9]  # from ^ commented out thing  # USE

        comp_mass = []
        comp_met = []
        comp_dust = []
        # START STACKING
        t1 = []
        draws = []
        n1 = 0
        n0 = 0
        for glxy in gals:
            file = pkls + glxy + '_extra_out.pkl'
            if os.path.exists(file):
                oute = git + 'out/' + out_folds[0]
                for outs in os.listdir(oute):
                    if outs.startswith(glxy):
                        if outs.endswith(".h5"):
                            stuff = gmd.printer(oute + outs)
                            mass = stuff[0]
                n1 += 1
                comp_mass.append(mass)
                comp_dust.append(stuff[1])
                comp_met.append(stuff[2])
                temp = one.randraw(file, mass)
                draws.append(temp[0])
                t1.append(temp[1])
            else:
                n0 += 1
                print(file)
        print(n1, n0)

        # stacker2(draws, t1)
        sig = 1  # what sigma error to show on plot
        perc1 = one.stacker(draws, sigma=sig)
        print(np.percentile(comp_mass, [16., 50., 84]))
        print(np.percentile(comp_dust, [16., 50., 84]))
        print(np.percentile(comp_met, [16., 50., 84]))

        smooth_perc = one.smooth(perc1)
        one.plot_sfhs(smooth_perc, t1[0], elist=eelgs, uvj_in=True, sigma=sig, priors=pri, tuniv=True)

'''
# COMP_04:
Mass [  9.91665936  10.11143303  10.37669609]
Dust [ 0.1015627   0.27033682  0.42009685]
Met [-1.80391939 -1.69474024 -1.36637911]
# COMP_10:
Mass [  9.56806931   9.9428606   10.26535469]
Dust [ 0.1366495   0.29065425  0.4745783 ]
Met [-1.85909033 -1.72481412 -1.57049279]
# COMP_14:
Mass [ 10.06390381  10.25686312  10.45410967]
Dust [ 0.1176697   0.23385095  0.43429242]
Met [-1.80541593 -1.58386713 -1.02036411]
'''
