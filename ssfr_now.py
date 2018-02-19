import os
import numpy as np
import pickle
import one_stack as one
import get_mass_dust as gmd


def get_ssfr(file):
    with open(file, 'rb') as exout:
        extra_output = pickle.load(exout)
        # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000

        ssfr = np.percentile(extra_output['extras']['ssfr'][0], [16., 50., 84.])[1]  # median ssfr
    return extra_output['bfit']['sfr_100'], ssfr


def get_sfh(file):
    with open(file, 'rb') as exout:
        extra_output = pickle.load(exout)
        # print(len(extra_output['extras']['sfh']), len(extra_output['extras']['sfh'][0]))  # 22, 2000

        sfh = extra_output['extras']['sfh']
        return sfh


if __name__ == "__main__":

    git = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
    vary = 0
    fifty = 1
    if vary:
        base = ['vary', 'vary']
        out_folds = ['out/out_evar/', 'out/out_nvary/']
        folders = ['pkl_evar/', 'pkl_nvary/']
        # mass = [10.15, 10.51]
        import eelg_varymet_params as e_params
        import eelg_varymet_params as n_params
    elif fifty:
        base = ['fifty', 'vary']
        out_folds = ['out/out_efifty/', 'out/out_nvary/']
        folders = ['pkl_efifty/', 'pkl_nvary/']
        import eelg_fifty_params as e_params
        import eelg_fifty_params as n_params

    import stellar_ages as sa
    eelgs, lbgs1 = sa.get_gal_lists(base)
    pkls = git + folders[0]
    l_pkls = git + folders[1]

    sfh = []
    ssfr = []
    sfr = []
    # BUILD COMP
    comp = []
    eelg_list = open('Comp_10.dat', 'r')
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
    eelgs1 = []
    count = 0
    for file in eelgs:
        for i in range(len(comp)):
            if file.startswith(comp[i]):
                if comp[i] == '20366_cdfs_vary':
                    pass
                else:
                    print(file)
                    eelgs1.append(file)
    eelgs = eelgs1
    # END checking COMP
    print(len(eelgs))

    second = []
    for glxy in eelgs:
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            get = get_ssfr(file)
            sfr.append(get[0])
            ssfr.append(get[1])
            sfh.append(get_sfh(file)[0])
            second.append(get_sfh(file)[4])
        else:
            print(file)

    sfh_l = []
    ssfr_l = []
    sfr_l = []
    for glxy in lbgs1:
        file = l_pkls + glxy + '_extra_out.pkl'

        if os.path.exists(file):
            get = get_ssfr(file)
            sfr_l.append(get[0])
            ssfr_l.append(get[1])
            sfh_l.append(get_sfh(file)[0])
        else:
            print(file)
    lbgs = []
    for i in range(len(lbgs1)):
        if lbgs1[i] == '22919_cdfs_vary':
            pass
        else:
            lbgs.append(lbgs1[i])

    import get_mass_dust as gmd
    out = git + out_folds[0]
    out_l = git + out_folds[1]
    eelgs1 = []
    lbgs1 = []
    for file in os.listdir(out):
        if file.endswith(".h5"):
            eelgs1.append(file)
    for file in os.listdir(out_l):
        if file.endswith(".h5"):
            lbgs1.append(file)

    print(len(eelgs))
    frac = np.zeros(shape=(len(eelgs), 10**3))
    frac2 = np.zeros(shape=(len(eelgs), 10**3))
    for i in range(len(eelgs)):
        print(i)
        get_e = gmd.printer(out + eelgs1[i], percs=False)
        for k in range(10**3):
            # print((get_e[np.random.randint(len(get_e))]))
            # print(sfh[i][np.random.randint(len(sfh[i]))])
            frac[i, k] = sfh[i][np.random.randint(len(sfh[i]))] * (10**8) / (10 ** get_e[np.random.randint(len(get_e))])
            frac2[i, k] = second[i][np.random.randint(len(second[i]))] * (10 ** 8) /\
                          (10 ** get_e[np.random.randint(len(get_e))])
    frac_l = np.zeros(shape=(len(lbgs), 10**3))
    for i in range(len(lbgs)):
        print(i)
        get_l = gmd.printer(out_l + lbgs1[i], percs=False)
        for k in range(10**3):
            # print((get_l[np.random.randint(len(get_l))]))
            # print(sfh_l[i][np.random.randint(len(sfh_l[i]))])
            frac_l[i, k] = sfh_l[i][np.random.randint(len(sfh_l[i]))] * (10**8) / (10**get_l[np.random.randint(len(get_l))])

    all_fracs2 = []
    all_fracs = []
    for i in range(len(frac)):
        all_fracs.append(np.percentile(frac[i], [16., 50., 84.]))
        all_fracs2.append(np.percentile(frac2[i], [16., 50., 84.]))
    all_fracs_l = []
    for i in range(len(frac_l)):
        all_fracs_l.append(np.percentile(frac_l[i], [16., 50., 84.]))
    print(np.percentile(all_fracs, [16., 50., 84.]))  # 86.7 - 26.1, 26.1 - 9.1
    print(np.percentile(all_fracs2, [16., 50., 84.]))  # 21.3 - 5.1, 5.1 - 0.91
    print(np.percentile(all_fracs_l, [16., 50., 84.]))  # 28.7% - 8.95%, 8.95% - 2.5%
    print('me')

    print(np.percentile(sfr, [16, 50, 84]) / (10 ** 9.84233))
    print(np.percentile(sfr_l, [16, 50, 84]) / (10 ** 10.55351))

    print(np.percentile(sfr, [16, 50, 84])[1] * 10 ** 8 / (10 ** 9.84233))
    print(np.percentile(sfr_l, [16, 50, 84])[1] * 10 ** 8 / (10 ** 10.55351))
    '''
    ssfr = []
    ssfr_l = []
    for i in range(len(sfr)):
        ssfr.append(sfr[i] / mass[i])
    for i in range(len(mass)):
        ssfr_l.append(sfr_l[i] / mass_l[i])
    '''
    print(np.percentile(ssfr, [16., 50., 84.]))
    print(np.percentile(ssfr_l, [16., 50., 84.]))


    '''
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
