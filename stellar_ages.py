import numpy as np
import os
import pickle
import random
import matplotlib.pyplot as plt


def randraw_sfr_perc(infile, num=1000):  # num=1000
    """
    For a given galaxy, randraw samples the posterior num times for each point in extra_output['extras']['sfh'][i]

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :param num: number of times to sample the galaxy posterior at each point in the SFH
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['sfh']), num))  # shape=(22, num)
    # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000
    # print(len(draw_from_sfh), len(draw_from_sfh[0]))  # 22, num

    for i in range(len(extra_output['extras']['sfh'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['sfh'][i][random.randint(0, num)]

    median_sfh = []
    for i in range(len(draw_from_sfh)):  # 22
        median_sfh.append(np.percentile(draw_from_sfh[i], [16., 50., 84.])[1])  # median value of draw_from_sfh[i]

    return median_sfh, extra_output['extras']['t_sfh']


def stacker(gal_draws):
    """
    stacker takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin
    gal_draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...]
    each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (i.e. comprised of the output from randraw) for a list of galaxies
    :return: perc = stored lists of the median and +/- sigma SFH values calculated from the gal_draws
    """

    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = num (1000)
    all_draws = np.zeros(shape=(len(gal_draws[0]), len(gal_draws[0][0]) * len(gal_draws)))
    for k in range(len(gal_draws)):
        # append the num=1000 values in each gal_draws[k] at each of the 22 points to all_draws:
        note = k * len(gal_draws[0][0])
        for i in range(len(gal_draws[k])):
            for j in range(len(gal_draws[k][i])):
                all_draws[i][note+j] += gal_draws[k][i][j]
    print(len(all_draws), len(all_draws[0]))  # 22, (number of galaxies in stack) * (num=1000)

    perc = np.zeros(shape=(len(gal_draws[0]), 3))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    for jj in xrange(len(gal_draws[0])):
        perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma

    return perc


def smooth(perc):
    """
    Takes the stacked sfh that is output from stacker and averages the sfr values in each bin, such that the sfr within
    each bin is flat, as is the case in the original extra_output['extras']['sfh'] output

    :param perc: stored lists of the median and +/1 1sigma SFH values calculated from the gal_draws
    :return: smoother = same shape as perc, but with all points within a bin averaged s.t. all points within a bin=flat
    """
    # from perc: bin1 0:2, bin2 3:6, bin3 7:10, bin4 11:14, bin5 15:18, bin6 19:22
    smoother = np.zeros(shape=(len(perc), 1))  # len(perc[0])))  # shape=(22, 3)
    yax = perc  # perc[:, j]
    for i in range(3):
        smoother[i] = (yax[0] + yax[1] + yax[2]) / 3
        smoother[i+19] = (yax[-1] + yax[-2] + yax[-3]) / 3
    for i in range(4):
        smoother[i+3] = (yax[3] + yax[4] + yax[5] + yax[6]) / 4
        smoother[i+7] = (yax[7] + yax[8] + yax[9] + yax[10]) / 4
        smoother[i+11] = (yax[11] + yax[12] + yax[13] + yax[14]) / 4
        smoother[i+15] = (yax[15] + yax[16] + yax[17] + yax[18]) / 4

    # print(smoother)
    return smoother


def stellar_age(sfr, agebins):
    # let's do it
    # from sfh: bin1 0:2, bin2 3:6, bin3 7:10, bin4 11:14, bin5 15:18, bin6 19:22
    sfr_per_bin = [sfr[1], sfr[4], sfr[8], sfr[13], sfr[17], sfr[21]]

    time_per_bin = np.diff(10**agebins, axis=-1)[:, 0]
    avg_bin_age = []  # average age in each bin
    for i in range(len(time_per_bin)):
        j = 0
        time = 0
        while j < i:
            time += time_per_bin[j]
            j += 1
        avg_bin_age.append(time + time_per_bin[i]/2)

    # SUM((SFR in each bin) * (time in each bin) * (average age of bin)) / SUM(SFR in each bin * time in each bin)
    num = []
    denom = []
    for i in range(len(time_per_bin)):
        num.append(sfr_per_bin[i] * time_per_bin[i] * avg_bin_age[i])
        denom.append(sfr_per_bin[i] * time_per_bin[i])
    numerator = np.sum(num)
    denominator = np.sum(denom)
    age = numerator / denominator
    return age


def get_gal_lists(base):
    eelg_list = open('eelg_specz_ids', 'r')
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

    return eelgs, lbgs


if __name__ == "__main__":

    vary = True
    if vary:
        folders = ['pkl_evar/', 'pkl_nvar/']
        base = ['vary', 'vary']
        import eelg_varymet_params as param
    else:
        folders = ['pkls/', 'nmpkls/']
        base = ['fixedmet', 'noelg']
        import eelg_fixedmet_params as param

    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
    l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]
    eelgs, lbgs = get_gal_lists(base)

    # START STACKING
    t1 = []
    e_draws = []
    boots = []
    num_e = 0
    e_sample = 0
    for glxy in eelgs:
        e_sample += 1  # number of galaxies in sample
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            num_e += 1  # number of galaxies for which we have output
            temp = randraw_sfr_perc(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
            e_draws.append(smooth(temp[0]))
            t1.append(temp[1])
        else:
            print(file)  # print galaxy if pkls don't exist for it
    print('enums', num_e, e_sample)

    l_draws = []
    num_l = 0
    l_sample = 0
    for glxy in lbgs:
        l_sample += 1  # number of galaxies in sample
        file = l_pkls + glxy + '_extra_out.pkl'

        if os.path.exists(file):
            num_l += 1  # number of galaxies for which we have output
            temp = randraw_sfr_perc(file)
            l_draws.append(smooth(temp[0]))
        else:
            print(file)  # print galaxy if pkls don't exist for it
    print('lnums', num_l, l_sample)

    model = param.load_model(objname='1824', field='cosmos')
    agebins = model.params['agebins']

    e_age = []
    for i in range(len(e_draws)):  # len(eelgs)
        e_age.append(stellar_age(e_draws[i], agebins) / 1e9)
    e_age_percs = np.percentile(e_age, [16., 50., 84.])

    l_age = []
    for i in range(len(l_draws)):  # len(eelgs)
        l_age.append(stellar_age(l_draws[i], agebins) / 1e9)
    l_age_percs = np.percentile(l_age, [16., 50., 84.])

    fig = plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    ax1.hist(e_age, bins=50, histtype="step", weights=[1./num_e]*len(e_age), normed=False, color='b', lw=2,
             label='EELGs')
    ax1.hist(l_age, bins=50, histtype="step", weights=[1./num_l]*len(l_age), normed=False, color='r', lw=2,
             label='LBGs')

    # plot median, +/-1sigma for both histograms
    ax1.axvline(x=e_age_percs[1], color='b', linestyle='--', lw=2)
    ax1.axvline(x=l_age_percs[1], color='r', linestyle='--', lw=2)

    # shade in +/-1sigma region
    ax1.axvspan(e_age_percs[0], e_age_percs[2], color='b', alpha=0.2)
    ax1.axvspan(l_age_percs[0], l_age_percs[2], color='r', alpha=0.2)
    # print(e_age_percs)
    print(e_age_percs[1] - e_age_percs[0], e_age_percs[1], e_age_percs[2] - e_age_percs[1])
    # print(l_age_percs)
    print(l_age_percs[1] - l_age_percs[0], l_age_percs[1], l_age_percs[2] - l_age_percs[1])

    ax1.set_ylim(0, 0.1)

    # figure labels
    fs = 20
    ax1.legend(numpoints=1, loc='upper left', prop={'size': fs})
    ax1.set_xlabel('Stellar ages [Gyr]', ha='center', fontsize=fs)
    ax1.set_ylabel(r'Fraction of galaxies', fontsize=fs)
    plt.show()
