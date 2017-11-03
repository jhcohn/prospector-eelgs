import smart_stack as stack
import os
import numpy as np
import matplotlib.pyplot as plt


def b_stacker(gal_draws, sigma=1):
    """
    stacker takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin
    gal_draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...]
    each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (i.e. comprised of the output from randraw) for a list of galaxies
    :param sigma: how many sigma of error we want to show in the plot
    :return: perc = stored lists of the median and +/- sigma SFH values calculated from the gal_draws
    """

    '''
    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = num (1000)
    # print(len(gal_draws[0])) = 1000
    # all_draws = np.zeros(shape=(len(gal_draws[0]), len(gal_draws[0][0]) * len(gal_draws)))
    all_draws = np.zeros(shape=(22, 1000))
    for k in range(len(gal_draws)):
        # append the num=1000 values in each gal_draws[k] at each of the 22 points to all_draws:
        for i in range(len(gal_draws[k])):  # range(num=1000)
            all_draws[i][k] += gal_draws[k][i]
    print(len(all_draws), len(all_draws[0]))  # 22, (number of galaxies in stack) * (num=1000)
    '''
    # print(len(gal_draws[0]))
    # print(len(gal_draws[0][0]))
    # print(len(gal_draws))
    all_draws = np.zeros(shape=(len(gal_draws[0]), len(gal_draws[0][0]) * len(gal_draws)))
    for k in range(len(gal_draws)):
        # append the num=1000 values in each gal_draws[k] at each of the 22 points to all_draws:
        note = k * len(gal_draws[0][0])
        for i in range(len(gal_draws[k])):
            for j in range(len(gal_draws[k][i])):
                all_draws[i][note+j] += gal_draws[k][i][j]
    # print(len(all_draws), len(all_draws[0]))  # 22, (number of galaxies in stack) * (num=1000)

    # print(len(gal_draws), len(gal_draws[0]), 'hi')
    perc = np.zeros(shape=(len(gal_draws[0]), 2*sigma + 1))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    for jj in xrange(len(gal_draws[0])):
        if sigma == 1:
            perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma
        elif sigma == 3:
            perc[jj, :] = np.percentile(all_draws[jj, :], [0.3, 2.4, 16.0, 50.0, 84.0, 97.6, 99.7])  # out to 3 sigma

    return perc


def plot_boots(draw_lists, t, n, visualize=True, lw=1, spec=True, sigma=1):
    """
    Plots SFH stacks for two different galaxy samples side-by-side

    :param draw_lists: list comprised of two lists of gal_draws, which have been repeated n times on each sample and
    appended together
    :param t: time vector output by randraw
    :param n: number of times we're repeating the bootstrap process
    :param lw: line width
    :param spec: if stacking specific SFR instead of plain SFR, spec=True
    :param sigma: how many sigma of error we want to show in the plot
    :return: plot
    """

    # INITIALIZE FIGURE
    fig = plt.figure()
    ax1 = plt.subplot(1, 2, 1)
    ax2 = plt.subplot(1, 2, 2)  # , sharey=ax1, sharex=ax1)  # don't need to share axis if plot same region & never zoom

    # INITIALIZE MATRICES FOR STORAGE OF SSFH MEDIANS & 1SIGMAS (INITIAL AND SMOOTHED)
    perc = np.zeros(shape=(n, len(draw_lists[0][0]), 2 * sig + 1))
    # len(perc)=N=number of bootstacks, len(perc[0])=len(draws[0])=22, len(perc[0][0])=2*sig+1=3
    perc2 = np.zeros(shape=(n, len(draw_lists[1][0]), 2 * sig + 1))
    smooth_perc = np.zeros(shape=(len(perc), len(perc[0]), len(perc[0][0])))  # shape=(N, 22, 3)
    smooth_perc2 = np.zeros(shape=(len(perc2), len(perc2[0]), len(perc2[0][0])))  # shape=(N, 22, 3)

    # CUT EACH DRAWS LIST EVERY 133 (OR 87) GALAXIES
    # THEN FOR EACH OF THOSE SETS OF 133 (OR 87) GALAXIES, STACK SSFRs IN EACH BIN (using b_stacker), THEN SMOOTH
    init = 0
    init2 = 0  # inits (and cuts, below) used to select correct slice (set) of 133 (or 87) galaxies from draw lists
    for i in range(n):  # len(draws) = len(draw_lists[0 or 1]) = N*len(eelgs or lbgs)
        cut = (i+1) * len(draw_lists[0]) / n  # 133 or 87, depending on eelgs or lbgs
        cut2 = (i+1) * len(draw_lists[1]) / n  # 133 or 87, depending on eelgs or lbgs

        # hopefully unnecessary safety cut
        if cut >= len(draw_lists[0]):
            cut = len(draw_lists[0])
        if cut2 >= len(draw_lists[1]):
            cut2 = len(draw_lists[1])
        print(init, cut, 'init, cut')
        print(init2, cut2, 'init2, cut2')

        # CALCULATE PERCS FOR BOOTSTRAPPED STACK
        perc[i, :, :] = b_stacker(draw_lists[0][init:cut], sigma=sig)  # stack.stacker --> b_stacker
        perc2[i, :, :] = b_stacker(draw_lists[1][init2:cut2], sigma=sig)
        # stacker function prints (22, 1000*(number of galaxies=N*len(eelgs or lbgs)))

        # SMOOTH PERCS
        smooth_perc[i, :, :] = stack.smooth(perc[i])
        smooth_perc2[i, :, :] = stack.smooth(perc2[i])

        init = cut
        init2 = cut2

    # INIT NEW LIST OF SMOOTHED PERCS FOR ALL BOOTSTRAPS IN BOTH SAMPLES
    percs_list = [smooth_perc, smooth_perc2] # len(percs_list)=2, len(percs_list[i])=len(draws)=N*len(galaxies)
    # len(percs_list[i][j][k])=3, len(percs_list[i][j])=len(draws[j])=22
    alpha = 0.03  # use for alpha while plotting individual bootstraps as errorbars in each bin
    if sigma == 1:
        if visualize:
            # INSERT VERTICAL DOTTED LINES SHOWING BIN EDGES ON EACH AXIS
            ax1.axvline(10**-1, ls='--')
            ax1.axvline(3*10**-1, ls='--')
            ax1.axvline(1, ls='--')
            ax1.axvline(1.3, ls='--')
            ax1.axvline(1.5, ls='--')
            ax2.axvline(10**-1, ls='--')
            ax2.axvline(3*10**-1, ls='--')
            ax2.axvline(1, ls='--')
            ax2.axvline(1.3, ls='--')
            ax2.axvline(1.5, ls='--')

            # initialize bootstrap net median & upper/lower (+/-)1 sigma
            usig1 = np.zeros(shape=(N, 6))
            usig2 = np.zeros(shape=(N, 6))
            lsig1 = np.zeros(shape=(N, 6))
            lsig2 = np.zeros(shape=(N, 6))
            med1 = np.zeros(shape=(N, 6))
            med2 = np.zeros(shape=(N, 6))
            glx = 0

            # for each perc'd stack (i.e. for each bootstrap that has been stacked):
            for i in range(len(percs_list[0])):  # range(N) (FOR EELGs)
                glx += 1  # index of stack; max index = n = num of stacks
                # CREATE new_t: ONLY WANT ONE POINT IN EACH BIN, SO DON'T NEED ORIGINAL FULL t VECTOR
                new_t = [10**-2, 10**-1, 3*10**-1, 1, 1.3, 1.5]  # left sides of bins
                delta = [0.02, 0.015, 0.01, 0.001, 0.001, 0.001]  # OFFSET EACH POINT IN new_t WITH EACH ITERATION
                idx = 0  # index indicating which bin (6 bins --> idx runs 0 to 5)
                # Num of points in each bin (22 total): bin1: 4, bin2: 4, bin3: 4: bin4: 4, bin5: 4, bin6: 2
                for j in (0, 4, 8, 12, 16, 20):  # at 1 point in each bin (formerly each of the 22 points)
                    new_t[idx] += new_t[idx] * delta[idx] * glx  # create offset in x-axis (add ~1% * index)
                    print(new_t[idx], 't')
                    '''
                    # print(glx/len(percs_list[0]), len(percs_list[0]), glx/2, 'hi', glx, i)  # python printing (1/2)=0?
                    # print(new_t[idx], i, (1 + i * (1 / len(percs_list[0]))))  # python what the fuck
                    '''
                    yerr = [[percs_list[0][i][j, 0]], [percs_list[0][i][j, 2]]]  # -1sigma, +1 sigma
                    lsig1[i, idx] = percs_list[0][i][j, 0]  # -1sigma
                    usig1[i, idx] = percs_list[0][i][j, 2]  # +1sigma
                    med1[i, idx] = percs_list[0][i][j, 1]  # median
                    # PRINT MEDIAN, (+/-)1sigma as point with errorbars
                    ax1.errorbar(new_t[idx], percs_list[0][i][j, 1], yerr, color='gray', alpha=0.5, fmt='o', lw=lw)
                    # 1 pt at a time
                    idx += 1

            # REPEAT FOR LBGs
            glx2 = 0
            for i in range(len(percs_list[1])):
                glx2 += 1
                new_t2 = [10**-2, 10**-1, 3*10**-1, 1, 1.3, 1.5]  # left sides of bins
                delta2 = [0.02, 0.015, 0.01, 0.001, 0.001, 0.001]  # 0.05
                idx2 = 0
                for j in (0, 4, 8, 12, 16, 20):
                    new_t2[idx2] += new_t2[idx2] * delta2[idx2] * glx2  # create offset in x-axis (add ~1% * index)
                    yerr2 = [[percs_list[1][i][j, 0]], [percs_list[1][i][j, 2]]]  # -1sigma, +1 sigma
                    lsig2[i, idx2] = percs_list[1][i][j, 0]  # -1sigma
                    usig2[i, idx2] = percs_list[1][i][j, 2]  # +1sigma
                    med2[i, idx2] = percs_list[1][i][j, 1]  # median
                    # PRINT MEDIAN, (+/-)1sigma as point with errorbars
                    ax2.errorbar(new_t2[idx2], percs_list[1][i][j, 1], yerr2, color='gray', alpha=0.5, fmt='o', lw=lw)
                    idx2 += 1

            # CALCULATE AND PLOT PERCS OF ALL BOOTSTRAPS IN EACH BIN
            new_t_med = [0.09, 0.25, 0.9, 1.2, 1.4, 2.10]  # right sides of bins

            # initializing new set of percs, comprised from bootstraps
            things1 = np.zeros(shape=(6, 3))
            things2 = np.zeros(shape=(6, 3))
            uthings1 = np.zeros(shape=(6, 3))
            uthings2 = np.zeros(shape=(6, 3))
            lthings1 = np.zeros(shape=(6, 3))
            lthings2 = np.zeros(shape=(6, 3))

            # FILL THINGS WITH PERCS FROM BOOTSTRAPS
            for i in range(len(new_t_med)):  # for each bin:
                # perc for median for EELG, LBG bootstrapped stacks
                things1[i, :] = np.percentile(med1[:, i], [16.0, 50.0, 84.0])
                things2[i, :] = np.percentile(med2[:, i], [16.0, 50.0, 84.0])
                # perc for +1sigma for EELG, LBG bootstrapped stacks
                uthings1[i, :] = np.percentile(usig1[:, i], [16.0, 50.0, 84.0])
                uthings2[i, :] = np.percentile(usig2[:, i], [16.0, 50.0, 84.0])
                # perc for -1sigma for EELG, LBG bootstrapped stacks
                lthings1[i, :] = np.percentile(lsig1[:, i], [16.0, 50.0, 84.0])
                lthings2[i, :] = np.percentile(lsig2[:, i], [16.0, 50.0, 84.0])

            # PLOT PERCS BUILT FROM ALL BOOTSTRAP STACKS IN EACH BIN
            for i in range(len(things1)):  # for each bin:
                print(i)
                # in each bin, plot median of medians, then plot medians of (+/-)1sigmas as errorbars
                ax1.errorbar(new_t_med[i], things1[i, 1], yerr=[[lthings1[i, 1]], [uthings1[i, 1]]], color='k', fmt='o',
                             lw=lw)
                ax2.errorbar(new_t_med[i], things2[i, 1], yerr=[[lthings2[i, 1]], [uthings2[i, 1]]], color='k', fmt='o',
                             lw=lw)

        else:
            for i in range(len(percs_list[0])):  # range(N)
                ax1.plot(t, percs_list[0][i][:, 1], '-', color='k', lw=lw)  # median
                ax1.fill_between(t, percs_list[0][i][:, 0], percs_list[0][i][:, 2], color='k', alpha=alpha)  # fill between +/- 1sigma
                ax1.plot(t, percs_list[0][i][:, 0], '-', color='k', alpha=alpha, lw=lw)  # -1sigma
                ax1.plot(t, percs_list[0][i][:, 2], '-', color='k', alpha=alpha, lw=lw)  # +1sigma
            for i in range(len(percs_list[1])):
                ax2.plot(t, percs_list[1][i][:, 1], '-', color='k', lw=lw)  # median
                ax2.fill_between(t, percs_list[1][i][:, 0], percs_list[1][i][:, 2], color='k', alpha=alpha)  # fill between +/- 1sigma
                ax2.plot(t, percs_list[1][i][:, 0], '-', color='k', alpha=alpha, lw=lw)  # -1sigma
                ax2.plot(t, percs_list[1][i][:, 2], '-', color='k', alpha=alpha, lw=lw)  # +1sigma

    elif sigma == 3:
        if spec:
            ymin, ymax = 1e-13, 1e-6
        else:
            ymin, ymax = 1e-3, 1e4
        for i in range(len(percs_list[0])):  # range(N)
            ax1.plot(t, percs_list[0][i][:, 3], '-', color='k', lw=lw)  # median
            ax1.fill_between(t, percs_list[0][i][:, 2], percs_list[0][i][:, 4], color='k', alpha=0.3)  # fill between +/- 1sigma
            ax1.fill_between(t, percs_list[0][i][:, 1], percs_list[0][i][:, 5], color='k', alpha=0.3)  # fill between +/- 2sigma
            ax1.fill_between(t, percs_list[0][i][:, 0], percs_list[0][i][:, 6], color='k', alpha=0.3)  # fill between +/- 3sigma
            ax1.plot(t, percs_list[0][i][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -3sigma
            ax1.plot(t, percs_list[0][i][:, 6], '-', color='k', alpha=0.3, lw=lw)  # +3sigma
        for i in range(len(percs_list[1])):
            ax2.plot(t, percs_list[1][i][:, 3], '-', color='k', lw=lw)  # median
            ax2.fill_between(t, percs_list[1][i][:, 2], percs_list[1][i][:, 4], color='k', alpha=0.3)  # fill between +/- 1sigma
            ax2.fill_between(t, percs_list[1][i][:, 1], percs_list[1][i][:, 5], color='k', alpha=0.3)  # fill between +/- 2sigma
            ax2.fill_between(t, percs_list[1][i][:, 0], percs_list[1][i][:, 6], color='k', alpha=0.3)  # fill between +/- 3sigma
            ax2.plot(t, percs_list[1][i][:, 0], '-', color='k', alpha=0.3, lw=lw)  # -3sigma
            ax2.plot(t, percs_list[1][i][:, 6], '-', color='k', alpha=0.3, lw=lw)  # +3sigma

    # TECHNICAL FIG STUFF
    ymin, ymax = 1e-11, 1e-6
    label = r'Stacked sSFH [yr$^{-1}$]'

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim(ymin, 10**-7)  # , ymax
    ax1.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    ax1.set_ylabel(label, fontsize=30)  # 30
    # ax1.text(4, 10**2.5, 'EELGs', fontsize=30)
    # ax1.text(1, 4*10**-8, 'EELGs', fontsize=30)  # use if not uvj_in
    ax1.text(2*10**-2, 4*10**-8, 'EELGs', fontsize=30)  # use if uvj_in

    ax2.set_yscale("log")
    ax2.set_xscale("log")
    ax2.set_ylim(ymin, 10**-7)  # , ymax
    ax2.set_xlim(10**-2, 2.5)  # (0, 2.5)  # (10**-2, 13.6)
    # ax2.text(4, 10**2.5, 'LBGs', fontsize=30)
    # ax2.text(1, 4*10**-8, 'LBGs', fontsize=30)  # use if not uvj_in
    ax2.text(2*10**-2, 4*10**-8, 'LBGs', fontsize=30)  # use if uvj_in

    plt.setp(ax2.get_yticklabels(), visible=False)  # hide y-axis labels on right-hand subplot to prevent overlap
    plt.subplots_adjust(wspace=0.05)  # vertical whitespace (i.e. the width) between the two subplots
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)
    plt.rcParams.update({'font.size': 22})
    fig.text(0.5, 0.04, 'Lookback time [Gyr]', ha='center', fontsize=30)  # 30
    plt.tight_layout()
    plt.show()

others = False
if others:
    base = ['otherbins', 'nother']  # use for otherbins
    folders = ['opkls/', 'nopkls/']
else:
    base = ['fixedmet', 'noelg']  # use for fixedmet
    folders = ['pkls/', 'nmpkls/']

eelg_list = open('eelg_specz_ids', 'r')
pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
eelgs = []
for line in eelg_list:
    if line[0] == '#':
        pass
    else:
        cols = line.split()
        eelgs.append(cols[1] + '_' + cols[0] + '_' + base[0])  # base[0] = fixedmet (or otherbins)
eelg_list.close()

# NOTE: NEED TO RUN OUTPUT ON 7624_cdfs_fixedmet!!!!!!!
i = 0
while i < len(eelgs):
    if eelgs[i] == '7624_cdfs_fixedmet':
        eelgs.pop(i)
    i += 1
# print(eelgs)

lbg_list = open('lbg_ids', 'r')
flist = {}
lbgs = []
l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]
for line in lbg_list:
    if int(line) - 200000 > 0:
        flist[str(int(line) - 200000)] = 'uds'
        lbgs.append(str(int(line) - 200000) + '_uds_' + base[1])  # base[1] = noelg (or nother)
    elif int(line) - 100000 > 0:
        flist[str(int(line) - 100000)] = 'cosmos'
        lbgs.append(str(int(line) - 100000) + '_cosmos_' + base[1])
    else:
        flist[str(int(line))] = 'cdfs'
        lbgs.append(str(int(line)) + '_cdfs_' + base[1])
lbg_list.close()

# cut = 'thirty_full_max'  # 'thirty_30'
N = 500  # 200
# newfile = open('booties/bootlists-' + str(cut) + '.txt', 'w+')
boot_array = np.zeros(shape=(N, 2))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
sig = 1

draws = []
draws2 = []
for num in range(N):
    print('loop:', num)

    # BOOTSTRAP POPULATIONS
    eelgs = stack.bootstrap(np.asarray(eelgs))
    lbgs = stack.bootstrap(np.asarray(lbgs))

    '''
    # WRITE GALAXIES USED IN BOOTSTRAP
    newfile.write(str(num) + '\n \n [ ')
    for i in range(len(eelgs)):
        newfile.write(str(eelgs[i]) + ' ')
    newfile.write('] \n \n')
    newfile.write('[ ')
    for i in range(len(lbgs)):
        newfile.write(str(lbgs[i]) + ' ')
    newfile.write('] \n \n')
    '''

    t1 = []
    # draws = []
    boots = []
    nummy = 0
    c = 0
    for glxy in eelgs:
        c += 1
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            nummy += 1
            temp = stack.randraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
            draws.append(temp[0])
            t1.append(temp[1])
        else:
            print(file)

    # stacker2(draws, t1)
    # sig = 1  # what sigma error to show on plot (now done earlier in the code)
    # perc1 = stack.stacker(draws, sigma=sig)  # BUCKET

    # draws2 = []
    numl = 0
    cl = 0
    # t2 = []
    for glxy in lbgs:
        cl += 1
        file = l_pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            numl += 1
            temp = stack.randraw(file)
            draws2.append(temp[0])
        else:
            print(file)

    print(nummy, c, 'nume')
    print(numl, cl, 'numl')

    # perc2 = stack.stacker(draws2, sigma=sig)  # BUCKET
    print(str(num))
#    boot_array[num, 0] = perc1
#    boot_array[num, 1] = perc2
print(len(draws[0][0]), len(draws[0]), len(draws), 'draw')  # 1000, 22, N*len(eelgs)
print(len(draws2[0][0]), len(draws2[0]), len(draws2), 'draw2')  # 1000, 22, N*len(lbgs)
plot_boots([draws, draws2], t1[0], n=N)

'''
# For N~100, loops take ~13 mins; rest takes ~8 mins

# Note: for ~100 loops, the stacker step takes a minute or two
perc1 = stack.stacker(draws, sigma=sig)  # stacker also prints (22, 1000*(number of galaxies=N*len(eelgs)))
perc2 = stack.stacker(draws2, sigma=sig)  # stacker also prints (22, 1000*(number of galaxies=N*len(lbgs)))
smooth_percs = [stack.smooth(perc1), stack.smooth(perc2)]  # len(gal_draws[0])=22=len(perc)=22, len(perc[0])=3

plot_boots(smooth_percs)
'''

'''
# full_array = np.zeros(shape=(2, N))  # len(full_array) = 2
for i in range(len(boot_array)):  # = len(N)
    new_perc1 = stack.stacker(boot_array[:, 0])
    new_perc2 = stack.stacker(boot_array[:, 1])
    # full_array[0, i] = boot_array
    smooth_percs = [stack.smooth(new_perc1), stack.smooth(new_perc2)]
    stack.plot_sfhs(smooth_percs, t1[0], sigma=sig, save=False)
# newfile.close()
'''
'''





# START STACKING
t1 = []
draws = []
boots = []
nummy = 0
c = 0
for glxy in eelgs:
    c += 1
    # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
    file = pkls + glxy + '_extra_out.pkl'
    if os.path.exists(file):
        nummy += 1
        temp = randraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
        # temp = bootdraw(file)  # temp[0] lists the num=1000 random posterior samples; temp[1] = time vector
        draws.append(temp[0])
        # boots.append(bootstrap(temp[0]))
        t1.append(temp[1])
    else:
        print(file)

# stacker2(draws, t1)
sig = 1  # what sigma error to show on plot
perc1 = stacker(draws, sigma=sig)

draws2 = []
numl = 0
cl = 0
# t2 = []
for glxy in lbgs:
    cl += 1
    # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
    file = l_pkls + glxy + '_extra_out.pkl'
    if os.path.exists(file):
        numl += 1
        temp = randraw(file)
        # temp = bootdraw(file)
        draws2.append(temp[0])
        # t2.append(temp[1])
    else:
        print(file)

perc2 = stacker(draws2, sigma=sig)

# smooth_percs = perc1, perc2
print(nummy, c, 'nume')
print(numl, cl, 'numl')
smooth_percs = [smooth(perc1), smooth(perc2)]
plot_sfhs(smooth_percs, t1[0], sigma=sig)
'''
