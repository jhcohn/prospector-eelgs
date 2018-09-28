import pickle
import matplotlib.pyplot as plt
import numpy as np
import argparse
import random
import os
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import uvj
from astropy.cosmology import WMAP9
import scipy
from matplotlib.ticker import MultipleLocator, NullFormatter, MaxNLocator
from matplotlib import rc
from matplotlib import ticker
import matplotlib.colors as colors
import matplotlib.gridspec as gridspec
import scipy.stats as stats
import get_mass_dust as gmd


def draw_sfrs(objs, flds, new_logmass=None, ndraw=1e4, alpha_sfh=1.0, pfile=None, show=True):
    # pfile=e_params or n_params

    # let's do it
    ndraw = int(ndraw)
    zred = np.array([3.5])
    logmass = np.array([10.])
    smass_factor = 0.8  # add in a not-unreasonable stellar mass <--> total mass conversion
    minssfr, maxssfr = 1e-14, 1e-5

    # new redshift, new model
    if new_logmass is None:
        new_logmass = 10.17
    new_SFRs = []
    for i in range(len(objs)):
        new_model = pfile.load_model(objname=objs[i], field=flds[i])
        new_zred = new_model.params['zred']
        tuniv = WMAP9.age(new_zred).value
        new_SFRs.append(new_logmass / tuniv)

    return new_SFRs


def tuniv_ssfr_prior(objs, flds, pfile=None):
    # NEW TUNIV PRIOR
    ssfrs = []
    for i in range(len(objs)):
        model = pfile.load_model(objname=objs[i], field=flds[i])  # load_model prints zred, tuniv
        tuniv = WMAP9.age(model.params['zred']).value
        ssfrs.append(1 / tuniv)
        print(i)
    perc = np.percentile(ssfrs, [16., 50., 84.])
    print(perc)

    return perc


def bootstrap(X, n=None, X_err=None):
    """
    Bootstrap resample an array_like
    Parameters
    :param X: array-like data to resample
    :param n: int, optional length of resampled array, equal to len(X) if n == None
    :param X_err:
    :return: X_resamples
    """
    if n is None:
        n = len(X)

    resample_i = np.floor(np.random.rand(n) * len(X)).astype(int)
    # creates len(n) array filled with random numbers from floor([0 to 1) * len(X))
    # i.e. it's a len(X) array filled with random integers from 0 to len(X)

    # resample_i=np.random.randint(low=0, high=len(X)-1, size=len(X))
    # print(resample_i)
    X_resample = X[resample_i]  # take X and use indices chosen randomly above
    if X_err != None:
        X_err_resample = X_err[resample_i]
        return X_resample, X_err_resample
    else:
        return X_resample, resample_i


def randraw(infile, logmass, num=10**4):  # num=1000
    """
    For a given galaxy, randraw samples the posterior num times for each point in extra_output['extras']['sfh'][i]

    :param infile: ID_field_base_extra_out.py, where extra_output is stored using output.py
    :param num: number of times to sample the galaxy posterior at each point in the SFH
    :return: draw_from_sfh = 22 x num, lists the num random posterior samples, t_sfh = time vector associated with SFH
    """
    with open(infile, 'rb') as exout:
        extra_output = pickle.load(exout)

#    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['ssfr']), num))  # shape=(22, num)  # BUCKET MASS WRONG?
    draw_from_sfh = np.zeros(shape=(len(extra_output['extras']['sfh']), num))  # shape=(22, num)
    # print(len(extra_output['extras']['ssfr']), len(extra_output['extras']['ssfr'][0]))  # 22, 2000
    print(len(draw_from_sfh), len(draw_from_sfh[0]))  # 22, num

    '''  # BUCKET MASS WRONG?
    for i in range(len(extra_output['extras']['ssfr'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['ssfr'][i][random.randint(0, num)] / 0.8  # BUCKET 1/0.8 here?
    '''
    for i in range(len(extra_output['extras']['sfh'])):  # at each of these 22 points
        for j in range(num):  # randomly draw from the ssfr posterior num times
            draw_from_sfh[i][j] = extra_output['extras']['sfh'][i][random.randint(0, 1999)] / (10**logmass) / 0.8

    return draw_from_sfh, extra_output['extras']['t_sfh']


def smooth(perc):
    """
    Takes the stacked sfh that is output from stacker and averages the sfr values in each bin, such that the sfr within
    each bin is flat, as is the case in the original extra_output['extras']['sfh'] output

    :param perc: stored lists of the median and +/1 1sigma SFH values calculated from the gal_draws
    :return: smoother = same shape as perc, but with all points within a bin averaged s.t. all points within a bin=flat
    """
    # from perc: bin1 0:2, bin2 3:6, bin3 7:10, bin4 11:14, bin5 15:18, bin6 19:22
    smoother = np.zeros(shape=(len(perc), len(perc[0])))  # shape=(22, 3)
    for j in range(len(perc[0])):  # 3
        yax = perc[:, j]
        # print(yax)  # columns from perc: j=0 --> -1sigma, j=1 --> median, j=2 --> +1sigma
        for i in range(3):
            smoother[i][j] = (yax[0] + yax[1] + yax[2]) / 3
            smoother[i+19][j] = (yax[-1] + yax[-2] + yax[-3]) / 3
        for i in range(4):
            smoother[i+3][j] = (yax[3] + yax[4] + yax[5] + yax[6]) / 4
            smoother[i+7][j] = (yax[7] + yax[8] + yax[9] + yax[10]) / 4
            smoother[i+11][j] = (yax[11] + yax[12] + yax[13] + yax[14]) / 4
            smoother[i+15][j] = (yax[15] + yax[16] + yax[17] + yax[18]) / 4

    # print(smoother)
    return smoother


def stacker(gal_draws, sigma=1):
    """
    stacker takes input of random points drawn from a list of galaxies' SFH posteriors, concatenates them within each
    bin, and then calculates the median and 1 sigma errors in each bin
    gal_draws should be in format draws = [draw_from_sfh1, draw_from_sfh2, ...]
    each draw_from_sfh has shape=(22,num)

    :param gal_draws: list comprised of draw_from_sfh (i.e. comprised of the output from randraw) for a list of galaxies
    :param sigma: how many sigma of error we want to show in the plot
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

    perc = np.zeros(shape=(len(gal_draws[0]), 2*sigma + 1))  # len(gal_draws[0])=22=len(t); len(perc)=22, len(perc[0])=3
    for jj in xrange(len(gal_draws[0])):
        if sigma == 1:
            perc[jj, :] = np.percentile(all_draws[jj, :], [16.0, 50.0, 84.0])  # median, +/- 34% = +/- 1sigma
        elif sigma == 3:
            perc[jj, :] = np.percentile(all_draws[jj, :], [0.3, 2.4, 16.0, 50.0, 84.0, 97.6, 99.7])  # out to 3 sigma

    return all_draws


# Define a function to make the ellipses
def ellipse(ra, rb, ang, x0, y0, Nb=100):
    xpos, ypos = x0, y0
    radm, radn = ra, rb
    an = ang
    co, si = np.cos(an), np.sin(an)
    the = np.linspace(0, 2 * np.pi, Nb)
    X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
    Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
    return X, Y


def percellipse(percs_a, percs_b, ang, Nb=100):
    xpos, ypos = percs_a[1], percs_b[1]
    radm, radn = percs_a[2] - percs_a[0], percs_b[2] - percs_b[0]
    an = ang
    co, si = np.cos(an), np.sin(an)
    the = np.linspace(0, 2 * np.pi, Nb)
    X = radm * np.cos(the) * co - si * radn * np.sin(the) + xpos
    Y = radm * np.cos(the) * si + co * radn * np.sin(the) + ypos
    return X, Y


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def density_estimation(m1, m2, xs=[-1.5,2.5], ys=[-1,2.5], num=100):  # 100j
    # X, Y = np.mgrid[xs[0]:xs[1]:num, ys[0]:ys[1]:num]
    X = np.logspace(xs[0], xs[1], num)
    Y = np.logspace(ys[0], ys[1], num)
    X, Y = np.meshgrid(X, Y)
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = stats.gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)
    return X, Y, Z


def simpler(recentx, secondy, logit=False, later=False, n_samps=10**4):
    if logit:
        x = [np.log10(rx) for rx in recentx[0]]
        y = [np.log10(sy) for sy in secondy[0]]
        x2 = [np.log10(rx2) for rx2 in recentx[1]]
        y2 = [np.log10(sy2) for sy2 in secondy[1]]

        xlims = [-2, 1.5]  # [10**-3, 13.] #[0., 13.]  # 17.
        ylims = [-2, 1.5]  # [10**-3, 13.]#[0., 8.]#13.]  # 17.
    else:
        x = recentx[0]
        y = secondy[0]
        x2 = recentx[1]
        y2 = secondy[1]

        xlims = [6*10**-2, 2*10**1.]  # 4*10**-2, 3*10**1.
        ylims = [6*10**-2, 2*10**1.]  # 4*10**-2, 3*10**1.

    # start with a Figure
    fig1 = plt.figure(1, figsize=(12, 12))
    cb = False
    '''
    if cb:
        gs = gridspec.GridSpec(1, 3, width_ratios=[12, 1, 1])
        ax1 = plt.subplot(gs[0])  # plt.subplot(111)
    else:
        ax1 = plt.subplot(111)
    '''
    if later:
        gs = gridspec.GridSpec(5, 5)
        gs.update(wspace=0, hspace=0)
        ax1 = plt.subplot(gs[1:, :-1])
        axHisty = plt.subplot(gs[1:, -1])
        axHistx = plt.subplot(gs[0, :-1])
        axHistx.set_xlim(xlims[0], xlims[1])
        axHisty.set_ylim(ylims[0], ylims[1])
        axHistx.set_xscale('log')
        axHisty.set_yscale('log')
        print('hi later')
    else:
        #ax1 = plt.subplot(121)
        #ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
        ax1 = plt.subplot(111)

    print(xlims, ylims)
    ax1.set_xlim(xlims[0], xlims[1])
    ax1.set_ylim(ylims[0], ylims[1])
    percx2 = np.percentile(x2, [16., 50., 84.])
    percy2 = np.percentile(y2, [16., 50., 84.])

    print(percx2, percy2)
    print(np.percentile(x, [16., 50., 84.]), np.percentile(y, [16., 50., 84.]))

    ang = 45
    xcenter = np.median(x2)
    ycenter = np.median(y2)
    rb = (percx2[2] - percx2[0])/2 # np.std(x2)
    ra = (percy2[2] - percy2[0])/2 # np.std(y2)
    '''
    X, Y = ellipse(ra, rb, ang, xcenter, ycenter)
    ax1.plot(X, Y, "b--", ms=1, linewidth=2.0)
    '''
    perc2sigx2 = np.percentile(x2, [2.5, 50., 97.5])

    not_diff = 0
    in_ell = [0, 0, 0]  # 1, 2, 3 sigma
    sames = 0
    for i in range(len(x)):
        if percx2[0] < x[i] < percx2[2]:  # and percy2[0] < y[i] < percy2[2]:
            not_diff += 1.
        if ((xcenter - x[i])**2 / rb**2) + ((ycenter - y[i])**2 / ra**2) <= 1.:
            in_ell[0] += 1.
        if perc2sigx2[0] < x[i] < perc2sigx2[2]:
            sames += 1.
    print(not_diff, not_diff / len(x))  # consistently ~24% to 25% are not different --> 75% are distinct
    print(in_ell)
    print(sames, sames/len(x))  # 15978, 84%

    # ax1.axvline(x=perc2sigx2[2], color='k')
    # ax1.axvline(x=percx2[2], color='k')

    # get colormaps, plot 2dhists
    cmap = plt.get_cmap('Purples')
    new_cmap = truncate_colormap(cmap, 0.4, 1.5)
    cmap2 = plt.get_cmap('Blues')
    new_cmap2 = truncate_colormap(cmap2, 0.4, 1.5)
    wantlog = True
    if wantlog:
        xbins = np.logspace(-2, 2, 100)  #10**np.linspace(-2, 2, 1000)
        ybins = np.logspace(-2, 2, 100)  # 10**np.linspace(-2, 2, 1000)
        xbins2 = np.logspace(-2, 2, 100*int(167/19))  # 10**np.linspace(-2, 2, 1000)
        ybins2 = np.logspace(-2, 2, 100*int(167/19))  # 10**np.linspace(-2, 2, 1000)

        # x_bins = np.logspace(np.log10(min(x)), np.log10(max(x)), np.sqrt(100))
        # y_bins = np.logspace(np.log10(min(y)), np.log10(max(y)), np.sqrt(100))
        x_bins = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), np.sqrt(len(x))/2)
        y_bins = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), np.sqrt(len(x))/2)
        H1, xedges, yedges = np.histogram2d(x, y, bins=[x_bins, y_bins])
        H1 = np.ma.masked_array(H1, H1 < len(x) / (1.5*n_samps))  # 20.
        # H1 /= 19
        # ax1.pcolormesh(xedges, yedges, H1.T, cmap='Purples')

        # x_bins2 = np.logspace(np.log10(min(x2)), np.log10(max(x2)), np.sqrt(100))
        # y_bins2 = np.logspace(np.log10(min(y2)), np.log10(max(y2)), np.sqrt(100))
        x_bins2 = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), np.sqrt(len(x))/2)
        y_bins2 = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), np.sqrt(len(x))/2)
        H2, xedges2, yedges2 = np.histogram2d(x2, y2, bins=[x_bins2, y_bins2])
        H2 = np.ma.masked_array(H2, H2 < len(x2) / (2*n_samps))  # 20.
        # H2 /= 167.

        im2 = ax1.pcolormesh(xedges2, yedges2, H2.T, cmap=new_cmap2)
        im1 = ax1.pcolormesh(xedges, yedges, H1.T, cmap=new_cmap)#, vmin=2, vmax=100)

        if cb:
            ax2 = plt.subplot(gs[1])
            ax3 = plt.subplot(gs[2])
            cbar1 = plt.colorbar(im1, cax=ax2)
            cbar2 = plt.colorbar(im2, cax=ax3)

            cbar1.set_ticks([20., 40., 60., 80., 100.])
            cbar2.set_ticks([200., 250., 300., 350., 400.])
            cbar1.set_ticklabels([r'20', r'40', r'60', r'80', r'100'])
            cbar2.set_ticklabels([r'200', r'250', r'300', r'350', r'400'])
            cbar1.ax.tick_params(labelsize=20)
            cbar2.ax.tick_params(labelsize=20)

            plt.subplots_adjust(wspace=0.1)  # , hspace=0.)

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax2.set_xscale('log')
        #ax2.set_yscale('log')
    elif logit:
        print("hi")
        X2, Y2, Z2 = density_estimation(x2, y2, xs=[-2., 2.], ys=[-2., 2.], num=10)
        levels2 = np.arange(0.5, np.amax(Z2), 0.02) + 0.02
        print("hi")
        plt.contourf(X2, Y2, Z2, locator=ticker.LogLocator(), levels=levels2, cmap='Blues', alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]
        X, Y, Z = density_estimation(x, y, xs=[-2., 2.], ys=[-2., 2.], num=10)
        levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
        plt.contourf(X, Y, Z, locator=ticker.LogLocator(), levels=levels, cmap='Purples', alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]
        plt.plot(np.percentile(x2, [10., 50., 90.])[1], np.percentile(y2, [10., 50., 90.])[1], 'k*')
        plt.plot(np.percentile(x, [10., 50., 90.])[1], np.percentile(y, [10., 50., 90.])[1], 'ko')
        '''
        xbins = np.logspace(-2, 1.5, 100)  #10**np.linspace(-2, 2, 1000)
        ybins = np.logspace(-2, 1.5, 100)  # 10**np.linspace(-2, 2, 1000)
        xbins2 = np.logspace(-2, 1.5, 50*127/19)  #10**np.linspace(-2, 2, 1000)
        ybins2 = np.logspace(-2, 1.5, 50*127/19)  # 10**np.linspace(-2, 2, 1000)
        ax1.hist2d(x2, y2, bins=[xbins2, ybins2], cmap=new_cmap2, cmin=5)  #, legend='SFGs')  # *(127/19), ..., 19
        ax1.hist2d(x, y, bins=[xbins, ybins], cmap=new_cmap, cmin=19)  # , legend='EELGs')  # 19
        '''
    else:
        ax1.hist2d(x2, y2, bins=100*167/19, cmap=new_cmap2, cmin=19)  #, legend='SFGs')  # *(127/19), ..., 19
        # bins with counts < cmin won't be displayed! cmin=len(x2)/1000
        ax1.hist2d(x, y, bins=100, cmap=new_cmap, cmin=19)#, legend='EELGs')  # 19
        # bins with counts < cmin won't be displayed! cmin=len(x)/1000

    ax1.plot(xlims, ylims, ls='--', color='k')#, lw=2)  # [0., 17.]
    #ax2.plot(xlims, ylims, ls='--', color='k')  # [0., 17.]

    # Set up your x and y labels
    if later:
        xlabel = r'$<$SFR$_{50-100}>$/M$_{\rm tot}$ [Gyr$^{-1}$]'
        # r'$<$SSFR$>_{0-50}$ [Gyr$^{-1}$]'  # r'SSFR, most recent bin [Gyr$^{-1}$]'
        ylabel = r'$<$SFR$_{100-1000}>$/M$_{\rm tot}$ [Gyr$^{-1}$]'
        # r'$<$SSFR$>_{50-100}$ [Gyr$^{-1}$]'  # r'SSFR, second most recent bin [Gyr$^{-1}$]'
        '''
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = left_h = left + width + 0.02

        # Set up the geometry of the three plots
        rect_temperature = [left, bottom, width, height]  # dimensions of temp plot
        rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
        rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram
        # Make the three plots
        axTemperature = plt.axes(rect_temperature)  # temperature plot
        axHistx = plt.axes(rect_histx)  # x histogram
        axHisty = plt.axes(rect_histy)  # y histogram
        '''
        print('hi later')
        from matplotlib.ticker import NullFormatter, MaxNLocator
        nullfmt = NullFormatter()
        axHistx.xaxis.set_major_formatter(nullfmt)
        axHistx.yaxis.set_major_formatter(nullfmt)
        axHisty.xaxis.set_major_formatter(nullfmt)
        axHisty.yaxis.set_major_formatter(nullfmt)

        H1x = H1.sum(axis=-1)
        H1y = H1.sum(axis=0)

        H2x = H2.sum(axis=-1)
        H2y = H2.sum(axis=0)

        print(np.amax(H1x))
        print(len(H1x))
        print(len(x_bins))
        x_binsx = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), len(H1x))
        x_binsx2 = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), len(H2x))
        y_binsy = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), len(H1y))
        y_binsy2 = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), len(H2y))
        print(len(x_binsx), len(x_binsx2))
        print(H1x)

        # axHistx.hist(H1x, bins=x_bins, color='purple', normed=True, alpha=0.7)
        axHistx.plot(x_binsx, H1x / np.amax(H1x), color='purple', lw=2)
        # axHistx.plot(range(len(H2x)), H2x / np.amax(H2x), color='b')
        axHistx.plot(x_binsx2, H2x / np.amax(H2x), 'b--', lw=2)
        axHisty.plot(H1y / np.amax(H1y), y_binsy, color='purple', lw=2)
        axHisty.plot(H2y / np.amax(H2y), y_binsy2, 'b--', lw=2)
        # axHisty.hist(H2.sum(axis=-1), bins=y_bins2, orientation='horizontal', facecolor='none', lw=1.5, edgecolor='b', hatch='/',
                     # normed=True)  # , range=[min(y2) * 10., max(y2)]

        # Set up the histogram limits
        axHisty.xaxis.set_major_locator(MaxNLocator(4))
        axHistx.yaxis.set_major_locator(MaxNLocator(4))
        # plt.draw()
    else:
        xlabel = r'$<$SFR$_{0-50}>$/M$_{\rm tot}$ [Gyr$^{-1}$]'
        # r'$<$SSFR$>_{0-50}$ [Gyr$^{-1}$]'  # r'SSFR, most recent bin [Gyr$^{-1}$]'
        ylabel = r'$<$SFR$_{50-100}>$/M$_{\rm tot}$ [Gyr$^{-1}$]'
        # r'$<$SSFR$>_{50-100}$ [Gyr$^{-1}$]'  # r'SSFR, second most recent bin [Gyr$^{-1}$]'
    ax1.set_xlabel(xlabel, fontsize=30)
    ax1.set_ylabel(ylabel, fontsize=30)
    #ax2.set_xlabel(xlabel, fontsize=30)

    # Make the tickmarks pretty
    fs = 20
    fs_ticks = 25
    log = True
    if logit:
        xtick = [-3, -2, -1, 0, 1, 2, 3]
        ytick = [-3, -2, -1, 0, 1, 2, 3]

        ax1.set_xticklabels([r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'],
                            size=fs_ticks)
        ax1.set_yticklabels([r'$10^{-3}$', r'$10^{-2}$', r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'],
                            size=fs_ticks)
    elif log:
        ax1.set_xscale('log')
        # ax1.set_xscale('log')
        xtick = [10**-1, 10**0, 10**1]  # 10**-2
        ytick = [10 ** -1, 10 ** 0, 10 ** 1]  # 10**-2
        ax1.set_xticklabels([r'$10^{-1}$', r'$10^0$', r'$10^1$'], size=fs_ticks)  # r'$10^{-2}$',
        ax1.set_yticklabels([r'$10^{-1}$', r'$10^0$', r'$10^1$'], size=fs_ticks)  # r'$10^{-2}$',
        ax1.tick_params(axis='x', which='major', pad=10)
    else:
        xtick = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.]  # , 15.]  # , 2e-8]  # used if log=0
        ytick = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.]  # , 15.]  # , 2e-8]  # used if log=0
        # ytick = [0., 1., 2., 3., 4., 5., 6., 7., 8.]  # , 10.]#, 15.]  # , 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
        ax1.set_xticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$', r'$9$', r'$10$',
                             r'$11$', r'$12$', r'$13$'], size=fs_ticks)  # , r'$15$'
        # ax1.set_yticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$'], size=fs_ticks)
        # r'$10$', r'$15$'
        ax1.set_xticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$', r'$9$', r'$10$',
                             r'$11$', r'$12$', r'$13$'], size=fs_ticks)  # , r'$15$'
    ax1.xaxis.set_ticks(xtick)
    ax1.yaxis.set_ticks(ytick)
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)

    from matplotlib.patches import Patch
    legend_elements = [Patch(edgecolor='purple', facecolor='purple', label='EELGs'),
                       Patch(edgecolor='b', facecolor='b', label='SFGs')]
    ax1.legend(handles=legend_elements, loc='upper left', fontsize=30)

    #plt.axes().set_aspect('equal')
    plt.show()


def simpler_mass(recentx, secondy, logit=False, later=False):
    if logit:
        x = [np.log10(rx) for rx in recentx[0]]
        y = [np.log10(sy) for sy in secondy[0]]
        x2 = [np.log10(rx2) for rx2 in recentx[1]]
        y2 = [np.log10(sy2) for sy2 in secondy[1]]

        xlims = [-2, 1.5]  # [10**-3, 13.] #[0., 13.]  # 17.
        ylims = [-2, 1.5]  # [10**-3, 13.]#[0., 8.]#13.]  # 17.
    else:
        x = recentx[0]
        y = secondy[0]
        x2 = recentx[1]
        y2 = secondy[1]

        xlims = [3*10**-3, 1.]  #[0., 13.]  # 17.
        ylims = [3*10**-3, 1.]  #[0., 8.]#13.]  # 17.

    # start with a rectangular Figure
    fig1 = plt.figure(1, figsize=(12,12))#(16, 16)) #(19.5, 12))
    if later:
        gs = gridspec.GridSpec(3, 3)
        ax1 = plt.subplot(gs[1:, :-1])
        axHisty = plt.subplot(gs[1:, -1])
        axHistx = plt.subplot(gs[0, :-1])
        print('hi later')
    else:
        #ax1 = plt.subplot(121)
        #ax2 = plt.subplot(122, sharex=ax1, sharey=ax1)
        ax1 = plt.subplot(111)
    #plt.subplots_adjust(wspace=0., hspace=0.)

    print(xlims, ylims)
    ax1.set_xlim(xlims[0], xlims[1])
    ax1.set_ylim(ylims[0], ylims[1])
    #ax2.set_xlim(xlims[0], xlims[1])
    #ax2.set_ylim(ylims[0], ylims[1])
    percx2 = np.percentile(x2, [16., 50., 84.])
    percy2 = np.percentile(y2, [16., 50., 84.])

    print(percx2, percy2)
    # (array([ 0.33044154,  0.94145865,  2.54337235]), array([ 0.40491095,  1.14943934,  5.03257121]))
    print(np.percentile(x, [16., 50., 84.]), np.percentile(y, [16., 50., 84.]))

    ang = 45
    xcenter = np.median(x2)
    ycenter = np.median(y2)
    rb = (percx2[2] - percx2[0])/2 # np.std(x2)
    ra = (percy2[2] - percy2[0])/2 # np.std(y2)
    '''
    X, Y = ellipse(ra, rb, ang, xcenter, ycenter)
    ax1.plot(X, Y, "b--", ms=1, linewidth=2.0)
    '''
    perc2sigx2 = np.percentile(x2, [2.5, 50., 97.5])

    not_diff = 0
    in_ell = [0, 0, 0]  # 1, 2, 3 sigma
    sames = 0
    for i in range(len(x)):
        if percx2[0] < x[i] < percx2[2]:  # and percy2[0] < y[i] < percy2[2]:
            not_diff += 1.
        if ((xcenter - x[i])**2 / rb**2) + ((ycenter - y[i])**2 / ra**2) <= 1.:
            in_ell[0] += 1.
        if perc2sigx2[0] < x[i] < perc2sigx2[2]:
            sames += 1.
    print(not_diff, not_diff / len(x))  # consistently ~24% to 25% are not different --> 75% are distinct
    print(in_ell)
    print(sames, sames/len(x))  # 15978, 84%

    # ax1.axvline(x=perc2sigx2[2], color='k')
    # ax1.axvline(x=percx2[2], color='k')

    # get colormaps, plot 2dhists
    cmap = plt.get_cmap('Purples')
    new_cmap = truncate_colormap(cmap, 0.4, 1.5)
    cmap2 = plt.get_cmap('Blues')
    new_cmap2 = truncate_colormap(cmap2, 0.4, 1.5)
    wantlog = True
    if wantlog:
        xbins = np.logspace(-2, 2, 100)  #10**np.linspace(-2, 2, 1000)
        ybins = np.logspace(-2, 2, 100)  # 10**np.linspace(-2, 2, 1000)
        xbins2 = np.logspace(-2, 2, 100*int(167/19))  # 10**np.linspace(-2, 2, 1000)
        ybins2 = np.logspace(-2, 2, 100*int(167/19))  # 10**np.linspace(-2, 2, 1000)

        # x_bins = np.logspace(np.log10(min(x)), np.log10(max(x)), np.sqrt(100))
        # y_bins = np.logspace(np.log10(min(y)), np.log10(max(y)), np.sqrt(100))
        x_bins = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), np.sqrt(len(x))/2)
        y_bins = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), np.sqrt(len(x))/2)
        H1, xedges, yedges = np.histogram2d(x, y, bins=[x_bins, y_bins])
        H1 = np.ma.masked_array(H1, H1 < 20.)
        # ax1.pcolormesh(xedges, yedges, H1.T, cmap='Purples')

        # x_bins2 = np.logspace(np.log10(min(x2)), np.log10(max(x2)), np.sqrt(100))
        # y_bins2 = np.logspace(np.log10(min(y2)), np.log10(max(y2)), np.sqrt(100))
        x_bins2 = np.logspace(np.log10(xlims[0]), np.log10(xlims[1]), np.sqrt(len(x))/2)
        y_bins2 = np.logspace(np.log10(ylims[0]), np.log10(ylims[1]), np.sqrt(len(x))/2)
        H2, xedges2, yedges2 = np.histogram2d(x2, y2, bins=[x_bins2, y_bins2])
        H2 = np.ma.masked_array(H2, H2 < 30.)
        ax1.pcolormesh(xedges2, yedges2, H2.T, cmap=new_cmap2)
        ax1.pcolormesh(xedges, yedges, H1.T, cmap=new_cmap)#, vmin=2, vmax=100)

        ax1.set_xscale('log')
        ax1.set_yscale('log')
        #ax2.set_xscale('log')
        #ax2.set_yscale('log')

    ax1.plot(xlims, ylims, ls='--', color='k')  # [0., 17.]
    #ax2.plot(xlims, ylims, ls='--', color='k')  # [0., 17.]

    # Set up your x and y labels
    if later:
        xlabel = r'M$_{50-100}$/M$_{\rm tot}$'
        ylabel = r'M$_{100-1000}$/M$_{\rm tot}$'

        '''
        left, width = 0.12, 0.55
        bottom, height = 0.12, 0.55
        bottom_h = left_h = left + width + 0.02

        # Set up the geometry of the three plots
        rect_temperature = [left, bottom, width, height]  # dimensions of temp plot
        rect_histx = [left, bottom_h, width, 0.25]  # dimensions of x-histogram
        rect_histy = [left_h, bottom, 0.25, height]  # dimensions of y-histogram
        # Make the three plots
        axTemperature = plt.axes(rect_temperature)  # temperature plot
        axHistx = plt.axes(rect_histx)  # x histogram
        axHisty = plt.axes(rect_histy)  # y histogram
        '''
        print('hi later')
        from matplotlib.ticker import NullFormatter, MaxNLocator
        nullfmt = NullFormatter()
        #axHistx.xaxis.set_major_formatter(nullfmt)
        #axHisty.yaxis.set_major_formatter(nullfmt)

        axHistx.hist(x, bins=x_bins, color='purple')
        axHistx.hist(x, bins=x_bins2, color='b')
        axHisty.hist(y, bins=y_bins, orientation='horizontal', color='purple')
        axHisty.hist(y, bins=y_bins2, orientation='horizontal', color='b')

        # Set up the histogram limits
        axHistx.set_xlim(min(x), max(x))
        axHisty.set_ylim(min(y), max(y))
        #axHisty.xaxis.set_major_locator(MaxNLocator(4))
        #axHistx.yaxis.set_major_locator(MaxNLocator(4))
        #plt.draw()

    else:
        xlabel = r'M$_{0-50}$/M$_{\rm tot}$'
        # r'$<$SSFR$>_{0-50}$ [Gyr$^{-1}$]'  # r'SSFR, most recent bin [Gyr$^{-1}$]'
        ylabel = r'M$_{50-100}$/M$_{\rm tot}$'
    # r'$<$SSFR$>_{50-100}$ [Gyr$^{-1}$]'  # r'SSFR, second most recent bin [Gyr$^{-1}$]'
    ax1.set_xlabel(xlabel, fontsize=30)
    ax1.set_ylabel(ylabel, fontsize=30)
    #ax2.set_xlabel(xlabel, fontsize=30)

    # Make the tickmarks pretty
    fs = 20
    fs_ticks = 25
    log = True
    if log:
        #plt.xscale('log')
        #plt.yscale('log')
        #plt.gca().set_yscale('log')
        #plt.gca().set_xscale('log')
        ax1.set_xscale('log')
        ax1.set_xscale('log')
        #ax2.set_xscale('log')
        #ax2.set_yscale('log')
        xtick = [10**-2, 10**-1, 10**0]  # 10**-2
        ytick = [10 ** -2, 10 ** -1, 10 ** 0]  # 10**-2
        ax1.set_xticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)  # r'$10^{-2}$',
        ax1.set_yticklabels([r'$10^{-2}$', r'$10^{-1}$', r'$10^0$'], size=fs_ticks)  # r'$10^{-2}$',
        ax1.tick_params(axis='x', which='major', pad=10)
        #ax2.set_xticklabels([r'$10^{-1}$', r'$10^0$', r'$10^1$'], size=fs_ticks)  # r'$10^{-2}$',
    else:
        xtick = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.]  # , 15.]  # , 2e-8]  # used if log=0
        ytick = [0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.]  # , 15.]  # , 2e-8]  # used if log=0
        # ytick = [0., 1., 2., 3., 4., 5., 6., 7., 8.]  # , 10.]#, 15.]  # , 0.45, 0.5, 0.6, 0.7, 0.8, 0.9]
        ax1.set_xticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$', r'$9$', r'$10$',
                             r'$11$', r'$12$', r'$13$'], size=fs_ticks)  # , r'$15$'
        # ax1.set_yticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$'], size=fs_ticks)
        # r'$10$', r'$15$'
        ax1.set_xticklabels([r'$0$', r'$1$', r'$2$', r'$3$', r'$4$', r'$5$', r'$6$', r'$7$', r'$8$', r'$9$', r'$10$',
                             r'$11$', r'$12$', r'$13$'], size=fs_ticks)  # , r'$15$'
    ax1.xaxis.set_ticks(xtick)
    ax1.yaxis.set_ticks(ytick)
    #ax2.xaxis.set_ticks(xtick)
    #ax2.yaxis.set_ticks(ytick)
    ax1.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    ax1.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    #ax2.tick_params('x', length=3, width=1, which='both', labelsize=fs)
    #ax2.tick_params('y', length=3, width=0.5, which='both', labelsize=fs)
    #plt.setp(ax2.get_yticklabels(), visible=False)
    plt.rc('xtick', labelsize=20)
    plt.rc('ytick', labelsize=20)

    from matplotlib.patches import Patch
    legend_elements = [Patch(edgecolor='purple', facecolor='purple', label='EELGs'),
                       Patch(edgecolor='b', facecolor='b', label='SFGs')]
    ax1.legend(handles=legend_elements, loc='lower left', fontsize=30)

    #plt.axes().set_aspect('equal')
    # plt.savefig('/home/jonathan/newfigs/butts.png',dpi=200)
    plt.show()


if __name__ == "__main__":

    boot = 0

    corr = 0
    fico = 1
    newsfg = 0

    if corr:
        base = ['corr', 'corr']
        folders = ['pkl_ecorr/', 'pkl_ncorr/']
        mass = [9.48, 10.12]
        import eelg_varymet_params as e_params
        import eelg_varymet_params as n_params
        normal = True
    elif fico:
        base = ['fico', 'fico']
        folders = ['pkl_efico/', 'pkl_nfico/']
        mass = [9.42, 10.13]
        import eelg_fifty_params as e_params
        import eelg_fifty_params as n_params
        normal = True
    elif newsfg:
        base = ['fico', 'newsfg']
        folders = ['pkl_efico/', 'pkl_nnewsfg/']
        mass = [9.42, 9.84]
        import eelg_fifty_params as e_params
        import eelg_fifty_params as n_params
        normal = False

    import stellar_ages as sa
    e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base, objlists=True, normal=normal)
    eelgs, lbgs = sa.get_gal_lists(base, objlists=False, normal=normal)

    pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[0]
    l_pkls = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/' + folders[1]
    out = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + 'out_efico/'
    l_out = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/out/' + 'out_nfico/'

    # START STACKING
    n_samps = 10**3  # 10**4

    t1 = []
    draws = []
    boots = []
    nummy = 0
    c = 0

    draws2 = []
    numl = 0
    cl = 0
    # t2 = []

    eelgs1 = []
    for file in os.listdir(out):
        if file.endswith(".h5"):
            eelgs1.append(file)
    sfgs1 = []
    for file in os.listdir(l_out):
        if file.endswith(".h5"):
            sfgs1.append(file)
    for glxy in eelgs:
        c += 1
        print(c)
        # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        file = pkls + glxy + '_extra_out.pkl'
        if os.path.exists(file):
            for mfile in eelgs1:
                if mfile.startswith(glxy):
                    use_mfile = mfile
            get_e = gmd.printer(out + use_mfile, percs=False, quiet=True)
            nummy += 1
            # temp = randraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            temp = randraw(file, get_e[np.random.randint(len(get_e))], num=n_samps)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            # temp = bootdraw(file)  # temp[0] lists num=1000 random posterior samples; temp[1] = time vector
            draws.append(temp[0])  # append random draw of ssfr
            # boots.append(bootstrap(temp[0]))
            t1.append(temp[1])
            print(len(temp[0]))  # 22, 1000
        else:
            print(file)

    lfile_loc = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfg_ssfrs'
    for glxy in lbgs:
        cl += 1
        print(cl)
        # file = glxy[0] + '_' + glxy[1] + '_' + glxy[2] + '_extra_out.pkl'
        file = l_pkls + glxy + '_extra_out.pkl'

        if os.path.exists(file):
            for mfile in sfgs1:
                if mfile.startswith(glxy):
                    use_mfile = mfile
            get_l = gmd.printer(l_out + use_mfile, percs=False, quiet=True)
            numl += 1
            # temp = randraw(file)
            temp = randraw(file, get_l[np.random.randint(len(get_l))], num=n_samps)
            # temp = bootdraw(file)
            draws2.append(temp[0])
            # t2.append(temp[1])
        else:
            print(file)

    # stacker2(draws, t1)
    sig = 1  # what sigma error to show on plot
    all1 = stacker(draws, sigma=sig)
    all2 = stacker(draws2, sigma=sig)

    means1 = np.zeros(shape=(n_samps))
    means2 = np.zeros(shape=(n_samps))
    # len(gal_draws) = number of galaxies in stack; len(gal_draws[0]) = 22, len(gal_draws[0][0]) = num (1000)
    # print(len(draws), len(draws[0]), len(draws[0][0]))  # 19, 22, 1000
    for nu in range(n_samps):
        each = []
        each2 = []
        for gal in range(len(draws)):
            each.append(draws[gal][0][nu])
        for gal2 in range(len(draws2)):
            each2.append(draws2[gal2][0][nu])
        means1[nu] = np.mean(each)
        means2[nu] = np.mean(each2)
    print('errors on mean:')
    print(np.percentile(means1, [16., 50., 84.]))
    print(np.percentile(means2, [16., 50., 84.]))

    rec1 = []
    sec1 = []
    thi1 = []
    mr1 = []
    ms1 = []
    mt1 = []
    for i in range(len(all1[0])):
        recent = []
        second = []
        third = []
        for j in (0, 1, 2):
            recent.append(all1[j][i])
        for k in (3, 4, 5, 6):
            second.append(all1[k][i])
        for l in (7, 8, 9, 10):
            third.append(all2[l][i])
        rec1.append(sum(recent) * 10**9 / 3)
        sec1.append(sum(second) * 10**9 / 4)
        thi1.append(sum(third) * 10**9 / 4)
        mr1.append(sum(recent) * 5*10**7 / 3)
        ms1.append(sum(second) * 5 * 10 ** 7 / 4)
        mt1.append(sum(third) * 5 * 10 ** 7 / 4)
    rec2 = []
    sec2 = []
    thi2 = []
    mr2 = []
    ms2 = []
    mt2 = []
    for i in range(len(all2[0])):
        recent = []
        second = []
        third = []
        for j in (0, 1, 2):
            recent.append(all2[j][i])
        for k in (3, 4, 5, 6):
            second.append(all2[k][i])
        for l in (7, 8, 9, 10):
            third.append(all2[l][i])
        rec2.append(sum(recent) * 10**9 / 3)  # convert Gyr to yr, and average all three samples within the bin
        sec2.append(sum(second) * 10**9 / 4)  # convert Gyr to yr, and average all four samples within the bin
        thi2.append(sum(third) * 10**9 / 4)
        mr2.append(sum(recent) * 5*10**7 / 3)
        ms2.append(sum(second) * 5 * 10 ** 7 / 4)
        mt2.append(sum(third) * 5 * 10 ** 7 / 4)

    randome1 = []
    randome2 = []
    idx = 4
    for i in range(len(draws[0][0])):
        recent = []
        for j in (0, 1, 2):
            recent.append(draws[idx][j][i])
        randome1.append(sum(recent) * 10**9 / 3)
    for i in range(len(draws[0][0])):
        second = []
        for k in (3, 4, 5, 6):
            second.append(draws[idx][k][i])
        randome2.append(sum(second) * 10**9 / 4)

    randoms1 = []
    randoms2 = []
    idx2 = 0
    for i in range(len(draws2[0][0])):
        recent = []
        for j in (0, 1, 2):
            recent.append(draws2[idx2][j][i])
        randoms1.append(sum(recent) * 10**9 / 3)
    for i in range(len(draws2[0][0])):
        second = []
        for k in (3, 4, 5, 6):
            second.append(draws2[idx2][k][i])
        randoms2.append(sum(second) * 10**9 / 4)

    # twodplot(randoms1, randoms2, 'b')
    # twodplot(randome1, randome2, 'purple')

    # smooth_percs = perc1, perc2
    print(nummy, c, 'nume')
    print(numl, cl, 'numl')
    # smooth_percs = [smooth(perc1), smooth(perc2)]

    rc('font', **{'family': 'serif', 'serif': ['Times']})
    rc('text', usetex=True)

    recents = [rec1, rec2]
    seconds = [sec1, sec2]
    thirds = [thi1, thi2]
    massrec = [mr1, mr2]
    masssec = [ms1, ms2]

    simpler(seconds, thirds, later=True, n_samps=n_samps)

    simpler(recents, seconds, n_samps=n_samps)
    # simpler(rec2, sec2, col='b')

    simpler_mass(massrec, masssec)
