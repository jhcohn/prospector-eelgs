# -*- coding: utf-8 -*-
"""
UVJ plotter
Fluxes from rest frame catalog
Use flag from main catalog
Photometric redshifts from zout catalog
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.colors as colors


def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap


def uvj_plot(objname, field, objlist=None, title=True, labels=True, lims=False, size=20, show=True, col='b',
             legend=None, hist=False):
    # UVJ plotter
    # Choose correct catalogs based on field name
    if objlist is not None:
        # print(type(col))
        # if type(col) == 'str':
        #     col = [col] * len(objlist)
        objs = []
        fs = []
        for obj in objlist:
            count = 0
            ob = ''
            f = ''
            for i in range(len(obj)):
                if obj[i] == '_':
                    count += 1
                elif count == 0:
                    ob += obj[i]
                elif count == 1:
                    f += obj[i]
                else:
                    pass
            objs.append(int(ob))
            fs.append(f)
        print('lens', len(objs), len(fs))

    if field == 'cdfs':
        rest = ['/home/jonathan/cdfs/cdfs.v1.6.9.rest.v0.9.cat']
        cat = ['/home/jonathan/cdfs/cdfs.v1.6.11.cat']
        zname = ['/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout']
    elif field == 'cosmos':
        rest = ['/home/jonathan/cosmos/cosmos.v1.3.6.rest.v0.9.cat']
        cat = ['/home/jonathan/cosmos/cosmos.v1.3.8.cat']  # main catalog
        zname = ['/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout']  # redshift catalog
    elif field == 'uds':
        rest = ['/home/jonathan/uds/uds.v1.5.8.rest.v0.9.cat']
        cat = ['/home/jonathan/uds/uds.v1.5.10.cat']
        zname = ['/home/jonathan/uds/uds.v1.5.8.zout']
    elif field == 'all':
        rest = ['/home/jonathan/cdfs/cdfs.v1.6.9.rest.v0.9.cat', '/home/jonathan/cosmos/cosmos.v1.3.6.rest.v0.9.cat',
                '/home/jonathan/uds/uds.v1.5.8.rest.v0.9.cat']
        cat = ['/home/jonathan/cdfs/cdfs.v1.6.11.cat', '/home/jonathan/cosmos/cosmos.v1.3.8.cat',
               '/home/jonathan/uds/uds.v1.5.10.cat']
        zname = ['/home/jonathan/cdfs/cdfs.v1.6.9.awk.zout', '/home/jonathan/cosmos/cosmos.v1.3.6.awk.zout',
                 '/home/jonathan/uds/uds.v1.5.8.awk.zout']
    fields = ['cdfs', 'cosmos', 'uds']
    special_uv, special_vj = [], []
    nummy = 0

    for j in range(len(rest)):
        # load catalogs
        fld = fields[j]
        s = np.loadtxt(rest[j])  # rest frame catalog
        main = np.loadtxt(cat[j])  # main catalog
        redshift = np.loadtxt(zname[j])  # redshift catalog

        # initialize everything we'll need
        objuse = []
        zcut = []
        UV = []
        VJ = []
        Ksnr = []
        table = []
        ax_uv = []
        ax_vj = []
        star = False  # for use in plotting a specific object as a blue star
        for i in range(len(main)):
            # USE FLAG AND Z_PHOT for each obj; note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec
            objuse.append(main[i][-4])  # index good for all main cats: last 4 elements are use, snr, use_nosnr, z_spec
            if main[i][-1] >= 0:  # if z_spec exists for this galaxy
                zcut.append(main[i][-1])  # index good for all main cats
            else:
                zcut.append(redshift[i][17])  # using zout catalog, z_peak = z_phot; index good for all zout cats

            # CREATE UVJ AXES
            UV.append(-2.5 * np.log10(s[i][11] / s[i][15]))  # indices good for all rest frame cats
            VJ.append(-2.5 * np.log10(s[i][15] / s[i][17]))  # indices good for all rest frame cats

            if objlist is not None:
                for j in range(len(objs)):
                    if fs[j] == fld and i == objs[j] - 1:
                        # print(fs[j], objs[j])  # showed LOTS (including COSMOS 1824!)
                        nummy += 1
                        special_uv.append(-2.5 * np.log10(s[i][11] / s[i][15]))
                        special_vj.append(-2.5 * np.log10(s[i][15] / s[i][17]))

            elif i == int(objname) - 1:
                star = True
                special_uv = -2.5 * np.log10(s[i][11] / s[i][15])
                special_vj = -2.5 * np.log10(s[i][15] / s[i][17])
                print(main[i][0], int(objname), special_vj, special_uv)
                print(zcut[i], main[i][-1], redshift[i][17])

            # SNR IN K BAND
            Ksnr.append(main[i][21] / main[i][22])  # f_Ksall / e_Ksall, indices good for all main cats

            # UVJ CORNER CUT, USE CUT, K-BAND SNR CUT, AND REDSHIFT CUT (zcut[i] < 2.369 if include F160W; else < 2.560)
            if objuse[i] == 1 and Ksnr[i] >= 10:  # 1 < zcut[i] < 2:  # Ksnr[i] >= 10:  # >=20
                ax_uv.append(UV[i])
                ax_vj.append(VJ[i])
                table.append(main[i][0])  # ID in catalog = main[i][0]

        # PLOT HIST OF REMAINING POINTS
        # plt.scatter(ax_vj, ax_uv, color='0.5', alpha=0.1, marker=".")  # x, y
        cmap = plt.get_cmap('binary')
        new_cmap = truncate_colormap(cmap, 0.0, 2.0)
        plt.hist2d(ax_vj, ax_uv, bins=100, cmap=new_cmap, cmin=1)  # 'binary')  # x, y

        '''
        if objlist is not None:
            for i in range(len(special_uv)):
                print(fld, i + 1)
                plt.scatter(special_uv[i], special_vj[i], color='b', marker="*", s=20)
        '''
        if star:
            plt.scatter(special_vj, special_uv, color=col, marker="*", s=100)

        if title:
            plt.title(field + '-' + objname)
        if labels:
            plt.text(-0.4, 1.35, 'Quiescent', fontsize=16)  # label quiescent region
            plt.text(-0.4, 1.1, 'Star-forming', fontsize=16)  # label star-forming region
        plt.plot([-0.5, 0.8], [1.3, 1.3], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2]) x2=0.9 or 0.8
        plt.plot([1.6, 1.6], [2.5, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2])
        plt.plot([0.8, 1.6], [1.3, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2]) x1=0.9 or 0.8

        if lims:
            plt.xlim(-1.5, 2.5)
            plt.xticks([-1, 0, 1, 2], size=size)
            plt.ylim(-1., 2.5)  # 3
            plt.yticks([-1., 0., 1., 2.], size=size)  # , 3.])

        plt.xlabel(r'$V - J$', fontsize=size)  # (Rest)
        plt.ylabel(r'$U - V$', fontsize=size)  # (Rest)
        # if show:
        #    if objlist is None:
        #        plt.show()

    if objlist is not None:
        print(nummy)
        print(len(special_vj), len(special_uv))
        if hist:
            cmap2 = plt.get_cmap('Blues')
            new_cmap2 = truncate_colormap(cmap2, 0.2, 2.0)
            plt.hist2d(special_vj, special_uv, bins=10, cmap=new_cmap2, cmin=2)  # 'binary')  # x, y
            # plt.hist2d(special_vj, special_uv, bins=10, cmap='Blues', cmin=2)#, alpha=0.75)  # x, y
            '''
            from scipy import stats

            def density_estimation(m1, m2):
                X, Y = np.mgrid[-1.5:2.5:100j, -1.:2.5:100j]
                positions = np.vstack([X.ravel(), Y.ravel()])
                values = np.vstack([m1, m2])
                kernel = stats.gaussian_kde(values)
                Z = np.reshape(kernel(positions).T, X.shape)
                return X, Y, Z

            X, Y, Z = density_estimation(special_vj, special_uv)
            levels = np.arange(0.5, np.amax(Z), 0.02) + 0.02
            plt.contourf(X, Y, Z, levels=levels, cmap='Blues', alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]
            # plt.contour(X, Y, Z, linewidths=3, cmap='Blues', alpha=0.5)  # levels=levels  # [0.1, 0.2, 0.5, 1., 25.]
            '''
            plt.xlim(-1.5, 2.5)
            plt.xticks([-1, 0, 1, 2], size=size)
            plt.ylim(-1., 2.5)  # 3
            plt.yticks([-1., 0., 1., 2.], size=size)  # , 3.])
        else:
            for i in range(len(special_vj)):
                if legend is None:
                    plt.scatter(special_vj[i], special_uv[i], color=col[i], marker="*", s=200)
                else:
                    plt.scatter(special_vj[i], special_uv[i], color=col[i], marker="*", s=200, label=legend[i])
                    plt.legend(scatterpoints=1, loc='upper left', prop={'size': 15})
    elif show:
        plt.show()

    '''
    # load catalogs
    s = np.loadtxt(rest)  # rest frame catalog
    main = np.loadtxt(cat)  # main catalog
    redshift = np.loadtxt(zname)  # redshift catalog

    # initialize everything we'll need
    objuse = []
    zcut = []
    UV = []
    VJ = []
    Ksnr = []
    table = []
    ax_uv = []
    ax_vj = []
    star = False  # for use in plotting a specific object as a blue star
    for i in range(len(main)):
        # USE FLAG AND Z_PHOT for each obj) (note: len(main[0] = 156; elements 153:155 = use, snr, use_nosnr, z_spec)
        objuse.append(main[i][-4])  # index good for all main cats: last 4 elements are use, snr, use_nosnr, z_spec
        if main[i][-1] >= 0:  # if z_spec exists for this galaxy
            zcut.append(main[i][-1])  # index good for all main cats
        else:
            zcut.append(redshift[i][17])  # using zout catalog, z_peak = z_phot; index good for all zout cats

        # CREATE UVJ AXES
        UV.append(-2.5 * np.log10(s[i][11] / s[i][15]))  # indices good for all rest frame cats
        VJ.append(-2.5 * np.log10(s[i][15] / s[i][17]))  # indices good for all rest frame cats
        if i == int(objname) - 1:
            star = True
            special_uv = -2.5 * np.log10(s[i][11] / s[i][15])
            special_vj = -2.5 * np.log10(s[i][15] / s[i][17])
            print(main[i][0], int(objname), special_vj, special_uv)
            print(zcut[i], main[i][-1], redshift[i][17])

        # SNR IN K BAND
        Ksnr.append(main[i][21] / main[i][22])  # f_Ksall / e_Ksall, indices good for all main cats

        # UVJ CORNER CUT, USE CUT, K-BAND SNR CUT, AND REDSHIFT CUT (zcut[i] < 2.369 if include F160W; else < 2.560)
        if objuse[i] == 1 and Ksnr[i] >= 10:  # 1 < zcut[i] < 2:  # Ksnr[i] >= 10:  # >=20
            ax_uv.append(UV[i])
            ax_vj.append(VJ[i])
            table.append(main[i][0])  # ID in catalog = main[i][0]

    # PLOT SCATTER OF REMAINING POINTS
    plt.scatter(ax_vj, ax_uv, color='0.5', alpha=0.1, marker=".")  # x, y
    if star:
        plt.scatter(special_vj, special_uv, color='b', marker="*", s=100)

    if title:
        plt.title(field + '-' + objname)
    if labels:
        plt.text(-0.4, 1.35, 'Quiescent', fontsize=16)  # label quiescent region
        plt.text(-0.4, 1.1, 'Star-forming', fontsize=16)  # label star-forming region
    plt.plot([-0.5, 0.8], [1.3, 1.3], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2]) x2=0.9 or 0.8?
    plt.plot([1.6, 1.6], [2.5, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2])
    plt.plot([0.8, 1.6], [1.3, 2.0], 'k')  # plot a line to show quiescent region ([x1,x2], [y1, y2]) x1=0.9 or 0.8

    if lims:
        plt.xlim(-1.5, 2.5)
        plt.xticks([-1, 0, 1, 2])
        plt.ylim(-1., 2.5)  # 3
        plt.yticks([-1., 0., 1., 2.])  # , 3.])

    plt.xlabel(r'$V - J$ (Rest)', fontsize=size)
    plt.ylabel(r'$U - V$ (Rest)', fontsize=size)
    if show:
        plt.show()
    '''


if __name__ == "__main__":
    # don't create keyword if not passed in!
    parser = argparse.ArgumentParser(argument_default=argparse.SUPPRESS)
    parser.add_argument('--obj')  # to print without blue star highlighting a given object, do --obj=-1
    parser.add_argument('--field')

    args = vars(parser.parse_args())
    kwargs = {}
    for key in args.keys():
        kwargs[key] = args[key]

    uvj_plot(kwargs['obj'], kwargs['field'])

'''
RUN WITH:
python uvj.py --obj=1824 --field=cosmos
'''
