import os
import cosmology_calculator as cc
import numpy as np

path = '/home/jonathan/'
eels = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/eelg_specz_ids'
z_comps = ['C_04_rshift.txt', 'C_10_rshift.txt', 'C_12_rshift.txt', 'C_14_rshift.txt']
all_r = []
for file in z_comps:
    print(file)
    ids = []
    zs = []
    zp = []
    r_pix = []
    fields = []
    r_phys = []

    info = open(path + file, 'r')
    for line in info:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            if int(cols[0]) - 200000 > 0:
                ids.append(int(cols[0]) - 200000)
                fields.append('uds')
            elif int(cols[0]) - 100000 > 0:
                ids.append(int(cols[0]) - 100000)
                fields.append('cosmos')
            else:
                ids.append(int(cols[0]))
                fields.append('cdfs')
            # ids.append(cols[0])
            zs.append(float(cols[1]))
            zp.append(float(cols[2]))
            r_pix.append(float(cols[3]))
    info.close()
    print(len(ids))

    eelgs = open(eels, 'r')
    e_objs = []
    e_fields = []
    for line in eelgs:
        if line[0] == '#':
            pass
        else:
            cols = line.split()
            e_objs.append(cols[1])
            e_fields.append(cols[0])
    eelgs.close()
    print(len(e_objs))

    county = 0
    for i in range(len(ids)):
        it_is = 1
        '''
        for j in range(len(e_objs)):
            if (ids[i] == int(e_objs[j])) and (fields[i] == e_fields[j]) and (zs[i] > 0):
                print('hi')
                it_is = 1  # True!
        '''
        if it_is:
            # print(fields[i], ids[i])
            county += 1
            z = zs[i]
            if z < 0:
                z = zp[i]
            # print(z)
            cmrd = cc.calc(z=z, H0=70, WM=0.3, WV=0.7)[2]
            # print(cmrd)
            r_phys.append(r_pix[i] * 0.149994 * 1.7*10**-6 * cmrd * 10**3 / (1 + z))
            # = r_pix * arcsec/pix * rad/arsec * comoving radial distance [Mpc] * 10**3 [kpc/Mpc] / (1+z)

    all_r += r_phys
    print(county)
    print(np.percentile(r_phys, [16., 50., 84.]))

print(len(all_r))
print(np.percentile(all_r, [16., 50., 84.]))
# [ 0.72900932  0.90181655  1.09167149]
# 0.90 +0.19/-0.17 kpc
