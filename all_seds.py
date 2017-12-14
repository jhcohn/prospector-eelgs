import os
import make_all_plots as map
import pickle
import numpy as np
import print_sfh

home = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
eelgs = 0  # True=1
sfh = True  # if not sfh, then seds
folders = ['pkl_evar', 'pkl_nvar']  # ['pkl_emask', 'pkl_nmask']
file = 'vary'  # 'newmask'

objs = []
fields = []
count = 0
if eelgs:
    folder = folders[0]
    galxs = open(home + 'eelg_specz_ids', 'r')  # open(base + 'specz_ids', 'r')
    for line in galxs:
        skip = False
        if line[0] == '#':  # ignore commented lines
            skip = True
        space = 0  # counts which column we're in, using python indexing (space = 0 --> 1st column)
        cols = line.split()
        field = cols[0]
        obj = cols[1]
        if not skip:
            count += 1
            objs.append(obj)
            fields.append(field)
    galxs.close()
else:
    folder = folders[1]
    galxs = open(home + 'lbg_ids', 'r')
    for line in galxs:
        count += 1
        if int(line) - 200000 > 0:
            objs.append(str(int(line) - 200000))
            fields.append('uds')
        elif int(line) - 100000 > 0:
            objs.append(str(int(line) - 100000))
            fields.append('cosmos')
        else:
            objs.append(str(int(line)))
            fields.append('cdfs')
    galxs.close()

max = []
for i in range(len(objs)):
    obj = objs[i]
    field = fields[i]
    pre = folder + '/' + obj + '_' + field + '_' + file + '_'
    base = '_out.pkl'
    extra = pre + 'extra' + base  # includes SFH *AND* rest of extra_output, so call it extra and not sfh
    res = pre + 'res' + base
    sed = pre + 'sed' + base
    restwave = pre + 'restwave' + base
    spec = pre + 'spec' + base
    spswave = pre + 'spswave' + base
    chisq = pre + 'chisq' + base
    justchi = pre + 'justchi' + base

    files = [extra, res, sed, restwave, spec, spswave, chisq, justchi]
    # map.all_plots(files, obj, field, file, loc='upper left', sfh=False, curves=False, sep_uvj=False)

    if os.path.exists(files[1]) and os.path.exists(files[7]):
        if sfh:
            print(files[0])
            print_sfh.plotter(files[0], specific=True)
        else:
            with open(files[1], 'rb') as res:
                results = pickle.load(res)
            mask = results['obs']['phot_mask']  # mask out -99.0 values
            with open(files[7], 'rb') as justchi:
                chi = pickle.load(justchi)
            max.append(abs(chi[mask]).max())
            print(1)

if not sfh:
    print(np.percentile(max, [16, 50, 84]))
