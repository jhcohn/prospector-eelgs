import random
import stellar_ages as sa

copyfrom = 'eelg_sfhtest_params.py'
m_pri = [8.5, 11.]  # [8.8, 10.]  # [8., 11.]
d_pri = [0., 2.]  # [0.1, 0.6]  # [0.0, 2.0]
z_pri = [-2., 0.]  # [-2., -1.]  # [-2., 0.19]
s1_pri = [0.01, 0.99]  # [0.05, 0.8]  # [0., 1.]
s2_pri = [0.01, 0.99]  # [0.05, 0.8]  # [0., 1.] (with min still)

em_pri = [9.0, 9.8]  # [8.8, 10.]  # [8., 11.]
ed_pri = [0.1, 0.5]  # [0.1, 0.6]  # [0.0, 2.0]
ez_pri = [-2., -1.5]  # [-2., -1.]  # [-2., 0.19]
es1_pri = [0.3, 0.7]  # [0.31, 0.68]  # [0.05, 0.8]  # [0., 1.]
es2_pri = [0.05, 0.15]  # [0.06, 0.15]  # [0.05, 0.8]  # [0., 1.] (with min still)
red_pri = [2.5, 4.0]
base = 'fico'
e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base, objlists=True, normal=True)

'''
for x in range(100):  # 200
'''
for x in range(2,50):
    print(x)
    '''
    m = random.uniform(m_pri[0], m_pri[1])
    d = random.uniform(d_pri[0], d_pri[1])
    z = random.uniform(z_pri[0], z_pri[1])
    s1 = random.uniform(s1_pri[0], s1_pri[1])
    s2 = random.uniform(s2_pri[0], min([1-s1, s2_pri[1]]))  # total fractions cannot add to more than 1!
    red = random.uniform(red_pri[0], red_pri[1])
    # '''

    '''
    m = 9.4
    d = 0.32
    z = -1.7
    s1 = 0.51
    s2 = 0.09
    red = 3.5
    # Identity: cdfs_20366
    # '''

    m = random.uniform(em_pri[0], em_pri[1])
    d = random.uniform(ed_pri[0], ed_pri[1])
    z = random.uniform(ez_pri[0], ez_pri[1])
    s1 = random.uniform(es1_pri[0], es1_pri[1])
    s2 = random.uniform(es2_pri[0], min([1-s1, es2_pri[1]]))  # total fractions cannot add to more than 1!
    red = random.uniform(red_pri[0], red_pri[1])

    idx = random.randint(0, len(e_objs)-1)
    obj = e_objs[idx]
    field = e_fields[idx]
    # '''

    mstr = str(m)
    dstr = str(d)
    zstr = str(z)
    s1str = str(s1)
    s2str = str(s2)

    '''
    newpar = 'ctest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    '''
    newpar = 'eehometest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    # newpar = 'eetest2/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    # newpar = 'newtest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    # e.g. newpar = testpars/sfhtest_0_params.py
    orig = open(copyfrom, 'r')
    writenew = open(newpar, 'w+')
    for line in orig.readlines():  # for each line in original param file
        towrite = line  # copy the line over, unless I need to change it:
        if line.startswith("              'field':"):
            towrite = "              'field': '" + field + "', \n"
        elif line.startswith("              'objname':"):
            towrite = "              'objname': '" + str(obj) + "',\n"
        elif line.startswith("    model_params[n.index('logmass')]['init']"):
            towrite = "    model_params[n.index('logmass')]['init'] = " + mstr + '\n'
        elif line.startswith("    model_params[n.index('dust2')]['init']"):
            towrite = "    model_params[n.index('dust2')]['init'] = " + dstr + '\n'
        elif line.startswith("    model_params[n.index('logzsol')]['init']"):
            towrite = "    model_params[n.index('logzsol')]['init'] = " + zstr + '\n'
        elif line.startswith("    eelg_frac1 = "):
            towrite = "    eelg_frac1 = " + s1str + '\n'
        elif line.startswith("    eelg_frac2 = "):
            towrite = "    eelg_frac2 = " + s2str + '\n'
        elif line.startswith("    zred = "):
            towrite = "    zred = " + str(red) + "\n"
        writenew.write(towrite)  # + '\n')
    writenew.write('# Codename: ' + mstr + '_' + dstr + '_' + zstr + '_' + s1str + '_' + s2str + '\n')
    writenew.write("# Identity: " + field + "_" + str(obj) + "_" + str(red))
    orig.close()
    writenew.close()

print('done!')
