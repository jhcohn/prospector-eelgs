import random
import stellar_ages as sa

copyfrom = 'eelg_sfhtest_params.py'
m_pri = [8.8, 10.]
d_pri = [0.1, 0.6]
z_pri = [-2., -1.]
s1_pri = [0.05, 0.8]
s2_pri = [0.05, 0.8]
base = 'fico'
e_objs, e_fields, l_objs, l_fields = sa.get_gal_lists(base, objlists=True, normal=True)

for x in range(200):
    m = random.uniform(m_pri[0], m_pri[1])
    d = random.uniform(d_pri[0], d_pri[1])
    z = random.uniform(z_pri[0], z_pri[1])
    s1 = random.uniform(s1_pri[0], s1_pri[1])
    s2 = random.uniform(s2_pri[0], min([1-s1, s2_pri[1]]))  # total fractions cannot add to more than 1!

    idx = random.randint(0, len(e_objs)-1)
    obj = e_objs[idx]
    field = e_fields[idx]

    mstr = str(m)
    dstr = str(d)
    zstr = str(z)
    s1str = str(s1)
    s2str = str(s2)

    newpar = 'besttest/sfhtest_' + str(x) + '_params.py'  # 'bettertest/sfhtest_' ...
    # e.g. newpar = testpars/sfhtest_0_params.py
    orig = open(copyfrom, 'r')
    writenew = open(newpar, 'w+')
    for line in orig.readlines():  # for each line in original param file
        towrite = line  # copy the line over, unless I need to change it:
        if line.startswith("              'field':"):
            towrite = "              'field': " + field + '\n'
        elif line.startswith("              'objname':"):
            towrite = "              'objname': " + obj + '\n'
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
        writenew.write(towrite)  # + '\n')
    writenew.write('# Codename: ' + mstr + '_' + dstr + '_' + zstr + '_' + s1str + '_' + s2str + '\n')
    writenew.write('# Identity: ' + field + '_' + obj)
    orig.close()
    writenew.close()

print('done!')
