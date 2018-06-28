import random

copyfrom = 'eelg_sfhtest_params.py'
m_pri = [8.8, 10.]
d_pri = [0.1, 0.6]
z_pri = [-2., -1.]
s1_pri = [0.05, 0.8]
s2_pri = [0.05, 0.2]

for x in range(200):
    m = random.uniform(m_pri[0], m_pri[1])
    d = random.uniform(d_pri[0], d_pri[1])
    z = random.uniform(z_pri[0], z_pri[1])
    s1 = random.uniform(s1_pri[0], s1_pri[1])
    s2 = random.uniform(s2_pri[0], s2_pri[1])

    mstr = str(m)
    dstr = str(d)
    zstr = str(z)
    s1str = str(s1)
    s2str = str(s2)

    newpar = 'bettertest/sfhtest_' + str(x) + '_params.py'
    # e.g. newpar = testpars/sfhtest_0_params.py
    orig = open(copyfrom, 'r')
    writenew = open(newpar, 'w+')
    for line in orig.readlines():  # for each line in original param file
        towrite = line  # copy the line over, unless I need to change it:
        if line.startswith("    model_params[n.index('logmass')]['init']"):
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
    writenew.write('# Codename: ' + mstr + '_' + dstr + '_' + zstr + '_' + s1str + '_' + s2str)
    orig.close()
    writenew.close()

print('done!')
