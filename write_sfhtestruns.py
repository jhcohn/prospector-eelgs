import numpy as np

# obj = str(11063)
# field = 'cosmos'
this = 'sfhtest'

time = '25:00'  # '30:00'

base = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/'
# newfile = open(base + 'lbg_otherruns.sh', 'w+')  # (base + 'lbg_runs.sh', 'w+')
# newfile = open(base + 'eelg_otherruns.sh', 'w+')  # (base + 'eelg_runs.sh', 'w+')
# newfile = open(base + 'dsfg_runs.lsf', 'w+')  # (base + 'dsfg_otherruns.lsf', 'w+')

line1 = '# NECESSARY JOB SPECIFICATIONS\n'
line2 = '#BSUB -J '  # add newname name + '\n' to this
line3 = '#BSUB -L /bin/bash\n'
line4 = '#BSUB -W ' + time + '\n'  # 40:00  # each galaxy expected to take 6 SUs, allow up to 10-15 per galaxy?
line5 = '#BSUB -n 4\n'
line6 = '#BSUB -R "span[ptile=1]"\n'
line7 = '#BSUB -R "rusage[mem=2560]"\n'
line8 = '#BSUB -M 2560\n'
line9 = '#BSUB -u joncohn@tamu.edu\n'
line10 = ''  # '#BSUB -o yepthirty_out.%J \n \n' # REPLACE "yepthirty" with newnames name
line11 = 'ml purge\n'
line12 = 'source /home/joncohn/myfsps.sh\n\n'

lines = [line1, line2, line3, line4, line5, line6, line7, line8, line9, line10, line11, line12]

# e_nums = [str(int(num)) for num in np.arange(1, 100)]  # ['1', ..., '10']
# newnames = ['eelg_' + this + 'runs' + n + '.lsf' for n in e_nums]
# NOTE: ALL GENERATED PHOT WAS BASED ON COSMOS-11063

idx = 0
'''
for z in range(50):  # 100
'''
for z in range(100):  # 100  # 1, 25 <--for eetest for now
    name = 'eelg_' + this + 'runs' + str(z) + '.lsf'
    parfile = 'eelg_test' + str(z) + '_params.py'  # modified from 'eelg_fifty_params.py' on the cluster
    '''
    parfile2 = 'eelg_test' + str(z+50) + '_params.py'  # modified from 'eelg_fifty_params.py' on the cluster
    # parfile2 = 'eelg_test' + str(z+100) + '_params.py'  # modified from 'eelg_fifty_params.py' on the cluster
    genfile = 'ctest/sfhtest_' + str(z) + '_params.py'
    genfile2 = 'ctest/sfhtest_' + str(z + 50) + '_params.py'
    # genfile = 'eetest/sfhtest_' + str(z) + '_params.py'
    # genfile2 = 'eetest/sfhtest_' + str(z+10) + '_params.py'
    # genfile = 'newtest/sfhtest_' + str(z) + '_params.py'
    # genfile2 = 'newtest/sfhtest_' + str(z+100) + '_params.py'
    '''
    genfile = 'eetest2/sfhtest_' + str(z) + '_params.py'
    run_name = 'mpirun -n 4 python prospector.py --param_file=' + parfile + ' --niter=2500 --outfile=out_sfhtest/'
    '''
    run_name2 = 'mpirun -n 4 python prospector.py --param_file=' + parfile2 + ' --niter=2500 --outfile=out_sfhtest/'
    '''

    field = ''
    obj = ''
    red = ''
    with open(genfile, 'r') as pfile:
        for line in pfile:
            if line.startswith('# Codename: '):
                pset = line[12:]
            elif line.startswith('# Identity: '):
                counterl = 0
                for l in line:
                    if l == ' ' or l == '_':
                        counterl += 1
                    elif counterl == 2:
                        field += l
                    elif counterl == 3:
                        obj += l
                    elif counterl == 4:
                        red += l

    field2 = ''
    obj2 = ''
    red2 = ''
    '''
    with open(genfile2, 'r') as pfile:
        for line in pfile:
            if line.startswith('# Codename: '):
                pset = line[12:]
            elif line.startswith('# Identity: '):
                counterl = 0
                for l in line:
                    if l == ' ' or l == '_':
                        counterl += 1
                    elif counterl == 2:
                        field2 += l
                    elif counterl == 3:
                        obj2 += l
                    elif counterl == 4:
                        red2 += l

    newfile = open(base + 'copy_stuff/' + name, 'w+')
    '''
    newfile = open(base + 'copy_eestuff2/' + name, 'w+')
    lines[1] = '#BSUB -J ' + name + '\n'  # lines[1] = line2
    lines[9] = '#BSUB -o ' + name[:-4] + '_out.%J\n\n'
    for line in lines:
        newfile.write(line)

    newfile.write(run_name + obj + '_' + field + '_' + this + '_' + str(z) + ' --field=' + field + ' --objname=' +
                  str(z) + '\n' + '\n')

    '''
    # newfile.write(run_name2 + obj2 + '_' + field2 + '_' + this + '_' + str(z+100) + ' --field=' + field2
    #               + ' --objname=' + str(z+100) + '\n' + '\n')
    newfile.write(run_name2 + obj2 + '_' + field2 + '_' + this + '_' + str(z + 50) + ' --field=' + field2
                  + ' --objname=' + str(z + 50) + '\n' + '\n')
    '''
    newfile.close()

    '''
    newpfile = open(base + 'copy_stuff/' + 'eelg_test' + str(z) + '_params.py', 'w+')
    newpfile2 = open(base + 'copy_stuff/' + 'eelg_test' + str(z+50) + '_params.py', 'w+')
    # newpfile2 = open(base + 'copy_stuff/' + 'eelg_test' + str(z+100) + '_params.py', 'w+')
    '''
    newpfile = open(base + 'copy_eestuff2/' + 'eelg_test' + str(z) + '_params.py', 'w+')
    newpfile2 = open(base + 'deleteme', 'w+')
    with open('eelg_test_params.py', 'r') as etp:
        for line_etp in etp:
            if line_etp.startswith("                if cols[0] == "):
                newpfile.write("                if cols[0] == str(" + str(z) + "):\n")
                newpfile2.write("                if cols[0] == str(" + str(z+50) + "):\n")
                # newpfile2.write("                if cols[0] == str(" + str(z+100) + "):\n")
            elif line_etp.startswith("        photname = '/home/jonathan/"):
                    newpfile.write("        photname = '/scratch/user/joncohn/cosmos.v1.3.8.cat'\n")
                    newpfile2.write("        photname = '/scratch/user/joncohn/cosmos.v1.3.8.cat'\n")
            elif line_etp.startswith("        testphotname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_cosmos'"):
                '''
                newpfile.write("        testphotname = '/scratch/user/joncohn/tphot_cosmos'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/tphot_cosmos'\n")
                '''
                # newpfile.write("        testphotname = '/scratch/user/joncohn/sfhphot_cosmos'\n")
                # newpfile2.write("        testphotname = '/scratch/user/joncohn/sfhphot_cosmos'\n")
                newpfile.write("        testphotname = '/scratch/user/joncohn/eephot2_cosmos'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/eephot2_cosmos'\n")
            elif line_etp.startswith("        testphotname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_cdfs'"):
                '''
                newpfile.write("        testphotname = '/scratch/user/joncohn/tphot_cdfs'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/tphot_cdfs'\n")
                '''
                # newpfile.write("        testphotname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/eephot_cdfs'\n")
                # newpfile2.write("        testphotname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/eephot_cdfs'\n")
                # newpfile.write("        testphotname = '/scratch/user/joncohn/sfhphot_cdfs'\n")
                # newpfile2.write("        testphotname = '/scratch/user/joncohn/sfhphot_cdfs'\n")
                newpfile.write("        testphotname = '/scratch/user/joncohn/eephot2_cdfs'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/eephot2_cdfs'\n")
            elif line_etp.startswith("        testphotname = '/home/jonathan/.conda/envs/snowflakes/lib/python2.7/site-packages/prospector/git/sfhphot_uds'"):
                '''
                newpfile.write("        testphotname = '/scratch/user/joncohn/tphot_uds'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/tphot_uds'\n")
                '''
                # newpfile.write("        testphotname = '/scratch/user/joncohn/sfhphot_uds'\n")
                # newpfile2.write("        testphotname = '/scratch/user/joncohn/sfhphot_uds'\n")
                newpfile.write("        testphotname = '/scratch/user/joncohn/eephot2_uds'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/eephot2_uds'\n")
            elif line_etp.startswith("        zname = '/home/jonathan/"):
                newpfile.write("        zname = '/scratch/user/joncohn/cosmos.v1.3.6.awk.zout'\n")
                newpfile2.write("        zname = '/scratch/user/joncohn/cosmos.v1.3.6.awk.zout'\n")
            elif line_etp.startswith("              'field':"):
                newpfile.write("              'field': '" + field + "',\n")
                newpfile2.write("              'field': '" + field2 + "',\n")
            elif line_etp.startswith("              'objname':"):
                newpfile.write("              'objname': '" + obj + "',\n")
                newpfile2.write("              'objname': '" + obj2 + "',\n")
            elif line_etp.startswith("    zred = "):
                newpfile.write("    zred = " + red + "\n")
                newpfile2.write("    zred = " + red2 + "\n")
            else:
                newpfile.write(line_etp)
                newpfile2.write(line_etp)
    newpfile.close()
    newpfile2.close()

    print(z, z+50)
    # print(z, z+100)

'''
elif line_etp.startswith("    zred_idx = (zout['id']"):
    newpfile.write("    zred_idx = (zout['id'] == str("+obj+"))\n")
    newpfile2.write("    zred_idx = (zout['id'] == str("+obj2+"))\n")
elif line_etp.startswith("    idx = dat['id'] == "):
    newpfile.write("    idx = dat['id'] == str("+obj+")\n")
    newpfile2.write("    idx = dat['id'] == str("+obj2+")\n")
'''