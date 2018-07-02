import numpy as np
import math

obj = str(11063)
field = 'cosmos'
this = 'sfhtest'

time = '25:00'

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

e_nums = [str(int(num)) for num in np.arange(1, 100)]  # ['1', ..., '10']
newnames = ['eelg_' + this + 'runs' + n + '.lsf' for n in e_nums]
# NOTE: ALL GENERATED PHOT WAS BASED ON COSMOS-11063

idx = 0
for z in range(100):
    name = 'eelg_' + this + 'runs' + str(z) + '.lsf'
    parfile = 'eelg_test' + str(z) + '_params.py'  # modified from 'eelg_fifty_params.py' on the cluster
    parfile2 = 'eelg_test' + str(z+100) + '_params.py'
    run_name = 'mpirun -n 4 python prospector.py --param_file=' + parfile + ' --niter=2500 --outfile=out_sfhtest/'
    run_name2 = 'mpirun -n 4 python prospector.py --param_file=' + parfile2 + ' --niter=2500 --outfile=out_sfhtest/'

    newfile = open(base + 'copy_stuff/' + name, 'w+')
    lines[1] = '#BSUB -J ' + name + '\n'  # lines[1] = line2
    lines[9] = '#BSUB -o ' + name[:-4] + '_out.%J\n\n'
    for line in lines:
        newfile.write(line)

    newfile.write(run_name + obj + '_' + field + '_' + this + '_' + str(z) + ' --field=' + field + ' --objname=' +
                  str(z) + '\n' + '\n')

    newfile.write(run_name2 + obj + '_' + field + '_' + this + '_' + str(z+100) + ' --field=' + field + ' --objname=' +
                  str(z+100) + '\n' + '\n')
    newfile.close()

    newpfile = open(base + 'copy_stuff/' + 'eelg_test' + str(z) + '_params.py', 'w+')
    newpfile2 = open(base + 'copy_stuff/' + 'eelg_test' + str(z+100) + '_params.py', 'w+')
    with open('eelg_test_params.py', 'r') as etp:
        for line_etp in etp:
            if line_etp.startswith("                if cols[0] == "):
                newpfile.write("                if cols[0] == str(" + str(z) + "):\n")
                newpfile2.write("                if cols[0] == str(" + str(z+100) + "):\n")
            elif line_etp.startswith("        photname = '/home/jonathan/"):
                    newpfile.write("        photname = '/scratch/user/joncohn/cosmos.v1.3.8.cat'\n")
                    newpfile2.write("        photname = '/scratch/user/joncohn/cosmos.v1.3.8.cat'\n")
            elif line_etp.startswith("        testphotname = '/home/jonathan/"):
                newpfile.write("        testphotname = '/scratch/user/joncohn/sfhphot_all'\n")
                newpfile2.write("        testphotname = '/scratch/user/joncohn/sfhphot_all'\n")
            elif line_etp.startswith("        zname = '/home/jonathan/"):
                newpfile.write("        zname = '/scratch/user/joncohn/cosmos.v1.3.6.awk.zout'\n")
                newpfile2.write("        zname = '/scratch/user/joncohn/cosmos.v1.3.6.awk.zout'\n")
            else:
                newpfile.write(line_etp)
                newpfile2.write(line_etp)
    newpfile.close()
    newpfile2.close()

    print(z, z+100)

