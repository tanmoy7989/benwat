#!/usr/bin/env python

import sys, os

data_dir = os.path.abspath('./data/gromacs')
pyMain = 'main_Gromacs.py'
sJobIn = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(jobname)s

#$ -m be
#$ -M tanmoy.7989@gmail.com

date
python %(pyMain)s %(NB)d %(NW)d %(Ncores)d
'''

NB = [250]
NW = [250]
Ncores = 8

Nruns = len(NB)
for i in range(Nruns):
    Prefix = 'NB%dNW%d' % (NB[i], NW[i])
    Dir = os.path.join(data_dir, Prefix)
    jobScriptName = os.path.join(Dir, Prefix+'.sh')
    
    if not os.path.isdir(Dir):  os.mkdir(Dir)
    file(jobScriptName, 'w').write(sJobIn % {'NB': NB[i], 'NW': NW[i],'jobname': Prefix, 'Ncores': Ncores, 'pyMain': pyMain})
    
    os.system('cp %s %s' % (pyMain, Dir))
    os.system('chmod 777 ' + jobScriptName)
    
