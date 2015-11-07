#!/usr/bin/env python

import sys, os


pyMain = os.path.join('data', 'gromacs', 'main.py')
sJobIn = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(jobname)s
#$ -pe ompi %(Ncores)d

#$ -m be
#$ -M tanmoy.7989@gmail.com

date
python main.py %(NB)d %(NW)d
'''

NB = [10, 20]
NW = [500, 100]


Nruns = len(NB)
for i in range(Nruns):
    Prefix = 'NB%dNW%d' % (NB[i], NW[i])
    jobScriptName =os.path.join('data', 'gromacs', Prefix+'.sh')
    file(jobScriptName, 'w').write(sJobIn % {'NB': NB[i], 'NW': NW[i],'jobname': Prefix, 'Ncores': 8})
    os.system('chmod 777 ' + jobScriptName)
    
