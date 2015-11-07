#!/usr/bin/env python

import sys, os

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
python makeGromacs.py %(NB)d %(NW)d %(Temp)g %(Ncores)d
'''

NB = [10, 20]
NW = [500, 100]
xB = []
Temp = 300
Ncores = 4

Nruns = len(NB)
for i in range(Nruns):
    Prefix = 'NB%dNW%d' % (NB[i], NW[i])
    jobScriptName =Prefix+'.sh'
    file(jobScriptName, 'w').write(sJobIn % {'NB': NB[i], 'NW': NW[i], 'Temp': Temp, 'Ncores': Ncores, 'jobname': Prefix})
    os.system('chmod 777 ' + jobScriptName)
    
