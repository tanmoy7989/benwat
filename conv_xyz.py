#usr/bin/env python

import mysim

import os, sys
import numpy as np
import sim, pickleTraj

DelTempFiles = True

srcDir = os.path.expanduser('~/benwat/data/Christine_Peter_data/AA_reftraj')
tarDir = os.path.expanduser('~/benwat/data/gromacs')
reftable = np.loadtxt(os.path.join(srcDir, 'table.txt'))


## VMD tcl script
tclscript = '''
mol new %(xyz)s waitfor all
animate write lammpstrj %(lammpstrj)s waitfor all
quit
'''

## user-input
conc = sys.argv[1]
for i in range(len(reftable)):
	if reftable[i,0] == float(conc):
		NB = reftable[i,1]
		NW = reftable[i,2]
		BoxL = reftable[i,3]

# parse filenames
zipTrj = os.path.join(srcDir, 'conc_%s' % conc, 'CG_trajectory_%s.xyz.bz2' % conc)
unzipTrj = os.path.join(tarDir, 'NB%dNW%d' % (NB,NW), 'CG_trajectory_%s.xyz' % conc)
LammpsTrj = os.path.join(tarDir, 'NB%dNW%d' % (NB,NW), 'NB%dNW%d_prod.lammpstrj' % (NB, NW))
TCLFile = os.path.join(tarDir, 'conv_xyz.tcl')

# unzip traj and copy it to required location
if not os.path.isdir(os.path.join(tarDir, 'NB%dNW%d' % (NB, NW))): os.mkdir(os.path.join(tarDir, 'NB%dNW%d' % (NB, NW)))
os.system('bunzip2 %s ; mv %s %s' % (zipTrj, zipTrj.split('.bz2')[0], unzipTrj))

# run vmd
VMDExec = '/home/cask0/home/tsanyal/software/tanmoy_vmd/vmd'
VMDParams = {'xyz' : unzipTrj, 'lammpstrj': LammpsTrj}
file(TCLFile, 'w').write(tclscript % VMDParams)
os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile))
os.system('gzip %s' % LammpsTrj)

# pickle final traj
pickleTraj(LammpsTrj+'.gz', Verbose = True)


if DelTempFiles: 
	[os.remove(x) for x in [TCLFile, unzipTrj]]
