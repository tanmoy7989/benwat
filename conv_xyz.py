#usr/bin/env python

import mysim

import os, sys
import numpy as np
import sim, pickleTraj

DelTempFiles = True

srcDir = os.path.expanduser('~/benwat/data/Christine_Peter_data/AA_reftraj')
tarDir = os.path.expanduser('~/benwat/data/gromacs')
reftable = np.loadtxt(os.path.join(srcDir, 'table.txt'))
curr_dir = os.getcwd()


## VMD tcl script
tclscript = '''
mol addfile %(xyz)s waitfor all
package require pbctools
pbc set {%(boxl)g %(boxl)g %(boxl)g} -all
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
xyzTrj_original = os.path.join(srcDir, 'conc_%s' % conc, 'CG_trajectory_%s.xyz' % conc)
xyzTrj = os.path.join(tarDir, 'NB%dNW%d' % (NB,NW), 'CG_trajectory_%s.xyz' % conc)
LammpsTrj = os.path.join(tarDir, 'NB%dNW%d' % (NB,NW), 'NB%dNW%d_prod.lammpstrj' % (NB, NW))
TCLFile = os.path.join(tarDir, 'NB%dNW%d' % (NB, NW), 'conv_xyz.tcl')

# copy traj to required location and unzip it
if not os.path.isdir(os.path.join(tarDir, 'NB%dNW%d' % (NB, NW))): os.mkdir(os.path.join(tarDir, 'NB%dNW%d' % (NB, NW)))
os.system('cp %s %s' % (xyzTrj_original, xyzTrj))

# run vmd
VMDExec = '/home/cask0/home/tsanyal/software/tanmoy_vmd/vmd'
VMDParams = {'xyz' : xyzTrj.split('/')[-1], 'lammpstrj': LammpsTrj.split('/')[-1], 'boxl': BoxL}
file(TCLFile, 'w').write(tclscript % VMDParams)
os.chdir(os.path.join(tarDir, 'NB%dNW%d' % (NB, NW)))
os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile))
os.system('gzip %s' % LammpsTrj)
os.chdir(curr_dir)

# pickle final traj
pickleTraj(LammpsTrj+'.gz', Verbose = True)


if DelTempFiles: 
	[os.remove(x) for x in [TCLFile, xyzTrj]]
