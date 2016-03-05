#usr/bin/env python

import os, sys
import numpy as np
import sim, pickleTraj

DelTempFiles = True

## VMD tcl script
tclscript = '''
mol addfile %(gro)s
package require topotools
topo writelammpsdata %(lammpsdata)s
mol addfile %(xtc)s waitfor all
animate write lammpstrj %(lammpstrj)s waitfor all
quit
'''

## user-input
GromacsStruct = sys.argv[1]
GromacsTrj = sys.argv[2]
NB = int(sys.argv[3])
NW = int(sys.argv[4])

# parse filenames
TrjDir = os.path.dirname(GromacsTrj)
TrjName = GromacsTrj.split('/')[-1]
LammpsData = os.path.join(TrjDir, 'temp.data')
LammpsTrj = GromacsTrj.split('.')[0] + '_unmapped.lammpstrj'
TCLFile = os.path.join(TrjDir, 'conv_trj.tcl')

# run vmd
VMDExec = '/home/cask0/home/tsanyal/software/tanmoy_vmd/vmd'
VMDParams = {'gro' : GromacsStruct, 'xtc': GromacsTrj, 'lammpsdata': LammpsData, 'lammpstrj': LammpsTrj}
file(TCLFile, 'w').write(tclscript % VMDParams)
os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile))
os.system('gzip %s' % LammpsTrj)

# map traj using sim
LammpsTrjIn = LammpsTrj + '.gz'
LammpsTrjOut = GromacsTrj.split('.')[0] + '.lammpstrj.gz'
Trj = pickleTraj(LammpsTrjIn, Verbose = False)
BoxL = Trj.FrameData['BoxL']
Map = sim.atommap.PosMap()
# vmd created lammpstraj places all benzene atoms at the beginning
for i in range(0, NB):  Map += [sim.atommap.AtomMap(Atoms1 = range(i*12, (i+1)*12), Atom2 = i)]
for i in range(0, NW):  Map += [sim.atommap.AtomMap(Atoms1 = NB*12 + 3*i, Atom2 = NB+i)]
AtomTypes = [1]*NB + [2]*NW
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes)
sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrjOut, Verbose = True)            

# pickle final traj
pickleTraj(LammpsTrjOut, Verbose = True)


if DelTempFiles:
    for x in [LammpsData, TCLFile]:
        os.remove(x)
