#usr/bin/env python

import mysim

import os, sys
import numpy as np
import sim, pickleTraj

DelTempFiles = True

## VMD tcl script (selects only C,H atoms from benzene and O atoms from water)
tclscript = '''
mol new %(gro)s
mol addfile %(xtc)s waitfor all
set selOxygen [atomselect top "not name HW1 and not name HW2"]
animate write lammpstrj %(lammpstrj)s waitfor all sel $selOxygen
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
VMDParams = {'gro' : GromacsStruct, 'xtc': GromacsTrj,'lammpstrj': LammpsTrj}
file(TCLFile, 'w').write(tclscript % VMDParams)
os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile))
os.system('gzip %s' % LammpsTrj)

# map traj using sim
LammpsTrjIn = LammpsTrj + '.gz'
LammpsTrjOut = GromacsTrj.split('.')[0] + '.lammpstrj.gz'
Trj = pickleTraj(LammpsTrjIn, Verbose = False)
BoxL = Trj.FrameData['BoxL']
Map = sim.atommap.PosMap()
# vmd created lammpstraj places all benzene atoms at the beginning in (C1,H1,C2,H2,...) sequence
MassList_B = [12.011, 1.008] * 6
for i in range(0, NB):  Map += [sim.atommap.AtomMap(Atoms1 = range(i*12, (i+1)*12), Atom2 = i, Mass1 = MassList_B)]
for i in range(0, NW):  Map += [sim.atommap.AtomMap(Atoms1 = NB*12 + i, Atom2 = NB+i)]
AtomTypes = [1]*NB + [2]*NW
MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)
sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, LammpsTrjOut, Verbose = True)            

# pickle final traj
pickleTraj(LammpsTrjOut, Verbose = True)


if DelTempFiles: os.remove(TCLFile)
