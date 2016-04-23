#!/usr/bin/env python
import os, sys
import pickleTraj

sys.path.append('/home/cask0/home/tsanyal/benwat')
import cgmodel as cg

MultiSrel = bool(sys.argv[1])
NB = int(sys.argv[2])
NW = int(sys.argv[3])
LammpsTraj = sys.argv[4]
LDCutBB = float(sys.argv[5])
LDCutBW = float(sys.argv[6])
Temp = 300

cg.NB = NB
cg.NW = NW
cg.TempSet = Temp
cg.LammpsTraj = LammpsTraj
cg.Prefix = 'NB%dNW%d' % (NB, NW)
cg.LDCutBB = LDCutBB
cg.LDCutBW = LDCutBW

trj = pickleTraj(LammpsTraj, Verbose = True)
cg.BoxL = trj.FrameData['BoxL']
Sys = cg.makeSys()


if MultiSrel: 
	TrajList = []
	for NB in [50,100,150,200,250,300,350,400,450,500]:
		NW = 500-NB
		Traj = os.path.join('/home/cask0/home/tsanyal/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod.lammpstrj.gz'
							% (NB, NW))
		TrajList.append(Traj)

	cg.MultiLammpsTraj = TrajList
	cg.runMultiSrel(Sys = Sys)

else: 
	cg.runSrel(Sys = Sys)
