#!/usr/bin/env python
import os, sys
import pickleTraj

sys.path.append('/home/cask0/home/tsanyal/benwat')
import cgmodel as cg

#TODO: populate these lists with distinct values for the different lammpstraj
# or with same value (need to enquire with Scott)
LDCutBB = []
LDCutBW = []
Temp = 300

cg.TempSet = Temp
cg.LammpsTraj = LammpsTraj
cg.Prefix = 'NB%dNW%d_multi' % (NB, NW)
cg.LDCutBB = LDCutBB
cg.LDCutBW = LDCutBW

trj = pickleTraj(LammpsTraj, Verbose = True)
cg.BoxL = trj.FrameData['BoxL']
Sys = cg.makeSys()

TrajList = []
NMol = 500
NBList = [50,100,150,200,250,300,350,400,450,500]
NWList = [] ; [NWList.append(500-x) for x in NBList]
for i, NB in enumerate(NBList):
	NW = NWList[i]
	Traj = os.path.join('/home/cask0/home/tsanyal/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod.lammpstrj.gz'
							% (NB, NW))
	
	TrajList.append(Traj)

cg.MultiLammpsTraj = TrajList
cg.MultiNBList = NBList
cg.MultiNWList = NWList
cg.runMultiSrel()