#!/usr/bin/env python
import os, sys
import pickleTraj

sys.path.append('/home/cask0/home/tsanyal/benwat')
import cgmodel as cg

NB = int(sys.argv[1])
NW = int(sys.argv[2])
LammpsTraj = sys.argv[3]
LDCutBW = float(sys.argv[4])
Temp = 300

cg.NB = NB
cg.NW = NW
cg.TempSet = Temp
cg.LammpsTraj = LammpsTraj
cg.Prefix = 'NB%dNW%d' % (NB, NW)
cg.LDCutBW = LDCutBW

trj = pickleTraj(LammpsTraj, Verbose = True)
cg.BoxL = trj.FrameData['BoxL']
Sys = cg.makeSys()
cg.runSrel(Sys = Sys)
