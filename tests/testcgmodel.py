#/usr/bin/env python
import os, sys
import numpy as np
import pickleTraj

sys.path.append('../')
import cgmodel

cgmodel.doMinimize = False
cgmodel.NB = 50 ; cgmodel.NW = 450
cgmodel.LDCutBW = 7.8
cgmodel.LammpsTraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB50NW450/NB50NW450_npt.lammpstrj.gz'
cgmodel.Prefix = 'test'

trj = pickleTraj(cgmodel.LammpsTraj, Verbose = True)
cgmodel.BoxL = trj.FrameData['BoxL']
Sys = cgmodel.makeSys()
cgmodel.runSrel(Sys = Sys)

