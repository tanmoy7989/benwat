#!/usr/bin/env python
import os, sys
import cPickle as pickle

import cgmodel as cg

NB = int(sys.argv[1])
NW = int(sys.argv[2])
FFType = sys.argv[3]
if len(sys.argv) > 4:
    useParallel = bool(sys.argv[4])
    NCores = int(sys.argv[5]) if useParallel else 1
else:
    useParallel = False
    NCores = 1

LammpsTraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod_mapped.lammpstrj.gz' % (NB, NW, NB, NW)
FF_file = '/home/cask0/home/tsanyal/benwat/data/cgff/NB250NW250/control/NB250NW250_%s_ff.dat' % FFType
Temp = 300

cg.NB = NB
cg.NW = NW
cg.Prefix = 'NB%dNW%d_%s' % (NB, NW, FFType)
cg.TempSet = Temp
cg.LammpsTraj = LammpsTraj
cg.LDCutBB = 7.5
cg.LDCutWW = 3.5
cg.LDCutBW = 0.5 * (cg.LDCutBB + cg.LDCutWW)
cg.LDCutWB = cg.LDCutBW

cg.MinSteps = 1000
cg.EquilSteps = 2000000
cg.ProdSteps = 20000000
cg.StepFreq = 1000

Sys = cg.makeSys()
modtraj, modtrajfile = cg.runMD(Sys, ParamString = file(FF_file).read(),
                                useParallel = useParallel, NCores = NCores)
