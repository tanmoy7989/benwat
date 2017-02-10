#!/usr/bin/env python
import os, sys
import pickleTraj

import cgmodel as cg

# LDCuts:--
# BB: 1st shell or 2nd shell
# WW: 1st shell or 2nd shell (both must be same shell)
# BW and WB : (BB+WW)/2

NB = int(sys.argv[1])
NW = int(sys.argv[2])
LammpsTraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod_mapped.lammpstrj.gz' % (NB, NW, NB, NW)
LDCutBB = float(sys.argv[3])
LDCutWW = float(sys.argv[4])
LDCutBW = 0.5 * (LDCutBB + LDCutWW)
LDCutWB = LDCutBW
Temp = 300

cg.NB = NB
cg.NW = NW
cg.Prefix = 'NB%dNW%d' % (NB, NW)
cg.TempSet = Temp
cg.LammpsTraj = LammpsTraj
cg.LDCutBB = LDCutBB
cg.LDCutWW = LDCutWW
cg.LDCutBW = LDCutBW
cg.LDCutWB = LDCutWB

cg.MinSteps = 1000
cg.EquilSteps = 1000000
cg.ProdSteps = 2000000
cg.StepFreq = 100

Sys = cg.makeSys()

cg.useLammps = True


# unconstrain knots
for P in Sys.ForceField:
    P.EneSlopeInner = None

# refine SP-only
cg.OptStageNames = ['SP']
cg.OptStages = {'SP': ['SP_WW', 'SP_BB', 'SP_BW']}
cg.runSrel(Sys = Sys, ParamString = file('ff_SP.dat').read())

# refine SPLD_BB
cg.OptStageNames = ['SPLD_BB']
cg.OptStages = {'SPLD_BB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB']}
cg.runSrel(Sys = Sys, ParamString = file('ff_SPLD_BB.dat').read())

#refine SPLD_WW
cg.OptStageNames = ['SPLD_WW']
cg.OptStages = {'SPLD_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WW']}
cg.runSrel(Sys = Sys, ParamString = file('ff_SPLD_WW.dat').read())

#refine SPLD_BW
cg.OptStageNames = ['SPLD_BW']
cg.OptStages = {'SPLD_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BW']}
cg.runSrel(Sys = Sys, ParamString = file('ff_SPLD_BW.dat').read())

#refine all
cg.OptStageNames = []
cg.runSrel(Sys = Sys, ParamString = file('ff_all.dat').read())
