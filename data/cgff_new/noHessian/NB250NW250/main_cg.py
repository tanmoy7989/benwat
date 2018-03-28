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
FFType = sys.argv[5]
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

# better resolved WW LD potential
cg.RhoMin = 0.0
cg.RhoMax_WW = 6.0
cg.NLDWWKnots = 60
cg.useMoreLDWWKnots = True

cg.MinSteps = 1000
cg.EquilSteps = 500000#1000000
cg.ProdSteps = 1000000#2000000
cg.StepFreq = 200

Sys = cg.makeSys()

# seed with previously obtaind forcefields
#ff_file = os.path.abspath('../../init/NB%dNW%d_%s_ff.dat' % (NB, NW, FFType))

# NB250NW250_SP run was halted so restarting
ff_file = os.path.abspath('./NB250NW250_SP_ff_halted.dat')

ff = file(ff_file).read()
Sys.ForceField.SetParamString(ff)

# stages
LDStageDict = {'SP'    : ['SP_WW', 'SP_BB', 'SP_BW'],
              'SPLD_BB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB'],
              'SPLD_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WW'],
              'SPLD_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BW'],
              'SPLD_WB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WB'],
              'SPLD_all': [P.Name for P in Sys.ForceField if not P is None]
            }

# optimize ff for user provided fftype
cg.OptStageNames = [FFType]
cg.OptStages = {FFType: LDStageDict[FFType]}

# run without Hessian
#cg.StartWithHessian = False
#cg.NoHessianMaxIter = 100

# running with Hessian the halted NB250NW250_SP job
cg.StartWithHessian = True

cg.runSrel(Sys = Sys, ParamString = ff)
