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

# stages
ffSP= file('NB250NW250_SP_ff.dat').read()
LDStageDict = {
              'SPLD_BB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB'],
              'SPLD_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WW'],
              'SPLD_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BW'],
              'SPLD_WB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WB'],
              
              'SPLD_BB_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB', 'LD_WW'],
              'SPLD_BB_BW': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_BB', 'LD_BW'],
              'SPLD_BB_WB': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_BB', 'LD_WB'],
              'SPLD_WW_BW': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_WW', 'LD_BW'],
              'SPLD_WW_WB': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_WW', 'LD_WB'],
              'SPLD_BW_WB': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_BW', 'LD_WB'],
                
              'SPLD_BB_WW_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB', 'LD_WW', 'LD_BW'],
              'SPLD_BB_WW_WB': ['SP_BB', 'SP_WW', 'SP_BW', 'LD_BB', 'LD_WW', 'LD_WB'],
              
              'all': [P.Name for P in Sys.ForceField if not P is None]
            }

                  
# refine user selection
StageName = sys.argv[5]
cg.OptStageNames = [StageName]
cg.OptStages = {StageName: LDStageDict[StageName]}
cg.runSrel(Sys = Sys, ParamString = ffSP)
