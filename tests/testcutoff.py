#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import sim

sys.path.append('../')
import measure

NB = 50 ; NW = 450
LammpsTraj = '../data/gromacs/NB50NW450/NB50NW450_npt.lammpstrj.gz'
LDCuts = np.linspace(4.0, 10.0, 30)
FSWCut = 6.0

measure.rdf(NB = NB, NW = NW, LammpsTraj = LammpsTraj, Prefix = 'rdf_test')
measure.fsw(NB = NB, NW = NW, LammpsTraj = LammpsTraj, Prefix = 'fsw_test',
            LDCuts_BW = LDCuts, LDCuts_WB = LDCuts, FirstShellCut_BW = FSWCut,
            FirstShellCut_WB = FSWCut)
        
