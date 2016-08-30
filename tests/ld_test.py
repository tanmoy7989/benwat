#!/usr/bin/env python

import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
from selLDCut import structcorr as ld

LDType = 'BW'

#choices for LDCuts for NB250NW250
#BB: 1st shell = 7.66, quite good
#WW: 1st shell = 3.5, 2nd shell = 5.75, 3rd shell = 8
#BW: 1st shell = 3.8, 2nd shell = 6.25, 3rd shell = 9

LDTypes = ['BB', 'WW', 'BW', 'WB']
LDCuts = dict(BB = 7.66 , WW = 3.5 , BW = 17, WB = 9)
LDAtomTypes = dict(BB = (1,1), WW = (2,2), BW = (1,2), WB = (2,1))

ld.Normalize = True
ld.MeasureFreq = 50
ld.Normalize = True
ld.Nbins = 100
ld.AtomNames2Types = True
ld.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
ld.Prefix = 'NB250NW250_AA_%s' % LDType

ld.calcErrorBar = False

ld.LDCut = LDCuts[LDType]
ld.LDDelta = 1.0

ld.genFileNames()
ld.makeLD(LDAtomTypes[LDType][0], LDAtomTypes[LDType][1])

r,h,e,x,ld  = pickle.load(open(ld.ldpickle, 'r'))
plt.plot(r,h)
plt.show()

