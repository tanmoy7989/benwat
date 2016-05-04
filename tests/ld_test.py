#!/usr/bin/env python

import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
from selLDCut import structcorr as ld

ld.Normalize = True
ld.MeasureFreq = 10
ld.Normalize = True
ld.Nbins = 100
ld.AtomNames2Types = True
ld.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
ld.Prefix = 'NB250'

ld.calcErrorBar = False
ld.NBlocks = 5

ld.LDCut = 8.0
ld.LDDelta = 0.5

ld.genFileNames()
ld.makeLD(1,2)

r,h,e,x  = pickle.load(open(ld.ldpickle, 'r'))
plt.plot(r,h)
plt.show()

