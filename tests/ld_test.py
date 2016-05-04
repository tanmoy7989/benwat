#!/usr/bin/env python

import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
from selLDCut import structcorr as ld

ld.Normalize = True
ld.MeasureFreq = 50
ld.Normalize = True
ld.Nbins = 100
ld.AtomNames2Types = True
ld.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
ld.Prefix = 'NB250NW250_AA_BB'

ld.calcErrorBar = False

ld.LDCut = 16.0
ld.LDDelta = 1.0

ld.genFileNames()
ld.makeLD(2,2)

r,h,e,x,ld  = pickle.load(open(ld.ldpickle, 'r'))
plt.plot(r,h)
plt.show()

