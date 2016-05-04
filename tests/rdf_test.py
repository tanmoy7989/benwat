#!/usr/bin/env python

import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
from selLDCut import structcorr as rdf

rdf.Normalize = True
rdf.MeasureFreq = 10
rdf.Normalize = True
rdf.Nbins = 100
rdf.AtomNames2Types = True
rdf.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
rdf.Prefix = 'NB250'

rdf.calcErrorBar = True
rdf.NBlocks = 5

rdf.genFileNames()
rdf.makeRDF([1], [2])

r,g,e  = pickle.load(open(rdf.rdfpickle, 'r'))
plt.errorbar(r,g, yerr = e)
plt.show()

