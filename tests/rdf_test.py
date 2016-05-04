#!/usr/bin/env python

import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
from selLDCut import structcorr as rdf

rdf.Normalize = True
rdf.MeasureFreq = 50
rdf.Normalize = True
rdf.Nbins = 50
rdf.AtomNames2Types = True																																																																																																																																																																																																																																																																																																																																																																																																																																								
rdf.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250test/NB250NW250_prod.lammpstrj.gz')
rdf.Prefix = 'NB250_CG_BB'

rdf.calcErrorBar = False

rdf.genFileNames()
rdf.makeRDF([1], [1])

r,g,e,a  = pickle.load(open(rdf.rdfpickle, 'r'))
plt.errorbar(r,g, yerr = e)
plt.show()

