#!/usr/bin/env python

import numpy as np
import os, sys, pickle
import matplotlib.pyplot as plt

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim

def __prepLD():
	from selLDCut import structcorr as ld
	ld.TrjIter = [0,100,1]
	#ld.MeasureFreq = 50
	ld.Nbins = 50
	ld.AtomNames2Types = True
	ld.calcErrorBar = False
	return ld

def testInfCut():
	ld = __prepLD()
	ld.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
	ld.LDDelta = 1.0
	ld.LDCut = 1e3
	ld.Normalize = True

	ld.Prefix = 'infCut'
	ld.genFileNames()

	for i, ldtype in enumerate([(1,1), (2,2), (1,2), (2,1)]):
		print ldtype
		ld.makeLD(ldtype[0], ldtype[1])
		r,h,e,x,y  = pickle.load(open(ld.ldpickle, 'r'))
		ax = plt.subplot(2,2,i+1)
		ax.plot(r,h, label = 'C: %d, N: %d' % (ldtype[0], ldtype[1]))
		ax.legend()

	os.remove(ld.ldpickle)
	plt.show()


def testNCoord():
	
	rdfcut = 7.5

	ld = __prepLD()
	ld.LammpsTraj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
	ld.LDDelta = 1e-3
	ld.LDCut = rdfcut

	ld.Prefix = 'NCoord'
	ld.genFileNames()
	
	ld.Normalize = True ; ld.makeRDF([1], [1])
	r,g,e1,x1 = pickle.load(open(ld.rdfpickle, 'r'))
	ind = np.where(abs(r-rdfcut) < 0.1)[0][0]
	rho_bulk = (ld.NCentAtoms + ld.NNeighAtoms) / (np.mean(ld.BoxL)**3.)
	N_rdf = 4 * np.pi * rho_bulk * np.trapz(y = r[:ind+1]**2. * g[:ind+1], x = r[:ind+1])

	ld.Normalize = True ; ld.makeLD(1,1)
	rho, h, e2, x2, y2 = pickle.load(open(ld.ldpickle, 'r'))
	N_ld = np.trapz(y = h * rho, x = rho)
	#N_ld = np.mean(np.mean(y2, axis = 0))

	print '\n\n'
	print 'N_RDF = ', N_rdf
	print 'N_LD = ', N_ld

	ax = plt.subplot(1,2,1) ; ax.plot(r,g, label = r'$g(r)$') ; ax.legend()
	ax = plt.subplot(1,2,2) ; ax.plot(rho, h, label = r'$h(\rho)$') ; ax.legend()
	plt.show()

	os.remove(ld.rdfpickle)
	os.remove(ld.ldpickle)


#testInfCut()
testNCoord()
