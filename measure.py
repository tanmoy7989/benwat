#!/usr/bin/env python

import os, sys, pickle
import numpy as np

# dependencies
import sim
import pickleTraj
import parse_potential
import cgmodel as cg

sys.path.append('~/')
from selLDCut import rdf

# flags and other globals
doBlockAvg = True ; NBlocks = 5
MeasureFreq = 1
Normalize = True
useAtomTypes = True
Prefix = 'Measure'



def __AtomName2Type(AtomNames):
	AtomTypes = []
	[AtomTypes.append(int(float(x)) for x in AtomNames]
	return np.array(AtomTypes, np.int32)


def makeAllRDF(Traj, TrajType = 'CG', CentAtomTypes = None, NeighAtomTypes = None):
	if TrajType == 'CG':
		rdfPairs = {'BB': ([1], [1]), 'BW': ([1], [2]), 'WW': ([2], [2])}
	else:
		rdfPairs = {'BB': (range(1,13), range(1,13)), 'BW': (range(1,13), [13]), 'WW': ([13], [13])}
	
	rdf.LammpsTraj = Traj
	rdf.Nbins = 50
	for ext in ['BB', 'BW', 'WW']:
		print '\nCalculating rdf for %s\n\n' % ext
		rdf.Prefix = Prefix + '_%rdf_%s' % ext
		rdf.genFileNames()
		rdf.makeRDF(CentAtomType = rdfPairs[ext][0], NeighAtomTypes = rdfPairs[ext][1])


def makeCluster(Traj, Cut = None, ClustAtomType = 1):
	Trj = pickleTraj(Traj)
	FrameRange = range(0, len(Trj), MeasureFreq)
	NFrames = len(FrameRange)
	BoxL = Trj.FrameData['BoxL']

	AtomNames = Trj.AtomNames
	if useAtomTypes: Trj.AtomTypes = __AtomName2Type(Trj.AtomNames)
	Inds = np.where(Trj.AtomTypes == ClustAtomType)
	NAtom = len(Inds)

	bin_centers = range(1, NAtom+1)
	bin_val_measure = np.zeros([NFrames, NAtom], np.float64)
	bin_val_hist = np.zeros(NAtom)
	err = None

	pb = sim.utility.ProgressBar(Text = 'Calculating cluster dist', Steps = NFrames)
	count = 0
	for frame in FrameRange:
		Pos = Trj[frame][Inds]
		clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = BoxL, Cutoff = Cut)
		bin_val_measure[count, :] = np.array(clustdist)
		pb.Update(count)
		count += 1

	bin_val_hist = np.mean(bin_val_measure, axis = 0)
	if Normalize: bin_val_hist /= np.sum(bin_val_hist)

	if doBlockAvg:
		print '\n\n'
		err = np.zeros(NAtom)
		bin_val_block = np.zeros([NBlocks, NAtom], np.float64)
		BlockLen = int(NFrames/NBlocks)
		count = 0
		pb = sim.utility.ProgressBar(Text = 'Calculating error bars', Steps = NBlocks)
		for block in range(NBlocks):
			bin_val_block[block] = np.mean(bin_val_measure[count:count+BlockLen], axis = 0)
			if Normalize: bin_val_block[block] /= np.sum(bin_val_block[block])
			pb.Update(block)
			count += BlockLen

		err = np.std(bin_val_block, axis = 0)	

	PickleName = Prefix + '.pickle'
	pickle.dump((bin_val_measure, (bin_centers, bin_val_hist, err)), open(PickleName, 'w'))
