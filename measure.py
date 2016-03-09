#!/usr/bin/env python

import os, sys, pickle
import numpy as np

# dependencies
import sim
import pickleTraj
import parse_potential
import cgmodel as cg
sys.path.append('~/'); from selLDCut import *

# flags and other globals
doBlockAvg = True ; NBlocks = 5
MeasureFreq = 1
Normalize = True
useAtomTypes = True
Prefix = 'Measure'



def __AtomName2Type(AtomNames):
	AtomTypes = []
	[AtomTypes.append(int(x)) for x in AtomNames]
	return np.array(AtomTypes)


def getBoxL(Traj):
    trj = pickleTraj(Traj)
    return trj.FrameData['BoxL']


def gen_Srel_rc(LammpsTraj, LDCuts = []):
    # prepare the system
    cg.LammpsTraj = LammpsTraj
    cg.LDCutBW = LDCuts[0]
    cg.BoxL = getBoxL(LammpsTraj)
    Sys = cg.makeSys()
    
    # initialize the Bennett and FEP based Delta_Srel calculator
    SrelBennett.Sys = Sys
    SrelBennett.LammpsTraj = LammpsTraj
    SrelBennett.Prefix = Prefix
    SrelBennett.LDCuts = LDCuts
    SrelBennett.LDDelta = 1.2
    SrelBennett.genFileNames()
    
    # setup cg model and CG-MD simulations for different cutoffs
    SrelBennett.runCGMD()
    
    # start Bennett method calculations
    SrelBennett.Srel_rc()
    

def gen_fsw(LammpsTraj, LDCuts = [], NB = 250, NW = 250):
    # initialize the fsw calculator
    fsw.LammpsTraj = LammpsTraj
    fsw.NCentAtoms = NB ; fsw.NNeighAtoms = NW
    fsw.CentAtomType = 1 ; fsw.NeighAtomType = 2
    fsw.LDCuts = LDCuts ; fsw.LDDelta = 1.2
    fsw.Prefix = 'NB%dNW%d' % (NB, NW) ; fsw.genFileNames()

    #calculate rdf
    fsw.makeRDF()
    
    # calculate maximal correlation
    fsw.makeFSWCorrelation(rdfCut = 4.) 
    #based on NB250 B-W rdf (but unusual)
    fsw.calcCorrelation()
    

def makeCluster(Traj, Cut = None, Atom = 'B', NAtom = 250):
	Trj = pickleTraj(Traj)
	FrameRange = range(0, len(Trj), MeasureFreq)
	NFrames = len(FrameRange)
	BoxL = getBoxL(Traj)

	AtomNames = Trj.AtomNames
	if useAtomTypes: Trj.AtomTypes = __AtomName2Type(Trj.AtomNames)
	if Atom == 'B': ClustAtomType = 1
	else: ClustAtomType = 2
	Inds = np.where(Trj.AtomTypes == ClustAtomType)

	bin_centers = range(0, NAtom+1)
	bin_val_measure = np.zeros([NFrames, NAtom+1], np.float64)
	bin_val_hist = np.zeros(NAtom+1)
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
		err = np.zeros(NAtom + 1)
		bin_val_block = np.zeros([NBlocks, NAtom+1], np.float64)
		BlockLen = int(NFrames/NBlocks)
		count = 0
		pb = sim.utility.ProgressBar(Text = 'Calculating error', Steps = NBlocks)
		for block in range(NBlocks):
			bin_val_block[block] = np.mean(bin_val_measure[count:count+BlockLen], axis = 0)
			bin_val_block[block] /= np.sum(bin_val_block[block])
			pb.Update(block)
			count += BlockLen

		err = np.std(bin_val_block, axis = 0)	

	MeasurePickleName = Prefix + '.pickle'
	pickle.dump((bin_val_measure, (bin_centers, bin_val_hist, err)), open(MeasurePickleName, 'w'))