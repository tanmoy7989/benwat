#!/usr/bin/env python

import os, sys, pickle, copy
import numpy as np

# dependencies
sys.path.append(os.path.expanduser('~/benwat')); import mysim
import sim
import pickleTraj
import parse_potential
import cgmodel as cg
import measurelib as lib

sys.path.append('~/')
from selLDCut import structcorr as sc

# flags and other globals
calcErrorBar = True ; NBlocks = 5
MeasureFreq = 1
Normalize = True
AtomNames2Types = True
Traj = None
Prefix = 'Measure'
isMappedTrj = True

def __isComputed(filename):
	if os.path.isfile(filename): return True
	else: return False 


def __AtomName2Type(AtomNames):
	AtomTypes = []
	[AtomTypes.append(int(float(x))) for x in AtomNames]
	return np.array(AtomTypes, np.int32)


def makeAllRDF():
	if isMappedTrj:
		rdfPairs = {'BB': ([1], [1]), 'BW': ([1], [2]), 'WW': ([2], [2])}
	else:
		rdfPairs = {'BB': (range(1,13), range(1,13)), 'BW': (range(1,13), [13]), 'WW': ([13], [13])}
	
	sc.LammpsTraj = Traj
	sc.Nbins = 100
	sc.TrjIter[2] = MeasureFreq
	sc.Normalize = Normalize
	sc.AtomNames2Types = AtomNames2Types
	if not isMappedTrj: sc.AtomNames2Types = False #(hotfix for a vmd screwup)

	sc.calcErrorBar = calcErrorBar
	sc.NBlocks = NBlocks
	
	for ext in ['BB', 'BW', 'WW']:
		print '\nCalculating rdf for %s\n\n' % ext
		sc.Prefix = Prefix + '_%s' % ext
		sc.genFileNames()
		if __isComputed(sc.rdfpickle): continue
		sc.makeRDF(CentAtomType = rdfPairs[ext][0], NeighAtomType = rdfPairs[ext][1])


def makeAllLD(LDCuts, LDDelta = 0.5):
	LDPairs = {'BB': (1,1), 'BW' : (1,2), 'WB': (2,1), 'WW' : (2,2)}
	if not type(LDCuts) is dict:
		raise TypeError("LDCuts must be a dict of the form {'BB': <>, 'BW': <>, 'WB': <>, 'WW': <>}")

	sc.LammpsTraj = Traj
	sc.Nbins = 100
	sc.TrjIter[2] = MeasureFreq
	sc.LDDelta = LDDelta
	sc.Normalize = Normalize
	sc.AtomNames2Types = AtomNames2Types

	sc.calcErrorBar = calcErrorBar
	sc.NBlocks = NBlocks

	for ext in ['BB', 'BW', 'WB', 'WW']:
		sc.LDCut = LDCuts[ext]
		sc.Prefix = Prefix + '_%s' % ext
		sc.genFileNames()
		if __isComputed(sc.ldpickle): continue
		sc.makeLD(LDPairs[ext][0], LDPairs[ext][1])


def makeCluster(Cut = None, ClustAtomType = 1):
	PickleName = Prefix + '.pickle'
	if __isComputed(PickleName): return

	Trj = pickleTraj(Traj)
	FrameRange = range(0, len(Trj), MeasureFreq)
	NFrames = len(FrameRange)
	BoxL = Trj.FrameData['BoxL']

	AtomNames = Trj.AtomNames
	if AtomNames2Types: Trj.AtomTypes = __AtomName2Type(Trj.AtomNames)
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

	pickle.dump((bin_val_measure, (bin_centers, bin_val_hist, err)), open(PickleName, 'w'))


def calcWaterCylinder():
	pickleName = Prefix + '_water_cylinder.pickle'
	if __isComputed(pickleName): return

	Trj = pickleTraj(Traj)
	BoxL = Trj.FrameData['BoxL']
	AtomTypes = __AtomName2Type(Trj.AtomNames)
	FrameRange = range(0, len(Trj), MeasureFreq)
	NFrames = len(FrameRange)
	
	minPos = np.zeros([NFrames, 3], np.float64)
	maxPos = np.zeros([NFrames, 3], np.float64)

	pb = sim.utility.ProgressBar(Text = 'Calculating min and max water coordinates...', Steps = NFrames)
	for i, frame in enumerate(FrameRange):
		Pos = Trj[frame]
		minPos_frame = [BoxL.min() * 10] * 3 
		maxPos_frame = [-BoxL.max() * 10] * 3

		minPos_frame, maxPos_frame = lib.water_cylinder_frame(pos = Pos, atomtypes = AtomTypes, wateratomtype = 2, boxl = BoxL,
															  minpos = minPos_frame, maxpos = maxPos_frame)

		minPos[i,:] = minPos_frame
		maxPos[i,:] = maxPos_frame

		pb.Update(i)

	pickle.dump((minPos, maxPos), open(pickleName, 'w'))

	# identify the cylinder axis
	v = maxPos - minPos
	diff = np.zeros(3)
	axs = np.array([0,1,2])
	for ax in axs: diff[ax] = np.abs(np.mean(v[:,ax] - Trj.FrameData['BoxL'][ax]))
	ind = np.argsort(diff)
	ax1 = axs[np.where(ind==1)[0][0]] ; ax2 = axs[np.where(ind==2)[0][0]]
	dsq = v[:,ax1]**2. + v[:, ax2]**2.
	r = np.mean(0.5 * np.sqrt(dsq))

	print '\n\nMean cylinder radius = ', r


def makeAllKBI(runAvg = True, delRDFPickle = True):
	pickleName = Prefix + '_KBI.pickle'
	if __isComputed(pickleName): return

	Trj = pickleTraj(Traj)
	BoxVol = np.prod(np.array(Trj.FrameData['BoxL']))
	if AtomNames2Types: AtomTypes = __AtomName2Type(Trj.AtomNames)
	else: AtomTypes = Trj.AtomTypes
	NB = len(np.where(AtomTypes == 1)[0])
	NW = len(np.where(AtomTypes == 2)[0])
	x_B = NB/(NB+NW)
	x_W = 1-x_B

	makeAllRDF()
	
	KBItypes = ['BB', 'WW', 'BW']
	G = {'BB': (), 'WW': (), 'BW': ()}
	Delta_N = copy.copy(G)
	rho_bulk = {'B': NB/BoxVol, 'W': NW/BoxVol}
	
	for KBItype in KBItypes:
		r,g,e,x = pickle.load(open('%s_%s_rdf.pickle' % (Prefix, KBItype), 'r'))
		
		R = r
		G[KBItype] = np.zeros([len(R), 2], np.float64)
		Delta_N[KBItype] = copy.copy(G[KBItype])
		dr = r[1]-r[0]
		for i, R_ in enumerate(R):
			G[KBItype][i,0] = R_
			x = np.array(r[:i]) ; y  = (4.*np.pi*r[:i]**2) * (g[:i]-1)
			G[KBItype][i,1] = np.sum(y*x*dr)
			# running avg between 1 and 1.4 nm to control oscillations (suggested by Pritam Ganguly)
			if runAvg:
				if R_ >= 12: G[KBItype][i,1] = np.mean(G[KBItype][i-5:i-1, 1])

			Delta_N[KBItype][i,0] = R_
			Delta_N[KBItype][i,1] = rho_bulk[KBItype[-1]] * G[KBItype][i,1]

		if delRDFPickle: os.remove('%s_%s_rdf.pickle' % (Prefix, KBItype))

	Delta_BW = G['BB'] + G['WW'] - 2*G['BW']
	dmudx_B = 1./(x_B * (1 + rho_bulk['B'] *x_W * Delta_BW))
	dgammadx_B = - (rho_bulk['W'] * x_B * Delta_BW) / (1 + rho_bulk['W'] * x_B * Delta_BW)

	output = (G, Delta_N, Delta_BW, dmudx_B, dgammadx_B)
	pickle.dump(output, open(pickleName, 'w'))
