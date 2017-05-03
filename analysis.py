#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import sim, pickleTraj, measure, KBI as kbi, benwatlib

kB = 0.001987
TempSet = 300.
AtomDict = {'B': 1, 'W': 2}

##### MEASURE FUNCTIONS #####
def make_rdf(Traj, Prefix, AtomPair):
    print 'RDF ', AtomPair[0], AtomPair[1]
    measure.LammpsTraj = Traj
    i = AtomPair[0]
    j = AtomPair[1]
    CentAtomType = AtomDict[i]
    NeighAtomType = AtomDict[j]
    return measure.makeRDF(CentAtomType, NeighAtomType, Prefix = Prefix)

def make_ld(Traj, Prefix, AtomPair, LDCut, **kwargs):
    print 'LD distribution ', AtomPair[0], AtomPair[1]
    measure.LammpsTraj = Traj
    LDDelta = kwargs.get('LDDelta', 1.0)
    i = AtomPair[0]
    j = AtomPair[1]
    CentAtomType = AtomDict[i]
    NeighAtomType = AtomDict[j]
    return measure.makeLDHist(CentAtomType, NeighAtomType, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)

def make_clust(Traj, Prefix, Cut, ClustAtom):
    print 'Cluster size distribution: ', ClustAtom
    ClustAtomType = AtomDict[ClustAtom]
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    clustpickle = Prefix + '.pickle'
    if measure.__isComputed(clustpickle):
        return ( pickle.load(open(clustpickle, 'r')), clustpickle)
    
    Trj = measure.Trj
    BoxL = measure.BoxL
    NFrames = measure.NFrames
    NBlocks = measure.NBlocks
    FrameRange = measure.FrameRange
    inds = np.where(measure.AtomTypes == ClustAtomType)
    NAtoms = len(inds[0])
    bin_vals_block = np.zeros([NAtoms+1, NBlocks], np.float64)
    clust_frame = np.zeros([NFrames, NAtoms+1])
    
    pb = sim.utility.ProgressBar(Text = '', Steps = NFrames)
    for frame_Ind, frame in enumerate(FrameRange):
        Pos = Trj[frame][inds]
        clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = BoxL, Cutoff = Cut)
        clustdist = np.array(clustdist)
        clust_frame[frame_Ind, :] = clustdist / np.sum(clustdist)
        pb.Update(frame_Ind)

    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        bin_vals_block[:, b] = np.mean(clust_frame[b*BlockSize:(b+1)*BlockSize, :], axis = 0)
    bin_centers = range(NAtoms+1)
    bin_vals = np.mean(bin_vals_block, axis = 1)
    bin_errs = np.std(bin_vals_block, axis = 1, ddof = 1)
    pickle.dump( (bin_centers, bin_vals, bin_errs), open(clustpickle, 'w'))
    return (bin_centers, bin_vals, bin_errs), clustpickle


def make_KBI(Traj, Prefix, NB, NW):
    print 'KBI'
    KBIPickle = Prefix + '.pickle'
    if measure.__isComputed(KBIPickle):
        return pickle.load(open(KBIPickle, 'r')), KBIPickle    
    
    kbi.LammpsTraj = Traj
    kbi.NB = NB ; kbi.NW = NW
    kbi.Prefix = Prefix
    print 'RKBI'
    methods = ['UnCorrected', 'TailCorrected', 'GeomCorrectedExact', 'GeomCorrectedExtrapolated']
    ret0 = {}
    for method in methods:
        ret = kbi.make_RKBI(method)
        ret0[method] = ret
    print 'SKBI'
    ret1 = kbi.make_SKBI()
    pickle.dump( (ret0, ret1), open(KBIPickle, 'w'))
    return (ret0, ret1), KBIPickle


def make_DensityProfile(Traj, Prefix, NB, NW, StepFreq = 1, NSlice = 50):
    densityPickle = Prefix + '.pickle'
    measure.LammpsTraj = Traj
    measure.StepFreq = StepFreq
    measure.__parseFrameData()
    FrameRange = measure.FrameRange ; NFrames = measure.NFrames
    Trj = measure.Trj
    BoxL = np.array(measure.BoxL)
    AtomTypes = [1]*NB + [2]*NW
    
    rhoB = np.zeros(NSlice) ; rhoW = np.zeros(NSlice)
    SubBoxVol = np.prod(BoxL) / NSlice # A^3
    z = np.linspace(0, BoxL[2], NSlice)
    
    pb = sim.utility.ProgressBar('Caclulating density along z axis...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        rhoB, rhoW = benwatlib.slicedensity(pos = Pos, boxl = BoxL, atomtypes = AtomTypes, atomtype_b = 1, atomtype_w = 2, nslice = NSlice)
        pb.Update(i)
    
    # convert numbers to densities
    Mass_B = 78.11 ; Mass_W = 18.01 # molecular weights in g/mol
    factor = (10/6.023) # convert particle densities to g/cc
    rhoB = (rhoB * Mass_B / SubBoxVol) * factor
    rhoW = (rhoW * Mass_W / SubBoxVol) * factor
    
    pickle.dump( (z, rhoB, rhoW), open(densityPickle, 'w'))



### MAIN ###
if __name__ == '__main__':
    OutputDir = os.path.abspath(sys.argv[1])
    SysPrefix = sys.argv[2]
    FilePrefix = os.path.join(OutputDir, SysPrefix)
    LammpsTraj = sys.argv[3]
    NB = int(sys.argv[4])
    NW = int(sys.argv[5])
    AtomNames2Types = False 
    if len(sys.argv) > 6: AtomNames2Types = bool(sys.argv[6])

    measure.AtomNames2Types = AtomNames2Types ; kbi.AtomNames2Types = AtomNames2Types
    measure.StepFreq = 10
    measure.NBins = 50
    measure.NBlocks = 4
    measure.Normalize = True

    LDCut_BB = 7.5
    LDCut_WW = 3.5
    LDCut_BW = 0.5 * (LDCut_BB + LDCut_WW)
    LDCut_WB = LDCut_BW

    print '\nComputing properties for ', SysPrefix
    print '------------------------------------------'
    make_rdf(LammpsTraj, FilePrefix + '_rdf_BB', 'BB')
    make_rdf(LammpsTraj, FilePrefix + '_rdf_WW', 'WW')
    make_rdf(LammpsTraj, FilePrefix + '_rdf_BW', 'BW')
    make_ld(LammpsTraj, FilePrefix + '_ld_BB', 'BB', LDCut_BB)
    make_ld(LammpsTraj, FilePrefix + '_ld_WW', 'WW', LDCut_WW)
    make_ld(LammpsTraj, FilePrefix + '_ld_BW', 'BW', LDCut_BW)
    make_ld(LammpsTraj, FilePrefix + '_ld_WB', 'WB', LDCut_WB)
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 7.5, 'B')
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 3.5, 'W')
    make_KBI(LammpsTraj, FilePrefix + '_KBI', NB, NW)
