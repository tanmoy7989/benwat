#!/usr/bin/env python

import numpy as np
import os, sys
import cPickle as pickle
import measure
import cgmodel as cg

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

def makeClusterHist(Traj, Prefix, Cut, ClustAtomType):
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
    print 'KBI '
    KBIPickle = Prefix + '.pickle'
    if measure.__isComputed(KBIPickle):
        return pickle.load(open(KBIPickle, 'r')), KBIPickle    

    # calculate densities
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    BoxL = measure.BoxL
    if isinstance(BoxL, list): BoxL = np.array(BoxL)
    BoxVol = np.prod(BoxL)
    rho_B = float(NB) / BoxVol ; rho_W = float(NW) / BoxVol
    
    # calculate rdfs
    rdfPrefix = Prefix.split('_KBI')[0]
    hist_BB, rdfpickle_BB = make_rdf(Traj, rdfPrefix + '_rdf_BB', 'BB')
    hist_WW, rdfpickle_WW = make_rdf(Traj, rdfPrefix + '_rdf_WW', 'WW')
    hist_BW, rdfpickle_BW = make_rdf(Traj, rdfPrefix + '_rdf_BW', 'BW')
    r_BB, g_BB, err_BB = hist_BB
    r_WW, g_WW, err_WW = hist_WW
    r_BW, g_BW, err_WW = hist_BW

    # calculate KB integrals
    def func_G(r, gof, L_corr):
        N_corr = len(L_corr)
        G = np.zeros(N_corr)
        dx = r[1] - r[0]
        for i in range(N_corr):
            try:
                CutInd = [list(r).index(this_r) for this_r in list(r) if this_r >= L_corr[i]][0]
            except IndexError:
                print i, L_corr[i]
                print r
                exit()
            x = r[:CutInd]
            y = (gof[:CutInd] - 1) * x**2.
            G[i] = 4 * np.pi * dx * np.sum(y)
        return G

    N = 50
    Rmax = BoxL[0]/2
    R = np.linspace(2.0, Rmax, N) # based on Nico's paper
    G_BB = func_G(r_BB, g_BB, R)
    G_WW = func_G(r_WW, g_WW, R)
    G_BW = func_G(r_BW, g_BW, R)
    
    # calculate excess coordination numbers
    N_BB = rho_B * G_BB ; N_WW = rho_W * G_WW; N_BW = rho_W * G_BW
    # calculate solvation parameter delta_BW
    Delta_BW = G_BB + G_WW - 2 * G_BW ; Delta_BW = np.mean(Delta_BW[-10:]) # average over last 10 values
    # calculate dmudx
    x_B =  float(NB) / (float(NB + NW)) ; x_W = 1.0 - x_B
    dmudx  = (kB*TempSet) * (x_B * (1 + rho_B * x_W * Delta_BW)) ** (-1.0)
    # calculate dgammadx
    dgammadx = - (rho_W * x_B * Delta_BW) * (1 + rho_W * x_B * Delta_BW) ** (-1.0)
    # output structure
    ret = {'R': R, 'G_BB': G_BB, 'G_WW': 'G_WW', 'G_BW': G_BW, 'N_BB': N_BB, 'N_WW': N_WW, 'N_BW': N_BW,
         'Delta_BW': Delta_BW, 'dmudx': dmudx, 'dgammadx': dgammadx}
    pickle.dump(ret, open(KBIPickle, 'w'))
    return ret, KBIPickle


### MAIN ###
OutputDir = os.path.abspath(sys.argv[1])
SysPrefix = sys.argv[2]
FilePrefix = os.path.join(OutputDir, SysPrefix)
LammpsTraj = sys.argv[3]
NB = int(sys.argv[4])
NW = int(sys.argv[5])
AtomNames2Types = False 
if len(sys.argv) > 6: AtomNames2Types = bool(sys.argv[6])

measure.AtomNames2Types = AtomNames2Types
measure.StepFreq = 10
measure.NBins = 50
measure.NBlocks = 4

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
#make_KBI(LammpsTraj, FilePrefix + '_KBI', NB, NW)
