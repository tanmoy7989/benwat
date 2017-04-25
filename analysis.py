#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import sim, pickleTraj, measure, cgmodel as cg
import SKBIlib as skbilib

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


def make_RKBI(r, gof):
    N = len(gof)
    G = np.zeros(N)
    dx = r[1] - r[0]
    for i in range(N):
        x = r[:i]
        y = (gof[:i] - 1) * x**2.
        G[i] = 4 * np.pi * dx * np.sum(y)
    # averaging the running integrals between 10 and 12 A
    inds = [list(r).index(x) for x in list(r) if x >= 10.0 and x <= 12.0]
    G_inf = np.mean(G[inds])
    return G, G_inf

def make_SKBI(Traj, AtomPair):
    Atomi = AtomDict[AtomPair[0]]
    Atomj = AtomDict[AtomPair[1]]
    Trj = pickleTraj(Traj)
    BoxL = Trj.FrameData['BoxL']
    
    StepFreq = 1
    FrameRange = range(0, len(Trj), StepFreq)
    NFrames = len(FrameRange)
    AtomTypes = np.array(Trj.AtomTypes, np.int32)
    for i in range(len(AtomTypes)):
        if not isinstance(AtomTypes[i], int):
            AtomTypes[i] = int(float(Trj.AtomNames[i]))
    
    L_s = np.arange(0.05, min(BoxL)/2, 0.01)
    Ni = np.zeros([NFrames, len(L_s)], np.int32)
    Nj = np.zeros([NFrames, len(L_s)], np.int32)
    G = np.zeros(len(L_s))
    for m, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        ni, nj = skbilib.skbi(pos = Pos, boxl = BoxL, cuts = L_s, atomtypes = AtomTypes,
                              atomtype_i = Atomi, atomtype_j = Atomj)
        Ni[m, :] = ni
        Nj[m, :] =  nj
    
    mu_i = np.mean(Ni, axis = 0)
    mu_j = np.mean(Nj, axis = 0)
    cov_ij = np.mean(Ni*Nj, axis = 0) - mu_i * mu_j
    V_s = (4/3.) * np.pi*L_s**3
    G = V_s * (cov_ij / (mu_i * mu_j) - int(Atomi==Atomj)/mu_j)
    return L_s, G
    
    
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

    # calculate RKBI
    RKBI_BB, RKBI_BB_inf = make_RKBI(r_BB, g_BB)
    RKBI_WW, RKBI_WW_inf = make_RKBI(r_WW, g_WW)
    RKBI_BW, RKBI_BW_inf = make_RKBI(r_BW, g_BW)
    R = r_BB
    
    #TODO: calculate SKBI and sanity check to make sure we are using the right KBI values
    # calculate SKBI
    #SKBI_BB_inf = make_SKBI(Traj, 'BB')
    #SKBI_WW_inf = make_SKBI(Traj, 'WW')
    #SKBI_BW_inf = make_SKBI(Traj, 'BW')
    
    
    # final KBI values to calculate thermodynamics
    # USING RKBI values for now
    G_BB, G_BB_inf = RKBI_BB, RKBI_BB_inf
    G_WW, G_WW_inf = RKBI_WW, RKBI_WW_inf
    G_BW, G_BW_inf = RKBI_BW, RKBI_BW_inf
    
    # calculate excess coordination numbers
    N_BB = rho_B * G_BB ; N_WW = rho_W * G_WW; N_BW = rho_W * G_BW
    N_BB_inf = rho_B * G_BB_inf ; N_WW_inf = rho_W * G_WW_inf; N_BW_inf = rho_W * G_BW_inf
    # calculate solvation parameter delta_BW
    Delta_BW = G_BB_inf + G_WW_inf - 2 * G_BW_inf
    # calculate dmudx
    x_B =  float(NB) / (float(NB + NW)) ; x_W = 1.0 - x_B
    dmudx  = (kB*TempSet) * (x_B * (1 + rho_B * x_W * Delta_BW)) ** (-1.0)
    # calculate dgammadx
    dgammadx = - (rho_W * x_B * Delta_BW) * (1 + rho_W * x_B * Delta_BW) ** (-1.0)
    # output structure
    ret = {'R': R, 
           'G_BB': G_BB, 'G_WW': G_WW, 'G_BW': G_BW,
           'G_BB_inf': G_BB_inf, 'G_WW_inf': G_WW_inf, 'G_BW_inf': G_BW_inf,
           'N_BB': N_BB, 'N_WW': N_WW, 'N_BW': N_BW,
           'N_BB_inf': N_BB_inf, 'N_WW_inf': N_WW_inf, 'N_BW_inf': N_BW_inf,
           'Delta_BW': Delta_BW, 'dmudx': dmudx, 'dgammadx': dgammadx}
    pickle.dump(ret, open(KBIPickle, 'w'))
    return ret, KBIPickle


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

    measure.AtomNames2Types = AtomNames2Types
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
