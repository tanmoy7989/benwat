#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import sim, pickleTraj, measure, benwatlib, cgmodel as cg

kB = 0.001987
TempSet = 300.
AtomDict = {'B': 1, 'W': 2}

## NOTE: need to set the following settings directly using measure
## StepFreq, NBins, NBlocks, AtomNames2Types

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


def make_RKBI(Traj, Prefix, NB, NW, RDFPrefix = None, method = 'UnCorrected'):
    ''' Families of formulae that directly integrate the rdf from definition.
    Corrects the g(r) or its integration based on heuristic formulae.'''
    ## Note: error bar calculation not supported currently
    
    print 'RKBI: ', method
    RKBIPickle = Prefix + '.pickle'
    if measure.__isComputed(RKBIPickle):
        return pickle.load(open(RKBIPickle, 'r')), RKBIPickle
    
    # extract Traj information
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    BoxL = np.array(measure.BoxL)
    NBins = 1000 # fine resolution bins since need to integrate
    measure.NBins = NBins
    
    # extract fine grained RDFs
    if RDFPrefix is None: RDFPrefix = 'NB%dNW%d' % (NB, NW)
    rdf_BB, rdf_BB_pickle = measure.makeRDF(1, 1, Prefix = RDFPrefix + '_rdf_BB_fine_grained')
    rdf_WW, rdf_WW_pickle = measure.makeRDF(2, 2, Prefix = RDFPrefix + '_rdf_WW_fine_grained')
    rdf_BW, rdf_BW_pickle = measure.makeRDF(1, 2, Prefix = RDFPrefix + '_rdf_BW_fine_grained')
    r_BB = rdf_BB[0] ; g_BB = rdf_BB[1] ; dr_BB = r_BB[1] - r_BB[0]
    r_WW = rdf_WW[0] ; g_WW = rdf_WW[1] ; dr_WW = r_WW[1] - r_WW[0]
    r_BW = rdf_BW[0] ; g_BW = rdf_BW[1] ; dr_BW = r_BW[1] - r_BW[0]
    
    # initialize all arrays
    h_BB = np.zeros(NBins) ; h_WW = np.zeros(NBins) ; h_BW = np.zeros(NBins)
    w_BB = lambda i: 1.0 ; w_WW = lambda i: 1.0 ; w_BW = lambda i: 1.0
    G_BB = np.zeros(NBins) ; G_WW = np.zeros(NBins) ; G_BW = np.zeros(NBins)
    
    if method == 'UnCorrected':
        # from definition of KBI for open systems
        h_BB = g_BB - 1
        h_WW = g_WW - 1
        h_BW = g_BW - 1
    
    elif method == 'TailCorrected':
        # Ganguly, van der Vegt et. al, JCTC, 2013, 9, 1347-1355, Eq (5)
        rho_B = float(NB) / np.prod(BoxL) ; rho_W = float(NW) / np.prod(BoxL)
        DeltaN_BB = rho_B * np.array( [ dr_BB * np.sum(4*np.pi * r_BB[:i]**2. * (g_BB[:i]-1)) for i in range(NBins) ] )
        DeltaN_WW = rho_W * np.array( [ dr_WW * np.sum(4*np.pi * r_WW[:i]**2. * (g_WW[:i]-1)) for i in range(NBins) ] )
        DeltaN_BW = rho_W * np.array( [ dr_BW * np.sum(4*np.pi * r_BW[:i]**2. * (g_BW[:i]-1)) for i in range(NBins) ] )
        Bulk_B = NB - rho_B * (4/3.)*np.pi * r_BB**3.
        Bulk_W = NW - rho_W * (4/3.)*np.pi * r_WW**3.
        h_BB = g_BB * Bulk_B / (Bulk_B - DeltaN_BB - 1) - 1
        h_WW = g_WW * Bulk_W / (Bulk_W - DeltaN_WW - 1) - 1
        h_BW = g_BW * Bulk_W / (Bulk_W - DeltaN_BW) - 1
        
    elif method == 'GeomCorrectedExact':
        # Krueger, Schenll et. al, J.Phys.Chem.Lett, 2013, 4, 235-238, Eqn (6)
        x_BB = [np.array([0])] + [r_BB[:i] / r_BB[:i][-1] for i in range(1,NBins)]
        x_WW = [np.array([0])] + [r_WW[:i] / r_WW[:i][-1] for i in range(1,NBins)]
        x_BW = [np.array([0])] + [r_BW[:i] / r_BW[:i][-1] for i in range(1,NBins)]
        w_BB = lambda i: (1 - 1.5 * x_BB[i] + 0.5 * x_BB[i]**3.)
        w_WW = lambda i: (1 - 1.5 * x_WW[i] + 0.5 * x_WW[i]**3.)
        w_BW = lambda i: (1 - 1.5 * x_BW[i] + 0.5 * x_BW[i]**3.)
        h_BB = g_BB - 1
        h_WW = g_WW - 1
        h_BW = g_BW - 1
    
    elif method == 'GeomCorrectedExtrapolated':
        # Krueger, Schenll et. al, J.Phys.Chem.Lett, 2013, 4, 235-238, Eqn (7)
        x_BB = [np.array([0])] + [r_BB[:i] / r_BB[:i][-1] for i in range(1,NBins)]
        x_WW = [np.array([0])] + [r_WW[:i] / r_WW[:i][-1] for i in range(1,NBins)]
        x_BW = [np.array([0])] + [r_BW[:i] / r_BW[:i][-1] for i in range(1,NBins)]
        w_BB = lambda i: (1 - x_BB[i]**3.)
        w_WW = lambda i: (1 - x_WW[i]**3.)
        w_BW = lambda i: (1 - x_BW[i]**3.)
        h_BB = g_BB - 1
        h_WW = g_WW - 1
        h_BW = g_BW - 1
    
    # integrate
    for n in range(1, NBins):
        G_BB[n] = dr_BB * np.sum(4*np.pi * r_BB[:n]**2 * h_BB[:n] * w_BB(n))
        G_WW[n] = dr_WW * np.sum(4*np.pi * r_WW[:n]**2 * h_WW[:n] * w_WW(n))
        G_BW[n] = dr_BW * np.sum(4*np.pi * r_BW[:n]**2 * h_BW[:n] * w_BW(n))
    
    ret = (r_BB, G_BB), (r_WW, G_WW), (r_BW, G_BW)
    pickle.dump(ret, open(RKBIPickle, 'w'))
    return ret, RKBIPickle
    

def make_SKBI(Traj, Prefix, NB, NW, L_corr = 12.0, NSubVols = 50, RandIter = 5000):
    ''' Calculates finite sized KBIs based on particle fluctuations in subvolumes
    These can be extrpolated to find the true KBI using the method as outlined in
    Cortes-Huerto, Kremer et.al , JCP, 2016, 145, 141103, Eq (1)'''
    
    print 'SKBI'
    SKBIPickle = Prefix + '.pickle'
    if measure.__isComputed(SKBIPickle):
        return pickle.load(open(SKBIPickle, 'r')), SKBIPickle
    
    # extract frame looping data
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    FrameRange = measure.FrameRange ; NFrames = measure.NFrames ; NBlocks = measure.NBlocks
    Trj = measure.Trj 
    BoxL = np.array(Trj.FrameData['BoxL'])
    AtomTypes = measure.AtomTypes
    
    # choose a grid of subvolume sizes greater than the correlation length
    L_corr = L_corr # 12A works OK for small systems and good for large systems
    NSubVols = NSubVols
    L0 = np.linspace(L_corr, BoxL.min(), NSubVols)
    
    # initialize all arrays
    G_BB_frame = np.zeros([NFrames, NSubVols]) ; G_WW_frame = np.zeros([NFrames, NSubVols]) ; G_BW_frame = np.zeros([NFrames, NSubVols])
    G_BB_block = np.zeros([NSubVols, NBlocks]) ; G_WW_block = np.zeros([NSubVols, NBlocks]) ; G_BW_block = np.zeros([NSubVols, NBlocks])
    G_BB = np.zeros(NSubVols) ; G_WW = np.zeros(NSubVols) ; G_BW = np.zeros(NSubVols)
    G_BB_err = np.zeros(NSubVols) ; G_WW_err = np.zeros(NSubVols) ; G_BW_err = np.zeros(NSubVols)
    
    # frame iteration
    # don't separate into blocks at this stage as parsing Trj in block mode will be very slow
    pb = sim.utility.ProgressBar('Calculating SKBI per frame...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        for j, this_L0 in enumerate(L0):
            G_BB_frame[i,j], G_WW_frame[i,j], G_BW_frame[i,j] = benwatlib.skbi(pos = Pos, boxl = BoxL, atomtypes = AtomTypes, atomtype_b = 1, atomtype_w = 2, l0 = this_L0, randiter = RandIter)
        pb.Update(i)
        
    # block average
    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        if NBlocks > 1: print 'Block: ', b
        start = b * BlockSize ; stop = (b+1) * BlockSize
        G_BB_block[:, b] = np.mean(G_BB_frame[start:stop, :], axis = 0)
        G_WW_block[:, b] = np.mean(G_WW_frame[start:stop, :], axis = 0)
        G_BW_block[:, b] = np.mean(G_BW_frame[start:stop, :], axis = 0)
    
    # output structure
    G_BB = np.mean(G_BB_block, axis = 1); G_WW = np.mean(G_WW_block, axis = 1) ; G_BW = np.mean(G_BW_block, axis = 1)
    if NBlocks > 1:
        G_BB_err = np.std(G_BB_block, axis = 1, ddof = 1)
        G_WW_err = np.std(G_WW_block, axis = 1, ddof = 1)
        G_BW_err = np.std(G_BW_block, axis = 1, ddof = 1)
    
    ret = ( (L0, G_BB, G_BB_err), (L0, G_WW, G_WW_err), (L0, G_BW, G_BW_err) )
    pickle.dump(ret, open(SKBIPickle, 'w'))
    return ret, SKBIPickle
        

def make_DensityProfile(Traj, Prefix, NB, NW, NSlice = 50):
    print 'Density profile along longest(z) axis'
    densityPickle = Prefix + '.pickle'
    if measure.__isComputed(densityPickle):
        return pickle.load(open(densityPickle, 'r')), densityPickle
    
    # extract frame looping data
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    FrameRange = measure.FrameRange ; NFrames = measure.NFrames; NBlocks = measure.NBlocks
    Trj = measure.Trj
    BoxL = np.array(measure.BoxL)
    AtomTypes = measure.AtomTypes
    
    # initialize all arrays
    rhoB_frame = np.zeros([NFrames, NSlice], np.float64) ; rhoW_frame = np.zeros([NFrames, NSlice])
    rhoB_block = np.zeros([NSlice, NBlocks]) ; rhoW_block = np.zeros([NSlice, NBlocks])
    rhoB = np.zeros(NSlice) ; rhoW = np.zeros(NSlice)
    errB = np.zeros(NSlice) ; errW = np.zeros(NSlice)
    
    # frame iteration
    # don't separate into blocks at this stage as parsing Trj in block mode will be very slow
    SubBoxVol = np.prod(BoxL) / NSlice # A^3
    z = np.linspace(0, BoxL[2], NSlice)
    pb = sim.utility.ProgressBar('Caclulating density along z axis...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        rhoB_frame[i,:], rhoW_frame[i,:] = benwatlib.slicedensity(pos = Pos, boxl = BoxL, atomtypes = AtomTypes, atomtype_b = 1, atomtype_w = 2, nslice = NSlice)
        pb.Update(i)
    
    # convert numbers to densities
    Mass_B = 78.11 ; Mass_W = 18.01 # molecular weights in g/mol
    factor = (10/6.023) # convert particle densities to g/cc
    rhoB_frame = (rhoB_frame * Mass_B / SubBoxVol) * factor
    rhoW_frame = (rhoW_frame * Mass_W / SubBoxVol) * factor
    
    # block average
    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        if NBlocks > 1: print 'Block: ', b
        start = b * BlockSize ; stop = (b+1) * BlockSize
        rhoB_block[:, b] = np.mean(rhoB_frame[start:stop], axis = 0)
        rhoW_block[:, b] = np.mean(rhoW_frame[start:stop], axis = 0)
    
    # output structure
    rhoB = np.mean(rhoB_block, axis = 1); rhoW = np.mean(rhoW_block, axis = 1)
    if NBlocks > 1:
        errB = np.std(rhoB_block, axis = 1, ddof = 1)
        errW = np.std(rhoW_block, axis = 1, ddof = 1)
    
    ret = ( (z, rhoB, errB), (z, rhoW, errW) )
    pickle.dump(ret, open(densityPickle, 'w'))
    return ret, densityPickle


def make_Widom(Traj, Prefix, NB, NW, ff, StepFreq = 1, RandIter = 100):
    ''' calculate excess chemical potential by random particle insertion'''
    WidomPickle = Prefix + '.pickle'
    if measure.__isComputed(WidomPickle):
        return pickle.load(open(WidomPickle, 'r')), WidomPickle
        
    # make trial system
    cg.NB = NB ; cg.NW = NW
    cg.LDCutBB = 7.5 ; cg.LDCutWW = 3.5
    cg.LDCutBW = 0.5 * (7.5+3.5) ; cg.LDCutWB = cg.LDCutBW
    SysB = cg.makeSys() ; SysB.ForceField.SetParamString(ff)
    SysW = cg.makeSys() ; SysW.ForceField.SetParamString(ff)
    
    # insert phantom particle
    testB = SysB.World[0] ; testW = SysW.World[1]
    SysB += testB ; SysW += testW
    
    # extract details of traj
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    FrameRange = measure.FrameRange ; NFrames = measure.NFrames ; NBlocks = measure.NBlocks
    Trj = measure.Trj
    BoxL = measure.BoxL ; BoxHi = Trj.FrameData['BoxHi'] ; BoxLo = Trj.FrameData['BoxLo']
    
    # frame iteration to calculate exp(-beta * Delta U)
    beta = 1. / (0.001987 * 300.0)
    exptermB = np.zeros(NFrames * RandIter) ; exptermW = np.zeros(NFrames * RandIter)
    LowerBound = 1.-3
    n = 0
    pb = sim.utility.ProgressBar('Caclulating test particle energies...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        SysB.Arrays.Pos[0: -1] = Pos ; SysW.Arrays.Pos[0: -1] = Pos
        for j in range(RandIter):
            testBPos = np.array([ BoxLo[i] + BoxL[i]*np.random.random() for i in [0,1,2] ])
            SysB.ForceField.Eval(ANumList = [NB+NW+1])
            exptermB[n] = np.exp(-beta * SysB.PEnergy)
            if not LowerBound is None and (np.isnan(exptermB[n]) or exptermB[n] < LowerBound): exptermB[n] = 0
        
            testWPos = np.array([ BoxLo[i] + BoxL[i]*np.random.random() for i in [0,1,2] ] )
            SysW.ForceField.Eval(ANumList = [NB+NW+1])
            exptermW[n] = np.exp(-beta * SysW.PEnergy)
            if not LowerBound is None and (np.isnan(exptermW[n]) or exptermW[n] < LowerBound): exptermW[n] = 0
            
            n += 1
        
    # TODO: block average
    measure.Normalize = False
    histB = measure.makeHist(exptermB) ; histW = measure.makeHist(exptermW)
    muB = -beta * np.log( np.sum(histB[1]) / np.sum(histB[0]) )
    muW = -beta * np.log( np.sum(histW[1]) / np.sum(histW[0]) )
    
    ret = (histB, histW, muB, muW)
    pickle.dump(ret, open(WidomPickle, 'w'))
    
    


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
    if len(sys.argv) > 7:
        ff_file = os.path.abspath(sys.argv)
        ff = file(ff_file).read()

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
    for method in ['UnCorrected', 'TailCorrected', 'GeomCorrectedExact', 'GeomCorrectedExtrapolated']:
        make_RKBI(LammpsTraj, '%s_RKBI_%s' % (FilePrefix, method), NB, NW, method = method)
    
    if not LammpsTraj.__contains__('SP') or NB == 10:
        measure.NBlocks = 4
        measure.StepFreq = 50
        make_SKBI(LammpsTraj, FilePrefix+'_SKBI', NB, NW)

    
    
    ### Deprecated calculations (don't yield meaningful results)
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 7.5, 'B')
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 3.5, 'W')
    
