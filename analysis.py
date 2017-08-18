#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle, subprocess, copy
from scipy.spatial import cKDTree as Octree
import sim, pickleTraj, measure, benwatlib, cgmodel as cg, parse_potential as pp

kB = 0.001987
TempSet = 300.
AtomDict = {'B': 1, 'W': 2}
LammpsExec = 'lmp_mpich2'

## NOTE: need to set the following settings directly using measure
## StepFreq, NBins, NBlocks, AtomNames2Types


##### RADIAL DISTRIBUTION FUNCTIONS #####
def make_rdf(Traj, Prefix, AtomPair):
    print 'RDF ', AtomPair[0], AtomPair[1]
    measure.LammpsTraj = Traj
    i = AtomPair[0]
    j = AtomPair[1]
    CentAtomType = AtomDict[i]
    NeighAtomType = AtomDict[j]
    return measure.makeRDF(CentAtomType, NeighAtomType, Prefix = Prefix)



##### LOCAL DENSITY DISTRIBUTION #####
def make_ld(Traj, Prefix, AtomPair, LDCut, **kwargs):
    print 'LD distribution ', AtomPair[0], AtomPair[1]
    measure.LammpsTraj = Traj
    LDDelta = kwargs.get('LDDelta', 1.0)
    i = AtomPair[0]
    j = AtomPair[1]
    CentAtomType = AtomDict[i]
    NeighAtomType = AtomDict[j]
    return measure.makeLDHist(CentAtomType, NeighAtomType, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)



##### CLUSTER SIZE DISTRIBUTION ##### (!! DEPRECATED !!)
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



##### KIRKWOOD BUFF INTEGRALS BASED ON RADIAL DISTRIBUTION FUNCTIONS #####
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
    NBlocks = 4
    measure.NBins = NBins
    measure.NBlocks = NBlocks
    measure.ExposeBlocks = True
    
    # extract fine grained RDFs
    if RDFPrefix is None: RDFPrefix = 'NB%dNW%d' % (NB, NW)
    rdf_BB, rdf_BB_pickle = measure.makeRDF(1, 1, Prefix = RDFPrefix + '_rdf_BB_fine_grained')
    rdf_WW, rdf_WW_pickle = measure.makeRDF(2, 2, Prefix = RDFPrefix + '_rdf_WW_fine_grained')
    rdf_BW, rdf_BW_pickle = measure.makeRDF(1, 2, Prefix = RDFPrefix + '_rdf_BW_fine_grained')
    r_BB = rdf_BB[0] ; g_BB_block = rdf_BB[-1] ; dr_BB = r_BB[1] - r_BB[0]
    r_WW = rdf_WW[0] ; g_WW_block = rdf_WW[-1] ; dr_WW = r_WW[1] - r_WW[0]
    r_BW = rdf_BW[0] ; g_BW_block = rdf_BW[-1] ; dr_BW = r_BW[1] - r_BW[0]
    
    # initialize all arrays
    h_BB_block = np.zeros([NBins, NBlocks]) ; h_WW_block = np.zeros([NBins, NBlocks]) ; h_BW = np.zeros([NBins, NBlocks])
    w_BB = lambda i: 1.0 ; w_WW = lambda i: 1.0 ; w_BW = lambda i: 1.0
    G_BB_block = np.zeros([NBins, NBlocks]) ; G_WW_block = np.zeros([NBins, NBlocks]) ; G_BW_block = np.zeros([NBins, NBlocks])
    G_BB = np.zeros(NBins) ; G_WW = np.zeros(NBins) ; G_BW = np.zeros(NBins)
    err_BB = np.zeros(NBins) ; err_WW = np.zeros(NBins) ; err_BW = np.zeros(NBins)
    
    if method == 'UnCorrected':
        # from definition of KBI for open systems
        h_BB_block = g_BB_block - 1
        h_WW_block = g_WW_block - 1
        h_BW_block = g_BW_block - 1
    
    elif method == 'TailCorrected':
        print "Block averaged version not implemented"
        exit()
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
        print "Block averaged version not implemented"
        exit()
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
        print "Block averaged version not implemented"
        exit()
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
    for b in range(NBlocks):
        for n in range(1, NBins):
            G_BB_block[n, b] = dr_BB * np.sum(4*np.pi * r_BB[:n]**2 * h_BB_block[:n, b] * w_BB(n))
            G_WW_block[n, b] = dr_WW * np.sum(4*np.pi * r_WW[:n]**2 * h_WW_block[:n, b] * w_WW(n))
            G_BW_block[n, b] = dr_BW * np.sum(4*np.pi * r_BW[:n]**2 * h_BW_block[:n, b] * w_BW(n))
    
    G_BB = np.mean(G_BB_block, axis = 1) ; G_WW = np.mean(G_WW_block, axis = 1) ; G_BW = np.mean(G_BW_block, axis = 1)
    if NBlocks > 1:
        err_BB = np.std(G_BB_block, ddof = 1, axis = 1)
        err_WW = np.std(G_WW_block, ddof = 1, axis = 1)
        err_BW = np.std(G_BW_block, ddof = 1, axis = 1)
        
    ret = (r_BB, G_BB, err_BB), (r_WW, G_WW, err_WW), (r_BW, G_BW, err_BW)
    pickle.dump(ret, open(RKBIPickle, 'w'))
    return ret, RKBIPickle
    


##### KIRKWOOD BUFF INTEGRALS BASED ON SMALL SUBVOLUMES #####
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
        


##### BULK DENSITY PROFILE FOR LONG BOXES #####
def make_DensityProfile(Traj, Prefix, NB, NW, NBins = 50, ReimageCOM = True, COM = None):
    print 'Density profile along longest(z) axis'
    densityPickle = Prefix + '.pickle'
    if measure.__isComputed(densityPickle):
        return pickle.load(open(densityPickle, 'r')), densityPickle
    
    # extract frame looping data
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    FrameRange = measure.FrameRange
    NFrames = measure.NFrames
    NBlocks = 4
    Trj = measure.Trj
    AtomTypes = [1]*NB + [2]*NW
    
    # initialize all arrays
    rhoB_frame = np.zeros([NFrames, NBins], np.float64) ;rhoW_frame = np.zeros([NFrames, NBins])
    rhoB_block = np.zeros([NBins, NBlocks]) ; rhoW_block = np.zeros([NBins, NBlocks])
    rhoB = np.zeros(NBins) ; rhoW = np.zeros(NBins)
    errB = np.zeros(NBins) ; errW = np.zeros(NBins)
    
    # number density to mass density conversion
    Mass_B = 78.11 ; Mass_W = 18.01 # molecular weights in g/mol
    factorAv = (10/6.023) # convert particle densities to g/cc
    rhoB_avg = 0.876
    rhoW_avg = 1.0
    
    # box details and bins
    BoxL = Trj.FrameData['BoxL']
    dz = (BoxL[2] + measure.HistPadding) / float(NBins)
    z = np.array([0.0 + (i+0.5)*dz for i in range(NBins)])
    SubBoxVol = np.prod(Trj.FrameData['BoxL']) / NBins # A^3
    
    if ReimageCOM and COM is None:
        # compute COM
        comB = np.zeros([3])
        comW = np.zeros([3])
        pb = sim.utility.ProgressBar('Computing COMs...', Steps = NFrames)
        for i, frame in enumerate(FrameRange):
            Pos = Trj[frame]
            Pos = sim.geom.Minimage(Pos, BoxL)
            PosB = Pos[0 : NB] ; PosW = Pos[NB : (NB+NW)]
            comB += np.mean(PosB, axis = 0)
            comW += np.mean(PosW, axis = 0)      
            pb.Update(i)
    
        comB /= NFrames ; comW /= NFrames
    
    if not COM is None:
        print 'Using supplied COM'
        comB = COM[0]
        comW = COM[1]
    
    # frame iteration to compute density
    # don't separate into blocks at this stage as parsing Trj in block mode will be very slow
    pb = sim.utility.ProgressBar('Caclulating density along z axis...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        # separate into B and W pos
        Pos = Trj[frame]
        PosB = Pos[0 : NB] ; PosW = Pos[NB : (NB+NW)]
        if ReimageCOM:
            # Reimage to the respective COMs
            PosB = sim.geom.Reimage(PosB, comB, BoxL)
            PosW = sim.geom.Reimage(PosW, comW, BoxL)
        # calculate number densities
        rhoB_frame[i,:] = benwatlib.zdensity(pos = PosB, boxl = BoxL, nbins = NBins)
        rhoW_frame[i,:] = benwatlib.zdensity(pos = PosW, boxl = BoxL, nbins = NBins)
        # convert to mass densities (g/cc)
        rhoB_frame[i,:] = (rhoB_frame[i,:] * Mass_B / SubBoxVol) * factorAv
        rhoW_frame[i,:] = (rhoW_frame[i,:] * Mass_W / SubBoxVol) * factorAv
        pb.Update(i)
    
    # block average
    print '\n'
    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        if NBlocks > 1: print 'Block: ', b
        start = b * BlockSize ; stop = (b+1) * BlockSize
        rhoB_block[:, b] = np.mean(rhoB_frame[start:stop], axis = 0)
        rhoW_block[:, b] = np.mean(rhoW_frame[start:stop], axis = 0)
    
    # output structure
    rhoB = np.mean(rhoB_block, axis = 1)
    rhoW = np.mean(rhoW_block, axis = 1)
    if NBlocks > 1:
        errB = np.std(rhoB_block, axis = 1, ddof = 1)
        errW = np.std(rhoW_block, axis = 1, ddof = 1)
    
    rhoB = np.mean(rhoB_frame, axis = 0)
    rhoW = np.mean(rhoW_frame, axis = 0)
    COM_ret = np.array([list(comB), list(comW)])
    ret = ( (z, rhoB, errB), (z, rhoW, errW), COM_ret)
    pickle.dump(ret, open(densityPickle, 'w'))
    return ret, densityPickle



##### INTERFACIAL SURFACE TENSION FOR LONG BOX #####
def make_SurfTens(Traj, Prefix, FF_File, NB = 380, NW = 1000, DelTempFiles = True):
    ''' calculate interfacial surface tension
    '''
    
    gammaPickle = Prefix + '.pickle'
    PressFile = Prefix + '_press.dat'
    if measure.__isComputed(gammaPickle):
        return pickle.load(open(gammaPickle, 'r')), gammaPickle
    
    s = '''
#set up integration and thermostat
fix             timeintegration all nve
fix             thermostat all langevin  3.0000e+02  3.0000e+02  1.0000e+02 18684

# dispay and dump pressure tensor
thermo_style    custom step c_thermo_press[1] c_thermo_press[2] c_thermo_press[3]
fix             record all ave/time %(WRITEFREQ)d 1 %(WRITEFREQ)d c_thermo_press[1] c_thermo_press[2] c_thermo_press[3] file %(PREFIX)s_press.dat

#run production
rerun           %(TRAJ)s dump x y z
    '''
    
    # make dict for filling in template
    d = {'PREFIX': Prefix, 'TRAJ': Traj, 'WRITEFREQ': 2500}
    
    # extract details from traj
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    Trj = measure.Trj
    FrameRange = measure.FrameRange ; NFrames = len(FrameRange)
    NBlocks = measure.NBlocks
    Lz = Trj.FrameData['BoxL'][2]
        
    if not os.path.isfile(PressFile):
        # make system
        cg.Prefix = Prefix
        cg.NB = NB ; cg.NW = NW
        cg.LDCutBB = 7.5 ; cg.LDCutWW = 3.5
        cg.LDCutBW = cg.LDCutWB = 0.5 * (7.5 + 3.5)
        Sys = cg.makeSys()
        Sys.BoxL = Trj.FrameData['BoxL']
    
        # load forcefield
        Sys.ForceField.SetParamString(file(FF_File).read())
    
        # make Lammps fles
        LammpsFiles = sim.export.lammps.MakeLammps(Sys = Sys, Prefix = Prefix, LammpsCommands = s % d)
        InFile = LammpsFiles[0]
    
        # run Lammps    
        print 'Rerunning Lammps...'
        DEVNULL = open(os.devnull, "w")
        p = subprocess.Popen('%s < %s' % (LammpsExec, InFile), shell = True, stdout = DEVNULL   , stderr = subprocess.PIPE)
        retout, reterr = p.communicate()
    
    # compute surface tension in mN/m
    factor = 0.01 # convert from atm-A to mN/m
    data = np.loadtxt(PressFile)
    pxx = data[:,1] ; pyy = data[:,2] ; pzz = data[:,3]
    gamma_frame = Lz * (pzz - 0.5*(pxx + pyy) )
    gamma_frame *= factor
    
    # block average
    gamma_block = np.zeros(NBlocks)
    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        start = b * BlockSize
        stop = (b+1) * BlockSize
        if b == NBlocks - 1: stop = len(Trj)-1
        gamma_block[b] = np.mean(gamma_frame[start:stop])
    
    gamma = np.mean(gamma_block)
    err = 0.0
    if NBlocks > 1: err = np.std(gamma_block, ddof = 1)
    
    print 'Gamma = %g mN/m, err = %g mN/m' % (gamma, err)
    
    # del temp files
    if DelTempFiles: os.system('rm *lammps*')
            
    # output
    ret = (gamma, err)
    pickle.dump(ret, open(gammaPickle, 'w'))
    return ret, gammaPickle



##### HARD SPHERE CHEMICAL POTENTIAL #####
def make_Grid(Box, Nx = 20, Ny = 20, Nz = 20):
    ''' generate a 3D grid on minimum imaged atom positions
    across the box and an associated octree
    '''
    BoxLo, BoxHi, BoxL = Box
    iBoxL = 1.0 / BoxL
    Lmin = BoxLo ; Lmax = Lmin + BoxL + np.array([measure.HistPadding]*3)
    dx = (Lmax[0] - Lmin[0]) / float(Nx)
    dy = (Lmax[1] - Lmin[1]) / float(Ny)
    dz = (Lmax[2] - Lmin[2]) / float(Nz)
    x = np.array( [ Lmin[0] + (i+0.5)*dx for i in range(Nx) ] )
    y = np.array( [ Lmin[1] + (i+0.5)*dy for i in range(Ny) ] )
    z = np.array( [ Lmin[2] + (i+0.5)*dz for i in range(Nz) ] )
    # initialize octrees
    X, Y, Z = np.meshgrid(x,y,z)
    M, N, P = np.meshgrid(range(Nx), range(Ny), range(Nz))
    Tree = Octree( zip(X.ravel(), Y.ravel(), Z.ravel()) )
    # row major indexing from 3D meshgrid to Octree
    grid2tree = lambda xi, yj, zk : xi*Nx*Ny + yj*Nz + zk
    return (Tree, grid2tree), (Lmin, Lmax), (dx, dy, dz), (X, Y, Z)

def make_EVM(Traj, Prefix, NB, NW, N = 20, Rcav = None, storeMap = False):
    ''' generate the excluded volume map of the 
    system i.e. find cavities and their locations
    '''
    print '\nComputing cavities in box for Rcav = %g...' % Rcav
    # pickling large amounts of data takes a lot of extra time so no need
    if storeMap:
        EVMPickle = Prefix + '_EVM.pickle'
        if measure.__isComputed(EVMPickle):
            return (pickle.load(open(EVMPickle, 'r')), EVMPickle)
        
    # van-der Waals radii
    rVDW = {'B': 1.77, 'W': 1.52} # ref: A.Bondi, 'van der Waals Volumes and Radii JPC, 68(3), 1994
    
    # extract Traj info
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    Trj = measure.Trj
    FrameRange = measure.FrameRange
    NFrames = measure.NFrames
    NAtoms = Trj.FrameData['NAtom']
    
    # initialize grid
    Nx = Ny = Nz = N
    BoxL = Trj.FrameData['BoxL'] ; BoxLo = Trj.FrameData['BoxLo'] ; BoxHi = Trj.FrameData['BoxHi']
    iBoxL = 1.0 / BoxL
    Box = (BoxLo, BoxHi, BoxL)
    (Tree, grid2tree), (Lmin, Lmax), (dx, dy, dz), (X, Y, Z) = make_Grid(Box, Nx, Ny, Nz)
    
    # initalize Map (Map has same order of elements as Tree)
    Map = np.zeros([ NFrames, Nx*Ny*Nz ], np.uint8)
    
    # organize atomnames since there
    # are inconsistencies between AA and CG
    AtomNames = ['B']*NB + ['W']*NW
    
    # frame iteration
    print '%d X %d X %d min-imaged grid, dx = %g , dy = %g, dz = %g\n' % (Nx, Ny, Nz, dx, dy, dz)    
    pb = sim.utility.ProgressBar('Locating cavities..', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        # locate all sites occupied by atoms
        pB = [] ; pW = []
        for j in range(NAtoms):
            Posj = Pos[j]
            dPosj = Posj - Lmin
            dPosj -= BoxL * (np.around(dPosj * iBoxL - 0.5) + 0.5)
            xi = int(dPosj[0]/dx) ; yj = int(dPosj[1]/dy) ; zk = int(dPosj[2]/dz)
            ind = grid2tree(xi, yj, zk)
            # record B and W coordinates separately
            if AtomNames[j] == 'B': pB.append(Tree.data[ind])
            if AtomNames[j] == 'W': pW.append(Tree.data[ind])
            # mark occupied sites 
            Map[i, ind] = 1
        
        # locate excluded volumes using octree nearest neighbor search
        # for efficient search use tree-tree search
        Bneigh = []
        Wneigh = []
        if pB:
            BTree = Octree(pB)
            if not Rcav is None: r = rVDW['B'] + Rcav
            else: r = 2 * rVDW['B']
            Bneigh = BTree.query_ball_tree(Tree, r = r)      
        if pW:
            WTree = Octree(pW)
            if not Rcav is None: r = rVDW['W'] + Rcav
            else: r = 2 *rVDW['W']
            Wneigh = WTree.query_ball_tree(Tree, r = r)
        allneigh = []
        for neigh in Bneigh+Wneigh: allneigh.extend(neigh)
        allneigh = list(set(allneigh))   
        # mark excluded volume sites
        for j in allneigh: Map[i, j] = 1
        
        pb.Update(i)
        
    # output
    ret = (Tree.data, Map)
    if storeMap: 
        pickle.dump(ret, open(EVMPickle, 'w'))
        return ret, EVMPickle
    else:
        return ret

def make_HardSphereMu_1(Traj, Prefix, NB, NW, N = 20):
    ''' computes Hard sphere mu vs cavity size
    for a range of cavity sizes
    '''
    
    HSPickle = Prefix + '.pickle'
    if measure.__isComputed(HSPickle):
        return (pickle.load(open(HSPickle, 'r')), HSPickle)
        
    mu = []
    GridSize = N**3.
    RcavList = np.array([0, 0.5, 1., 1.5, 2, 2.5, 3, 3.5])
    for i, Rcav in enumerate(RcavList):
        # compute EVM
        thisPrefix = Prefix + '_Rcav_%d' % i
        ret = make_EVM(Traj, Prefix = thisPrefix, NB = NB, NW = NW, N = N, Rcav = Rcav)
        # compute HS mu
        print '\nComputing HS mu_ex...'
        Grid, Map = ret
        pfill = np.mean(np.sum(Map, axis = 1) / float(GridSize))
        pcav = 1.0 - pfill
        mu.append( -(kB*TempSet) * np.log(pcav) )
    mu = np.array(mu)
    
    ret = (RcavList, mu)
    pickle.dump(ret, open(HSPickle, 'w'))
    return ret, HSPickle

def make_HardSphereMu(Traj, Prefix, NB, NW, Algorithm = 'Random', RanIter = None, NBins = None, RanSeed = 84321, Rcavlist = None):    
    ''' hard sphere mu with random or grid
    based insertions
    '''
    
    HSPickle = '%s_%s.pickle' % (Prefix, Algorithm)
    if measure.__isComputed(HSPickle):
        return (pickle.load(open(HSPickle, 'r')), HSPickle)
    
    # set the random number seed
    np.random.seed(RanSeed)
    
    # get Traj box data
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    Trj = measure.Trj
    FrameRange = measure.FrameRange
    NFrames = len(FrameRange)
    BoxL = Trj.FrameData['BoxL']
    BoxHi = Trj.FrameData['BoxLo']
    BoxLo = Trj.FrameData['BoxHi']
    AtomTypes = [1]*NB + [2]*NW
    
    # number of insertion attemps and one time grid pos
    if Algorithm == 'Random':
        if RanIter is None: RanIter = 10000
        NIns = RanIter * NFrames
    if Algorithm == 'GridSearch':
        if NBins is None: NBins = 20
        NIns = (NBins**3) * NFrames
        dx = BoxL[0] / float(NBins)
        dy = BoxL[1] / float(NBins)
        dz = BoxL[2] / float(NBins)
        x = np.array([BoxLo[0] + (i+0.5)*dx for i in range(NBins)])
        y = np.array([BoxLo[1] + (i+0.5)*dy for i in range(NBins)])
        z = np.array([BoxLo[2] + (i+0.5)*dz for i in range(NBins)])
        X, Y, Z = np.meshgrid(x,y,z)
        Grid = np.array(zip(X.ravel(), Y.ravel(), Z.ravel()))
        
    # calculate insertion coordinates
    def genInsPos():
        if Algorithm == 'Random':
            x = np.random.uniform(BoxLo[0], BoxHi[0], RanIter)
            y = np.random.uniform(BoxLo[1], BoxHi[1], RanIter)
            z = np.random.uniform(BoxLo[2], BoxHi[2], RanIter)
            return np.array(zip(x,y,z))
        if Algorithm == 'GridSearch':
            return Grid
               
    # declare all arrays
    if Rcavlist is None:
        Rcavlist = np.array([0, 0.5, 1., 1.5, 2, 2.5, 3, 3.5])
    N = len(Rcavlist)
    mu = np.zeros(N)
    
    # iterate over cavity radii
    for i, Rcav in enumerate(Rcavlist):
        print '\nRcav = %g' % Rcav
        Ncav = 0.0
        # frame iteration
        pb = sim.utility.ProgressBar(Text = 'Trying insertions...', Steps = NFrames)
        for n, frame in enumerate(FrameRange):
            Pos = Trj[frame]
            InsPos = genInsPos()
            Ncav += benwatlib.gridinsert(pos = Pos, boxl = BoxL, inspos = InsPos, 
                                         atomtypes = AtomTypes, atomtype_b = 1, 
                                         atomtype_w = 2, rcav = Rcav)
            pb.Update(n)
        # calculate insertion probability
        pcav = float(Ncav) / float(NIns)
        # calculate mu for this radii
        mu[i] = -(kB * TempSet) * np.log(pcav)
    
    # output
    ret = (Rcavlist, mu)
    pickle.dump(ret, open(HSPickle, 'w'))
    return ret, HSPickle
        




                
###### MAIN ######
if __name__ == '__main__':
    
    # user input
    OutputDir = os.path.abspath(sys.argv[1])
    SysPrefix = sys.argv[2]
    FilePrefix = os.path.join(OutputDir, SysPrefix)
    LammpsTraj = sys.argv[3]
    NB = int(sys.argv[4])
    NW = int(sys.argv[5])
    AtomNames2Types = False 
    if len(sys.argv) > 6: AtomNames2Types = bool(sys.argv[6])
    if len(sys.argv) > 7:
        fftype = sys.argv[7]
        if not fftype == 'AA':
            ff_file = os.path.abspath('/home/cask0/home/tsanyal/benwat/data/cgff/NB250NW250/control/NB250NW250_%s_ff.dat' % fftype)
            ff = file(ff_file).read()
    
    # histogram settings
    measure.AtomNames2Types = AtomNames2Types
    measure.StepFreq = 10
    measure.NBins = 50
    measure.NBlocks = 4
    measure.Normalize = True
    
    # LD settings
    LDCut_BB = 7.5
    LDCut_WW = 3.5
    LDCut_BW = 0.5 * (LDCut_BB + LDCut_WW)
    LDCut_WB = LDCut_BW
    
    # RKBI methods
    RKBI_Methods = ['UnCorrected', 'TailCorrected', 'GeomCorrectedExact', 'GeomCorrectedExtrapolated']
    #RKBI_Methods = ['UnCorrected']
    
    # small systems
    print '\nComputing properties for ', SysPrefix
    print '------------------------------------------'
    make_rdf(LammpsTraj, FilePrefix + '_rdf_BB', 'BB')
    make_rdf(LammpsTraj, FilePrefix + '_rdf_WW', 'WW')
    make_rdf(LammpsTraj, FilePrefix + '_rdf_BW', 'BW')
    
    make_ld(LammpsTraj, FilePrefix + '_ld_BB', 'BB', LDCut_BB)
    make_ld(LammpsTraj, FilePrefix + '_ld_WW', 'WW', LDCut_WW)
    make_ld(LammpsTraj, FilePrefix + '_ld_BW', 'BW', LDCut_BW)
    make_ld(LammpsTraj, FilePrefix + '_ld_WB', 'WB', LDCut_WB)
    
    #for method in RKBI_Methods:
    #    make_RKBI(LammpsTraj, '%s_RKBI_%s' % (FilePrefix, method), NB, NW, method = method)
    
    make_HardSphereMu(LammpsTraj, FilePrefix + '_HS', NB, NW, N = 30)
    
    
    ### Deprecated calculations (don't yield meaningful results)
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 7.5, 'B')
    #make_clust(LammpsTraj, FilePrefix + '_clust_W', 3.5, 'W')
    
    #if not LammpsTraj.__contains__('SP') or [10, 19, 38, 57].__contains__(NB):
    #    measure.NBlocks = 4
    #    measure.StepFreq = 50
    #    make_SKBI(LammpsTraj, FilePrefix+'_SKBI', NB, NW)
