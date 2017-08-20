#!/ usr/bin/env python
import os, sys, numpy as np, cPickle as pickle, copy, subprocess
from scipy.spatial import cKDTree as Octree
import sim, parse_potential as pp, measure, pickleTraj
import cgmodel as cg

# user input (globals)
Traj = None
NB = None
NW = None
TPType = 'None'
Prefix = 'tpi'
NBins = (10, 10, 10)
measure.StepFreq = 10

# constants
Pad = 1e-3
kB = 0.001987
TempSet = 300.0
LammpsExec = 'lmp_mpich2'

# flags
Verbose = True
DelTempFiles = True

def make_Grid(BoxL, BoxLo, BoxHi, NBins = (10, 10, 10), visMap = False):
    ''' generate a 3D grid on minimum imaged atom positions
    across the box and an associated octree
    '''
    iBoxL = 1.0 / BoxL
    Lmin = BoxLo ; Lmax = Lmin + BoxL + np.array([Pad]*3)
    dx = (Lmax[0] - Lmin[0]) / float(NBins[0])
    dy = (Lmax[1] - Lmin[1]) / float(NBins[1])
    dz = (Lmax[2] - Lmin[2]) / float(NBins[2])
    x = np.array( [ Lmin[0] + (i+0.5)*dx for i in range(NBins[0]) ] )
    y = np.array( [ Lmin[1] + (i+0.5)*dy for i in range(NBins[1]) ] )
    z = np.array( [ Lmin[2] + (i+0.5)*dz for i in range(NBins[2]) ] )
    
    # initialize octrees
    X, Y, Z = np.meshgrid(x,y,z)
    M, N, P = np.meshgrid(range(NBins[0]), range(NBins[1]), range(NBins[2]))
    Tree = Octree( zip(X.ravel(), Y.ravel(), Z.ravel()) )
    
    # row major indexing from 3D meshgrid to Octree
    grid2tree = lambda xi, yj, zk : xi*NBins[1]*NBins[2] + yj*NBins[2] + zk
    
    return (Tree, grid2tree), (Lmin, Lmax), (dx, dy, dz), (X, Y, Z)


def make_EVM(visMap = False, Rcav = None, rVDW = None):
    ''' generate the excluded volume map of the 
    system i.e. find cavities and their locations
    '''
    global Traj, Prefix, TPType, NBins

    EVMPickle = Prefix + '_EVM.pickle'
    
    # van-der Waals radii
    if rVDW is None:
        rVDW = {'B': 1.77, 'W': 1.4} # ref: A.Bondi, 'van der Waals Volumes and Radii JPC, 68(3), 1994
    
    # extract Traj info
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    Trj = measure.Trj
    FrameRange = measure.FrameRange
    NFrames = measure.NFrames
    NAtoms = Trj.FrameData['NAtom']
    
    # initialize grid
    BoxL = Trj.FrameData['BoxL'] ; BoxLo = Trj.FrameData['BoxLo'] ; BoxHi = Trj.FrameData['BoxHi']
    iBoxL = 1.0 / BoxL
    # reset NBins to get an appropriate cavity size
    if NBins is None:
        NBins_each = int(BoxL/ (2*rVDW[TPType])) + 1
        NBins = (NBins_each, )*3
    
    (Tree, grid2tree), (Lmin, Lmax), (dx, dy, dz), (X, Y, Z) = make_Grid(BoxL, BoxLo, BoxHi, NBins = NBins)
    
    # initalize Map 
    # Map has same order of elements as Tree
    Map = np.zeros([ NFrames, NBins[0]*NBins[1]*NBins[2] ], np.uint8)
    
    if Verbose: print '%d X %d X %d min-imaged grid, dx = %g , dy = %g, dz = %g\n' \
           % (NBins[0], NBins[1], NBins[2], dx, dy, dz)
           
    # check if pickle already present
    if measure.__isComputed(EVMPickle): return
    
    # frame iteration
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
            if Trj.AtomNames[j] == 'B': pB.append(Tree.data[ind])
            if Trj.AtomNames[j] == 'W': pW.append(Tree.data[ind])
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
    pickle.dump(ret, open(EVMPickle, 'w'))


def make_TPITraj():
    ''' parses Traj to create new traj with test
    particles inserted
    '''
    global Traj, Prefix, TPType, NBins
    
    # check if traj already present
    if os.path.isfile(Prefix + '_tpi.lammpstrj.gz'): return
    
    # extract Traj info
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    Trj = measure.Trj
    BoxL = Trj.FrameData['BoxL'] ; BoxLo = Trj.FrameData['BoxLo'] ; BoxHi = Trj.FrameData['BoxHi']
    FrameRange = measure.FrameRange
    NFrames = measure.NFrames
    
    # get excluded volume map and tree
    ret = pickle.load(open(Prefix + '_EVM.pickle', 'r'))
    Grid, Map = ret
    
    # new traj properties
    NewTraj = Prefix + '_tpi.lammpstrj.gz'
    NewAtomName = TPType
    NewAtomType = 1 if TPType == 'B' else 2
    FrameData = copy.copy(Trj.FrameData)
    FrameData['AtomNames'] = Trj.AtomNames + [NewAtomName]
    FrameData['AtomTypes'] = np.append(Trj.AtomTypes, NewAtomType)
    FrameData['NAtom'] += 1
    NewTrj = sim.traj.LammpsWrite(NewTraj, MinImage = True)
    
    # frame iteration to insert test particles
    # (with minimaged coordinates)
    FrameCount = 0
    pb = sim.utility.ProgressBar('Inserting particles in cavities..', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        InsPos = Grid[np.logical_not(Map[i])]
        for j in range(len(InsPos)):
            Pos = Trj[i]
            Pos = np.append(Pos, [InsPos[j]], axis = 0)
            FrameData['TimeStep'] = FrameCount
            NewTrj.WriteFrame(Pos = Pos, FrameData = FrameData)
            FrameCount += 1
        
        pb.Update(i)
            
    if Verbose:
        print '\n\n'        
        print 'Frames in traj: ', NFrames
        print 'Wrote %d frames to new traj' % FrameCount
    NewTrj.Close()    
    
    
def make_TPI_Lammps(FF_file):
    ''' calculate excess chemical potential by random particle insertion''' 
    global Traj, Prefix, NB, NW, TPType
      
    # template for Lammps test particle insertion
    s_tpi = '''
    
#Thermostat
fix             thermostat all langevin %(TEMP)11.4e %(TEMP)11.4e %(LANGEVINDAMP)11.4e %(RANSEED)d

#group for inserted test particle
group           x id %(NMOL)d

#Calculate test particle energy
compute         pertpe x pe/atom
compute         tpe x reduce sum c_pertpe

#Dump test particle energy to file and screen
thermo_style    custom step atoms pe c_tpe
fix             write2file all ave/time 1 1 1 c_tpe file %(DELTAUFILE)s

#Rerun
rerun           %(TRAJ)s dump x y z

#Clear computes
uncompute       pertpe
uncompute       tpe
unfix           write2file
unfix           thermostat
    '''
    
    DeltaUFile = Prefix + '_tpi.dat'
    
    # check if ene data file already present
    if os.path.isfile(DeltaUFile): return
    
    # get BoxL to pass to Sys
    Trj = pickleTraj(Traj)
    BoxL = Trj.FrameData['BoxL']
    
    # get NewTraj
    NewTraj = Prefix + '_tpi.lammpstrj.gz'
    
    # make trial system
    if TPType == 'B': NB += 1
    else : NW += 1
    cg.NB = NB ; cg.NW = NW
    cg.LDCutBB = 7.5 ; cg.LDCutWW = 3.5
    cg.LDCutBW = 0.5 * (7.5+3.5) ; cg.LDCutWB = cg.LDCutBW
    Sys = cg.makeSys()
    Sys.TempSet = TempSet
    Sys.BoxL = BoxL
    
    # load in the forcefield
    pp.loadParam(Sys, FF_file)
    
    # parameters for Lammps
    d = {'NMOL'     : NB+NW,
         'TRAJ': NewTraj,
         'TEMP': Sys.TempSet,
         'LANGEVINDAMP': 1. / Sys.Int.Methods.VVIntegrate.LangevinGamma,
         'RANSEED': 84321,
         'DELTAUFILE': DeltaUFile}
    
    # create Lammps files and template for input script
    LammpsFiles = sim.export.lammps.MakeLammps(Sys, Prefix = Prefix)
    InFile = LammpsFiles[0]
    s_in = file(InFile).read()
    os.remove(InFile)
    s_in += s_tpi % d
    NewInFile = Prefix + '_tpi.in'
    NewLogFile = Prefix + '_tpi.log'
    with open(NewInFile, 'w') as of: of.write(s_in)
    if Verbose: print 'Running Lammps rerun...'
    p = subprocess.Popen('%s -in %s -log %s' % (LammpsExec, NewInFile, NewLogFile), shell = True, 
                          stdout = subprocess.PIPE, stderr = subprocess.STDOUT)
    p.communicate()
    
    # delete temporary lammps scripts
    if DelTempFiles:
        for i in LammpsFiles + [NewInFile, NewLogFile]:
                if os.path.isfile(i): os.remove(i)

    
def calc_mu():    
    ''' calculate the excess chemical potential 
    from completed Lammps run and extended volume map
    '''
    global Prefix 
    
    # filter for excess tpi energies
    UMAX = 0.0
    accept = lambda u : u <= UMAX
    
    # calculate excess expterm applying max/min filter for u
    DeltaUFile = Prefix + '_tpi.dat'
    DeltaU_raw = np.loadtxt(DeltaUFile)[:,1]
    if Verbose: print '\nAveraging %d snapshots' % len(DeltaU_raw)
    DeltaU = np.array( [x for x in DeltaU_raw if accept(x)] )
    beta = 1. / (kB * TempSet)
    expterm = np.exp(-beta * DeltaU)
    
    # calculate HS cavity probability
    EVMPickle = Prefix + '_EVM.pickle'
    ret = pickle.load(open(EVMPickle, 'r'))
    Grid, Map = ret
    NMapFrames = len(Map)
    pcav = 0.0 
    for i in range(NMapFrames):
        nfill = np.sum(Map[i])
        pfill = float(nfill) / float(NBins[0]*NBins[1]*NBins[2])
        pcav += (1-pfill)
    pcav /= NMapFrames
    
    # calculate excess mu
    mu1 = -(kB * TempSet) * np.log(pcav)
    mu2 = -(kB * TempSet) * np.log(np.mean(expterm))
    mu = mu1 + mu2
    
    # output in kJ/mol
    print "pcav = %g"  % pcav
    print "mu_cav = %g kcal/mol" % mu1
    print "mu_in_cav = %g kcal/mol" % mu2
    print "mu = %g kcal/mol = %g kJ/mol" % (mu, mu*4.184)
    





######## MAIN ########
if __name__ == '__main__':
        Traj = sys.argv[1]
        FF_file = sys.argv[2]
        NB = int(sys.argv[3])
        NW = int(sys.argv[4])
        TBType = sys.argv[5]
        Prefix = sys.argv[6]
        if len(sys.argv) > 7: NBins_each = int(sys.argv[7])
        else: NBins_each = None
        
        if not NBins_each is None: NBins = (NBins_each, )*3
        make_EVM()
        make_TPITraj()
        make_TPI_Lammps(FF_file)
        calc_mu()




