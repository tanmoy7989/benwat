#usr/bin/env python

import os, sys
import numpy as np

import sim
import pickleTraj

# System description
Prefix = None
TempSet = 300
NB = 250
NW = 250
Name_W = 'W' ; Name_B = 'B'
Mass_W = 18.01 ; Mass_B = 78.11
Dia_W = 2.8 ; Dia_B = 5.3

# Pair potential settings
SPCutScale = 2.5
NSPKnots = 30
ManageInnerCore = True # True to capture inner core tails of corresponding histograms

# Local density potential settings
LDCutWW = None
LDCutBB = None
LDCutBW = None
LDCutWB = None
NLDKnots = 30
RhoMin = 0
RhoMax = 20
LD_Delta = 1.0

# MD settings
MinSteps = 1000
NPTSteps = 1000
EquilSteps = 2000000
ProdSteps = 2000000
StepFreq = 500

# Srel optimization settings
OptStages = None
OptStageNames = None
LammpsTraj = None
LangevinGamma = 0.01

# Lammps settings
sim.export.lammps.LammpsExec = os.path.expanduser('~/mysoftware/tanmoy_lammps/lammps-15May15/src/lmp_ZIN')
sim.export.lammps.InnerCutoff = 0.02
sim.srel.base.DiffEneFracTol = 0.1 #increase tolerance to prevent sim-Lammps mismatch blowing up the srel run

def makeSys():
    global NB, NW, BoxL
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW
    global LammpsTraj, Prefix

    # read trajectory
    if LammpsTraj: Trj = pickleTraj(LammpsTraj)

    # system chemistry
    if Prefix is None: Prefix = 'NB%dNW%d' % (NB, NW)
    AtomTypeW = sim.chem.AtomType(Name_W, Mass = Mass_W, Charge = 0.0, Color = (0,0,1))
    AtomTypeB = sim.chem.AtomType(Name_B, Mass = Mass_B, Charge = 0.0, Color = (1,1,0))
    MolTypeW = sim.chem.MolType(Name_W, [AtomTypeW])
    MolTypeB = sim.chem.MolType(Name_B, [AtomTypeB])
    World = sim.chem.World([MolTypeB, MolTypeW], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = Prefix)
    for i in range(NB): Sys += MolTypeB.New()
    for i in range(NW): Sys += MolTypeW.New()

    # box lengths and cutoffs
    SPCutWW = SPCutScale * Dia_W
    SPCutBB = SPCutScale * Dia_B
    SPCutBW = SPCutScale * 0.5 * (Dia_B + Dia_W)

    # atom selection filters
    FilterWW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW])
    FilterBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB])
    FilterBW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB])
    FilterWW_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW], Ordered = True)
    FilterBB_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], Ordered = True)
    FilterBW_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeW], Ordered = True)
    FilterWB_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB], Ordered = True)

    # forcefield design
    SP_WW = None ; SP_BB = None; SP_BW = None
    LD_WW = None; LD_BB = None; LD_BW = None; LD_WB = None
    SP = sim.potential.PairSpline
    LD = sim.potential.LocalDensity
    SP_BB = SP(Sys, Cut = SPCutBB, NKnot = NSPKnots, Filter = FilterBB, Label = "SP_BB")
    SP_WW = SP(Sys, Cut = SPCutWW, NKnot = NSPKnots, Filter = FilterWW, Label = "SP_WW")
    SP_BW = SP(Sys, Cut = SPCutBW, NKnot = NSPKnots, Filter = FilterBW, Label = "SP_BW")
    if not LDCutWW is None: LD_WW = LD(Sys, Cut = LDCutWW, LowerCut = LDCutWW - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WW", Filter = FilterWW_ordered)
    if not LDCutBB is None: LD_BB = LD(Sys, Cut = LDCutBB, LowerCut = LDCutBB - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BB", Filter = FilterBB_ordered)
    if not LDCutBW is None: LD_BW = LD(Sys, Cut = LDCutBW, LowerCut = LDCutBW - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BW", Filter = FilterBW_ordered)
    if not LDCutWB is None: LD_WB = LD(Sys, Cut = LDCutWB, LowerCut = LDCutWB - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WB", Filter = FilterWB_ordered)
    # only BB spline if dry benzene
    if not NW: SP_WW = LD_WW = LD_BW = LD_WB = None

    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_WW, LD_BW, LD_WB]:
        if P: Sys.ForceField.extend([P])
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)

    # system setup
    Sys.Load()
    if LammpsTraj:
        Sys.Arrays.Pos = Trj[0]
    else:
        sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    Int = Sys.Int

    # integrator setup
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = LangevinGamma
    Sys.TempSet = TempSet

    ## other settings based on previous trial runs
    if ManageInnerCore:
        SP_BB.EneInner = "50kT"
        SP_BB.EneSlopeInner = None
        #SP_BB.KnotMinHistFrac = 0.005
        #SP_BB.KnotMinHistFracInner = 0.0
        
        #SP_BW.EneInner = "20kT"
        #SP_BW.EneSlopeInner = None
        SP_BW.KnotMinHistFrac = 0.01
        SP_BB.KnotMinHistFracInner = 0.0
                    
    return Sys


def mapTrj(InTraj, OutTraj =  None):
	# this function has been intentionally kept free of global variables
	# and thus unhinged from the rest of the module, to enable different
	# kinds of mapping if need be
	
	global NB, NW
	Map = sim.atommap.PosMap()
	for i in range(0, NB):  Map += [sim.atommap.AtomMap(Atoms1 = range(i*12, (i+1)*12, 2), Atom2 = i)]
	for i in range(0, NW):  Map += [sim.atommap.AtomMap(Atoms1 = range(NB*12+i, NB*12+i+3, 3), Atom2 = NB+i)]
	AtomTypes = [1]*NB + [2]*NW
	Trj = pickleTraj(InTraj)
	BoxL = Trj.FrameData['BoxL']
	if OutTraj is None: OutTraj = InTraj.split('.lammpstrj.gz')[0] + '_mapped.lammpstrj.gz'
	MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)
	sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, OutTraj, Verbose = True)


def ResetForceField(Sys, OptStage):
    if not OptStage.startswith('SPLD'): return
    LDList = OptStage.split('_')[1:]
    for P in Sys.ForceField:
        if P is None: continue
        if P.Name.startswith('LD') and not LDList.__contains__(P.Name):
            P.SetParam(Knots = 0.0)


def runSrel(Sys, ParamString = None):
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW, BoxL
    global LammpsTraj, Prefix
    global MinSteps, EquilSteps, ProdSteps, StepFreq
    global OptStages, OptStageNames

    Sys = makeSys()
    Map = sim.atommap.PosMap() # note all AA trajectories are in COM format
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

    Trj = pickleTraj(LammpsTraj, Verbose = False)
    Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix, Verbose = True)
    Opt = sim.srel.UseLammps(Opt)
    sim.srel.optimizetraj.PlotFmt = 'svg'
    sim.srel.optimizetrajlammps.LammpsStepsMin  = MinSteps
    Opt.TempFileDir = os.getcwd()
    Opt.MinReweightFrac = 0.15

    # freeze all potentials
    if ParamString: Sys.ForceField.SetParamString(ParamString)
    for P in Sys.ForceField:
        if not P is None: P.FreezeParam()

    # optimization stages
    if OptStages is None:
        OptStages = {
                     # quick optimization stages
                     'SP': ['SP_WW', 'SP_BB', 'SP_BW'],
                     'LD_WW': ['LD_WW'],
                     'LD_BB': ['LD_BB'],
                     'LD_BW': ['LD_BW'],
                     'LD_WB': ['LD_WB'],
                     
                     # refinement/control stages 
                     'SPLD_BB': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB'],
                     'SPLD_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_WW'],
                     'SPLD_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BW'],
                     'SPLD_BB_WW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB', 'LD_WW'],
                     'SPLD_BB_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB', 'LD_BW'],
                     'SPLD_BB_WW_BW': ['SP_WW', 'SP_BB', 'SP_BW', 'LD_BB', 'LD_WW', 'LD_BW'],
                     'All': [P.Name for P in Sys.ForceField if not P is None]}
    if OptStageNames is None:
        OptStageNames = ['SP', 'LD_WW', 'LD_BB', 'LD_BW', 'LD_WB',
                         'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 
                         'SPLD_BB_WW', 'SPLD_BB_BW',
                         'SPLD_BB_WW_BW', 'All']

    # stagewise optimization
    print '\n\n'
    print len(OptStageNames), 'Stages:', ', '.join(OptStageNames)
    for k in OptStageNames:
        Opt.Reset()
        Opt.FilePrefix = Prefix + '_' + k
        print "\nOptimizing stage: %s\n " %  k
        
        FrozenList = [] ;UnfrozenList = []
        for P in Sys.ForceField:
            if P is None: continue
            if OptStages[k].__contains__(P.Name):
                #print 'Estimating potential: %s' % P.Name
                #if ParamString is None: P.Estimate()
                P.UnfreezeParam()
                UnfrozenList.append(P.Name)
            else:
                P.FreezeParam()
                FrozenList.append(P.Name)
        
        print "\nFreezing: %s" % (', '.join(FrozenList))
        print "Unfreezing: %s\n" % (', '.join(UnfrozenList))

        # forcefields are reset during refinement runs
        ResetForceField(Sys, k)
     
        print Sys.ForceField.ParamString()
        
        Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)



def runMD(Sys, ParamString, MDPrefix = None, useParallel = False, NCores = 1):
    global LammpsTraj, Prefix, TempSet
    global MinSteps, EquilSteps, ProdSteps, StepFreq
    
    Trj = pickleTraj(LammpsTraj, Verbose = False)
    BoxL = Trj.FrameData['BoxL']
    
    Sys.ForceField.SetParamString(ParamString)
    
    # sanity check of system parameters
    if not Sys.TempSet == TempSet: Sys.TempSet = TempSet
    Sys.BoxL = BoxL
    
    # set run prefix
    if MDPrefix is None: MDPrefix = Prefix
    
    # parallelization settings
    sim.export.lammps.useParallel = useParallel
    sim.export.lammps.NCores = NCores
    
    ModTraj, ModTrajFile = sim.export.lammps.MakeLammpsTraj(Sys, Prefix = MDPrefix, TrajFile = ".lammpstrj.gz",
                                                            NStepsMin = MinSteps, NStepsEquil = EquilSteps, 
                                                            NStepsProd = ProdSteps, WriteFreq = StepFreq,
                                                            Verbose = True)
    return ModTraj, ModTrajFile
 


#TODO: correct this routine
def genMultiOpt(Sys, Map):
	global BoxL, LammpsTraj

	if not MultiLammpsTraj or not MultiNBList or not MultiNWList:
		raise IOError('Need more information to generate multi Opt object')

	Opts = []
	for i, Traj in enumerate(MultiLammpsTraj):
		NB = MultiNBList[i]
		NW = MultiNWList[i]
		print '\nGenerating CG system for NB = %d, NW = %d from extended ensemble\n' % (NB, NW)
    	LammpsTraj = Traj
    	Trj = pickleTraj(LammpsTraj, Verbose = False)
    	BoxL = pickleTraj(LammpsTraj).FrameData['BoxL']
    	Sys.ScaleBox(BoxL)
    	Opt = sim.srel.OptimizeTrajLammpsClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
    	sim.srel.optimizetraj.PlotFmt = 'svg'
    	Opts.append(Opt)

	del Opt
	Weights = [1.]*len(Opts) #TODO
	MultiOpt = sim.srel.OptimizeMultiTrajClass(OptimizeTrajList = Opts, Weights = Weights,
											   FilePrefix = Prefix + '_multi')

	return MultiOpt
