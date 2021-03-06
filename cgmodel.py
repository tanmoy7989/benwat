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
SPCutBB = None
SPCutWW = None
SPCutBW = None
SPCutScale = 2.5
NSPKnots = 30
ManageInnerCore = True # True to capture inner core tails of corresponding histograms

# Local density potential settings
LDCutWW = None
LDCutBB = None
LDCutBW = None
LDCutWB = None
RhoMin = 0
RhoMax = RhoMax_BB = RhoMax_WW = RhoMax_BW = RhoMax_WB = 20
LD_Delta = 1.0
NLDKnots = 30
NLDWWKnots = 50
useMoreLDWWKnots = False

# Histogram settings
ReportNBin = 100

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
StartWithHessian = False
NoHessianMaxIter = 100

# Extended Ensemble approach settings
doMultiSrel = False
LammpsTrajList = []
NBList = []
NWList = []

# Lammps settings
sim.export.lammps.LammpsExec = 'lmp_mpich2'
sim.export.lammps.InnerCutoff = 0.02
sim.srel.base.DiffEneFracTol = 0.1 #increase tolerance to prevent sim-Lammps mismatch blowing up the srel run

def makeSys():
    global NB, NW, BoxL
    global SPCutBB, SPCutWW, SPCutBW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW
    global RhoMax_BB, RhoMax_WW, RhoMax_BW, RhoMax_WB
    global LammpsTraj, Prefix
    global ReportNBin

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
    if SPCutWW is None: SPCutWW = SPCutScale * Dia_W
    if SPCutBB is None: SPCutBB = SPCutScale * Dia_B
    if SPCutBW is None: SPCutBW = SPCutScale * 0.5 * (Dia_B + Dia_W)

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
    
    n = NLDKnots if not useMoreLDWWKnots else NLDWWKnots
    if not LDCutWW is None: LD_WW = LD(Sys, Cut = LDCutWW, LowerCut = LDCutWW - LD_Delta, NKnot = n,
                       				   RhoMin = RhoMin, RhoMax = RhoMax_WW, Label = "LD_WW", Filter = FilterWW_ordered)
    if not LDCutBB is None: LD_BB = LD(Sys, Cut = LDCutBB, LowerCut = LDCutBB - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax_BB, Label = "LD_BB", Filter = FilterBB_ordered)
    if not LDCutBW is None: LD_BW = LD(Sys, Cut = LDCutBW, LowerCut = LDCutBW - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax_BW, Label = "LD_BW", Filter = FilterBW_ordered)
    if not LDCutWB is None: LD_WB = LD(Sys, Cut = LDCutWB, LowerCut = LDCutWB - LD_Delta, NKnot = NLDKnots,
                       				   RhoMin = RhoMin, RhoMax = RhoMax_WB, Label = "LD_WB", Filter = FilterWB_ordered)
    
    # only BB spline and BB local density if dry benzene
    if NW == 0: SP_WW = SP_BW = LD_WW = LD_BW = LD_WB = None

    # only WW spline and WW local density if pure water
    if NB == 0: SP_BB = SP_BW = LD_BB = LD_BW = LD_WB = None

    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_WW, LD_BW, LD_WB]:
        if not P is None: Sys.ForceField.extend([P])
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = ReportNBin)

    # system setup
    Sys.Load()
    if LammpsTraj:
        Sys.BoxL = Trj.FrameData['BoxL']
        Sys.Arrays.Pos = Trj[0]
    else: sim.system.init.positions.CubicLatticeFill(Sys, L = 1000., Random = 0.1)
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    Int = Sys.Int

    # integrator setup
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = LangevinGamma
    Sys.TempSet = TempSet

    ## other settings based on previous trial runs
    if ManageInnerCore:
        if SP_BB:
            SP_BB.EneInner = "50kT"
            SP_BB.KnotMinHistFracInner = 0.0
            SP_BB.EneSlopeInner = None
        if SP_BW:
            SP_BW.KnotMinHistFrac = 0.01
                    
    
    for P in [LD_BB, LD_WW, LD_BW, LD_WB]:
        print P.Name, P.RhoMin, P.RhoMax, P.Cut, len(P.Knots)
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
	print BoxL
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
    # do not use this routine to CG pure benzene or pure water !
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW, BoxL
    global LammpsTraj, Prefix
    global MinSteps, EquilSteps, ProdSteps, StepFreq
    global OptStages, OptStageNames, StartWithHessian, NoHessianMaxIter

    Sys = makeSys()
    Map = sim.atommap.PosMap() # note all AA trajectories are in COM format
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]

    Trj = pickleTraj(LammpsTraj, Verbose = False)
    
    if doMultiSrel:
        Opt = genMultiOpt(Sys, Map)
    else:
        Opt = sim.srel.OptimizeTrajClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix, Verbose = True)
        Opt = sim.srel.UseLammps(Opt)
        sim.srel.optimizetraj.PlotFmt = 'svg'
        Opt.TempFileDir = os.getcwd()
        Opt.MinReweightFrac = 0.15
    
    sim.srel.optimizetrajlammps.LammpsStepsMin  = MinSteps
    
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
                     'SPLD_all': [P.Name for P in Sys.ForceField if not P is None]}
    if OptStageNames is None:
        OptStageNames = ['SP', 'LD_WW', 'LD_BB', 'LD_BW', 'LD_WB',
                         'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 
                         'SPLD_BB_WW', 'SPLD_BB_BW',
                         'SPLD_BB_WW_BW', 'SPLD_all']

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
        
        # Hessian on or off (only for SPLD cases)
        if StartWithHessian:
            # check if Hessian is requested from start
            print "Using Hessian descent from the start"
            Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)
        elif not StartWithHessian:
            # turn off hessian for requested number of steps
            print "Turning off Hessian descent for first %d steps" % NoHessianMaxIter
            Opt.UseHessian = False
            Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq,
                                     MaxIter = NoHessianMaxIter)
            print "Turning on Hessian descent"
            Opt.UseHessian = True
            Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)
        

def runMD(Sys, ParamString, MDPrefix = None, BoxL = None, useParallel = False, NCores = 1, autoSubmit = True):
    global LammpsTraj, Prefix, TempSet
    global MinSteps, EquilSteps, ProdSteps, StepFreq
    
    if LammpsTraj is None:
        if not BoxL is None:
            Sys.BoxL = BoxL
            sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
    
    Sys.ForceField.SetParamString(ParamString)
    
    # sanity check of system parameters
    if not Sys.TempSet == TempSet: Sys.TempSet = TempSet
    
    # set run prefix
    if MDPrefix is None: MDPrefix = Prefix
    
    # parallelization settings
    sim.export.lammps.useParallel = useParallel
    sim.export.lammps.NCores = NCores
    sim.export.lammps.autoSubmit = autoSubmit
    
    # increase time-step to 2 fs (default is 1 fs)
    Sys.Int.Method.TimeStep *= 2
    
    ret = sim.export.lammps.MakeLammpsTraj(Sys, Prefix = MDPrefix, TrajFile = ".lammpstrj.gz",
                                           NStepsMin = MinSteps, NStepsEquil = EquilSteps, 
                                           NStepsProd = ProdSteps, WriteFreq = StepFreq,
                                           Verbose = True)
    return 
  
def genMultiOpt(Sys, Map):
    global LammpsTrajList, NBList, NWList, Prefix
    if not (LammpsTrajList or not NBList or not NWList):
        raise IOError('Need more information to generate multi Opt object')
    
    Opts = []
    for i, Traj in enumerate(LammpsTrajList):
        NB = NBList[i] ; NW = NWList[i]
        print '\nGenerating optimizer object for NB = %d, NW = %d' % (NB, NW)
        Trj = pickleTraj(Traj, Verbose = False)
        BoxL = pickleTraj(Traj).FrameData['BoxL']
        Sys.ScaleBox(BoxL)
        Opt = sim.srel.OptimizeTrajClass(ModSys = Sys, Map = Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = '%s_NB%dNW%d' % (Prefix, NB, NW))
        Opt = sim.srel.UseLammps(Opt)
        sim.srel.optimizetraj.PlotFmt = 'svg'
        Opt.TempFileDir = os.getcwd()
        Opt.MinReweightFrac = 0.15
        Opts.append(Opt)
    
    Weights = [1.]*len(Opts) # for starters give equal weight to all concentrations
    MultiOpt = sim.srel.OptimizeMultiTrajClass(OptimizeTrajList = Opts, Weights = Weights, FilePrefix = Prefix)
    return MultiOpt
