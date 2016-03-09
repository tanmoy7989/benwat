#usr/bin/env python

import os, sys
import numpy as np
import sim
import pickleTraj

doMinimize = False
DEBUG = True

############################ INPUTS ############################################
# system description
NB = 250
NW = 250
Prefix = None
BoxL = None
LDCutBB = None
LDCutBW = None
LDCutWB = None
ParamString = None

# MD settings
MinSteps = 1000
NPTSteps = 1000
EquilSteps = 1000000
ProdSteps = 20000000
StepFreq = 1000

# Srel optimization settings
LammpsTraj = None
CG_Thermostat = None
LangevinGamma = 0.01

# Lammps settings
sim.export.lammps.LammpsExec = '/home/cask0/home/tsanyal/software/tanmoy_lammps/lammps-15May15/src/lmp_ZIN'
sim.export.lammps.InnerCutoff = 0.02
sim.srel.base.ErrorDiffEneFracTol = 1.0 #increase tolerance to prevent sim-Lammps mismatch blowing up the srel run 

################################################################################
																				
# system parameters
TempSet = 300
Name_W = 'W' ; Name_B = 'B'
Mass_W = 18.01 ; Mass_B = 78.11
Dia_W = 2.8 ; Dia_B = 5.3
SPCutScale = 2.5
RhoMin = 0 ; RhoMax = 50 ; LD_Delta = 0.5
NSPKnots = 40 ; NLDKnots = 50


def makeSysPrefix():
    if Prefix: return Prefix
    else: return 'NB%dNW%d' % (NB, NW)


def isImageOverlap(Cut = None):
    flag = 0
    HalfBoxL = 0.5 * np.array(BoxL)
    for x in HalfBoxL:
        if Cut > x: 
            flag = 1
            break
    return flag


def SetSPCutoff():
    SPCutWW = SPCutScale * Dia_W
    SPCutBB = SPCutScale * Dia_B
    SPCutBW = SPCutScale * 0.5 * (Dia_B + Dia_W)
    return (SPCutWW, SPCutBB, SPCutBW)
    

def makeSys():
    # system chemistry
    AtomTypeW = sim.chem.AtomType(Name_W, Mass = Mass_W, Charge = 0.0, Color = (0,0,1))
    AtomTypeB = sim.chem.AtomType(Name_B, Mass = Mass_B, Charge = 0.0, Color = (1,1,0))
    MolTypeW = sim.chem.MolType(Name_W, [AtomTypeW])
    MolTypeB = sim.chem.MolType(Name_B, [AtomTypeB])
    World = sim.chem.World([MolTypeB, MolTypeW], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = makeSysPrefix())
    for i in range(NB): Sys += MolTypeB.New()
    for i in range(NW): Sys += MolTypeW.New()
    
    # box length and cutoffs
    Sys.BoxL[:] = BoxL
    SPCutWW, SPCutBB, SPCutBW = SetSPCutoff()
    for cut in [SPCutWW, SPCutBB, SPCutBW, LDCutBW, LDCutBW]:
        if cut and isImageOverlap(cut): raise ValueError('Images are overlapping. Revise cutoffs')
        
    # atom selection filters
    FilterBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB])    
    FilterWW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW])
    FilterBW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB])
    FilterBW_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeW], Ordered = True)
    FilterWB_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB], Ordered = True)
    FilterBB_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], Ordered = True)
            
    # potential energy objects (only BB spline if dry benzene)
    SP_BB = None ; SP_WW = None; SP_BW = None
    LD_BW = None; LD_WB = None ; LD_BB = None
    SP = sim.potential.PairSpline
    LD = sim.potential.LocalDensity
    SP_BB = SP(Sys, Cut = SPCutBB, NKnot = NSPKnots, Filter = FilterBB, Label = "SP_BB")
    if NW:
        SP_WW = SP(Sys, Cut = SPCutWW, NKnot = NSPKnots, Filter = FilterWW, Label = "SP_WW")
        SP_BW = SP(Sys, Cut = SPCutBW, NKnot = NSPKnots, Filter = FilterBW, Label = "SP_BW")
        if LDCutBW:
            LD_BW = LD(Sys, Cut = LDCutBW, LowerCut = LDCutBW - LD_Delta, NKnot = NLDKnots, 
                       RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BW", Filter = FilterBW_ordered)
    
        if LDCutWB:
            LD_WB = LD(Sys, Cut = LDCutWB, LowerCut = LDCutWB - LD_Delta, NKnot = NLDKnots, 
                       RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WB", Filter = FilterWB_ordered)
        if LDCutBB:
            LD_BB = LD(Sys, Cut = LDCutBB, LowerCut = LDCutWB - LD_Delta, NKnot = NLDKnots, 
                       RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BB", Filter = FilterBB_ordered)                  
    
    # system forcefield
    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_BW, LD_WB]:
        if P: Sys.ForceField.extend([P])
      
    # set up the histograms, must be done for Srel to work properly
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    
    # spline treatment
    Sys.ForceField.SetSplineTreatment(NonbondEneSlope = 50., BondEneSlope = 10., AngleEneSlope = 60.)
    
    # compile and load the system
    Sys.Load()

    # set up initial positions and velocities
    sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    
    # reference the integrator
    Int = Sys.Int
    
    # time-step
    #For now run with default sim time-step
    
    # load the force field if supplied
    if ParamString: Sys.ForceField.SetParamString(ParamString)
    
    # preliminary energy minimization (not advisable, since sim serial routines are slow)
    if doMinimize:
    	Int.Method = Int.Methods.VVQuench
    	Int.Run(10000, "Minimizing")
    
    #change to MD
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = LangevinGamma
    Sys.TempSet = TempSet
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    
    if DEBUG: print Sys.ForceField.ParamString()
    return Sys


def restart(Sys, opts, ParamString = None, RestartFrom = 'SP'):
    if DEBUG: print "Restart point = ", RestartFrom
    x = opts.index(RestartFrom)
    opts = opts[x:]
    if ParamString: Sys.ForceField.SetParamString(ParamString)
    Sys.ForceField.Update
    return Sys, opts

    
def runSrel(Sys, ParamString = None, RestartFrom = "SP"):
    # make map (note: all AA trajectories are stored in mapped format and so 1:1 mapping here)
    Map = sim.atommap.PosMap()
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
    
    # optimizer class
    Trj = pickleTraj(LammpsTraj, Verbose = False)
    Opt = sim.srel.OptimizeTrajLammpsClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
    sim.srel.optimizetraj.PlotFmt = 'svg'
    Opt.StepsMin = MinSteps

    # Lammps log inconsistency debugging
    if DEBUG:
        Opt.TempFilePrefix = 'logtest'
        Opt.TempFileDir = os.getcwd()

    # output initial histograms if required
    if DEBUG:
        Opt.MakeModTraj(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)
        Opt.OutputModHistFile()
        Opt.OutputPlot()

    # freeze all potentials
    [P.FreezeParam() for P in Sys.ForceField]
    
    # relative entropy minimization
    Opt_cases = ["SP", "SPLD_BW", "SPLD_all"]
    for i, case in enumerate(Opt_cases):
        Opt.Reset()
        Opt.FilePrefix = Prefix + '_' + case
        print "\n\nOptimizing for the case: ", case
        
        if case == "SP":
            for P in Sys.ForceField:
                if ['SP_WW', 'SP_BB', 'SP_BW'].__contains__(P.Name): P.UnfreezeParam()
        
        if case == 'SPLD_BW':
            for P in Sys.ForceField:
                if P.Name == 'LD_BW': P.UnfreezeParam()
        
        if case == "SPLD_all":
            for P in Sys.ForceField:
                if P.Name == 'LD_BB': P.UnfreezeParam()
        
        # add more cases here if required
        
        Sys.ForceField.Update()
        Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)
       
    print "Srel minimization for different cases finished"
    del Opt


def runCGMD(Sys):
    Sys.ForceField.SetParamString(ParamString)
    # export to Lammps
    Trj, TrjFile = sim.export.lammps.MakeLammpsTraj(Sys = Sys, Prefix = Prefix, NStepsMin = MinSteps, 
                                                    NStepsEquil = EquilSteps, NStepsProd = ProdSteps, 
                                                    WriteFreq = StepFreq, TrajFile = '.lammpstrj.gz', Verbose = True)
    return Trj, TrjFile  

