#usr/bin/env python

import os, sys
import numpy as np

sys.path.append(os.path.expanduser('~/benwat'))
import mysim
import sim

import pickleTraj
import parse_potential as pp

doMinimize = False
DEBUG = False

############################ INPUTS ############################################
global NB, NW, LDCutBB, LDCutBW, LDCutWB, LDCutWW
global LammpsTraj

# system description
NB = 250
NW = 250
Prefix = None
BoxL = []
LDCutBB = None
LDCutBW = None
LDCutWB = None
LDCutWW = None
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

# Multi Srel optimization settings
doMultiSrel = False
MultiLammpsTraj = []
MultiNBList = []
MultiNWList = []

# Lammps settings
sim.export.lammps.LammpsExec = '/home/cask0/home/tsanyal/software/tanmoy_lammps/lammps-15May15/src/lmp_ZIN'
sim.export.lammps.InnerCutoff = 0.02
sim.srel.base.ErrorDiffEneFracTol = 0.2 #increase tolerance to prevent sim-Lammps mismatch blowing up the srel run 

################################################################################
																				
# system parameters
TempSet = 300
Name_W = 'W' ; Name_B = 'B'
Mass_W = 18.01 ; Mass_B = 78.11
Dia_W = 2.8 ; Dia_B = 5.3
SPCutScale = 2.5
RhoMin = 0 ; RhoMax = 50 ; LD_Delta = 1.2
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
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW, BoxL
    global LammpsTraj

    # system chemistry
    AtomTypeW = sim.chem.AtomType(Name_W, Mass = Mass_W, Charge = 0.0, Color = (0,0,1))
    AtomTypeB = sim.chem.AtomType(Name_B, Mass = Mass_B, Charge = 0.0, Color = (1,1,0))
    MolTypeW = sim.chem.MolType(Name_W, [AtomTypeW])
    MolTypeB = sim.chem.MolType(Name_B, [AtomTypeB])
    World = sim.chem.World([MolTypeB, MolTypeW], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = makeSysPrefix())
    for i in range(NB): Sys += MolTypeB.New()
    for i in range(NW): Sys += MolTypeW.New()
    
    # box lengths
    if doMultiSrel: BoxL = [100,100,100] # large dummy boxlength
    else: 
        if not BoxL: BoxL = pickleTraj(LammpsTraj).FrameData['BoxL']
    Sys.BoxL = BoxL
    
    # cutoffs
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
    FilterWW_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW], Ordered = True)
            
    # potential energy objects (only BB spline if dry benzene)
    SP_BB = None ; SP_WW = None; SP_BW = None
    LD_BB = None; LD_WW = None; LD_BW = None; LD_WB = None 
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
            LD_BB = LD(Sys, Cut = LDCutBB, LowerCut = LDCutBB - LD_Delta, NKnot = NLDKnots, 
                       RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BB", Filter = FilterBB_ordered)

        if LDCutWW:
        	LD_WW = LD(Sys, Cut = LDCutWW, LowerCut = LDCutWW - LD_Delta, NKnot = NLDKnots,
        		       RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WW", Filter = FilterWW_ordered)                 
    
    # system forcefield
    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_BW, LD_WB, LD_WW]:
        if P: Sys.ForceField.extend([P])
      
    # set up the histograms, must be done for Srel to work properly
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)
    
    # spline treatment
    Sys.ForceField.SetSplineTreatment(NonbondEneSlope = 40., BondEneSlope = 10., AngleEneSlope = 60.)
    
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
    
    if DEBUG: print Sys.ForceField.ParamString()
    return Sys

    
def runSrel(Sys, Opt_cases = None):
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW, BoxL
    global LammpsTraj

    Sys = makeSys()
    Map = sim.atommap.PosMap() # note all AA trajectories are in COM format
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
    
    if not doMultiSrel:
    	Trj = pickleTraj(LammpsTraj, Verbose = False)
    	Opt = sim.srel.OptimizeTrajLammpsClass(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
    else:
    	Opt = genMultiOpt(Sys, Map)
    
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
    
    # create separate lists for pair and local density potentials
    SPList = []
    LDList = []
    for P in Sys.ForceField:
    	if P.Name.__contains__('SP'): SPList.append(P)
    	if P.Name.__contains__('LD'): LDList.append(P)


    # relative entropy minimization
    if not Opt_cases: Opt_cases = ["SP", "SPLD_BB", "SPLD_WW", "SPLD_BW", "SPLD_WB", "SPLD_BB_WW", "SPLD_BB_BW", "SPLD_BB_WB", "SPLD_WW_BW", "SPLD_WW_WB", "SPLD_BB_WW_BW", "SPLD_BB_WW_WB"]
    for i, case in enumerate(Opt_cases):
        Opt.Reset()
        Opt.FilePrefix = Prefix + '_' + case
        print "\n\nOptimizing for the case: ", case
        
        if case == "SP":
            for P in Sys.ForceField:
            	if SPList.__contains__(P):
            		P.UnfreezeParam()
            	else:
            		P.SetParam(Knots = 0.)
            		P.FreezeParam() 
        

        if not LDList: print "No LD potentials present."

        
        if ["SPLD_BB", "SPLD_WW", "SPLD_BW", "SPLD_WB"].__contains__(case):
        	ptype = case[2:]
        	for P in Sys.ForceField:
        		if SPList.__contains__(P):
        			P.UnfreezeParam()
        			if os.path.isfile('%s_SP_sum.txt' % Prefix):
        				SPKnots = pp.parseParam(sumfile = '%s_SP_sum.txt' % Prefix, ptype = P.Name)['Knots']
        				P.SetParam(Knots = SPKnots)

        		if (not SPList.__contains__(P)) and (LDList.__contains__(P)) and (P.Name != ptype):
        			P.SetParam(Knots = 0.)
        			P.FreezeParam()

        		if P.Name == ptype: P.UnfreezeParam()

        
        if ["SPLD_BB_WW", "SPLD_BB_BW", "SPLD_BB_WB", "SPLD_WW_BW", "SPLD_WW_WB"].__contains__(case):
        	ptype1 = 'LD_' + case.split('_')[1]
        	ptype2 = 'LD_' + case.split('_')[-1]
        	for P in Sys.ForceField:
        		if SPList.__contains__(P):
        			P.UnfreezeParam()
        			if os.path.isfile('%s_SP_sum.txt' % Prefix):
        				SPKnots = pp.parseParam(sumfile = '%s_SP_sum.txt' % Prefix, ptype = P.Name)['Knots']
        				P.SetParam(Knots = SPKnots)

        		if (not SPList.__contains__(P)) and (LDList.__contains__(P)) and (not [ptype1, ptype2].__contains__(P.Name)):
        			P.SetParam(Knots = 0.)
        			P.FreezeParam()

        		for ptype in [ptype1, ptype2]:
        			if P.Name == ptype:
        				P.UnfreezeParam()
        				if os.path.isfile('%s_SP%s_sum.txt' % (Prefix, ptype)):
        					LDKnots = pp.parseParam(sumfile = '%s_SP%s_sum.txt' % (Prefix, ptype), ptype = ptype)['Knots']
        					P.SetParam(Knots = LDKnots)

        
        if ["SPLD_BB_WW_BW", "SPLD_BB_WW_WB"].__contains__(case):
            ptype1 = 'LD_' + case.split('_')[1]
            ptype2 = 'LD_' + case.split('_')[2]
            ptype3 = 'LD_' + case.split('_')[3]
            for P in Sys.ForceField:
                if SPList.__contains__(P):
                    P.UnfreezeParam()
                    if os.path.isfile('%s_SP_sum.txt' % Prefix):
                        SPKnots = pp.parseParam(sumfile = '%s_SP_sum.txt' % Prefix, ptype = P.Name)['Knots']
                        P.SetParam(Knots = SPKnots)

                if (not SPList.__contains__(P)) and (LDList.__contains__(P)) and (not [ptype1, ptype2, ptype3].__contains__(P.Name)):
                    P.SetParam(Knots = 0.)
                    P.FreezeParam()

                for ptype in [ptype1, ptype2, ptype3]:
                    if P.Name == ptype:
                        P.UnfreezeParam()
                        if os.path.isfile('%s_SP%s_sum.txt' % (Prefix, ptype)):
                            LDKnots = pp.parseParam(sumfile = '%s_SP%s_sum.txt' % (Prefix, ptype), ptype = ptype)['Knots']
                            P.SetParam(Knots = LDKnots)
        
        Sys.ForceField.Update()

        if DEBUG:
            print Sys.ForceField.ParamString()
            for P in Sys.ForceField:
                count = 0
                for param in P.Param.Fixed:
                    if param: count += 1
                if count == len(P.Knots): print '\n%s Frozen' % P.Name

        Opt.RunConjugateGradient(StepsEquil = EquilSteps, StepsProd = ProdSteps, StepsStride = StepFreq)
       
    print "Srel minimization for different cases finished"
    del Opt


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


def runCGMD(Sys, forcefield):
    pp.loadParam(Sys = Sys, sumfile = forcefield)
    # export to Lammps
    Trj, TrjFile = sim.export.lammps.MakeLammpsTraj(Sys = Sys, Prefix = Prefix, NStepsMin = MinSteps, 
                                                    NStepsEquil = EquilSteps, NStepsProd = ProdSteps, 
                                                    WriteFreq = StepFreq, TrajFile = '.lammpstrj.gz', Verbose = True)
    return Trj, TrjFile  

