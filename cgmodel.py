#usr/bin/env python

import os, sys
import numpy as np

sys.path.append(os.path.expanduser('~/software')) ; import mysim ; import sim

import pickleTraj
import parse_potential as pp

DEBUG = False

# system description
Prefix = None
NB = 250
NW = 250
BoxL = None
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
LangevinGamma = 0.01
useLammps = True
isInit = False

# Lammps settings
sim.export.lammps.LammpsExec = os.path.expanduser('~/software/tanmoy_lammps/lammps-15May15/src/lmp_ZIN')
sim.export.lammps.InnerCutoff = 0.02
sim.srel.base.DiffEneFracTol = 0.1 #increase tolerance to prevent sim-Lammps mismatch blowing up the srel run 

# Parallelization settings
useParallel = False
NCores = 2
																				
# System parameters
TempSet = 300
Name_W = 'W' ; Name_B = 'B'
Mass_W = 18.01 ; Mass_B = 78.11
Dia_W = 2.8 ; Dia_B = 5.3
SPCutScale = 2.5
RhoMin = 0 ; RhoMax = 50 ; LD_Delta = 1.0
NSPKnots = 30 ; NLDKnots = 30


def isImageOverlap(Cut = None):
    if Cut is None: return
    global BoxL
    flag = 0
    HalfBoxL = 0.5 * np.array(BoxL)
    for x in HalfBoxL:
        if Cut > x: 
            flag = 1
            break
    return flag


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
    if BoxL is None: BoxL = Trj.FrameData['BoxL']
    Sys.BoxL = BoxL
    for cut in [SPCutWW, SPCutBB, SPCutBW, LDCutWW, LDCutBB, LDCutBW, LDCutWB]:
    	if not cut is None and isImageOverlap(cut): raise ValueError('Cutoff %g violates min image distance' % cut)

    FilterWW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW])
    FilterBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB])    
    FilterBW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB])
    FilterWW_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW], Ordered = True)
    FilterBB_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], Ordered = True)
    FilterBW_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeW], Ordered = True)
    FilterWB_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB], Ordered = True)
            
    # forcefield (only BB spline if dry benzene)
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
    
    if not NW: SP_WW = LD_WW = LD_BW = LD_WB = None
    for P in [SP_WW, SP_BB, SP_BW, LD_WW, LD_BB, LD_BW, LD_WB]:
        if P: Sys.ForceField.extend([P])
    Sys.ForceField.SetSplineTreatment(NonbondEneSlope = 40., BondEneSlope = 10., AngleEneSlope = 60.)
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)

    # system setup
    Sys.Load()
    if LammpsTraj: Sys.Arrays.Pos = Trj[0]
    else: sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    Int = Sys.Int
    for Method in Int.Methods:
    	if hasattr(Method, 'TimeStep'): Method.TimeStep = 2.0
    
    if not useLammps:
    	Int.Method = Int.Methods.VVQuench
    	Int.Run(1000, "Minimizing")
    
    #change to MD
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = LangevinGamma
    Sys.TempSet = TempSet

    return Sys


def mapTrj(InTraj, OutTraj =  None):
	# this function has been intentionally kept free of global variables
	# and thus unhinged from the rest of the module, to enable different 
	# kinds of mapping if need be
	Map = sim.atommap.PosMap()
	for i in range(0, NB):  Map += [sim.atommap.AtomMap(Atoms1 = range(i*12, (i+1)*12, 2), Atom2 = i)]
	for i in range(0, NW):  Map += [sim.atommap.AtomMap(Atoms1 = range(NB*12+i, NB*12+i+3, 3), Atom2 = NB+i)]
	AtomTypes = [1]*NB + [2]*NW
	Trj = pickleTraj(InTraj)
	BoxL = Trj.FrameData['BoxL']
	if OutTraj is None: OutTraj = InTraj.split('.lammpstrj.gz')[0] + '_mapped.lammpstrj.gz'
	MappedTrj = sim.traj.Mapped(Trj, Map, AtomNames = AtomTypes, BoxL = BoxL)
	sim.traj.Convert(MappedTrj, sim.traj.LammpsWrite, OutTraj, Verbose = True)            


def runSrel(Sys, ParamString = None, Opt_cases = None):
    # Lammps export setting
    global useLammps, useParallel, NCores
    sim.export.lammps.useParallel = useParallel
    sim.export.lammps.NCores = NCores
    
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW, BoxL
    global LammpsTraj, Prefix

    Sys = makeSys()
    Map = sim.atommap.PosMap() # note all AA trajectories are in COM format
    for (i, a) in enumerate(Sys.Atom): Map += [sim.atommap.AtomMap(Atoms1 = i, Atom2 = a)]
    
    Trj = pickleTraj(LammpsTraj, Verbose = False)
    OptObj = sim.srel.OptimizeTrajLammpsClass if useLammps else sim.srel.OptimizeTrajClass
    Opt = OptObj(Sys, Map, Traj = Trj, SaveLoadArgData = True, FilePrefix = Prefix)
    sim.srel.optimizetraj.PlotFmt = 'svg'
    Opt.StepsMin = MinSteps
    Opt.TempFileDir = os.getcwd()
    Opt.MinReweightFrac = 0.2

    # freeze all potentials
    if ParamString: Sys.ForceField.SetParamString(ParamString)
    [P.FreezeParam() for P in Sys.ForceField]
    
    # create separate lists for pair and local density potentials
    SPList = []
    LDList = []
    for P in Sys.ForceField:
    	if P.Name.__contains__('SP'): SPList.append(P)
    	if P.Name.__contains__('LD'): LDList.append(P)

    def checkInit(Forcefield, ptypes):
    	if not isInit: return
    	for potential in Forcefield:
    		if SPList.__contains__(potential): potential.FreezeParam()
    		if LDList.__contains__(potential) and not ptypes.__contains__(potential.Name): potential.FreezeParam()

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
        		checkInit(Sys.ForceField, [ptype])

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

        		checkInit(Sys.ForceField, [ptype1, ptype2])
        
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

                checkInit(Sys.ForceField, [ptype1, ptype2, ptype3])

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

