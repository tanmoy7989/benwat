#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle, copy
import sim, pickleTraj, parse_potential as pp

# user input
NB = None
NW = None
AATraj = None
Prefix = None
FF_File = None

# System conditions
TempSet = 300
Name_W = 'W' ; Name_B = 'B'
Mass_W = 18.01 ; Mass_B = 78.11

# Pair potential settings
Dia_W = 2.8 ; Dia_B = 5.3
SPCutScale = 2.5
NSPKnots = 30

# Local density potential settings
LDCutWW = 3.5
LDCutBB = 7.5
LDCutBW = LDCutWB = 0.5 * (LDCutBB + LDCutWW)
RhoMin = 0
RhoMax = 20
LD_Delta = 1.0
NLDKnots = 30

# MD settings
MinSteps = 10000 # since test particle is placed at origin --> greater repulsion
EquilSteps = 1000000
ProdSteps = 2000000
StepFreq = 500
AutoSubmit = False

# Lammps settings
LammpsExec = 'lmp_mpich2'
sim.export.lammps.InnerCutoff = 0.02

# TI settings
lambda_factors = np.linspace(0.0, 1.0, 11)

# make TI Sys
def makeTISys(TPType = 'W', lambda_factor = 1.0):
    print 'Making TI Sys object at lambda = %g' % lambda_factor
    global NB, NW
    global LDCutBB, LDCutBW, LDCutWB, LDCutWW
    global AATraj, Prefix

    # system chemistry
    if Prefix is None: Prefix = 'NB%dNW%d' % (NB, NW)
    AtomTypeW = sim.chem.AtomType(Name_W, Mass = Mass_W, Charge = 0.0, Color = (0,0,1))
    AtomTypeB = sim.chem.AtomType(Name_B, Mass = Mass_B, Charge = 0.0, Color = (1,1,0))
    # insert test particle 
    if TPType == 'B': AtomTypeX = sim.chem.AtomType('X', Mass = Mass_B, Charge = 0.0)
    elif TPType == 'W': AtomTypeX = sim.chem.AtomType('X', Mass = Mass_W, Charge = 0.0)
    MolTypeW = sim.chem.MolType(Name_W, [AtomTypeW])
    MolTypeB = sim.chem.MolType(Name_B, [AtomTypeB])
    MolTypeX = sim.chem.MolType('X', [AtomTypeX])
    World = sim.chem.World([MolTypeB, MolTypeW, MolTypeX], Dim = 3, Units = sim.units.AtomicUnits)
    Sys = sim.system.System(World, Name = Prefix)
    # add 1 less particle of whichever is greater
    if NB > NW:
        thisNB = NB - 1
        thisNW = NW
    elif NB < NW:
        thisNB = NB 
        thisNW = NW - 1
    for i in range(thisNB): Sys += MolTypeB.New()
    for i in range(thisNW): Sys += MolTypeW.New()
    Sys += MolTypeX.New()

    # pair potential cutoffs
    SPCutWW = SPCutScale * Dia_W
    SPCutBB = SPCutScale * Dia_B
    SPCutBW = SPCutScale * 0.5 * (Dia_B + Dia_W)

    # atom selection filters for system
    FilterWW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW])
    FilterBB = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB])
    FilterBW = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB])
    FilterWW_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeW], Ordered = True)
    FilterBB_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeB], Ordered = True)
    FilterBW_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeW], Ordered = True)
    FilterWB_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeB], Ordered = True)
    
    # atom selection filters for test particle
    # no X-X filters required since only 1 particle is inserted
    FilterXW = sim.atomselect.PolyFilter([AtomTypeX, AtomTypeW])
    FilterXB = sim.atomselect.PolyFilter([AtomTypeX, AtomTypeB])
    FilterXW_ordered = sim.atomselect.PolyFilter([AtomTypeX, AtomTypeW], Ordered = True)
    FilterWX_ordered = sim.atomselect.PolyFilter([AtomTypeW, AtomTypeX], Ordered = True)
    FilterXB_ordered = sim.atomselect.PolyFilter([AtomTypeX, AtomTypeB], Ordered = True)
    FilterBX_ordered = sim.atomselect.PolyFilter([AtomTypeB, AtomTypeX], Ordered = True)

    # system forcefield
    SP_BB = SP_WW = SP_BW = LD_BB = LD_WW = LD_BW = LD_WB = None
    SP = sim.potential.PairSpline
    LD = sim.potential.LocalDensity
    if NB > 0:
        SP_BB = SP(Sys, Cut = SPCutBB, NKnot = NSPKnots, Filter = FilterBB, Label = "SP_BB")
    if NW > 0:
        SP_WW = SP(Sys, Cut = SPCutWW, NKnot = NSPKnots, Filter = FilterWW, Label = "SP_WW")
    if NB > 0 and NW > 0:
        SP_BW = SP(Sys, Cut = SPCutBW, NKnot = NSPKnots, Filter = FilterBW, Label = "SP_BW")
    if NW > 0:
        LD_WW = LD(Sys, Cut = LDCutWW, LowerCut = LDCutWW - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WW", Filter = FilterWW_ordered)
    if NB > 0:
        LD_BB = LD(Sys, Cut = LDCutBB, LowerCut = LDCutBB - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BB", Filter = FilterBB_ordered)
    if NB > 0 and NW > 0:
        LD_BW = LD(Sys, Cut = LDCutBW, LowerCut = LDCutBW - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BW", Filter = FilterBW_ordered)
        LD_WB = LD(Sys, Cut = LDCutWB, LowerCut = LDCutWB - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WB", Filter = FilterWB_ordered)
    
    # test particle forcefield
    if TPType == 'B':
        SPCutXB = SPCutBB
        SPCutXW = SPCutBW
        LDCutXB = LDCutBX = LDCutBB
        LDCutXW = LDCutWX = LDCutBW
    elif TPType == 'W':
        SPCutXB = SPCutBW
        SPCutXW = SPCutWW
        LDCutXB = LDCutBX = LDCutBW
        LDCutXW = LDCutWX = LDCutWW
    SP_XB = SP_XW = LD_XB = LD_BX = LD_XW = LD_WX = None
    if NB > 0:
        SP_XB = SP(Sys, Cut = SPCutXB, NKnot = NSPKnots, Filter = FilterXB, Label = "SP_XB")
        LD_XB = LD(Sys, Cut = LDCutXB, LowerCut = LDCutXB - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_XB", Filter = FilterXB_ordered)
        LD_BX = LD(Sys, Cut = LDCutBX, LowerCut = LDCutBX - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_BX", Filter = FilterBX_ordered)
    if NW > 0:
        SP_XW = SP(Sys, Cut = SPCutXW, NKnot = NSPKnots, Filter = FilterXW, Label = "SP_XW")
        LD_XW = LD(Sys, Cut = LDCutXW, LowerCut = LDCutXW - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_XW", Filter = FilterXW_ordered)
        LD_WX = LD(Sys, Cut = LDCutWX, LowerCut = LDCutWX - LD_Delta, NKnot = NLDKnots, RhoMin = RhoMin, RhoMax = RhoMax, Label = "LD_WX", Filter = FilterWX_ordered)
    
    # load forcefield parameters for system
    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_WW, LD_BW, LD_WB]:
        if not P is None:
            Knots = pp.parseParamDict(FF_File, P.Name)['Knots']
            P.SetParam(Knots = Knots)
    
    # load forcefied parameters for test particle
    if TPType == 'B':
        pmap = {SP_XB: 'SP_BB', SP_XW: 'SP_BW', 
                LD_XB: 'LD_BB', LD_BX: 'LD_BB', LD_XW: 'LD_BW', LD_WX: 'LD_WB'}
    elif TPType == 'W':
        pmap = {SP_XB: 'SP_BW', SP_XW: 'SP_WW',
                LD_XB: 'LD_WB', LD_BX: 'LD_BW', LD_XW: 'LD_WW', LD_WX: 'LD_WW'}
    for P in [SP_XB, SP_XW, LD_XB, LD_BX, LD_XW, LD_WX]:
        if (not P is None):
            Knots = pp.parseParamDict(FF_File, pmap[P])['Knots']
            P.SetParam(Knots = lambda_factor * Knots)
    
    # add only relevant potentials to the forcefield
    for P in [SP_BB, SP_WW, SP_BW, LD_BB, LD_WW, LD_BW, LD_WB, SP_XB, SP_XW, LD_XB, LD_BX, LD_XW, LD_WX]:
        if not P is None: Sys.ForceField.append(P)
    for P in Sys.ForceField: P.Arg.SetupHist(NBin = 10000, ReportNBin = 100)

    # system setup
    Sys.Load()
    if not AATraj is None:
        Trj = pickleTraj(AATraj)
        Sys.BoxL = Trj.FrameData['BoxL']
        Sys.Arrays.Pos = Trj[0]
        # place test particle at origin
        #Sys.Arrays.Pos[-1, :] = [0.0, 0.0, 0.0]
    else: 
        Sys.BoxL = 0.0
        sim.system.init.positions.CubicLatticeFill(Sys, L = 1000., Random = 0.1)
    sim.system.init.velocities.Canonical(Sys, Temp = TempSet)
    Int = Sys.Int

    # integrator setup
    Int.Method = Int.Methods.VVIntegrate
    Int.Method.Thermostat = Int.Method.ThermostatLangevin
    Int.Method.LangevinGamma = 0.01
    Sys.TempSet = TempSet
    
    return Sys


def write2Lammps(Sys, lambda_index = 0):
    global Prefix
    global MinSteps, EquilSteps, ProdSteps, StepFreq
    global AutoSubmit
    
    MDPrefix = '%s_lambda%d' % (Prefix, lambda_index)
    
    #add compute for inserted particle energy
    s_before = '''
group           X id %(XID)d
compute         pertpe all pe/atom
compute         tpe X reduce sum c_pertpe
fix             write2file all ave/time %(STEPFREQ)d 1 %(STEPFREQ)d c_tpe file %(DELTAUFILE)s
'''
    d_before = {'XID': NB+NW, 'STEPFREQ': StepFreq, 'DELTAUFILE': MDPrefix + '_tpe.dat'}
    
    s_after = '''
unfix write2file
uncompute tpe
uncompute pertpe
''' 
    
    # write lammps input script
    LammpsFiles, TrajFile = sim.export.lammps_tsanyal.MakeLammpsMD(Sys, NStepsMin = MinSteps, NStepsEquil = EquilSteps,
                                                           NStepsProd = ProdSteps, WriteFreq = StepFreq, 
                                                           Prefix = MDPrefix, TrajFile = ".lammpstrj.gz",
                                                           LammpsCommandsBefore = s_before % d_before, 
                                                           LammpsCommandsAfter = s_after)
    
    # write jobscript and submit
    s_job = '''
#!/bin/bash
#
#$ -V
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -N %(JOBNAME)s
date
%(LAMMPSEXEC)s -in %(INFILE)s -log %(LOGFILE)s
'''
    d = {'JOBNAME': MDPrefix, 'LAMMPSEXEC': LammpsExec, 'INFILE': LammpsFiles[0], 'LOGFILE': MDPrefix + '.log'}
    JobFile = MDPrefix + '.sh'
    file(JobFile, 'w').write(s_job % d)
    if AutoSubmit: os.system('qsub %s' % JobFile) # need to submit from head node

    
def compute_Mu():
    global Prefix, lambda_factors
    
    # initialize arrays
    Nlambda = len(lambda_factors)
    dudl = []
    var_dudl = []
    mu= 0.0
    err = 0.0
    
    # reject data recorded during equilbration
    start = int(ProdSteps / StepFreq)
    
    for i, lambda_factor in enumerate(lambda_factors):
        if lambda_factor == 0.0:
            dudl.append(0.0)
            var_dudl.append(0.0)
        else:
            DeltaUFile = '%s_lambda%d_tpe.dat' % (Prefix, i)
            u_lambda = np.loadtxt(DeltaUFile)[-start:, 1]
            dudl.append(np.mean(u_lambda))
            var_dudl.append(np.std(u_lambda, ddof = 1))
    
    # calculate mu and err using trapezoidal rule
    dudl = np.array(dudl)
    var_dudl = np.array(var_dudl)
    print dudl, var_dudl
    for i in range(Nlambda):
        if i == 0 or i == Nlambda - 1:
            mu += 0.5 * dudl[i]
            err += 0.25 * var_dudl[i]
        else:
            mu += dudl[i]
            err += var_dudl[i]
    err /= np.sqrt(err)
    
    ret = (mu, err)
    OutPickle = Prefix + '_mu.pickle'
    pickle.dump(ret, open(OutPickle, 'w'))
    
    print 'mu = %g kcal/mol' % mu
    print 'err = %g kcal/mol' % mu
    

def main(Mode = 'test'):
    global NB, NW, AATraj, TPType, FF_File, Prefix
    global AutoSubmit
    
    if Mode == 'test' or Mode == 'lammps':
        AutoSubmit = (Mode == 'lammps')
        for i, l in enumerate(lambda_factors):
            Sys = makeTISys(TPType = TPType, lambda_factor = l)
            write2Lammps(Sys, lambda_index = i)
    
    if Mode == 'mu':
        compute_Mu()





######## MAIN ########
if __name__ == '__main__':
    NB = int(sys.argv[1])
    NW = int(sys.argv[2])
    TPType = sys.argv[3]
    FFType = sys.argv[4]
    Mode = sys.argv[5]
    if len(sys.argv) > 6: refNB = int(sys.argv[6])
    else: refNB = 250
    
    AATraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod_mapped.lammpstrj.gz' % (NB, NW, NB, NW)
    FF_File = '/home/cask0/home/tsanyal/benwat/data/cgff/ff_ref/NB%d/NB%dNW%d_%s_ff.dat' % (refNB, refNB, 500-refNB, FFType)
    
    Prefix = 'ti_%s' % FFType
    main(Mode = Mode)
