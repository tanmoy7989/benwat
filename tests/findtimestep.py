#/usr/bin/env python
import os, sys
import numpy as np

sys.path.append(os.path.expanduser('~/software')) ; import mysim ; import sim
sys.path.append('../') ; import cgmodel as cg
import pickleTraj
import parse_potential as pp

# build test system of 50 benzene and 450 waters
print 'Making system...'
cg.NB = 50 
cg.NW = 450
cg.LammpsTraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB50NW450/NB50NW450_prod.lammpstrj.gz'
cg.Prefix = 'findtimestep'
Trj = pickleTraj(cg.LammpsTraj, Verbose = True)
cg.BoxL = Trj.FrameData['BoxL']
Sys = cg.makeSys()

# load in converged pair splines
sumfile = '/home/cask0/home/tsanyal/benwat/data/cgff/NB50NW450_init/NB50NW450_SP_sum.txt'
SP_BB = Sys.ForceField[0] ; SP_WW = Sys.ForceField[1] ; SP_BW = Sys.ForceField[2]
for P in [SP_BB, SP_WW, SP_BW]: P.SetParam(**pp.parseParam(sumfile, P.Name))
Sys.ForceField.Update()

# other params
Sys.TempSet = 300.0
sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = 298.0)
Int = Sys.Int
Sys.Measures.VerboseOutput(StepFreq = 100)

print Sys.ForceField.ParamString()

# find timestep

sim.export.lammps.LammpsExec = 'lmp_tsanyal'
sim.integrate.velverlet.FindTimeStepLammps(Sys, NSteps = 1e6, GuessTimeStep = 2.e-3) # guess delta_t from C.Peter's data
