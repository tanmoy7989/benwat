#/usr/bin/env python
import os, sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.expanduser('~/software')) ; import mysim ; import sim
sys.path.append('../') ; import cgmodel as cg
import pickleTraj
import parse_potential as pp

# test system of 5 benzene and 10 waters
cg.doMinimize = False
cg.NB = 5 
cg.NW = 10
cg.LammpsTraj = '/home/cask0/home/tsanyal/benwat/data/gromacs/NB50NW450/NB50NW450_prod.lammpstrj.gz'
cg.Prefix = 'matchEne'

# local density potential params
cg.LDCutBW = 7.5  # can change to BB, WW, BW etc to experiment
cg.LD_Delta = 1.2
cg.RhoMin = 0
cg.RhoMax = cg.NB
testLDKnots = 2 * (1 - 0.3 * np.linspace(0, cg.NB, cg.NLDKnots))

# box length
LammpsTrj = pickleTraj(cg.LammpsTraj, Verbose = True)
cg.BoxL = LammpsTrj.FrameData['BoxL']

print 'Making system...'
Sys = cg.makeSys()
sumfile = '/home/cask0/home/tsanyal/benwat/data/cgff/NB50NW450_init/NB50NW450_SP_sum.txt' # load in converged pair splines or 0 pair potential from NB50NW450 and test LD Knots
SP_BB = Sys.ForceField[0] ; SP_WW = Sys.ForceField[1] ; SP_BW = Sys.ForceField[2]
LD = Sys.ForceField[-1]
for P in [SP_BB, SP_WW, SP_BW]: P.SetParam(Knots = [0.] * cg.NSPKnots) #P.SetParam(**pp.parseParam(sumfile, P.Name))
LD.SetParam(Knots = testLDKnots)
Sys.ForceField.Update()
Sys.TempSet = 300.0
sim.system.init.positions.CubicLatticeFill(Sys, Random = 0.1)
sim.system.init.velocities.Canonical(Sys, Temp = 298.0)
Int = Sys.Int
Sys.Measures.VerboseOutput(StepFreq = 100)

# export to lammps
sim.export.lammps.LammpsExec = 'lmp_tsanyal'
Trj, TrjFile = sim.export.lammps.MakeLammpsTraj(Sys = Sys, Prefix = cg.Prefix, NStepsMin = 1000, NStepsEquil = 10000, NStepsProd = 10000, WriteFreq = 100)

# match energies
print 'Computing energies...'
sim_pe = []
lmp_pe = []
diff = []
fracdiff = []
for (i,Pos) in enumerate(Trj):
	Sys.Arrays.Pos = Pos
	Sys.ForceField.Eval()
	sim_pe.append(Sys.PEnergy)
	lmp_pe.append(Trj.FrameData['PEnergy'])
	diff.append(np.abs(sim_pe[-1] - lmp_pe[-1]))
	fracdiff.append(diff[-1]/sim_pe[-1])

N = len(sim_pe)
fig = plt.figure()
ax1 = fig.add_subplot(2,1,1) ; ax2 = fig.add_subplot(2,1,2)
ax1.plot(range(N), sim_pe, linestyle = 'solid', lw = 2, label = 'sim')
ax1.plot(range(N), lmp_pe, linestyle = 'none', marker = 'o', markersize = 5, label = 'lmp')
ax1.set_xlabel('timestep') ; ax1.set_ylabel('total pe')

ax2.plot(range(N), diff)
#ax2.set_ylim([-1.e-6, 1.e-4])
ax2.set_xlabel('timestep') ; ax2.set_ylabel('relative error')

plt.show()







