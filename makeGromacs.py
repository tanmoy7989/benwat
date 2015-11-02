#!/usr/bin/env python
import os, sys
import numpy as np


DelTempFiles = True
NB = 10
NW = 500
TempSet = 300

# iteration steps
CALCSTEPS = 5
MINSTEPS = 1000
NPTSTEPS = 1000
NVTSTEPS = 1000
PRODSTEPS = 1000
STEPFREQ = 100
RESTARTFREQ = 100

# source benzene gro file (Christine Peter)
benzene_gro = '''Benzene
12
    1PHE    CD1    1   3.394   2.743   3.898  0.7515  0.0011 -0.1349
    1PHE    HD1    2   3.451   2.814   3.959  2.7392 -0.3604 -1.4893
    1PHE    CE1    3   3.472   2.668   3.810  0.1297  0.2876 -0.9463
    1PHE    HE1    4   3.563   2.712   3.770  0.8472 -1.1164 -0.9008
    1PHE     CZ    5   3.426   2.549   3.755  0.0697  0.2074 -0.7270
    1PHE     HZ    6   3.495   2.492   3.693  1.2402 -0.8774  1.4960
    1PHE    CE2    7   3.294   2.515   3.780  0.2128  0.1484 -0.0516
    1PHE    HE2    8   3.242   2.435   3.727 -2.2415  1.4410  0.2523
    1PHE    CD2    9   3.218   2.588   3.871  0.5985 -0.1622  0.5283
    1PHE    HD2   10   3.112   2.562   3.877  0.7848 -1.5420 -1.3155
    1PHE     CG   11   3.265   2.705   3.930  0.7087  0.1364 -0.1418
    1PHE     HG   12   3.197   2.781   3.968  1.4071 -0.7699  3.2040
'''

# source water pdb file generated spc216.gro
water_gro = '''H20
3 
    1SOL     OW    1    .230    .628    .113
    1SOL    HW1    2    .137    .626    .150
    1SOL    HW2    3    .231    .589    .021
'''

# gromacs topology file (Christine Peter)
topology = '''
; Include forcefield parameters
#include "gromos53a6.ff/forcefield.itp"

[ moleculetype ]
; Name nrexcl
Protein     3

[ atoms ]
;   nr      type  resnr resid  atom  cgnr   charge     mass    total_charge
     1         C     1  PHE     CD1     1   -0.14    12.011     ; -0.140
     2        HC     1  PHE     HD1     1    0.14     1.008     ;  0.000
     3         C     1  PHE     CE1     2   -0.14    12.011     ; -0.140   
     4        HC     1  PHE     HE1     2    0.14     1.008     ;  0.000   
     5         C     1  PHE     CZ      3   -0.14    12.011     ; -0.140   
     6        HC     1  PHE     HZ      3    0.14     1.008     ;  0.000   
     7         C     1  PHE     CE2     4   -0.14    12.011     ; -0.140   
     8        HC     1  PHE     HE2     4    0.14     1.008     ;  0.000   
     9         C     1  PHE     CD2     5   -0.14    12.011     ; -0.140   
    10        HC     1  PHE     HD2     5    0.14     1.008     ;  0.000   
    11         C     1  PHE      CG     6   -0.14    12.011     ; -0.140   
    12        HC     1  PHE      HG     6    0.14     1.008     ;  0.000   

[ bonds ]
; ai  aj funct    c0        c1
   1   2   2    0.1090  1.2300e+07 ; gb_3  (C-CH)
   1   3   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
   1  11   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
   3   4   2    0.1090  1.2300e+07 ; gb_3  (C-CH)
   3   5   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
   5   6   2    0.1090  1.2300e+07 ; gb_3  (C-CH)
   5   7   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
   7   8   2    0.1090  1.2300e+07 ; gb_3  (C-CH)
   7   9   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
   9  10   2    0.1090  1.2300e+07 ; gb_3  (C-CH)
   9  11   2    0.1390  1.0800e+07 ; gb_16 (CR1-CR2, 6-ring)
  11  12   2    0.1090  1.2300e+07 ; gb_3  (C-CH)

[ pairs ]
; ai  aj funct
   1   6   1
   1   7   1   
   1  10   1    
   2   4   1   
   2   5   1   
   2   9   1  
   2  12   1   
   3   8   1   
   3   9   1   
   3  12   1   
   4   6   1   
   4   7   1   
   4  11   1   
   5  10   1   
   5  11   1   
   6   8   1   
   6   9   1   
   7  12   1   
   8  10   1  
   8  11   1   
  10  12   1   

[ angles ]
; ai  aj  ak funct   c0      c1 
   2   1   3   2    120.0   505.0 ; ga_25 (HC - 6-ring)   
   2   1  11   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   3   1  11   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H)
   1   3   4   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   1   3   5   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H)
   4   3   5   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   3   5   6   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   3   5   7   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H)
   6   5   7   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   5   7   8   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   5   7   9   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H)
   8   7   9   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   7   9  10   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   7   9  11   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H)
  10   9  11   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   1  11   9   2    120.0   560.0 ; ga_27 (CR1 - 6-ring, no H) 
   1  11  12   2    120.0   505.0 ; ga_25 (HC - 6-ring)
   9  11  12   2    120.0   505.0 ; ga_25 (HC - 6-ring)

[ dihedrals ]
; ai  aj  ak  al funct   angle  fc
   1  11   3   2   2      0.0  167.42309 ; gi_1 (planar group)   
   3   1   5   4   2      0.0  167.42309 ; gi_1 (planar group)
   5   3   7   6   2      0.0  167.42309 ; gi_1 (planar group)
   7   5   9   8   2      0.0  167.42309 ; gi_1 (planar group)
   9   7  11  10   2      0.0  167.42309 ; gi_1 (planar group)
  11   1   9  12   2      0.0  167.42309 ; gi_1 (planar group)
   1   3   5   7   2      0.0  167.42309 ; gi_1 (planar group)
   3   5   7   9   2      0.0  167.42309 ; gi_1 (planar group)
   5   7   9  11   2      0.0  167.42309 ; gi_1 (planar group)
   7   9  11   1   2      0.0  167.42309 ; gi_1 (planar group)
   9  11   1   3   2      0.0  167.42309 ; gi_1 (planar group)
  11   1   3   5   2      0.0  167.42309 ; gi_1 (planar group)

; Include Position restraint file
#ifdef POSRES
#include "posre.itp"
#endif

; Include water topology
#include "spce.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include generic topology for ions
#include "ions.itp"

[system]
Benzene in water

[molecules]
Protein          %(NB)d
SOL              %(NW)d
'''

#packmol input script
packmol_script = '''
tolerance 2.0
filetype pdb
output %(Prefix)s.pdb
add_amber_ter

structure benzene.pdb
  number %(NB)d
  resnumbers 2
  centerofmass
  inside box -%(Len)d -%(Len)d -%(Len)d %(Len)d %(Len)d %(Len)d
end structure

structure water.pdb
  number %(NW)d
  resnumbers 2
  inside box -%(Len)d -%(Len)d -%(Len)d %(Len)d %(Len)d %(Len)d
end structure
''' 

# common params (Christine Peter)
common_params = '''
; neighbor searching
nstcalclr = %(calcsteps)d
nstlist = %(calcsteps)d
nstcomm = %(calcsteps)d
comm-mode = Linear
ns_type = grid
cutoff-scheme = group
verlet-buffer-drift  = 0
rlist = 1
rlistlong = 1.4

; electrostatics
coulombtype = PME
coulomb-modifier = None
rcoulomb-switch = 0
rcoulomb = 1.0
epsilon-r = 1
epsilon-rf = 1

; ewald summation
fourierspacing = 0
fourier-nx = 54
fourier-ny = 54
fourier-nz = 54
pme-order = 4
ewald-rtol  = 1e-05
ewald-geometry = 3d
optimize-fft = yes

; van der Waals
vdwtype = Cut-off
vdw-modifier = None
rvdw-switch = 0
rvdw = 1.4

; constraints
constraint-algorithm = LINCS
shake-tol = 0.0001
lincs-order = 4
lincs-iter = 1
lincs-warnangle = 30

; box 
pbc = xyz

; time-stepping
dt = 0.002

; output frequency settings
nstcalcenergy = %(calcsteps)d
nstenergy = %(stepfreq)d
nstlog = %(stepfreq)d
nstxout = %(stepfreq)d
nstvout = %(stepfreq)d
nstfout = 0
nstxtcout = %(stepfreq)d
'''

# energy minimization
minim_mdp = '''
integrator = steep
emtol = 10
emstep = 0.01
nstcgsteep = 1000
nsteps = %(minsteps)d
'''

# nvt equlibration
nvt_mdp = '''
integrator = md
nsteps %(nvtsteps)d

nsttcouple = %(calcsteps)d
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = %(TempSet)d %(TempSet)d

pcoupl = no
DispCorr = no

continuation = no
gen_vel = yes
gen_temp = %(TempSet)g
gen_seed = -1
'''

# npt equlibration
npt_mdp = '''
;define = -DPOSRES

integrator = md
nsteps %(nptsteps)d

nsttcouple = %(calcsteps)d
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = %(TempSet)d %(TempSet)d

nstpcouple = -1
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 1.0 1.0
ref-p = 0.0 0.0
compressibility = 4.5e-5
refcoord_scaling = no
DispCorr = no

continuation = yes
gen_vel = no
'''

# production
prod_mdp = '''
integrator = md
nsteps %(prodsteps)d

nsttcouple = %(calcsteps)d
tcoupl = v-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = %(TempSet)d %(TempSet)d

nstpcouple = -1
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 1.0 1.0
ref-p = 0.0 0.0
compressibility = 4.5e-5
refcoord_scaling = no
DispCorr = no

continuation = yes
gen_vel = no
'''


def makePrefix():
    Prefix = 'NB%dNW%d' % (NB, NW)
    return Prefix

def getBoxL():
    #specific vol of water
    N_A = 6.023e23
    v = (18./N_A) * (1e-2 * 1e9)**3.
    BoxVol = v * NW
    BoxL = BoxVol ** (1./3.)
    return BoxL

def genData(BoxL = 3., Prefix = 'benwat'):
    BoxL *= 10.
    file('benzene.gro', 'w').write(benzene_gro)
    file('water.gro', 'w').write(water_gro)
    file('solvate.inp', 'w').write(packmol_script % {'Prefix': Prefix, 'NB': NB, 'NW': NW, 'Len': BoxL/2.})
    file('%s.top' % Prefix, 'w'). write(topology % {'NB': NB, 'NW': NW})
    
    cmdstring = '''
editconf -f benzene.gro -o benzene.pdb
editconf -f water.gro -o water.pdb
packmol < solvate.inp
editconf -f %s.pdb -o %s.gro -d 1.0 -bt cubic
''' % (Prefix, Prefix)

    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['benzene.gro', 'water.gro', 'solvate.inp', 'benzene.pdb', 'water.pdb', '%s.pdb' % Prefix]
        [os.remove(x) for x in filenames]


def doEneMin(Prefix = 'benwat'):
    s =(common_params + minim_mdp) % {'calcsteps': CALCSTEPS, 'minsteps': MINSTEPS, 'stepfreq': STEPFREQ}
    file('%s_minim.mdp' % Prefix, 'w').write(s)
    cmdstring = '''
grompp -f %s_minim.mdp -c %s.gro -p %s.top -o %s_minim.tpr
mdrun -v -nt 1 -deffnm %s_minim
''' % ((Prefix,) * 5)
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%s_minim.mdp' % Prefix]
        [os.remove(x) for x in filenames]


def doNVT(Prefix = 'benwat'):
    s = (common_params + nvt_mdp) % {'calcsteps': CALCSTEPS, 'nvtsteps': NVTSTEPS, 'stepfreq': STEPFREQ, 'TempSet': TempSet}
    file('%s_nvt.mdp' % Prefix, 'w').write(s)
    cmdstring = '''
grompp -f %s_nvt.mdp -c %s_minim.gro -p %s.top -o %s_nvt.tpr
mdrun -v -nt 1 -deffnm %s_nvt
''' % ((Prefix, ) * 5)
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%s_nvt.mdp' % Prefix]
        [os.remove(x) for x in filenames]


def doNPT(Prefix = 'benwat'):
    s = (common_params + npt_mdp) % {'calcsteps': CALCSTEPS, 'nptsteps': NPTSTEPS, 'stepfreq': STEPFREQ, 'TempSet': TempSet}
    file('%s_npt.mdp' % Prefix, 'w').write(s)
    cmdstring = '''
grompp -f %s_npt.mdp -c %s_nvt.gro -p %s.top -o %s_npt.tpr
mdrun -v -nt 1 -deffnm %s_npt
''' % ((Prefix, ) * 5)
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%s_npt.mdp' % Prefix]
        [os.remove(x) for x in filenames]


def doProd(Prefix = 'benwat'):
    s = (common_params + prod_mdp) % {'calcsteps': CALCSTEPS, 'prodsteps': PRODSTEPS, 'stepfreq': STEPFREQ, 'TempSet': TempSet}
    file('%s_prod.mdp' % Prefix, 'w').write(s)
    cmdstring = '''
grompp -f %s_prod.mdp -c %s_npt.gro -p %s.top -o %s_prod.tpr
mdrun -v -nt 1 -deffnm %s_prod
''' % ((Prefix, ) * 5)
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%s_prod.mdp' % Prefix]
        [os.remove(x) for x in filenames]
                                                                                               
    

Prefix = makePrefix()
BoxL = getBoxL(); print BoxL
genData(BoxL = BoxL, Prefix = Prefix)
doEneMin(Prefix = Prefix)
doNVT(Prefix = Prefix)
doNPT(Prefix = Prefix)
doProd(Prefix = Prefix)
