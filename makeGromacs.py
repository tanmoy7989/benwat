#!/usr/bin/env python
import os, sys
import numpy as np


DelTempFiles = True

# constants
N_A = 6.023e23
Mass_B = 78.11  #g/mol
Mass_W = 18.0   #g/mol
rho_B = 0.874   #g/cc
rho_W = 1.00    #g/cc

# system parameters
NB = 20
NW = 50
TempSet = 300
Packmol_tol = 4.0

# iteration steps
# (dummy values and must be set by calling script)
TIMESTEP = 0.002 
NEIGHCALCSTEPS = 5
MINSTEPS = 100
NPTSTEPS = 10000
EQUILSTEPS = 1000
PRODSTEPS = 1000
STEPFREQ = 10
RESTART_TIME_MINS = 30

# executables and parallelization
os.system('export PATH=$PATH:/share/apps/x86_64/gromacs/bin')
os.system('source /share/apps/x86_64/gromacs/bin/GMXRC.bash')
Ncores = 8

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
tolerance %(Packmol_tol)g
discale 1.5
filetype pdb
output %(Prefix)s.pdb
add_amber_ter

structure benzene.pdb
  number %(NB)d
  resnumbers 2
  centerofmass
  inside box -%(Packmol_halfboxL)g -%(Packmol_halfboxL)g -%(Packmol_halfboxL)g %(Packmol_halfboxL)g %(Packmol_halfboxL)g %(Packmol_halfboxL)g 
end structure

structure water.pdb
  number %(NW)d
  resnumbers 2
  inside box -%(Packmol_halfboxL)g -%(Packmol_halfboxL)g -%(Packmol_halfboxL)g %(Packmol_halfboxL)g %(Packmol_halfboxL)g %(Packmol_halfboxL)g 
end structure
''' 

# common params (Christine Peter)
common_params = '''
;neighbor searching
nstcalclr = %(calcsteps)d
nstlist = %(calcsteps)d
nstcomm = %(calcsteps)d
comm-mode = Linear
ns_type = grid
cutoff-scheme = group
verlet-buffer-drift  = 0
rlist = 1
rlistlong = 1.0

;electrostatics
coulombtype = PME
coulomb-modifier = None
rcoulomb-switch = 0
rcoulomb = 1.0
epsilon-r = 1
epsilon-rf = 1

;ewald summation
fourierspacing = 0
fourier-nx = 54
fourier-ny = 54
fourier-nz = 54
pme-order = 4
ewald-rtol  = 1e-05
ewald-geometry = 3d
optimize-fft = yes

;van der Waals
vdwtype = Cut-off
vdw-modifier = None
rvdw-switch = 0
rvdw = 1.0

;constraints
constraints = all-bonds
constraint-algorithm = LINCS
shake-tol = 0.0001
lincs-order = 4
lincs-iter = 4 # increasing this leads to more minimum energy configs
lincs-warnangle = 30

;box 
pbc = xyz

; time-stepping
dt = %(timestep)g

;output frequency settings
nstcalcenergy = %(calcsteps)d
nstenergy = %(stepfreq)d
nstlog = %(stepfreq)d
;nstxout = %(stepfreq)d
;nstvout = %(stepfreq)d
;nstfout = 0
nstxtcout = %(stepfreq)d
'''

# initial energy minimization
minim1_mdp = '''
integrator = steep
emtol = 1
emstep = 0.01
nstcgsteep = 20000
nsteps = %(minsteps)d
'''

# npt
npt_mdp = '''
; v-rescale temp. coupling and berendsen press. coupling used 
; for quick and dirty equilbration

integrator = md
nsteps = %(nptsteps)d

;temp-coupling params (nose-hoover used in Nico's paper)
nsttcouple = -1
tcoupl = nose-hoover
tc-grps = Protein Non-Protein
tau-t = 0.5 0.5
ref-t = %(TempSet)g %(TempSet)g

;press-coupling params (Parrinello-Rahman used in Nico's paper)
nstpcouple = -1
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 3.0 3.0
ref-p = %(PressSet)g %(PressSet)g
compressibility = 4.5e-5
refcoord_scaling = no
DispCorr = EnerPres

;startup 
continuation = no
gen_vel = yes
gen_temp = %(TempSet)g 
gen_seed = -1
'''

# equilibration 1 (post npt)
equil1_mdp = '''
integrator = md
nsteps = %(equilsteps)d

nsttcouple = -1
tcoupl = nose-hoover
tc-grps = Protein SOL
tau-t = 0.5 0.5
ref-t = %(TempSet)g %(TempSet)g

pcoupl = no

continuation = yes
gen_vel = no
'''

# energy minimization 2 (right after rescaling box)
minim2_mdp = '''
integrator = steep
emtol = 1
emstep = 0.01
nstcgsteep = 20000
nsteps = %(minsteps)d
'''

# equilibration 2 (pre prod NVT)
equil2_mdp = '''
integrator = md
nsteps = %(equilsteps)d

nsttcouple = -1
tcoupl = nose-hoover
tc-grps = Protein SOL
tau-t = 0.5 0.5
ref-t = %(TempSet)g %(TempSet)g

pcoupl = no

continuation = yes
gen_vel = no
'''

# production
prod_mdp = '''
integrator = md
nsteps = %(prodsteps)d

nsttcouple = -1
tcoupl = nose-hoover
tc-grps = Protein SOL
tau-t = 0.5 0.5
ref-t = %(TempSet)g %(TempSet)g

pcoupl = no

continuation = yes
gen_vel = no
'''

def makePrefix():
    Prefix = 'NB%dNW%d' % (NB, NW)
    return Prefix

def makeBoxL(MolWt = Mass_W, rho = rho_W):
    #specific vol of species
    v = (MolWt/(rho*N_A)) * (1e-2 * 1e9)**3.
    BoxVol = v * NW
    BoxL = BoxVol**(1./3.) #boxlength based on packing of specific molecule (Water or benzene)
    return BoxL + 2.8      #2.8 nm is 2*cutoff sent by Christine 

def makeParamDict(BoxL, Prefix = 'benwat'):
    d = {'minsteps': MINSTEPS, 'nptsteps': NPTSTEPS, 'equilsteps': EQUILSTEPS, 'prodsteps': PRODSTEPS, 
         'calcsteps': NEIGHCALCSTEPS,'stepfreq': STEPFREQ, 'restart_time_mins': RESTART_TIME_MINS, 'timestep': TIMESTEP, 
         'Prefix': Prefix, 'BoxL': BoxL, 'Packmol_halfboxL': 0.5 *(10.*BoxL-2.), 'Packmol_tol': Packmol_tol,
         'TempSet': TempSet, 'PressSet': 1.0, 'NB': NB, 'NW': NW, 'Ncores': Ncores}

    return d
    
    
def makeData(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    file('benzene.gro', 'w').write(benzene_gro)
    file('water.gro', 'w').write(water_gro)
    file('solvate.inp', 'w').write(packmol_script % paramdict)
    file('%(Prefix)s.top' % paramdict, 'w'). write(topology % paramdict)
    
    cmdstring = '''
editconf -f benzene.gro -o benzene.pdb
editconf -f water.gro -o water.pdb
packmol < solvate.inp
editconf -f %(Prefix)s.pdb -o %(Prefix)s.gro -bt cubic -box %(BoxL)g %(BoxL)g %(BoxL)g
''' % paramdict

    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['benzene.gro', 'water.gro', 'solvate.inp', 'benzene.pdb', 'water.pdb', '%(Prefix)s.pdb' % paramdict]
        [os.remove(x) for x in filenames]


def doEneMin1(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s =(common_params + minim1_mdp) % paramdict
    file('%(Prefix)s_minim1.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_minim1.mdp -c %(Prefix)s.gro -p %(Prefix)s.top -o %(Prefix)s_minim1.tpr
mdrun -ntmpi %(Ncores)d -npme 2 -dlb yes -deffnm %(Prefix)s_minim1
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_minim1.mdp' % paramdict]
        [os.remove(x) for x in filenames]



def doNPT(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s = (npt_mdp + common_params) % paramdict
    file('%(Prefix)s_npt.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_npt.mdp -c %(Prefix)s_minim1.gro -p %(Prefix)s.top -o %(Prefix)s_npt.tpr -maxwarn 100
mdrun -nt %(Ncores)d -npme -1 -dlb no -cpt %(restart_time_mins)g -deffnm %(Prefix)s_npt
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_npt.mdp' % paramdict]
        [os.remove(x) for x in filenames]



def doEquil1(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s = (common_params + equil1_mdp) % paramdict
    file('%(Prefix)s_equil1.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_equil1.mdp -c %(Prefix)s_npt.gro -p %(Prefix)s.top -o %(Prefix)s_equil1.tpr -t %(Prefix)s_npt.cpt
mdrun -nt %(Ncores)d -npme -1 -dlb yes -cpt %(restart_time_mins)g -deffnm %(Prefix)s_equil1
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_equil1.mdp' % paramdict]
        [os.remove(x) for x in filenames]



def rescaleBox(paramdict = None):
    print 'Calculating avg. box size...'
    cmdstring = '''
g_traj -f %(Prefix)s_equil1.xtc -s %(Prefix)s_equil1.tpr -ob box.xvg << EOF
0
''' % paramdict
    os.system(cmdstring)
   
    # get average BoxL after equilibration
    BoxL = np.array([0., 0., 0.]) 
    of = open('box.xvg', 'r')
    nframes = 0
    for line in of:
        if line.startswith('#') or line.startswith('@'):
            continue
        else:
            l = line.split()
            BoxL[0] += float(l[1])
            BoxL[1] += float(l[2])
            BoxL[2] += float(l[3])
            nframes += 1
    of.close()
    BoxL /= nframes
    
    # get final BoxL after equilbration
    line = file('%(Prefix)s_equil1.gro' % paramdict, 'r').readlines()[-1].strip()
    BoxL0 = []
    [BoxL0.append(float(x)) for x in line.split()]
    
    scalefactor = (BoxL/BoxL0)
    Target_vol = np.prod(BoxL) * (1e-9) * 1e3 #liter
    Target_rho = ((Mass_B*NB + Mass_W*NW)/N_A) / Target_vol #g/L
    paramdict['scaleX'] = scalefactor[0]; paramdict['scaleY'] = scalefactor[1]; paramdict['scaleZ'] = scalefactor[2]
    paramdict['Target_rho'] = Target_rho
    cmdstring = 'editconf -f %(Prefix)s_equil1.gro -o %(Prefix)s_rescaled.gro -scale %(scaleX)g %(scaleY)g %(scaleZ)g' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles: os.remove('box.xvg')
    


def doEneMin2(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s =(common_params + minim2_mdp) % paramdict
    file('%(Prefix)s_minim2.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_minim2.mdp -c %(Prefix)s_rescaled.gro -p %(Prefix)s.top -o %(Prefix)s_minim2.tpr
mdrun -ntmpi %(Ncores)d -npme 2 -dlb yes -deffnm %(Prefix)s_minim2
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_minim2.mdp' % paramdict]
        [os.remove(x) for x in filenames]



def doEquil2(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s = (common_params + equil2_mdp) % paramdict
    file('%(Prefix)s_equil2.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_equil2.mdp -c %(Prefix)s_minim2.gro -p %(Prefix)s.top -o %(Prefix)s_equil2.tpr
mdrun -nt %(Ncores)d -npme -1 -dlb yes -cpt %(restart_time_mins)g -deffnm %(Prefix)s_equil2
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_equil2.mdp' % paramdict]
        [os.remove(x) for x in filenames]



def doProd(paramdict = None):
    if paramdict is None: raise TypeError('First populate param dict')
    s = (common_params + prod_mdp) % paramdict
    file('%(Prefix)s_prod.mdp' % paramdict, 'w').write(s)
    cmdstring = '''
grompp -f %(Prefix)s_prod.mdp -c %(Prefix)s_equil2.gro -p %(Prefix)s.top -o %(Prefix)s_prod.tpr -t %(Prefix)s_equil2.cpt
mdrun -nt %(Ncores)d -npme -1 -dlb yes -cpt %(restart_time_mins)g -deffnm %(Prefix)s_prod
''' % paramdict
    
    os.system(cmdstring)
    if DelTempFiles:
        filenames = ['mdout.mdp', '%(Prefix)s_prod.mdp' % paramdict]
        [os.remove(x) for x in filenames]
                                                                                               
    
# test
if __name__ == '__main__': 
    Prefix = 'gro_test'
    BoxL = makeBoxL()
    paramdict = makeParamDict(Prefix = Prefix, BoxL = BoxL)
    
    #makeData(paramdict)
    #doEneMin1(paramdict)
    #doNPT(paramdict)
    #doEquil1(paramdict)
    #rescaleBox(paramdict)
    #doEneMin2(paramdict)
    doEquil2(paramdict)
    doProd(paramdict)
