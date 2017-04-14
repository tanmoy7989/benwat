#usr/bin/env python

import os, sys
import numpy as np
import cgmodel as cg

convtype = sys.argv[1]

VMDExec = 'vmd'
TCLFile_xtc = 'conv_xtc.tcl'
TCLFile_xyz = 'conv_xyz.tcl'

def conv_xtc(NB, NW):
    # VMD tcl script (selects only O atoms from water and all atoms from benzene)
    tclscript = '''
mol new %(gro)s
mol addfile %(xtc)s waitfor all
set sel [atomselect top "not name HW1 and not name HW2"]
animate write lammpstrj %(lammpstrj)s waitfor all sel $sel
quit
    '''
    GromacsDir = '/home/cask0/home/tsanyal/benwat/data/gromacs'
    GromacsStruct = os.path.join(GromacsDir, 'NB%dNW%d' % (NB, NW), 'NB%dNW%d_equil2.gro' % (NB, NW))
    GromacsTrj = os.path.join(GromacsDir, 'NB%dNW%d' % (NB, NW), 'NB%dNW%d_prod.xtc' % (NB,NW))
    
    LammpsTraj = GromacsTrj.split('.')[0] + '.lammpstrj.gz'
    
    #VMDParams = {'gro' : GromacsStruct, 'xtc': GromacsTrj,'lammpstrj': LammpsTraj}
    #file(TCLFile_xtc, 'w').write(tclscript % VMDParams)
    #os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile_xtc))
    #os.system('gzip %s' % LammpsTraj)
    
    cg.NB = NB
    cg.NW = NW
    cg.mapTrj(LammpsTraj)


def conv_xyz(conc):
    ## VMD tcl script
    tclscript = '''
mol addfile %(xyz)s waitfor all
package require pbctools
pbc set {%(boxl)g %(boxl)g %(boxl)g} -all
animate write lammpstrj %(lammpstrj)s waitfor all
quit
    '''
    xyzDir = '/home/cask0/home/tsanyal/benwat/data/Christine_Peter_data/AA_reftraj'
    GromacsDir = '/home/cask0/home/tsanyal/benwat/data/gromacs'
    reftable = np.loadtxt(os.path.join(xyzDir, 'table.txt'))
    
    for i in range(len(reftable)):
	   if reftable[i,0] == float(conc):
	       NB = reftable[i,1]
	       NW = reftable[i,2]
	       BoxL = reftable[i,3]

    xyzTraj = os.path.join(xyzDir, 'conc_%s' % conc, 'CG_trajectory_%s.xyz' % conc)
    LammpsTraj = os.path.join(GromacsDir, 'NB%dNW%d' % (NB,NW), 'NB%dNW%d_prod.lammpstrj' % (NB, NW))
    if not os.path.isdir(os.path.join(GromacsDir, 'NB%dNW%d' % (NB, NW))): os.mkdir(os.path.join(GromacsDir, 'NB%dNW%d' % (NB, NW)))
    
    #VMDParams = {'xyz' : xyzTraj, 'lammpstrj': LammpsTraj, 'boxl': BoxL}
    #file(TCLFile_xyz, 'w').write(tclscript % VMDParams)
    #os.system('%s -dispdev text -e %s' % (VMDExec, TCLFile_xyz))
    #os.system('gzip %s' % LammpsTrj)
    
    cg.NB = NB
    cg.NW = NW
    cg.mapTrj(LammpsTraj+'.gz')

if convtype == 'xtc':
    for NB in [200, 300, 350, 400, 450, 500]:
        NW = 500 - NB
        print 'NB = %d, NW = %d\n' % (NB, NW)
        conv_xtc(NB, NW)

if convtype == 'xyz':
    conc = [0.0019, 0.0038, 0.0116] #[0.0019, 0.0038, 0.0057, 0.0076, 0.0095, 0.0076, 0.0116, 0.0526]
    for c in conc:
        print 'Concentration = %g' % c
        conv_xyz(c)
        
for i in [TCLFile_xtc, TCLFile_xyz]:
    if os.path.isdir(i): os.remove(i)
