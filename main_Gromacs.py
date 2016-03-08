#!/usr/bin/env python

import os, sys

sys.path.append('/home/cask0/home/tsanyal/benwat')
import makeGromacs as md

NB = int(sys.argv[1])
NW = int(sys.argv[2])
Ncores = int(sys.argv[3])
Temp = 300

md.NB = NB
md.NW = NW
md.TempSet = Temp

md.Ncores = Ncores

md.TIMESTEP = 0.002             #.002 ps i.e. 2 fs 
md.NEIGHCALCSTEPS = 5           #00.01 ps 
md.MINSTEPS = 50000             #200   ps
md.NPTSTEPS = 1000000           #02    ns
md.EQUILSTEPS = 10000000        #20    ns
md.PRODSTEPS = 20000000         #40    ns
md.STEPFREQ = 2000              #04    ps
md.RESTART_TIME_MINS = 60       #60    mins

Prefix = md.makePrefix()
BoxL = md.makeBoxL()
paramdict = md.makeParamDict(Prefix = Prefix, BoxL = BoxL)
print paramdict 

md.makeData(paramdict)
md.doEneMin1(paramdict)
md.doNPT(paramdict)
md.doEquil1(paramdict)
md.rescaleBox(paramdict)
md.doEneMin2(paramdict)
md.doEquil2(paramdict)
md.doProd(paramdict)
