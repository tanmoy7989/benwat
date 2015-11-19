#!/usr/bin/env python

import os, sys

sys.path.append(os.path.abspath('../../../'))
import makeGromacs as md

NB = int(sys.argv[1])
NW = int(sys.argv[2])
Ncores = int(sys.argv[3])
Temp = 300

md.NB = NB
md.NW = NW
md.TempSet = Temp

md.Ncores = Ncores

md.TIMESTEP = 0.002             #00.02 ps 
md.NEIGHCALCSTEPS = 5           #00.01 ps 
md.MINSTEPS = 50000             #100   ps
md.NPTSTEPS = 1000000           #02    ns
md.EQUILSTEPS = 1000000         #02    ns
md.PRODSTEPS = 10000000         #20    ns
md.STEPFREQ = 2000              #04    ps
md.RESTART_TIME_MINS = 30       #30    mins

Prefix = md.makePrefix()
BoxL = md.makeBoxL()
paramdict = md.makeParamDict(Prefix = Prefix, BoxL = BoxL)
print paramdict 
md.makeData(paramdict)
md.doEneMin(paramdict)
md.doNPT(paramdict)
md.doProd(paramdict)
