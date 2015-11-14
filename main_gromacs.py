#!/usr/bin/env python

import os, sys

sys.path.append('../../')
import makeGromacs as md

NB = int(sys.argv[1])
NW = int(sys.argv[2])
Ncores = int(sys.argv[3])
Temp = 300

md.NB = NB
md.NW = NW
md.TempSet = Temp
md.Ncores = Ncores

Prefix = md.makePrefix()
BoxL = md.makeBoxL()
paramdict = md.makeParamDict(Prefix = Prefix, BoxL = BoxL)
print paramdict 
md.makeData(paramdict)
md.doEneMin(paramdict)
md.doNPT(paramdict)
md.doProd(paramdict)
