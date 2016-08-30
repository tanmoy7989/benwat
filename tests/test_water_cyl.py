#!/usr/bin/env python

import os, sys
import numpy as np

sys.path.append(os.path.expanduser('~/benwat'))
import measure as m

m.calcErrorBar = False
m.Traj = os.path.expanduser('~/benwat/data/gromacs/NB250NW250/NB250NW250_prod.lammpstrj.gz')
m.Prefix = 'NB250'
m.calcWaterCylinder()
