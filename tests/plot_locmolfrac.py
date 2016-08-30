#!/usr/bin/env python
import numpy as np
import os, sys, pickle
import matplotlib.pyplot as plt

sys.path.append(os.path.expanduser('~/benwat'))
import measure as m

m.MeasureFreq = 10
m.Normalize = True

fmt = os.path.expanduser('~/benwat/data/gromacs/NB%dNW%d/NB%dNW%d_prod.lammpstrj.gz')

NB = int(sys.argv[1])
NW = int(sys.argv[2])

m.Traj = fmt % (NB, NW, NB, NW)
m.Prefix = 'NB%d' % NB
m.makeLocalMolFrac(Cutoff = 7.0)
data = pickle.load(open(m.Prefix+'_locmolfrac.pickle', 'r'))
x = data[0]; hb = data[1]['B'] ; hw = data[1]['W']
plt.plot(x,hb, 'k-', linewidth = 3, label = 'Benzene')
plt.plot(x,hw, 'b-', linewidth = 3, label = 'Water')
plt.legend()
plt.title('NB = %d' % NB)
plt.show()
