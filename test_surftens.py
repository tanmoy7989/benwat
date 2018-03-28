#!/usr/bin/env python

import os, sys, numpy as np
import matplotlib.pyplot as plt

refs = [50, 250, 450]
fftypes = ['SP', 'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 'SPLD_WB', 'SPLD_all']
ffclrs = ['red', 'gray', 'blue', 'yellow', 'cyan', 'green']

Lx = Ly = 25.7301
Lz = 129.6577
BoxVol = Lx*Ly*Lz
P_id = (1380 * 0.6) / BoxVol
conv_factor = (4.184*1e5) / 6.023 # convert press to atm
P_id *= conv_factor 

def getVals(presstensorfile):
    data = np.loadtxt(presstensorfile)
    pxx = data[:,1]
    pyy = data[:,2]
    pzz = data[:,3]
    p_tot = np.mean(pzz)
    #p_tot += P_id
    gamma = 0.01 * (Lz/2) * (np.mean(pzz) - 0.5*(np.mean(pxx) + np.mean(pyy)))
    return gamma, p_tot
    
fig = plt.figure(facecolor = 'w', edgecolor = 'w')
axs = []
ind = 1
for ref in refs:
    ax = fig.add_subplot(3,1,ind)
    axs.append(ax)
    xB = ref / float(500)
    gamma = []
    p_tot = []
    datadir = os.path.abspath('./data/analysis_new/ref_NB%d/longbox_final' % ref)
    for i, fftype in enumerate(fftypes):
        presstensorfile = os.path.join(datadir, 'NB380_gamma_%s_press.dat' % fftype)
        this_gamma, this_p_tot = getVals(presstensorfile)
        gamma.append(this_gamma)
        p_tot.append(this_p_tot)
        ax.plot(this_p_tot, this_gamma, 'k-', marker = 'o', markersize = 8, markerfacecolor = ffclrs[i], label = fftype)
    ax.set_title('refB = %g' % xB)
    ind += 1
    
axs[-1].set_xlabel('Pressure (atm)')
axs[1].set_ylabel('Surf tens (mN/m)')
axs[1].legend(prop = {'size': 8})
    
plt.tight_layout()
plt.show()
        
        
