#!/usr/bin/env python
import os, sys, numpy as np, cPickle as pickle
from scipy.stats import linregress
import matplotlib, matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

NB = int(sys.argv[1])
NW = int(sys.argv[2])
DataDir = os.path.abspath('./')


#Cut_RKBI = [10.0, 12.0]
Cut_RKBI = [12.0, 18.0]
CutDict_SKBI = {'AA': {50: 0.07, 100: 0.07, 150: 0.07, 200: 0.07, 250: 0.07, 300:0.07, 350: None, 400: None, 450: None, 10: 0.06}
               }


def linreg(Ls, G, cut):
    indices = np.where(1./Ls >= cut)[0]
    x = 1./Ls[indices] ; y = G[indices]
    slope, intercept, r, p, stderr = linregress(x,y)
    return slope, intercept

def getInfKBI(intuple):
    r, G = intuple
    ind = [list(r).index(x) for x in list(r) if x>= Cut_RKBI[0] and x<=Cut_RKBI[1]]
    G_inf = np.mean(G[ind[0]:ind[-1]])
    return G_inf


def compareKBI(NB, NW, TrajType = 'AA'):
    
    # extract data
    ret0, ret1 = pickle.load(open(pickleName, 'r'))
    Methods = ['UnCorrected', 'TailCorrected', 'GeomCorrectedExact', 'GeomCorrectedExtrapolated']
    clrs = ['red', 'blue', 'green', 'magenta']
    
    fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (10, 10))
    
    # RKBI
    ax_RKBI_BB = fig.add_subplot(3,2,1) ; ax_RKBI_WW = fig.add_subplot(3,2,3) ; ax_RKBI_BW = fig.add_subplot(3,2,5)
    for i, method in enumerate(Methods):
        RKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, TrajType, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        r_BB, G_BB = ret[0] ; r_WW, G_WW = ret[1] ; r_BW, G_BW = ret[2]
        ax_RKBI_BB.plot(r_BB, G_BB, lw = 3, color = clrs[i], label = method)
        ax_RKBI_WW.plot(r_WW, G_WW, lw = 3, color = clrs[i])
        ax_RKBI_BW.plot(r_BW, G_BW, lw = 3, color = clrs[i])
    
    # SKBI
    ax_SKBI_BB = fig.add_subplot(3,2,2) ; ax_SKBI_WW = fig.add_subplot(3,2,4) ; ax_SKBI_BW = fig.add_subplot(3,2,6)
    SKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_SKBI.pickle' % (NB, NW, TrajType))
    (L, G_BB), (L, G_WW), (L, G_BW) = pickle.load(open(SKBIPickle, 'r'))
    ax_SKBI_BB.plot(1./L, G_BB, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'BB')
    ax_SKBI_WW.plot(1./L, G_WW, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'WW')
    ax_SKBI_BW.plot(1./L, G_BW, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'BW')
    
    # estimate limiting value of SKBI from linear region
    m_BB, c_BB = linreg(L, G_BB, CutDict_SKBI[TrajType][NB])
    m_WW, c_WW = linreg(L, G_WW, CutDict_SKBI[TrajType][NB])
    m_BW, c_BW = linreg(L, G_BW, CutDict_SKBI[TrajType][NB])
    x = np.array([0] + list(1./L))
    ax_SKBI_BB.plot(x, m_BB*x + c_BB, lw = 2, color = 'k')
    ax_SKBI_WW.plot(x, m_WW*x + c_WW, lw = 2, color = 'k')
    ax_SKBI_BW.plot(x, m_BW*x + c_BW, lw = 2, color = 'k')
    
    ax_RKBI_BB.legend(loc = 'best', prop = {'size': 10})
    ax_RKBI_BB.set_ylabel(r'$G_BB (\AA^3)$', fontsize = 10)
    ax_RKBI_WW.set_ylabel(r'$G_WW (\AA^3)$', fontsize = 10)
    ax_RKBI_BW.set_ylabel(r'$G_BW (\AA^3)$', fontsize = 10)
    for ax in [ax_RKBI_BB, ax_RKBI_WW, ax_RKBI_BW]:  ax.set_xlabel(r'$R(\AA)$', fontsize = 10)
    
    for ax in [ax_SKBI_BB, ax_SKBI_WW, ax_SKBI_BW]:
        ax.set_xlabel(r'$1/R_s (\AA^{-1})$', fontsize = 10)
        ax.legend(loc = 'best', prop = {'size': 10})
        
    
    # print limiting values of various KBIs
    print '\nG_BB'
    print '-------'
    for method in Methods:
        RKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, TrajType, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        print '%s: %g' % (method, getInfKBI(ret[0]))
    print 'SKBI: %g' % (c_BB)
    
    print '\nG_WW'
    print '-------'
    for method in Methods:
        RKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, TrajType, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        print '%s: %g' % (method, getInfKBI(ret[1]))
    print 'SKBI: %g' % (c_WW)
    
    print '\nG_BW'
    print '-------'
    for method in Methods:
        RKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, TrajType, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        print '%s: %g' % (method, getInfKBI(ret[1]))
    print 'SKBI: %g' % (c_BW)
   
    
    
    fig.tight_layout()
    




compareKBI(NB, NW)
plt.show()    
