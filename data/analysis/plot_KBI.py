#!/usr/bin/env python
import os, sys, numpy as np, cPickle as pickle
from scipy.stats import linregress
import matplotlib, matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

NB = int(sys.argv[1])
NW = int(sys.argv[2])
DataDir = os.path.abspath('./')

Cut_RKBI = [15.0, 20.0]
CutDict_SKBI = {'AA': {50: 0.07, 100: 0.07, 150: 0.07, 200: 0.07, 250: 0.07, 300:0.07, 350: None, 400: None, 450: None, 10: 0.06}}

def linreg(Ls, G, cut):
    indices = np.where(1./Ls >= cut)[0]
    x = 1./Ls[indices] ; y = G[indices]
    slope, intercept, r, p, stderr = linregress(x,y)
    return slope, intercept

def getInfKBI(r, G):
    ind = [list(r).index(x) for x in list(r) if x>= Cut_RKBI[0] and x<=Cut_RKBI[1]]
    G_inf = np.mean(G[ind[0]:ind[-1]])
    return G_inf

def compareMethods(NB, NW, TrajType = 'AA'):
    Methods = ['UnCorrected', 'GeomCorrectedExact'] #'TailCorrected', 'GeomCorrectedExtrapolated']
    clrs = ['red', 'blue', 'green', 'magenta']
    
    # initialize inf value arrays
    G_BB_inf = np.zeros(len(Methods)+1) ; G_WW_inf = np.zeros(len(Methods)+1) ; G_BW_inf = np.zeros(len(Methods)+1)
    
    fig = plt.figure(facecolor = 'w', edgecolor = 'w')
    # RKBI
    ax_RKBI_BB = fig.add_subplot(3,2,1) ; ax_RKBI_WW = fig.add_subplot(3,2,3) ; ax_RKBI_BW = fig.add_subplot(3,2,5)
    for i, method in enumerate(Methods):
        RKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, TrajType, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        r_BB, G_BB = ret[0] ; r_WW, G_WW = ret[1] ; r_BW, G_BW = ret[2]
        ax_RKBI_BB.plot(r_BB, G_BB, lw = 3, color = clrs[i], label = method)
        ax_RKBI_WW.plot(r_WW, G_WW, lw = 3, color = clrs[i])
        ax_RKBI_BW.plot(r_BW, G_BW, lw = 3, color = clrs[i])
        # get limiting values
        G_BB_inf[i] = getInfKBI(r_BB, G_BB)
        G_WW_inf[i] = getInfKBI(r_WW, G_WW)
        G_BW_inf[i] = getInfKBI(r_BW, G_BW)
        
    # SKBI
    ax_SKBI_BB = fig.add_subplot(3,2,2) ; ax_SKBI_WW = fig.add_subplot(3,2,4) ; ax_SKBI_BW = fig.add_subplot(3,2,6)
    SKBIPickle = os.path.join(DataDir, TrajType, 'NB%dNW%d_%s_SKBI.pickle' % (NB, NW, TrajType))
    (L, G_BB, errBB), (L, G_WW, errWW), (L, G_BW, errBW) = pickle.load(open(SKBIPickle, 'r'))
    ax_SKBI_BB.errorbar(1./L, G_BB, yerr = errBB, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'BB')
    ax_SKBI_WW.errorbar(1./L, G_WW, yerr = errWW, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'WW')
    ax_SKBI_BW.errorbar(1./L, G_BW, yerr = errBW, lw = 1, marker = 'o', markersize = 6, color = 'black', label = 'BW')
    # estimate limiting value of SKBI from linear region
    m_BB, c_BB = linreg(L, G_BB, CutDict_SKBI[TrajType][NB])
    m_WW, c_WW = linreg(L, G_WW, CutDict_SKBI[TrajType][NB])
    m_BW, c_BW = linreg(L, G_BW, CutDict_SKBI[TrajType][NB])
    G_BB_inf[-1] = c_BB
    G_WW_inf[-1] = c_WW
    G_BW_inf[-1] = c_BW
    x = np.array([0] + list(1./L))
    ax_SKBI_BB.plot(x, m_BB*x + c_BB, lw = 2, color = 'k')
    ax_SKBI_WW.plot(x, m_WW*x + c_WW, lw = 2, color = 'k')
    ax_SKBI_BW.plot(x, m_BW*x + c_BW, lw = 2, color = 'k')
    
    # design
    ax_RKBI_BB.set_title('RKBI', fontsize = 10)
    ax_RKBI_BB.legend(loc = 'best', prop = {'size': 10})
    ax_RKBI_BB.set_ylabel(r'$G_BB (\AA^3)$', fontsize = 10)
    ax_RKBI_WW.set_ylabel(r'$G_WW (\AA^3)$', fontsize = 10)
    ax_RKBI_BW.set_ylabel(r'$G_BW (\AA^3)$', fontsize = 10)
    for ax in [ax_RKBI_BB, ax_RKBI_WW, ax_RKBI_BW]:  ax.set_xlabel(r'$R(\AA)$', fontsize = 10)
    
    ax_SKBI_BB.set_title('SKBI', fontsize = 10)
    for ax in [ax_SKBI_BB, ax_SKBI_WW, ax_SKBI_BW]:
        ax.set_xlabel(r'$1/R_s (\AA^{-1})$', fontsize = 10)
        ax.legend(loc = 'best', prop = {'size': 10})
        
    # print limiting values of various KBIs
    print '\nG_BB'
    print '-------'
    for i, method in enumerate(Methods): print '%s: %g' % (method, G_BB_inf[i])
    print 'SKBI: %g' % (G_BB_inf[-1])
    
    print '\nG_WW'
    print '-------'
    for i, method in enumerate(Methods): print '%s: %g' % (method, G_WW_inf[i])
    print 'SKBI: %g' % (G_WW_inf[-1])
    
    print '\nG_BW'
    print '-------'
    for i, method in enumerate(Methods): print '%s: %g' % (method, G_BW_inf[i])
    print 'SKBI: %g' % (G_BW_inf[-1])
   

def compareKBI(NB = 10, NW = 4990, RKBI_Method = 'GeomCorrectedExact'):
    CGTrajTypes = ['SP', 'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 'SPLD_BB_WW', 'SPLD_all'] #, 'SPLD_WB', 'SPLD_BB_WW_BW']
    TrajTypeLabels = ['SP', 'SPLD-BB', 'SPLD-WW', 'SPLD-BW', 'SPLD-BB-WW', 'SPLD-all'] #, 'SPLD-WB', 'SPLD-BB-WW-BW']
    clrs = ['red', 'blue', 'green', 'magenta', 'cyan', 'yellow']
    markers = ['^', '1', '*', '.', '|', '3']
    
    # initialize inf value arrays
    G_BB_inf = np.zeros([len(CGTrajTypes)+1, 2]) ; G_WW_inf = np.zeros([len(CGTrajTypes)+1, 2]) ; G_BW_inf = np.zeros([len(CGTrajTypes)+1, 2])
    
    fig = plt.figure(facecolor = 'w', edgecolor = 'w')
    # RKBI
    method = RKBI_Method
    ax_RKBI_BB = fig.add_subplot(3,2,1) ; ax_RKBI_WW = fig.add_subplot(3,2,3) ; ax_RKBI_BW = fig.add_subplot(3,2,5)
    for i, trajtype in enumerate(CGTrajTypes):
        if i == 0:
            RKBIPickle = os.path.join(DataDir, 'AA', 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, 'AA', method))
            ret  = pickle.load(open(RKBIPickle, 'r'))
            r_BB, G_BB = ret[0] ; r_WW, G_WW = ret[1] ; r_BW, G_BW = ret[2]
            ax_RKBI_BB.plot(r_BB, G_BB, color = 'black', marker = 'o', lw = 0, label = 'AA')
            ax_RKBI_WW.plot(r_BB, G_WW, color = 'black', marker = 'o', lw = 0)
            ax_RKBI_BW.plot(r_BB, G_BW, color = 'black', marker = 'o', lw = 0)
            # get limiting AA values
            G_BB_inf[0,0] = getInfKBI(r_BB, G_BB)
            G_WW_inf[0,0] = getInfKBI(r_WW, G_WW)
            G_BW_inf[0,0] = getInfKBI(r_BW, G_BW)
        
        RKBIPickle = os.path.join(DataDir, trajtype, 'NB%dNW%d_%s_RKBI_%s.pickle' % (NB, NW, trajtype, method))
        ret  = pickle.load(open(RKBIPickle, 'r'))
        r_BB, G_BB = ret[0] ; r_WW, G_WW = ret[1] ; r_BW, G_BW = ret[2]
        ax_RKBI_BB.plot(r_BB, G_BB, color = clrs[i], lw = 2, label = TrajTypeLabels[i])
        ax_RKBI_WW.plot(r_BB, G_WW, color = clrs[i], lw = 2)
        ax_RKBI_BW.plot(r_BB, G_BW, color = clrs[i], lw = 2)
        # get CG limiting values
        G_BB_inf[i+1, 0] = getInfKBI(r_BB, G_BB)
        G_WW_inf[i+1, 0] = getInfKBI(r_WW, G_WW)
        G_BW_inf[i+1, 0] = getInfKBI(r_BW, G_BW)
            
    # SKBI
    ax_SKBI_BB = fig.add_subplot(3,2,2) ; ax_SKBI_WW = fig.add_subplot(3,2,4) ; ax_SKBI_BW = fig.add_subplot(3,2,6)
    for i, trajtype in enumerate(CGTrajTypes):
        if i == 0:
            SKBIPickle = os.path.join(DataDir, 'AA', 'NB%dNW%d_%s_SKBI.pickle' % (NB, NW, 'AA'))
            (L, G_BB, errBB), (L, G_WW, errWW), (L, G_BW, errBW) = pickle.load(open(SKBIPickle, 'r'))
            ax_SKBI_BB.plot(1./L, G_BB, lw = 0, marker = 'o', color = 'black')
            ax_SKBI_WW.plot(1./L, G_WW, lw = 0, marker = 'o', color = 'black')
            ax_SKBI_BW.plot(1./L, G_BW, lw = 0, marker = 'o', color = 'black')
            # get limiting AA values
            G_BB_inf[0,1] = linreg(L, G_BB, CutDict_SKBI['AA'][NB])[1]
            G_WW_inf[0,1] = linreg(L, G_WW, CutDict_SKBI['AA'][NB])[1]
            G_BW_inf[0,1] = linreg(L, G_BW, CutDict_SKBI['AA'][NB])[1]
        
        SKBIPickle = os.path.join(DataDir, trajtype, 'NB%dNW%d_%s_SKBI.pickle' % (NB, NW, trajtype))
        (L, G_BB, errBB), (L, G_WW, errWW), (L, G_BW, errBW) = pickle.load(open(SKBIPickle, 'r'))
        ax_SKBI_BB.plot(1./L, G_BB, lw = 2, color = clrs[i])
        ax_SKBI_WW.plot(1./L, G_WW, lw = 2, color = clrs[i])
        ax_SKBI_BW.plot(1./L, G_BW, lw = 2, color = clrs[i])
        # get limiting CG values
        G_BB_inf[i+1, 1] = linreg(L, G_BB, CutDict_SKBI['AA'][NB])[1]
        G_WW_inf[i+1, 1] = linreg(L, G_WW, CutDict_SKBI['AA'][NB])[1]
        G_BW_inf[i+1, 1] = linreg(L, G_BW, CutDict_SKBI['AA'][NB])[1]
	
    # design
    ax_RKBI_BB.set_title('RKBI', fontsize = 10)
    ax_RKBI_BB.legend(loc = 'best', prop = {'size': 8})
    ax_RKBI_BB.set_ylabel(r'$G_BB (\AA^3)$', fontsize = 10)
    ax_RKBI_WW.set_ylabel(r'$G_WW (\AA^3)$', fontsize = 10)
    ax_RKBI_BW.set_ylabel(r'$G_BW (\AA^3)$', fontsize = 10)
    for ax in [ax_RKBI_BB, ax_RKBI_WW, ax_RKBI_BW]:  ax.set_xlabel(r'$R(\AA)$', fontsize = 10)
    
    ax_SKBI_BB.set_title('SKBI', fontsize = 10)
    for ax in [ax_SKBI_BB, ax_SKBI_WW, ax_SKBI_BW]:
        ax.set_xlabel(r'$1/R_s (\AA^{-1})$', fontsize = 10)
        ax.legend(loc = 'best', prop = {'size': 10})
 
    # print limiting values
    fmtstring_AA = '%s \t\t\t %g \t\t %g'
    fmtstring_CG = '%s \t\t %g \t\t %g'
    print '\nG_BB'
    print '------'
    print '\t\t\t RKBI \t\t\t SKBI'
    for i, trajtype in enumerate(['AA'] + CGTrajTypes):
        s = fmtstring_AA if trajtype == 'AA' or trajtype == 'SP' else fmtstring_CG
        print s % (trajtype, G_BB_inf[i,0], G_BB_inf[i,1])
    
    # print limiting values
    print '\nG_WW'
    print '------'
    print '\t\t\t RKBI \t\t\t SKBI'
    for i, trajtype in enumerate(['AA'] + CGTrajTypes):
        s = fmtstring_AA if trajtype == 'AA' or trajtype == 'SP' else fmtstring_CG
        print s % (trajtype, G_WW_inf[i,0], G_WW_inf[i,1])         

    # print limiting values
    print '\nG_BW'
    print '------'
    print '\t\t\t RKBI \t\t\t SKBI'
    for i, trajtype in enumerate(['AA'] + CGTrajTypes):
        s = fmtstring_AA if trajtype == 'AA' or trajtype == 'SP' else fmtstring_CG
        print s % (trajtype, G_BW_inf[i,0], G_BW_inf[i,1])


#### MAIN ####
compareMethods(NB = NB, NW = NW)
#compareKBI(NB = NB, NW = NW, RKBI_Method = 'GeomCorrectedExact')
plt.show()    
