#!/usr/bin/env python

import os, sys, numpy as np, cPickle as pickle
import sim, measure, pickleTraj, benwatlib

# User Input
LammpsTraj = None
NB = None
NW = None
Prefix = 'KBI_test'

# Traj frames
Trj = None
AtomNames2Types = False
StepFreq = 10
FrameRange = None
NFrames = None
BoxL = None

# Histogram settings
NBlocks = 1
NBins = 1000 # fine grained histogram

def calcRDF():
    '''calculate fine grained RDFS'''
    global LammpsTraj, NB, NW, Prefix
    global Trj, StepFreq, FrameRange, NFrames, BoxL
    
    print 'Calculating fine grained RDFs...'
    measure.LammpsTraj = LammpsTraj
    measure.StepFreq = StepFreq
    measure.AtomNames2Types = AtomNames2Types
    measure.NBlocks = NBlocks
    measure.NBins = NBins
    measure.Normalize = True
    measure.__parseFrameData()
    Trj = measure.Trj ; FrameRange = measure.FrameRange ; NFrames = measure.NFrames
    BoxL = np.array(Trj.FrameData['BoxL'])
    
    rdf_BB, rdf_BB_pickle = measure.makeRDF(1, 1, Prefix = Prefix + '_rdf_BB_fine_grained')
    rdf_WW, rdf_WW_pickle = measure.makeRDF(2, 2, Prefix = Prefix + '_rdf_WW_fine_grained')
    rdf_BW, rdf_BW_pickle = measure.makeRDF(1, 2, Prefix = Prefix + '_rdf_BW_fine_grained')
    

def make_RKBI(method = 'UnCorrected'):
    ''' Families of formulae that directly integrate the rdf from definition.
    Corrects the g(r) or its integration based on heuristic formulae.
    Does not save to file, since calculation is instantaneous'''
    global LammpsTraj, NB, NW, Prefix
    global Trj, StepFreq, FrameRange, NFrames, BoxL
    
    # initialize all arrays
    h_BB = np.zeros(NBins) ; h_WW = np.zeros(NBins) ; h_BW = np.zeros(NBins)
    w_BB = lambda i: 1.0 ; w_WW = lambda i: 1.0 ; w_BW = lambda i: 1.0
    G_BB = np.zeros(NBins) ; G_WW = np.zeros(NBins) ; G_BW = np.zeros(NBins)
    
    # extract fine grained RDFs
    calcRDF()
    rdf_BB = pickle.load(open(Prefix + '_rdf_BB_fine_grained.pickle', 'r'))
    rdf_WW = pickle.load(open(Prefix + '_rdf_WW_fine_grained.pickle', 'r'))
    rdf_BW = pickle.load(open(Prefix + '_rdf_BW_fine_grained.pickle', 'r'))
    r_BB = rdf_BB[0] ; g_BB_0 = rdf_BB[1] ; dr_BB = r_BB[1] - r_BB[0]
    r_WW = rdf_WW[0] ; g_WW_0 = rdf_WW[1] ; dr_WW = r_WW[1] - r_WW[0]
    r_BW = rdf_BW[0] ; g_BW_0 = rdf_BW[1] ; dr_BW = r_BW[1] - r_BW[0]
    
    tmp_r = r_BB, r_WW, r_BW
    tmp_g = g_BB_0, g_WW_0, g_BW_0
    
    if method == 'UnCorrected':
        # from definition of KBI for open systems
        h_BB = g_BB_0 - 1
        h_WW = g_WW_0 - 1
        h_BW = g_BW_0 - 1
    
    elif method == 'TailCorrected':
        # Ganguly, van der Vegt et. al, JCTC, 2013, 9, 1347-1355, Eq (5)
        rho_B = float(NB) / np.prod(BoxL) ; rho_W = float(NW) / np.prod(BoxL)
        DeltaN_BB = rho_B * np.array( [ dr_BB * np.sum(4*np.pi * r_BB[:i]**2. * (g_BB_0[:i]-1)) for i in range(NBins) ] )
        DeltaN_WW = rho_W * np.array( [ dr_WW * np.sum(4*np.pi * r_WW[:i]**2. * (g_WW_0[:i]-1)) for i in range(NBins) ] )
        DeltaN_BW = rho_W * np.array( [ dr_BW * np.sum(4*np.pi * r_BW[:i]**2. * (g_BW_0[:i]-1)) for i in range(NBins) ] )
        Bulk_B = NB - rho_B * (4/3.)*np.pi * r_BB**3.
        Bulk_W = NW - rho_W * (4/3.)*np.pi * r_WW**3.
        h_BB = g_BB_0 * Bulk_B / (Bulk_B - DeltaN_BB - 1) - 1
        h_WW = g_WW_0 * Bulk_W / (Bulk_W - DeltaN_WW - 1) - 1
        h_BW = g_BW_0 * Bulk_W / (Bulk_W - DeltaN_BW) - 1
        
    elif method == 'GeomCorrectedExact':
        # Krueger, Schenll et. al, J.Phys.Chem.Lett, 2013, 4, 235-238, Eqn (6)
        x_BB = [np.array([0])] + [r_BB[:i] / r_BB[:i][-1] for i in range(1,NBins)]
        x_WW = [np.array([0])] + [r_WW[:i] / r_WW[:i][-1] for i in range(1,NBins)]
        x_BW = [np.array([0])] + [r_BW[:i] / r_BW[:i][-1] for i in range(1,NBins)]
        w_BB = lambda i: (1 - 1.5 * x_BB[i] + 0.5 * x_BB[i]**3.)
        w_WW = lambda i: (1 - 1.5 * x_WW[i] + 0.5 * x_WW[i]**3.)
        w_BW = lambda i: (1 - 1.5 * x_BW[i] + 0.5 * x_BW[i]**3.)
        h_BB = g_BB_0 - 1
        h_WW = g_WW_0 - 1
        h_BW = g_BW_0 - 1
    
    elif method == 'GeomCorrectedExtrapolated':
        # Krueger, Schenll et. al, J.Phys.Chem.Lett, 2013, 4, 235-238, Eqn (7)
        x_BB = [np.array([0])] + [r_BB[:i] / r_BB[:i][-1] for i in range(1,NBins)]
        x_WW = [np.array([0])] + [r_WW[:i] / r_WW[:i][-1] for i in range(1,NBins)]
        x_BW = [np.array([0])] + [r_BW[:i] / r_BW[:i][-1] for i in range(1,NBins)]
        w_BB = lambda i: (1 - x_BB[i]**3.)
        w_WW = lambda i: (1 - x_WW[i]**3.)
        w_BW = lambda i: (1 - x_BW[i]**3.)
        h_BB = g_BB_0 - 1
        h_WW = g_WW_0 - 1
        h_BW = g_BW_0 - 1
    
    # integrate
    for n in range(1, NBins):
        G_BB[n] = dr_BB * np.sum(4*np.pi * r_BB[:n]**2 * h_BB[:n] * w_BB(n))
        G_WW[n] = dr_WW * np.sum(4*np.pi * r_WW[:n]**2 * h_WW[:n] * w_WW(n))
        G_BW[n] = dr_BW * np.sum(4*np.pi * r_BW[:n]**2 * h_BW[:n] * w_BW(n))
    
    return  (r_BB, G_BB), (r_WW, G_WW), (r_BW, G_BW)



def make_SKBI():
    ''' Calculates finite sized KBIs based on particle fluctuations in subvolumes
    These can be extrpolated to find the true KBI using the method as outlined in
    Cortes-Huerto, Kremer et.al , JCP, 2016, 145, 141103, Eq (1). Loops over frames
    so need to be saved'''
    global LammpsTraj, NB, NW, Prefix
    global Trj, StepFreq, FrameRange, NFrames, BoxL
    
    measure.LammpsTraj = LammpsTraj
    measure.StepFreq = StepFreq
    measure.__parseFrameData()
    Trj = measure.Trj ; FrameRange = measure.FrameRange ; NFrames = measure.NFrames
    BoxL = np.array(Trj.FrameData['BoxL'])
    AtomTypes = measure.AtomTypes
    
    # choose a grid of subvolume sizes greater than the correlation length
    L_corr = 12.0 #A
    NSubVols = 50
    L0 = np.linspace(L_corr, BoxL.min(), NSubVols)
    G_BB = np.zeros([NFrames, NSubVols]) ; G_WW = np.zeros([NFrames, NSubVols]) ; G_BW = np.zeros([NFrames, NSubVols])
    
    # frame iteration
    pb = sim.utility.ProgressBar('Calculating SKBI per frame...', Steps = NFrames)
    for i, frame in enumerate(FrameRange):
        Pos = Trj[frame]
        for j, this_L0 in enumerate(L0):
            G_BB[i,j], G_WW[i,j], G_BW[i,j] = benwatlib.skbi(pos = Pos, boxl = BoxL, atomtypes = AtomTypes, atomtype_b = 1, atomtype_w = 2, l0 = this_L0, randiter = 5000)
        pb.Update(i)
        
    # frame averaging
    G_BB = np.mean(G_BB, axis = 0)
    G_WW = np.mean(G_WW, axis = 0)
    G_BW = np.mean(G_BW, axis = 0)
    
    return (L0, G_BB, G_WW, G_BW)
                                
            



#### TEST ####
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    LammpsTraj = os.path.abspath(sys.argv[1])
    NB = int(sys.argv[2])
    NW = int(sys.argv[3])
    Prefix = os.path.abspath(sys.argv[4])
    
    # plot RKBI 
    methods = ['UnCorrected', 'TailCorrected', 'GeomCorrectedExact', 'GeomCorrectedExtrapolated']
    fig = plt.figure(facecolor = 'w', edgecolor = 'w')
    ax_BB = fig.add_subplot(3,1,1) ; ax_WW = fig.add_subplot(3,1,2) ; ax_BW = fig.add_subplot(3,1,3)
    for i, method in enumerate(methods):
        ret, ret_inf = make_RKBI(method)
        r_BB, G_BB = ret[0] ; G_BB_inf = ret_inf[0]
        r_WW, G_WW = ret[1] ; G_WW_inf = ret_inf[1]
        r_BW, G_BW = ret[2] ; G_BW_inf = ret_inf[2]
        ax_BB.plot(r_BB, G_BB, label = method) ; ax_BB.set_title('BB', fontsize = 15)
        ax_WW.plot(r_WW, G_WW) ; ax_WW.set_title('WW', fontsize = 15)
        ax_BW.plot(r_BW, G_BW) ; ax_BW.set_title('BW', fontsize = 15)
        ax_BB.legend(loc = 'best', prop = {'size': 10})
        for ax in [ax_BB, ax_WW, ax_BW]:
            ax.set_xlabel(r'$R(\AA)$', fontsize = 15)
            ax.set_ylabel(r'$G(\AA^3)$', fontsize = 15)
        print '\n\nMethod = %s' % method
        print '------------------------------'
        print 'G_BB_inf = %g' % G_BB_inf
        print 'G_WW_inf = %g' % G_WW_inf
        print 'G_BW_inf = %g' % G_BW_inf    
        print '------------------------------'
        
    # plot SKBI
    l, G_BB, G_WW, G_BW = make_SKBI()
    fig = plt.figure(facecolor = 'w', edgecolor = 'w')
    ax_BB = fig.add_subplot(3,1,1) ; ax_WW = fig.add_subplot(3,1,2) ; ax_BW = fig.add_subplot(3,1,3)
    ax_BB.plot(1./l, G_BB, 'k-', lw = 2, markersize = 10, marker = 'o') ; ax_BB.set_title('BB', fontsize = 15)
    ax_WW.plot(1./l, G_WW, 'k-', lw = 2, markersize = 10, marker = 'o') ; ax_WW.set_title('WW', fontsize = 15) 
    ax_BW.plot(1./l, G_BW, 'k-', lw = 2, markersize = 10, marker = 'o') ; ax_BW.set_title('BW', fontsize = 15)
    for ax in [ax_BB, ax_WW, ax_BW]:
            ax.set_xlabel(r'$1/L_0 (\AA^{-1})$', fontsize = 15)
            ax.set_ylabel(r'$G$', fontsize = 15)
    
    plt.show()    
