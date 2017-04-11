#!/usr/bin/env python

import numpy as np
import os, pickle, shelve
import matplotlib.pyplot as plt

# design
from matplotlib.ticker import MaxNLocator

import sim, pickleTraj, measure
import cgmodel as cg

kB = 0.001987
TempSet = 300.

##### MEASURE FUNCTIONS #####
def make_rdf_BB(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    return measure.makeRDF(1,1, Prefix = Prefix)

def make_rdf_WW(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    return measure.makeRDF(2,2, Prefix = Prefix)
    
def make_rdf_BW(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    return measure.makeRDF(1,2, Prefix = Prefix)

def make_ld_BB(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    LDCut = kwargs['LDCutBB']
    LDDelta = kwargs.get('LDDelta', 1.0)
    return measure.makeLDHist(1,1, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)

def make_ld_WW(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    LDCut = kwargs['LDCutWW']
    LDDelta = kwargs.get('LDDelta', 1.0)
    return measure.makeLDHist(2,2, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)

def make_ld_BW(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    LDCut = kwargs['LDCutBW']
    LDDelta = kwargs.get('LDDelta', 1.0)
    return measure.makeLDHist(1,2, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)
    
def make_ld_WB(Traj, Prefix, **kwargs):
    measure.LammpsTraj = Traj
    LDCut = kwargs['LDCutWB']
    LDDelta = kwargs.get('LDDelta', 1.0)
    return measure.makeLDHist(2,1, Prefix = Prefix, LDCut = LDCut, LDDelta = LDDelta)

def makeClusterHist(Traj, Cut, ClustAtomType, Prefix = 'clust', Normalize = True):
    global StepFreq, Trj, BoxL, AtomTypes, FrameRange, NFrames
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    
    clustpickle = Prefix + '.pickle'
    if measure.__isComputed(clustpickle):
        return ( pickle.load(open(clustpickle, 'r')), clustpickle)
    
    Trj = measure.Trj
    BoxL = measure.BoxL
    NFrames = measure.NFrames
    NBlocks = measure.NBlocks
    FrameRange = measure.FrameRange
    
    inds = np.where(measure.AtomTypes == ClustAtomType)
    NAtoms = len(inds[0])
    hist = np.zeros([NAtoms, 3])
    bin_vals_block = np.zeros([NAtoms+1, measure.NBlocks], np.float64)
    clust_frame = np.zeros([NFrames, NAtoms+1])
    
    pb = sim.utility.ProgressBar(Text = '', Steps = NFrames)
    for frame_Ind, frame in enumerate(FrameRange):
        Pos = Trj[frame][inds]
        clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = BoxL, Cutoff = Cut)
        clustdist = np.array(clustdist)
        clust_frame[frame_Ind, :] = clustdist / np.sum(clustdist)
        pb.Update(frame_Ind)

    BlockSize = int(NFrames / NBlocks)
    for b in range(NBlocks):
        bin_vals_block[:, b] = np.mean(clust_frame[b*BlockSize:(b+1)*BlockSize, :], axis = 0)
    
    bin_centers = range(NAtoms+1)
    bin_vals = np.mean(bin_vals_block, axis = 1)
    bin_errs = np.std(bin_vals_block, axis = 1, ddof = 1)
   
    pickle.dump( (bin_centers, bin_vals, bin_errs), open(clustpickle, 'w'))
    return hist, clustpickle
    
def make_clust_BB(Traj, Prefix, **kwargs):
    Cut = kwargs['ClustCutBB']
    return makeClusterHist(Traj = Traj, Prefix = Prefix, ClustAtomType = 1, Cut = Cut)

def make_clust_WW(Traj, Prefix, **kwargs):
    Cut = kwargs['ClustCutWW']
    return makeClusterHist(Traj = Traj, Prefix = Prefix, ClustAtomType = 2, Cut = Cut)


def make_KBI(Traj, Prefix = 'KBI', **kwargs):
    KBIPickle = Prefix + '.pickle'
    if measure.__isComputed(KBIPickle):
        return pickle.load(open(KBIPickle, 'r')), KBIPickle    

    # calculate densities
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    BoxL = measure.BoxL
    if isinstance(BoxL, list): BoxL = np.array(BoxL)
    BoxVol = np.prod(BoxL)
    NB = kwargs.get('NB') ; NW = kwargs.get('NW')
    rho_B = float(NB) / BoxVol ; rho_W = float(NW) / BoxVol
    
    # calculate rdfs
    hist_BB, rdfpickle_BB = make_rdf_BB(Traj, Prefix.split('KBI')[0] + 'rdf_BB')
    hist_WW, rdfpickle_WW = make_rdf_WW(Traj, Prefix.split('KBI')[0] + 'rdf_WW')
    hist_BW, rdfpickle_BW = make_rdf_BW(Traj, Prefix.split('KBI')[0] + 'rdf_BW')
    r_BB, g_BB, err_BB = hist_BB
    r_WW, g_WW, err_WW = hist_WW
    r_BW, g_BW, err_WW = hist_BW

    # calculate KB integrals
    def func_G(r, gof, L_corr):
        N_corr = len(L_corr)
        G = np.zeros(N_corr)
        dx = r[1] - r[0]
        for i in range(N_corr):
            try:
                CutInd = [list(r).index(this_r) for this_r in list(r) if this_r >= L_corr[i]][0]
            except IndexError:
                print i, L_corr[i]
                print r
                exit()
            x = r[:CutInd]
            y = (gof[:CutInd] - 1) * x**2.
            G[i] = 4 * np.pi * dx * np.sum(y)
        return G

    N = 50
    Rmax = min([np.max(r_BB), np.max(r_WW), np.max(r_BW)])
    R = np.linspace(2.0, Rmax, N) # based on Nico's paper
    G_BB = func_G(r_BB, g_BB, R)
    G_WW = func_G(r_WW, g_WW, R)
    G_BW = func_G(r_BW, g_BW, R)
    
    # calculate excess coordination numbers
    N_BB = rho_B * G_BB ; N_WW = rho_W * G_WW; N_BW = rho_W * G_BW

    # calculate solvation parameter delta_BW
    Delta_BW = G_BB + G_WW - 2 * G_BW ; Delta_BW = np.mean(Delta_BW[-10:]) # average over last 10 values
    
    # calculate dmudx
    x_B =  float(NB) / (float(NB + NW)) ; x_W = 1.0 - x_B
    dmudx  = (kB*TempSet) * (x_B * (1 + rho_B * x_W * Delta_BW)) ** (-1.0)

    # calculate dgammadx
    dgammadx = - (rho_W * x_B * Delta_BW) * (1 + rho_W * x_B * Delta_BW) ** (-1.0)
    
    # output structure
    ret = {'R': R, 'G_BB': G_BB, 'G_WW': 'G_WW', 'G_BW': G_BW, 'N_BB': N_BB, 'N_WW': N_WW, 'N_BW': N_BW,
         'Delta_BW': Delta_BW, 'dmudx': dmudx, 'dgammadx': dgammadx}
    pickle.dump(ret, open(KBIPickle, 'w'))
    return ret, KBIPickle


##### GLOBALS #####
AtomTypes = {'B': 1, 'W': 2}

# analysis parameters
TrajTypes = ['AA', 'SP', 'all']
NB, NW =  [50, 100, 150, 200, 250, 300, 350, 400, 450], [450, 400, 350, 300, 250, 200, 150, 100, 50]
MeasureFuncs = {'rdf_BB'   : make_rdf_BB,
                'rdf_WW'   : make_rdf_WW,
                'rdf_BW'   : make_rdf_BW,
                'ld_BB'    : make_ld_BB,
                'ld_WW'    : make_ld_WW, 
                'ld_BW'    : make_ld_BW,
                'ld_WB'    : make_ld_WB,
                'clust_BB' : make_clust_BB,
                'clust_WW' : make_clust_WW,
                'KBI'	   : make_KBI
                } 

# global settings for measure
measure.AtomNames2Types = True
measure.StepFreq = 10
measure.NBins = 50
measure.NBlocks = 1

# plot parameters
ref_styles = {'AA': 'ro-', 'SP': 'bx-:', 'all': 'k-'}
ref_lbls = {'AA': 'AA', 'SP': 'SP', 'all': 'CG'}
transf_styles = {'AA': 'ro-', 'SP': 'bx-', 'all': 'k-'}
transf_lbls = {'AA': 'AA', 'SP': 'AA-SP', 'all': 'CG'}
units = {'rdf_BB'   : (r'$r(\AA)$', r'$g(r)$'),
         'rdf_WW'   : (r'$r(\AA)$', r'$g(r)$'), 
         'rdf_BW'   : (r'$r(\AA)$', r'$g(r)$'), 
         'ld_BB'    : (r'$\rho_{BB}$', r'$P(\rho)$'),
         'ld_WW'    : (r'$\rho_{WW}$', r'$P(\rho)$'),
         'ld_BW'    : (r'$\rho_{BW}$', r'$P(\rho)$'), 
         'ld_WB'    : (r'$\rho_{WB}$', r'$P(\rho)$')
        }

# paths
CGTrajDir = '/home/cask0/home/tsanyal/benwat/data/modtraj/trans_test'
AATrajDir = '/home/cask0/home/tsanyal/benwat/data/gromacs'
OutputDir = os.getcwd()
RawDir = os.path.join(OutputDir, 'raw')



##### TRANSFERABILITY CLASS #####
class Transferability:
    def __init__(self, AATraj = None, CGTraj = None, Prefix = 'trans_test', Measures = None):
        global TrajTypes, MeasureFuncs, NB, NW
        global AATrajDir, CGTrajDir, OutputDir, RawDataDir
        
        self.Prefix = Prefix
        self.OutputDir = OutputDir
        self.RawDir = RawDir
        if not os.path.isdir(self.RawDir): os.mkdir(self.RawDir)
        self.Shelf = os.path.join(self.OutputDir, self.Prefix+'.shelf')
        
        if AATraj is None:
            AATraj = [os.path.join(AATrajDir, 'NB%dNW%d' % (NB[i], NW[i]), 'NB%dNW%d_prod_mapped.lammpstrj.gz' % (NB[i], NW[i]))
                      for i in range(len(NB))]
        if CGTraj is None:
            CGTraj = {}
            for i in TrajTypes[1:]:
                CGTraj[i] = [os.path.join(CGTrajDir, 'NB%dNW%d_cgmd_%s.lammpstrj.gz' % (NB[j], NW[j], i)) 
                             for j in range(len(NB))]
        self.Traj = {}
        self.Traj['AA'] = AATraj
        for trajtype, trajlist in CGTraj.iteritems(): self.Traj[trajtype] = trajlist
        self.TrajTypes = self.Traj.keys()
        self.Measures = MeasureFuncs.keys() if Measures is None else Measures
        
    def genKey(self, measuretype, trajtype, nb, nw):
        s1 = 'NB%dNW%d_%s_%s' % (nb, nw, trajtype, measuretype)
        s2 = os.path.join(self.RawDir, s1)
        return s1, s2
        
    def Compute(self, **kwargs):
        global NB, NW, MeasureFuncs
        mDict = shelve.open(self.Shelf)
        for i in self.Measures:
            print '\n\nProperty: ', i
            for j in self.TrajTypes:
                print ' Computing properties for traj type: ', j
                trajlist = self.Traj[j]
                for k in range(len(NB)):
                    print '  NB = %d, NW = %d\n' % (NB[k], NW[k])
                    mKey, mPrefix = self.genKey(i, j, NB[k], NW[k])
                    f = MeasureFuncs[i]
                    if i == 'KBI': ret, retPickle = f(Traj = trajlist[k], Prefix = mPrefix, NB = NB[k], NW = NW[k])
                    else: ret, retPickle = f(Traj = trajlist[k], Prefix = mPrefix, **kwargs)
                    mDict[mKey] = ret
        mDict.close()
    
            
    def PlotRef(self, Measure, RefNB = None, RefNW = None, showLegend = True):
        mDict = shelve.open(self.Shelf)
        if RefNB is None: RefNB = 250
        if RefNW is None: RefNW = 250
        
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        ax.set_title('NB%dNW%d' % (RefNB, RefNW))
        ax.set_xlabel(units[Measure][0], fontsize = 15)
        ax.set_ylabel(units[Measure][1], fontsize = 15)
        for i in self.TrajTypes:
            key = self.genKey(Measure, i, RefNB, RefNW)[0]
            x, y, err = mDict[key]
            ax.plot(x, y, ref_styles[i], lw = 3, markersize = 6, label = ref_lbls[i])
            ax.hold(True)
        if showLegend: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
        
                
    def PlotHist(self, Measure, RefNB = None, RefNW = None):
        global NB, NW
        mDict = shelve.open(self.Shelf)
        if RefNB is None: RefNB = 250
        if RefNW is None: RefNW = 250
        
        nrows = len(self.TrajTypes) - 1
        ncols = len(NB)
        fig = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (10, 3))
        
	axs = []
        for r in range(nrows):
            for c in range(ncols):
                ind = 1 + r*len(NB) + c
                key_AA = self.genKey(Measure, 'AA', NB[c], NW[c])[0]
                key_CG = self.genKey(Measure, self.TrajTypes[r+1], NB[c], NW[c])[0]
		pickleAA = os.path.join(self.RawDir, key_AA + '.pickle')
		pickleCG = os.path.join(self.RawDir, key_CG + '.pickle')                
		x_AA, y_AA, err_AA = pickle.load(open(pickleAA, 'r')) #mDict[key_AA]
                x_CG, y_CG, err_CG = pickle.load(open(pickleCG, 'r')) #mDict[key_CG]
                
                ax = fig.add_subplot(nrows, ncols, ind)
                ax.plot(x_AA, y_AA, transf_styles['AA'], lw = 0, markersize = 6, label = ref_lbls['AA'])
                ax.plot(x_CG, y_CG, transf_styles[self.TrajTypes[r+1]], lw = 2, label = ref_lbls[self.TrajTypes[r+1]])
		axs.append(ax)
                
                if r == 0: 
                    if NB[c] == RefNB: s = r'$x_{B} = $' + '%g (ref)' % ( float(NB[c]) / (NB[c] + NW[c]) )
                    else: s = s = r'$x_{B} = $' + '%g' % ( float(NB[c]) / (NB[c] + NW[c]) )
                    ax.set_title(s, fontsize = 20)
                #if c == 0: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
	
	# design
	for ax in axs:
		ax.xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
		ax.yaxis.set_major_locator(MaxNLocator(nbins = 5, prine = 'both'))
	

##### Tests #####
def dummy_test():
    global NB, NW
    NB, NW = [200, 300], [300, 200]
    CGTraj = {'all': ['/home/cask0/home/tsanyal/benwat/data/gromacs/NB150NW350/NB150NW350_prod_mapped.lammpstrj.gz',
                      '/home/cask0/home/tsanyal/benwat/data/gromacs/NB450NW50/NB450NW50_prod_mapped.lammpstrj.gz'],
                     
              'SP': ['/home/cask0/home/tsanyal/benwat/data/gromacs/NB300NW200/NB300NW200_prod_mapped.lammpstrj.gz',
                     '/home/cask0/home/tsanyal/benwat/data/gromacs/NB100NW400/NB100NW400_prod_mapped.lammpstrj.gz']}

    t = Transferability(CGTraj = CGTraj, Prefix = 'test', Measures = ['ld_WW', 'ld_BB'])
    t.Compute(LDCutBB = 7.5, LDCutWW = 3.5)
    t.Plot('ld_WW', RefNB = 200, RefNW = 300)
    plt.show()


def NB250test():
    global NB, NW
    NB, NW = [250], [250]
    CGTraj = {'all': ['/home/cask0/home/tsanyal/benwat/data/modtraj/NB250NW250_cgmd_all.lammpstrj.gz']} #add SP traj
    t = Transferability(CGTraj = CGTraj, Prefix = 'NB250test')
    t.Compute
    for m in t.Measures: t.Plot(m, RefNB = 250, RefNW = 250)
    plt.show()


if __name__ == '__main__': dummy_test()
                
