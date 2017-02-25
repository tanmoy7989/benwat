#!/usr/bin/env python

import numpy as np
import os, pickle, shelve
import matplotlib.pyplot as plt

import sim, pickleTraj, measure
import cgmodel as cg



##### MEASURE FUNCTIONS #####
def make_rdf_BB(Traj, Prefix):
    measure.LammpsTraj = Traj
    return measure.makeRDF(1,1, Prefix = Prefix)

def make_rdf_WW(Traj, Prefix):
    measure.LammpsTraj = Traj
    return measure.makeRDF(2,2, Prefix = Prefix)
    
def make_rdf_BW(Traj, Prefix):
    measure.LammpsTraj = Traj
    return measure.makeRDF(1,2, Prefix = Prefix)

def make_ld_BB(Traj, Prefix, LDCutBB, LDDelta):
    measure.LammpsTraj = Traj
    return measure.makeLDHist(1,1, LDCut = LDCutBB, LDDelta = LDDelta)

def make_ld_WW(Traj, Prefix, LDCutWW, LDDelta):
    measure.LammpsTraj = Traj
    return measure.makeLDHist(2,2, LDCut = LDCutWW, LDDelta = LDDelta)

def make_ld_BW(Traj, Prefix, LDCutBW, LDDelta):
    measure.LammpsTraj = Traj
    return measure.makeLDHist(1,2, LDCut = LDCutBW, LDDelta = LDDelta)
    
def make_ld_WB(Traj, Prefix, LDCutWB, LDDelta):
    measure.LammpsTraj = Traj
    return measure.makeLDHist(2,1, LDCut = LDCutWB, LDDelta = LDDelta)

def makeClusterHist(Traj, Cut, ClustAtomType, Prefix = 'clust', Normalize = True):
    global StepFreq, Trj, BoxL, AtomTypes, FrameRange, NFrames
    measure.LammpsTraj = Traj
    measure.__parseFrameData()
    
    clustpickle = Prefix + '.pickle'
    if __isComputed(clustpickle):
        return ( pickle.load(open(clustpickle, 'r')), clustpickle)
    
    inds = np.where(measure.AtomTypes == ClustAtomType)
    NAtoms = len(inds[0])
    hist = np.zeros([NAtoms, 3])
    bin_vals_block = np.zeros([NAtoms, measure.NBlocks], np.float64)
    clust_frame = np.zeros([measure.NFrames, NAtoms])
    
    pb = sim.utility.ProgressBar(Text = '', Steps = measure.NFrames)
    for frame_Ind, frame in enumerate(measure.FrameRange):
        Pos = Trj[frame][inds]
        clustdist, clustgroups = sim.geom.ClusterStats(Pos = Pos, BoxL = BoxL, Cutoff = Cut)
        if Normalize: clustdist /= np.sum(np.array(clustdist))
        clust_frame[frame_Ind, :] = np.array(clustdist)
        pb.Update(frame_Ind)

    BlockSize = int(measure.NFrames / measure.NBlocks)
    for b in range(measure.NBlocks):
        bin_vals_block[:, b] = np.mean(clust_frame[b*BlockSize:(b+1)*BlockSize, :], axis = 0)
    
    bin_centers = range(1, NAtoms+1)
    bin_vals = np.mean(bin_vals_block, axis = 1)
    bin_errs = np.std(bin_vals_block, axis = 1, ddof = 1)
   
    pickle.dump( (bin_centers, bin_vals, bin_errs), open(clustpickle, 'w'))
    return hist, clustpickle
    
def make_clust_BB(Traj, Prefix, ClustCutBB):
    return makeClusterHist(Traj = Traj, Prefix = Prefix, ClustAtomType = 1, Cut = ClustCutBB)

def make_clust_WW(Traj, Prefix, ClustCutWW):
    return makeClusterHist(Traj = Traj, Prefix = Prefix, ClustAtomType = 2, Cut = ClustCutWW)



##### GLOBALS #####
DelTempFiles = True
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
                'clust_WW' : make_clust_WW
                } 

# global settings for measure
measure.AtomNames2Types = True
measure.StepFreq = 10
measure.NBins = 50
measure.NBlocks = 1

# plot parameters
ref_styles = {'AA': 'ro-', 'SP': 'k:', 'all': 'k-'}
ref_lbls = {'AA': 'AA', 'SP': 'SP', 'all': 'all'}
transf_styles = {'AA': 'k-', 'SP': 'r-', 'all': 'b-'}
transf_lbls = {'SP': 'AA-SP', 'all': 'AA-all'}
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
                    ret, retPickle = f(Traj = trajlist[k], Prefix = mPrefix, **kwargs)
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
            ax.plot(x, y, ref_styles[i], label = ref_lbls[i])
            ax.hold(True)
        if showLegend: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
        plt.show()
        
                
    def Plot(self, Measure, RefNB = None, RefNW = None):
        global NB, NW
        mDict = shelve.open(self.Shelf)
        if RefNB is None: RefNB = 250
        if RefNW is None: RefNW = 250
        
        nrows = len(self.TrajTypes)-1
        ncols = len(NB)
        fig = plt.figure(facecolor = 'w', edgecolor = 'w')
        
        for r in range(nrows):
            for c in range(ncols):
                ind = 1 + r*len(NB) + c
                key_AA = self.genKey(Measure, 'AA', NB[c], NW[c])[0]
                key_CG = self.genKey(Measure, self.TrajTypes[r+1], NB[c], NW[c])[0]
                x_AA, y_AA, err_AA = mDict[key_AA]
                x_CG, y_CG, err_CG = mDict[key_CG]
                
                ax = fig.add_subplot(nrows, ncols, ind)
                ax.plot(x_AA, y_AA, transf_styles['AA'])
                ax.plot(x_CG, y_CG, transf_styles[self.TrajTypes[r+1]], label = transf_lbls[self.TrajTypes[r+1]])
                
                if r == 0: ax.set_title('NB = %d, NW = %d' % (NB[c], NW[c]), fontsize = 10)
                if c == 0: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
        plt.show()


## Tests
if __name__ == '__main__':
    NB, NW = [200, 300], [300, 200]
    CGTraj = {'all': ['/home/cask0/home/tsanyal/benwat/data/gromacs/NB150NW350/NB150NW350_prod_mapped.lammpstrj.gz',
                      '/home/cask0/home/tsanyal/benwat/data/gromacs/NB450NW50/NB450NW50_prod_mapped.lammpstrj.gz'],
                     
              'SP': ['/home/cask0/home/tsanyal/benwat/data/gromacs/NB300NW200/NB300NW200_prod_mapped.lammpstrj.gz',
                     '/home/cask0/home/tsanyal/benwat/data/gromacs/NB100NW400/NB100NW400_prod_mapped.lammpstrj.gz']}

    t = Transferability(CGTraj = CGTraj, Prefix = 'test', Measures = ['rdf_BB', 'rdf_BW'])
    t.Compute()
    t.Plot('rdf_BB', RefNB = 200, RefNW = 300)
    
