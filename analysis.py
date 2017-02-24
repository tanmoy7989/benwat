#!/usr/bin/env python

import numpy as np
import os, pickle, shelve
import matplotlib.pyplot as plt

import pickleTraj
import cgmodel as cg
import measure

DelTempFiles = True
AtomTypes = {'B': 1, 'W': 2}

# analysis parameters
TrajTypes = ['AA', 'SP', 'all']
MeasureTypes = ['rdf', 'ld']
NB, NW =  [50, 100, 150, 200, 250, 300, 350, 400, 450], [450, 400, 350, 300, 250, 200, 150, 100, 50]
MeasureFuncs = {'rdf_BB'   :,
                'rdf_WW'   :,
                'rdf_BW'   :,
                'ld_BB'    :,
                'ld_WW'    :, 
                'ld_BW'    :,
                'ld_WB'    :,
                } 

# global settings for measure
measure.AtomNames2Types = True
measure.StepFreq = 10
measure.NBins = 50
measure.NBlocks = 1

# plot parameters
ref_styles = {'AA': 'ro', 'SP': 'k:', 'all': 'k-'}
lbls = {'AA': 'AA', 'SP': 'SP', 'all': 'all'}
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

class Transferability:
    def __init__(self, AATraj = None, CGTraj = None, Prefix = 'trans_test', Measures = None, **kwargs):
        global TrajTypes, MeasureTypes, NB, NW
        global AATrajDir, CGTrajDir, OutputDir, RawDataDir
        
        self.Prefix = Prefix
        self.OutputDir = OutputDir
        self.RawDir = RawDir
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
        self.Measures = MeasureTypes if Measures is None else Measures
        
    def genKey(measuretype, trajtype, nb, nw):
        s1 = 'NB%dNW%d_%s_%s' % (nb, nw, trajtype, measuretype)
        s2 = os.path.join(self.RawDir, s1)
        return s1, s2
        
    def Compute(self, **kwargs):
        global NB, NW, MeasureFuncs
        mDict = shelve.open(self.Shelf)
        for i in self.Measures:
            print 'Measure: ', i
            for j in self.TrajTypes:
                print ' Computing properties for traj type: ', j
                trajlist = self.Traj[j]
                for k in range(len(NB)):
                    print '  NB = %d, NW = %d\n' % (NB, NW)
                    mKey, mPrefix = genKey(i, j, NB[k], NW[k])
                    f = MeasureFunc[i]
                    ret, retPickle = f(Traj = trajlist[k], Prefix = mPrefix, **kwargs)
                    mDict[mKey] = ret
        mDict.close()
    
            
    def PlotRef(self, Measure, RefNB = None, RefNW = None, showLegend = True):
        mDict = shelve.open(self.Shelf)
        if not Measure in mDict.keys():
            raise IOError('Measure %s has not yet been computed' % Measure)
        if RefNB is None: RefNB = 250
        if RefNW is None: RefNW = 250
        
        fig = plt.figure(figsize = (5,5), facecolor = 'w', edgecolor = 'w')
        ax = fig.add_subplot(1,1,1)
        ax.set_title('NB%dNW%d' % (RefNB, RefNW))
        ax.set_xlabel(units[Measure][0], fontsize = 15)
        ax.set_ylabel(units[Measure][1], fontsize = 15)
        for i in self.TrajTypes():
            key = genKey(Measure, i, RefNB, RefNW)[0]
            x, y, err = mDict[key]
            ax.plot(x, y, refstyles[i], label = reflbls[i])
            ax.hold(True)
        if showLegend: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
        plt.show()
        
                
    def Plot(self, Measure, RefNB = None, RefNW = None):
        global NB, NW
        mDict = shelve.open(self.Shelf)
        if not Measure in mDict.keys():
            raise IOError('Measure %s has not yet been computed' % Measure)
        if RefNB is None: RefNB = 250
        if RefNW is None: RefNW = 250
        
        nrows = len(self.TrajTypes[1:])
        ncols = len(self.Traj['AA'])
        fig = plt.figure(figsize = (5*nrows, 5*ncols), facecolor = 'w', edgecolor = 'w')
        
        for r in range(self.TrajTypes[1:]):
            for c in range(len(NB)):
                ind = 1 + r*len(NB) + c
                key_AA = genKey(Measure, 'AA', NB[c], NW[c])[0]
                key_CG = genKey(Measure, self.TrajTypes[r], NB[c], NW[c])[0]
                x_AA, y_AA, err_AA = mDict[key_AA]
                x_CG, y_CG, err_CG = mDict[key_CG]
                
                ax = fig.add_subplot(nrow, ncols, ind)
                ax.plot(x_AA, y_AA, transf_styles['AA'])
                ax.plot(x_CG, y_CG, transf_styles[self.TrajType[r]], label = trans_lbls[self.TrajType[r]])
                
                if r == 0: ax.set_title('NB = %d, NW = %d' % (NB[c], NW[c]), fontsize = 10)
                if c == 0: ax.legend(loc = 'best', prop = {'size': 15})
        mDict.close()
        plt.show()


##### MEASURE FUNCTIONS #####
