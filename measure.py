#!/usr/bin/env python

import os, sys, pickle
import numpy as np

# dependencies
import sim
import pickleTraj
import parse_potential
import cgmodel as cgsys.path.append('~/'); from selLDCut import *

def getBoxL(LammpsTraj):
    trj = pickleTraj(LammpsTraj)
    return trj.FrameData['BoxL']

def gen_Srel_rc(LammpsTraj, LDCuts = []):
    # prepare the system
    cg.LammpsTraj = LammpsTraj
    Prefix = cg.makeSysPrefix()
    cg.LDCutBW = LDCuts[0]
    cg.BoxL = getBoxL(LammpsTraj)
    Sys = cg.makeSys()
    
    # initialize the Bennett and FEP based Delta_Srel calculator
    SrelBennett.Sys = Sys
    SrelBennett.LammpsTraj = LammpsTraj
    SrelBennett.Prefix = Prefix
    SrelBennett.LDCuts = LDCuts
    SrelBennett.LDDelta = 1.2
    SrelBennett.genFileNames()
    
    # setup cg model and CG-MD simulations for different cutoffs
    SrelBennett.runCGMD()
    
    # start Bennett method calculations
    SrelBennett.Srel_rc()
    

def gen_fsw(LammpsTraj, LDCuts = [], NB = 250, NW = 250):
    # initialize the fsw calculator
    fsw.LammpsTraj = LammpsTraj
    fsw.NCentAtoms = NB ; fsw.NNeighAtoms = NW
    fsw.CentAtomType = 1 ; fsw.NeighAtomType = 1
    fsw.LDCuts = LDCuts ; fsw.LDDelta = 1.2
    fsw.Prefix = 'NB%dNW%d' % (NB, NW) ; fsw.genFileNames()

    # calculate rdf
    fsw.makeRDF()
    
    # calculate maximal correlation
    #fsw.makeFSWCorrelation(rdfCut = ?)
    #fsw.calcCorrelation()
    
