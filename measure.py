#!/usr/bin/env python

import os, sys, pickle
import numpy as np

# dependencies
import sim
import pickleTraj
import parse_potential
import cgmodel as cg
import measurelib


#### FREE ENERGY ESTIMATION ROUTINES ####

#1) Free Energy Pertturbation

def FEP(BE11, BE22, BE12, BE21, Dir = 1, Verbose = True):
    '''
    calculates forward or avg of forward and reverse 
    free energy differences using simple free energy
    perturbation
    '''
    
    DeltaBE = BE11 - BE12
    FE1 = -np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) - DeltaBE.max()
    DeltaBE = BE22 - BE21
    FE2 = np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) + DeltaBE.max()

    if Dir == 1:
        FE = FE1 
    elif Dir == 2:
        FE = 0.5*(FE1 + FE2)
        
    if Verbose:
        print "\t---> %dD Free Energy Perturbation gives: %f" % (Dir, FE)
        
    return FE



#2) Bennett's Method

def BennettFE(BE11, BE21, BE22, BE12, Verbose = False, 
              Tol = 1.e-10, MaxIter = 10000):
    '''
    Returns -log(Q2/Q1) = beta2 F2 - beta1 F1 using the Bennett method.  
    BEij is the list of beta*E values for system i evaluated in system j.
    '''
    
    fd = lambda x,w : 1./(np.exp(w) + np.exp(x))
    
    
    #check array sizes
    N1 = len(BE11)
    if len(BE12) != len(BE11):
        raise ValueError("BE11 and BE12 arrays are not the same length.")
    N2 = len(BE21)
    if len(BE22) != len(BE21):
        raise ValueError("BE21 and BE22 arrays are not the same length.")        
    if not (np.all(np.isfinite(BE11)) and np.all(np.isfinite(BE12)) 
            and np.all(np.isfinite(BE21)) and np.all(np.isfinite(BE22))):
        if DEBUG: 
            print "BE11 BE12 BE22 BE21"
            for i in xrange(n):
                print i, BE11[i], BE12[i], BE22[i], BE21[i]
        raise ValueError("Found non-finite value in BE11, BE12, BE22, or BE21")
    #fe perturbation for initial guesses
    DeltaBE = BE11 - BE12
    FE1 = -np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) - DeltaBE.max()
    DeltaBE = BE22 - BE21
    FE2 = np.log(np.mean(np.exp(DeltaBE - DeltaBE.max()))) + DeltaBE.max()
    if Verbose:
        print "Bennett initial guesses: %f, %f" % (FE1, FE2)
    M = np.log(float(N2)/float(N1))
    #setup first points: x is beta2 F2 - beta1 F1, y should be zero
    xo = FE1
    NVals = BE21 - BE22 + xo + M
    DVals = BE12 - BE11 - xo - M
    maxterm = -max(NVals.max(), DVals.max())
    Num = np.sum(fd(NVals+maxterm, maxterm))
    Dem = np.sum(fd(DVals+maxterm, maxterm))
    yo = np.log(Num / Dem)
    x = FE2
    y = 1.
    #loop until tolerance is met
    Iter = 0
    while np.abs(y) > Tol and Iter < MaxIter:
        Iter += 1
        #evaluate y for current x
        NVals = BE21 - BE22 + x + M
        DVals = BE12 - BE11 - x - M
        maxterm = -max(NVals.max(), DVals.max())
        Num = np.sum(fd(NVals+maxterm, maxterm))
        Dem = np.sum(fd(DVals+maxterm, maxterm))
        y = np.log(Num / Dem)
        #predict new x
        xn = (y*xo-yo*x)/(y-yo)
        xo = x
        x = xn
        yo = y
        #print messages
        if Verbose:
            print "Bennett iteration %d: current error is %.3e" % (Iter, np.abs(y))
    #now compute the estimated error
    FE = xo
    FEerr = BennettErr(BE11, BE21, BE22, BE12, FE)
    if Verbose:
        print "Bennett final free energy: %f +- %f" % (FE, FEerr)
    return FE, FEerr
    


#3) Error estimate in Bennett's method

def BennettErr(BE11, BE21, BE22, BE12, FE):
    '''
    Computes the error in the Bennett calculation, sqrt(var(ln(q2/q1))).
    Based on Shirts et al, PRL 91, 140601 (2003)
    '''
    #compute all "work" measurements
    W = np.concatenate((BE12 - BE11, BE22 - BE21))
    #add free energy (ln (Q2/Q1))
    W = W + FE
    n = len(W)
    c = np.max(np.abs(W))
    terms = 1. / (2.*np.exp(-c) + np.exp(W-c) + np.exp(-W-c))
    err = (np.exp(c) / np.mean(terms) - 4.) / n
    err = np.sqrt(err)
    return err    


#### RELATIVE ENTROPY VS. CUTOFF SPACE GENERATOR ####

def Srel_rc(NB, NW, LammpsTraj, data_dir = os.getcwd()):
    # flags
    DEBUG = False
    RefreshLogs = True

    def EvalSysEne(Sys, sampleTrj, TempSet = 300.0, Iter = (0,-1,1)):
        '''
        sampleTrj refer to the state (AA or CG) from which the position
        co-ordintates are being sampled. Sys is a system with the forcefield
        of the state at which energies are to be evaluated
        returns beta * Penergy of system
        '''
        start = Iter[0]; stop = Iter[1]; freq = Iter[2]
        if stop == -1: stop = len(sampleTrj)  
        FrameRange = range(start, stop, freq)
        Nframes = len(FrameRange)
        Ene = np.zeros([Nframes], np.float64)
        beta = 1./(TempSet * Sys.Units.kB)
    
        # for each frame, minimage and calculate energy
        Ene_Ind = 0
        pb = sim.utility.ProgressBar(Text = 'Processing frames...',  Steps = Nframes)
        for frame in FrameRange:
            Pos = sampleTrj[frame]
            Sys.Arrays.Pos = Pos
            Sys.ForceField.Eval()
            Ene[Ene_Ind] = Sys.PEnergy
            if DEBUG:
                print frame, Ene_Ind, Ene[Ene_Ind]
                raw_input()    
            pb.Update(Ene_Ind)
            Ene_Ind += 1
        if DEBUG: 
            print '(sample, eval)', Sys.sample_trajtype, Sys.eval_trajtype
            print len(Ene), len(sampleTrj)
            raw_input()
        
        return beta * Ene
    
    ########### MAIN #################        
    Prefix = 'NB%dNW%d' % (NB, NW)
    
    TempSet = 300.0
    Delta = 1.2
    base_fmt = {'cg_ff': '%s_%d_SPLD_sum.txt' % Prefix, 'cg_trj': '%s_%d_MD.lammpstrj.gz' % Prefix}

    # write output headers
    LDCuts = np.arange(4.5, 10.0, 0.1)
    delta_Srel = 0.
    datalog = {'Bennett': '%s_Srel_Bennett.dat' % datalog_prefix, 'Fep': '%s_Srel_Fep.dat' % datalog_prefix }
    for key in datalog.keys():
        this_log = datalog[key]
        if not os.path.isfile(this_log) or RefreshLogs:
            of = open(this_log, 'w')
            of.write('#LDCut \t Srel \t Err\n')
            of.write('%g \t %f \t %f\n' % (LDCuts[0], delta_Srel, 0.0))
            of.close()

    # extract all necessary AA data
    Trj_AA = pickleTraj(LammpsTraj)
    BoxL = Trj_AA.FrameData['BoxL']
    
    # loop over different LDCuts
    for i, LDCut in enumerate(LDCuts[0:-1]):
        print '\n\nLDCUTs = (%g, %g)\n' % (LDCuts[i], LDCuts[i+1])
        sumfile1 = base_fmt['cg_ff'] % (fftype, i); sumfile2 = base_fmt['cg_ff'] % (fftype, (i+1))
        trajfile1 = base_fmt['cg_trj'] % (fftype, i) ; trajfile2 = base_fmt['cg_trj'] % (fftype, (i+1))
        ParamString1 = parse_potential.parseParamString(sumfile1); ParamString2 = parse_potential.parseParamString(sumfile2)
        Trj_CG_1 = pickleTraj(trajfile1); Trj_CG_2 = pickleTraj(trajfile2)
    
    # create systems with forcefields of LDCuts i and i+1
    cg.NB = NB
    cg.NW = NW
    cg.BoxL = BoxL
    cg.LDCutBW = LDCuts[i] ; cg.ParamString = ParamString1 ; Sys1 = cg.makeSys()
    cg.LDCutBW = LDCuts[i+1] ; cg.ParamString = ParamString2 ; Sys2 = cg.makeSys()
                       
    # calculate Energy difference in AA ensemble
    print '--> Calculating energies in AA ensemble'
    Sys1.__setattr__('sample_trajtype', 'AA') ; Sys2.__setattr__('sample_trajtype', 'AA')
    Sys1.sample_trajtype = 'AA'; Sys2.sample_trajtype = 'AA'
    BE1_AA = EvalSysEne(Sys = Sys1, sampleTrj = Trj_AA)
    BE2_AA = EvalSysEne(Sys = Sys2, sampleTrj = Trj_AA)
    
    # reprocess trajectories to calculate free energy differences
    print '--> Reprocessing trajectories'
    Sys1.sample_trajtype = 'CG'; Sys2.sample_trajtype = 'CG'
    
    BE11 = EvalSysEne(Sys = Sys1, sampleTrj = Trj_CG_1)
    BE21 = EvalSysEne(Sys = Sys1, sampleTrj = Trj_CG_2)
    BE22 = EvalSysEne(Sys = Sys2, sampleTrj = Trj_CG_2)
    BE12 = EvalSysEne(Sys = Sys2, sampleTrj = Trj_CG_1)
    Nframes1 = len(BE11); Nframes2 = len(BE22)
    
    # running Bennett's algorithm
    print "--> Running Bennett's method..."
    FE1, Err1 = BennettFE(BE11 = BE11, BE22 = BE22, BE12 = BE12, BE21 = BE21)
    delta_Srel1 = (np.mean(BE2_AA) - np.mean(BE1_AA)) - FE1
    
    # running Free Energy perturbation
    print "--> Running 1D Free Energy Perturbation..."
    FE2 = FEP(BE11 = BE11, BE22 = BE22, BE12 = BE12, BE21 = BE21)
    delta_Srel2 = (np.mean(BE2_AA) - np.mean(BE1_AA)) - FE2
    Err2 = 0.0
    
    # logging data
    for key in datalog.keys():
        if key == 'Bennett':
            delta_Srel = delta_Srel1
            Err = Err1
        elif key == 'Fep':
            delta_Srel = delta_Srel2
            Err = Err2
        
        this_log = datalog[key]
        of = open(this_log, 'a')
        of.write('%g \t %f \t %f\n' %(LDCuts[i+1], delta_Srel, Err))
        of.close()
        
        
        
#### CORRELATION CALCULATOR WITH HYDRATION SHELL WATERS ####

#1) radial distribution function

def rdf(NB, NW, LammpsTraj, Prefix = None):
    Trj = pickleTraj(LammpsTraj, Verbose = True)
    BoxL = Trj.FrameData['BoxL'][0]
	Cut = 0.5 * BoxL
	if BoxL == 0: BoxL = float(raw_input('BoxL = 0, found. Enter nonperiodic boxlength: '))
	
	FrameRange = range(0, len(Trj), 1)
	NFrames = len(FrameRange)
	Bin_min = 1.00
	Bin_max = Cut
	Bin_delta = (Bin_max - Bin_min)/float(Nbins)
	Bin_centers = np.zeros([Nbins], np.float64)
	for i in range(Nbins): Bin_centers[i] = Bin_min + (i+0.5)*Bin_delta	
	g = np.zeros([Nbins, 3], np.float64) #1-BB #2-WW #3-BW
	 
	# frame stepping
	pb = sim.utility.ProgressBar(Text = 'Processing frame by frame...', Steps = NFrames)
	for frame in FrameRange:
		Pos = Trj[frame]
		g = measurelib.rdf(g = g, bin_centers = Bin_centers, bin_delta = Bin_delta, 
					       pos = Pos, atomtypes = Trj.AtomTypes, boxl = BoxL)										
		pb.Update(int(frame/stepfreq))
		
	# normalize the rdf
	rdf = {'bin_centers': Bin_centers, 'BB': '', 'WW': '', 'BW': ''}
	print 'Normalizing g(r)...'
	BoxVol = BoxL**3.
	g /= NFrames
	gofr = g * BoxVol
	for j in range(Nbins):
		r = Bin_centers[j] - 0.5*Bin_delta
		next_r = Bin_centers[j] + 0.5*Bin_delta
		gofr[j,:] /= (4.*np.pi/3.) * (next_r**3. - r**3.)

	rdf['BB'] = gofr[:,0]/(NB * (NB - 1)/2.)
	rdf['WW'] = gofr[:,1]/(NW * (NW - 1)/2.)
	rdf['BW'] = gofr[:,2]/(NB * NW)
		
	# dumping data
	if not Prefix: Prefix = 'NB%dNW%d_rdf' % (NB, NW)
	pickleName = os.path.join(Prefix + '.pickle')
	pickle.dump(rdf, open(pickleName, 'w'))	

    

#2) first shell waters

def fsw():
    pass

                                                                
