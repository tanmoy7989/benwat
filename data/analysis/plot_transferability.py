#!/usr/bin/env python
import os, numpy as np, cPickle as pickle
import matplotlib, matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

lbldict = {'rdf': {'xlabel': r'$r$', 'ylabel': r'$g(r)$'},
	  'ld':  {'xlabel': r'$\rho_{LD}$', 'ylabel': 'distribution'},
	   'delta_BW': {'xlabel' :r'$x_B$', 'ylabel': r'$\Delta_{BW}$'},
	   'dmudx': {'xlabel': r'$x_B$', 'ylabel': r'$\frac{d\mu}{dx_B}$'},
	   'dgammadx': {'xlabel': r'$x_B$', 'ylabel': r'$\frac{d ln \gamma}{dx_B}$'}
	  }
for k in ['rdf_BB', 'rdf_WW', 'rdf_BW']: lbldict[k] = lbldict['rdf']
for k in ['ld_BB', 'ld_BW', 'ld_WB', 'ld_WW']: lbldict[k] = lbldict['ld']
AADataDir = os.path.abspath('./')
CGDataDir = os.path.abspath('./ref_NB50')

CGTrajTypes = ['SP', 'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 'SPLD_WB', 'SPLD_BB_WW', 'SPLD_BB_WW_BW', 'SPLD_all']
#CGTrajTypes = ['SP', 'SPLD_BB', 'SPLD_WW', 'SPLD_BW', 'SPLD_WB']  # set 1
#CGTrajTypes = ['SP', 'SPLD_BB_WW', 'SPLD_BB_WW_BW', 'SPLD_all']    # set 2


clrdict = {} ; clrs = [name for name, hex in matplotlib.colors.cnames.items()]
for i, trajtype in enumerate(CGTrajTypes):
	clrdict[trajtype] = clrs[i]
#{'AA': 'black', 'SP': 'red', 'SPLD_BB': 'green', 'SPLD_BW': 'blue', 'SPLD_WB': 'magenta'}
NBlist = [50, 100, 150, 200, 250, 300, 350, 400, 450]
NWlist = [450, 400, 350, 300, 250, 200, 150, 100, 50]
#NBlist = [10] ; NWlist = [4990]

xlim = None
ylim = None

def plot_1Dhist(ax, x, y, OrderParam, trajtype, NB, NW, flagdict):
	global xlim, ylim	
	if trajtype == 'AA':
		ax.plot(x,y, linestyle = 'solid', lw = 1, marker = 'o', markersize = 4, color = 'black', label = 'AA')
	else:	
		ax.plot(x, y, linestyle = 'solid', lw = 1, color = clrdict[trajtype], label = trajtype)
	
	if not xlim is None: ax.set_xlim(xlim)
	if not ylim is None: ax.set_ylim(ylim)	
	if flagdict['has_legend']: ax.legend(loc = 'best', prop = {'size': 10})
	if flagdict['has_xlabel']: ax.set_xlabel(lbldict[OrderParam]['xlabel'], fontsize = 15)
	if flagdict['has_ylabel']: ax.set_ylabel(lbldict[OrderParam]['ylabel'], fontsize = 15)
	ax.xaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	ax.yaxis.set_major_locator(MaxNLocator(nbins = 5, prune = 'both'))
	if flagdict['has_title']:	
		x_B = float(NB) / (NB+NW)
		title = r'$x_B = $' + ' %g' % x_B
		if flagdict['isRef']: title += ' (ref)'
		ax.set_title(title)


def plotpanel_1Dhist(OrderParam):	
	nrows = len(CGTrajTypes)
	ncols = len(NBlist)
	KB_div = np.zeros([nrows, ncols])
	x_B = np.zeros(ncols)
	fig1 = plt.figure(figsize = (20, 10), facecolor = 'w', edgecolor = 'w')
	
	axdict = {}
	for i, trajtype in enumerate(CGTrajTypes):
		for j in range(len(NBlist)):
			# get data
			AApickle = os.path.join(AADataDir, 'AA', 'NB%dNW%d_AA_%s.pickle' % (NBlist[j], NWlist[j], OrderParam))
			CGpickle = os.path.join(CGDataDir, trajtype, 'NB%dNW%d_%s_%s.pickle' % (NBlist[j], NWlist[j], trajtype, OrderParam))
			x_AA, y_AA, err_AA = pickle.load(open(AApickle, 'r'))
			x_CG, y_CG, err_CG = pickle.load(open(CGpickle, 'r'))
			
			# calculate kullback-liebler divergence
			x_B[j] = float(NBlist[j] ) / (NBlist[j] + NWlist[j])
			dx_AA = x_AA[1] - x_AA[0]
			logterm = y_AA / y_CG
			logterm = logterm[~np.isnan(logterm)] ; p_AA = y_AA[~np.isnan(logterm)]
			logterm = logterm[~np.isinf(logterm)] ; p_AA = p_AA[~np.isinf(logterm)]
			KB_div[i,j] = dx_AA * np.sum(p_AA * logterm)
			
			# set up axes properties
			ind = i*ncols + j + 1
			flagdict = {'has_xlabel': (i == nrows-1 and j == ncols/2),
				    'has_ylabel': (j == 0 and i == nrows/2),
				    'has_legend': (j == 0),
				    'has_title': (i == 0),
				    'isRef': NBlist[j] == 250
				   }
			ax = fig1.add_subplot(nrows, ncols, ind)
			axdict[(i,j)] = ax
			
			# plot it			
			plot_1Dhist(ax, x_AA, y_AA, OrderParam, 'AA', NBlist[j], NWlist[j], flagdict)
			plot_1Dhist(ax, x_CG, y_CG, OrderParam, trajtype, NBlist[j], NWlist[j], flagdict)
	
	# plot kullback liebler divergence
	fig2 = plt.figure(facecolor = 'w', edgecolor = 'w', figsize = (5,5))
	ax = fig2.add_subplot(1,1,1)
	for i, trajtype in enumerate(CGTrajTypes):
		ax.plot(x_B, KB_div[i], lw = 1, markersize = 6, marker = 'o', color = clrdict[trajtype], label = trajtype)
	ax.set_xlabel(r'$x_B$', fontsize = 15)
	ax.set_ylabel(r'$D_{KL}$', fontsize = 15)
	ax.legend(loc = 'best', prop = {'size': 15})	
	
	# panel space settings	
	fig1.tight_layout()
	fig2.tight_layout()


def plot_KBI():
	dmudx = np.zeros([len(CGTrajTypes)+1, len(NBlist)])
	dgammadx = np.zeros([len(CGTrajTypes)+1, len(NBlist)])
	delta_BW = np.zeros([len(CGTrajTypes)+1, len(NBlist)])
	x_B_list = np.zeros(len(NBlist))
	
	fig = plt.figure(figsize = (15, 5), facecolor = 'w', edgecolor = 'w')
	ax1 = fig.add_subplot(1,3,1) ; ax2 = fig.add_subplot(1,3,2); ax3 = fig.add_subplot(1,3,3)	
	for i, trajtype in enumerate(CGTrajTypes):
		for j in range(len(NBlist)):
			# get data
			AApickle = os.path.join(DataDir, 'AA', 'NB%dNW%d_AA_KBI_scalar.pickle' % (NBlist[j], NWlist[j]))
			CGpickle = os.path.join(DataDir, trajtype, 'NB%dNW%d_%s_KBI_scalar.pickle' % (NBlist[j], NWlist[j], trajtype))
			dAA = pickle.load(open(AApickle, 'r'))
			dCG = pickle.load(open(CGpickle, 'r'))
			
			if i == 0:
				dmudx_AA = dAA['dmudx']
				dgammadx_AA = dAA['dgammadx']
				delta_BW_AA = dAA['Delta_BW'] ; 
				x_B = float(NBlist[j]) / (NBlist[j] + NWlist[j])
				dmudx[0, j] = dmudx_AA 
				x_B_list[j] = x_B
				dgammadx[0, j] = dgammadx_AA
				delta_BW[0, j] = delta_BW_AA
				

			dmudx_CG = dCG['dmudx']
			dgammadx_CG = dCG['dgammadx']
			delta_BW_CG = dCG['Delta_BW']
			
			# collect into arrays
			dmudx[i+1, j] = dmudx_CG
			dgammadx[i+1, j] = dgammadx_CG
			delta_BW[i+1, j] = delta_BW_CG

		# plot it
		if i == 0: ax1.plot(x_B_list, dmudx[0,:], lw = 1, marker = 'o', markersize = 10, color = 'black', label = 'AA')
		ax1.plot(x_B_list, dmudx[i+1,:], lw = 1, marker = 'o', markersize = 10, color = clrdict[trajtype], label = trajtype)
		
		if i == 0: ax2.plot(x_B_list, dgammadx[0,:], lw = 1, marker = 'o', markersize = 10, color = 'black', label = 'AA')
		ax2.plot(x_B_list, dgammadx[i+1,:], lw = 1, marker = 'o', markersize = 10, color = clrdict[trajtype], label = trajtype)

		if i == 0: ax3.plot(x_B_list, delta_BW[0,:], lw = 1, marker = 'o', markersize = 10, color = 'black', label = 'AA')
		ax3.plot(x_B_list, delta_BW[i+1,:], lw = 1, marker = 'o', markersize = 10, color = clrdict[trajtype], label = trajtype)
	
	ax1.set_xlabel(lbldict['dmudx']['xlabel'], fontsize = 15) ; ax1.set_ylabel(lbldict['dmudx']['ylabel'], fontsize = 15)
	ax2.set_xlabel(lbldict['dgammadx']['xlabel'], fontsize = 15) ; ax2.set_ylabel(lbldict['dgammadx']['ylabel'], fontsize = 15)
	ax3.set_xlabel(lbldict['delta_BW']['xlabel'], fontsize = 15) ; ax3.set_ylabel(lbldict['delta_BW']['ylabel'], fontsize = 15)
	ax1.legend(loc = 'best', prop = {'size': 8})
	
	# panel space settings
	fig.tight_layout()
			
			
#### MAIN
#xlim = [2, 10]
#ylim = None
plotpanel_1Dhist('ld_BW')
plt.show()
