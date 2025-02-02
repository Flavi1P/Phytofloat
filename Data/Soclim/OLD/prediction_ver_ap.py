# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from datetime import datetime, timedelta
import matplotlib.dates as mdates
from F490 import F490
from netCDF4 import Dataset
import pandas as pd
from plot_tools import *

plt.interactive(True)
plt.close('all')

#Define a moving average function
def moving(x,window):
    res = np.empty(x.shape[0]-window)
    for i in np.arange(x.shape[0]-window):
        res[i] = np.nanmean(x[i:i+window])
    return np.hstack([np.repeat(np.nan,window/2),res,np.repeat(np.nan,window/2)])

#define nonan
def nonan(data):
	return data[~np.isnan(data)]

#Import the Garcia and Gordon oxygen solubility
from O2sol import O2sol

def scale(data):  
	return data   # here the scaling function is the identity function -> no modification of the prediction!
	# might be usefull for furture application  : retrieve absolute abundance from relative
	
#==============================================================================
# Perform PLS
#==============================================================================
execfile('PLS.py')

#==============================================================================
# Declare the matrices
#==============================================================================
floats = ['049b','036b','037c','107c','104c']

plankton = ['Bact','Pico','Nano','Micro']

# Structures to store float data
time = []
lon = []
lat = []
temp = []
sal = []
ox = []
ox_sat = []
aou = []
sigma = []
chl = []
bbp = []
cp = []
par = []
ed = []
bbp_spikes = []
chl_spikes = []
chl_bbp = []
chl_cp = []
bbp_cp = []

bact = []
pico = []
nano = []
diat = []
bact_norm = []
pico_norm = []
nano_norm = []
diat_norm = []

ctot = []
cphyto = []

chl_mld = []
temp_mld = []
mld = []
chl_int  = []
bbp_mld = []
par_mld = []

F_all = []

threshold = 0.03 # Density crierion (De Boyer Montégut et al., 2004)
depth = np.arange(1001)
dep_sel = (depth<250)


#==============================================================================
# Open the data and calculate F490
#==============================================================================
#Start the loop for each float
for f in np.arange(len(floats)): #For each float
	#Open the data
	datadir = '/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/data/'+floats[f]+'/'
##datadir = '/home/apogon/Documents/SOCLIM/SOCLIM_float/data/output_lov/'+floats[f]+'/'

##  datadir = '/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/output_lov/'+floats[f]+'/'
	time.append(np.load(datadir+'TIME.npy'))
	lon.append(np.load(datadir+'LON.npy'))
	lat.append(np.load(datadir+'LAT.npy'))
	temp.append(np.load(datadir+'TEMP.npy'))
	sal.append(np.load(datadir+'SAL.npy'))
	ox.append(np.load(datadir+'OX.npy'))
	sigma.append(np.load(datadir+'SIGMA.npy'))
	chl.append(np.load(datadir+'CHLA.npy'))
	bbp.append(np.load(datadir+'BBP.npy'))
	cp.append(np.load(datadir+'CP.npy'))
	par.append(np.load(datadir+'PAR.npy'))
	ed.append(np.load(datadir+'ED490.npy'))
	bbp_spikes.append(np.load(datadir+'BBP_SPIKES.npy'))
	chl_spikes.append(np.load(datadir+'CHL_SPIKES.npy'))
	bbp_spikes[f][bbp_spikes[f]>0.048] = np.nan # threshold for acceptable data -> See Nathan for a justification
	chl_spikes[f][chl_spikes[f]>0.05] = np.nan # threshold for acceptable data -> See Nathan for a justification
		
	mld.append(np.repeat(np.nan,sigma[f].shape[1])) #Empy array for MLD
	ox_sat.append(np.empty(temp[f].shape)) ; ox_sat[f][:] = np.nan
	aou.append(np.copy(ox_sat[f]))
	F_all.append(np.repeat(np.nan,bbp[f].shape[1]))
	
	#Start the loop for each  profile
	for p in np.arange(sigma[f].shape[1]): # For each profile
		
		#Calculate oxygen saturation
		dens = (sigma[f][:,p]+1000)/1000
		ox_sat[f][:,p] = ox[f][:,p] / O2sol(temp[f][:,p],sal[f][:,p]*dens) *100 # Not used in the prediction !
		aou[f][:,p] = (O2sol(temp[f][:,p],sal[f][:,p])*dens) - ox[f][:,p]
		
		#Calculate MLD
		sigma_ref = np.nanmean(sigma[f][19:21,p])
		dep_val = depth[np.where(sigma[f][:,p]> sigma_ref + threshold)]
		if any(dep_val):
			mld[f][p] = min(dep_val)
		else:
			mld[f][p] = 50 # a minmum summer value if daily superficial temperature create a fake MLD (should not be necessary)
	
		#Unquench chlorophyll using the Xing et al., 2012 method
		mld_sel = (depth <= mld[f][p])
		chl_mld_values = chl[f][mld_sel,p]
		dep_unquench = min(depth[np.where(chl_mld_values == np.nanmax(chl_mld_values))])
		chl[f][:dep_unquench,p] = np.nanmean(chl[f][dep_unquench-1:dep_unquench+1,p])
		
		#Calculate F490 unsing Xing et al., 2011
		F_all[f][p] = F490(depth,chl[f][:,p],ed[f][:,p],100) # much deeper than Ze
		
		#Set bio-optical values to 0 in the deep		
		chl[f][:,p] = chl[f][:,p] - np.nanmin(chl[f][:,p])
		bbp[f][:,p] = bbp[f][:,p] - np.nanmin(bbp[f][:,p])
		cp[f][:,p] = cp[f][:,p] - np.nanmin(cp[f][:,p])
	#Apply a median F490 to each float
	F_all[f][F_all[f]<0.05] = np.nan # delete erroneous calculations that lead to extreme F values. Be careful with locations other than SO	
	F_all[f][F_all[f]>0.8] = np.nan # delete erroneous calculations lead to extreme F values. Be careful with locations other than SO !
	#F_all[f] = F_all[f]*0.7 # At this step you might want to deliberately bias F490. I don't think it is usefull (see Roesler et al. 2017).
	F_all[f] = nonan(F_all[f])
	chl[f] = chl[f] * (np.nanmedian(F_all[f]))
	print floats[f], 'F490=',np.nanmedian(F_all[f])

#==============================================================================
# Calcualte ratios and predict plankton groups
#==============================================================================
for f in np.arange(len(floats)): #For each float
	
	#Declare empty arrays for mld and chl_mld
	chl_int.append(np.repeat(np.nan,chl[f].shape[1]))
	chl_mld.append(np.repeat(np.nan,chl[f].shape[1]))
	temp_mld.append(np.repeat(np.nan,chl[f].shape[1]))
	bbp_mld.append(np.repeat(np.nan,bbp[f].shape[1]))
	par_mld.append(np.repeat(np.nan,bbp[f].shape[1]))
	
	#Calculate the ratio of bio-optical data
	chl_bbp.append(chl[f]/bbp[f])
	chl_bbp[f][(chl_bbp[f]<0) | (chl_bbp[f]>2e3)] = np.nan # cleanup extreme values
	chl_cp.append(chl[f]/cp[f])
	chl_cp[f][(chl_cp[f]<0) | (chl_cp[f]>1e3)] = np.nan  # cleanup extreme values
	bbp_cp.append(bbp[f]/cp[f])
	bbp_cp[f][(bbp_cp[f]<0) | (bbp_cp[f]>1)] = np.nan  # cleanup extreme values
		
	#Declare empty vectors to store predictions
	bact.append(np.empty(chl[f].shape))
	bact[f][:,:] = np.nan
	pico.append(np.copy(bact[f]))
	nano.append(np.copy(bact[f]))
	diat.append(np.copy(bact[f]))
	
	# Predict phyto classes from the surface data--------------------------------------
	for p in np.arange(sigma[f].shape[1]): # For each profile

		#calculate values within mld
		chl_int[f][p]  = np.nansum(chl[f][mld_sel,p])
		chl_mld[f][p] = np.nanmean(chl[f][mld_sel,p])
		temp_mld[f][p] = np.nanmean(temp[f][mld_sel,p])
		bbp_mld[f][p] = np.nanmean(bbp[f][mld_sel,p])
		par_mld[f][p] = np.nanmean(par[f][mld_sel,p])
		
		X = np.vstack([depth[dep_sel],
						temp[f][dep_sel,p],
						sal[f][dep_sel,p],
						ox_sat[f][dep_sel,p],
						chl[f][dep_sel,p],
						bbp[f][dep_sel,p],
						cp[f][dep_sel,p],
						chl_bbp[f][dep_sel,p],
						chl_cp[f][dep_sel,p], 
						bbp_cp[f][dep_sel,p]]).T
		if np.isfinite(X).all():
			bact[f][dep_sel,p] = pls.predict(X[:,btl_sel])[:,0]
			pico[f][dep_sel,p] = pls.predict(X[:,btl_sel])[:,1]
			nano[f][dep_sel,p] = pls.predict(X[:,btl_sel])[:,2]
			diat[f][dep_sel,p] = pls.predict(X[:,btl_sel])[:,3]
		else:
			bact[f][dep_sel,p] = np.nan
			pico[f][dep_sel,p] = np.nan
			nano[f][dep_sel,p] = np.nan
			diat[f][dep_sel,p] = np.nan

	#Scale data
	bact[f] = scale(bact[f])
	pico[f] = scale(pico[f])
	nano[f] = scale(nano[f])
	diat[f] = scale(diat[f])
	
	# calculate normalized carbon
	ctot.append(np.nansum(np.dstack((bact[f],pico[f],nano[f],diat[f])),2))
	bact_norm.append(bact[f]/ctot[f])
	pico_norm.append(pico[f]/ctot[f])
	nano_norm.append(nano[f]/ctot[f])
	diat_norm.append(diat[f]/ctot[f])



lab_fig =['a Chl','b Bact','c Pico','d Nano','e Micro']
#==============================================================================
# Plot the sections of plankton relative abundance
#==============================================================================
def contour_plot(ax,x,y,z):
	z = np.copy(z)
	z[z<0] = 0 ; z[z>1] = 1
	contours = np.arange(0,1.05,.05)
	p = ax.contourf(x,y,z,contours,cmap = 'RdYlBu_r',vmin=0,vmax=1)
	return p

def contour_plot_chl(ax,x,y,z):
	z = np.copy(z)
	z[z<0] = 0 ; z[z>1.5] = 1.5
	contours = np.arange(0,1.55,.05)
	p = ax.contourf(x,y,z,contours,cmap = 'PuBuGn',vmin=0,vmax=1.5)
	return p

for f in np.arange(len(floats)):
	T = mdates.date2num(time[f])
	fig, ax = plt.subplots(5,1,figsize=(5,8))
	plt.subplots_adjust(left = 0.15,right=0.8,top=0.95, bottom=0.07)
	ax[0].set_title(floats[f])
		
	# Plot Chlorophyll
	ctot_plot = contour_plot_chl(ax[0],T,depth,chl[f])
	pos = ax[0].get_position()
	cbar_ax = fig.add_axes([pos.x1 +0.02, pos.y0,0.02, pos.height])
	cbar = plt.colorbar(ctot_plot,cax = cbar_ax,ticks = [0,.5,1,1.5])
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label('$\mathdefault{Chl\ (mg\ m^{-3})}$')

	#Plot Bacteria
	contour_plot(ax[1],T,depth,bact[f])
	contour_plot(ax[2],T,depth,pico[f])
	contour_plot(ax[3],T,depth,nano[f])
	crel_plot = contour_plot(ax[4],T,depth,diat[f])

	# Set cosmetics
	for i in [0,1,2,3,4]:
		nice_axes(ax[i])
		ax[i].plot(T,mld[f],'-k')
		ax[i].set_ylim([250,0])
		ax[i].set_xticks([])
		ax[i].set_ylabel('Depth (m)')
	pos = ax[1].get_position()
	cbar_ax = fig.add_axes([pos.x1 +0.02, pos.y0,0.02, pos.height])
	cbar = plt.colorbar(crel_plot,cax = cbar_ax,ticks = [0,.2,.4,.6,.8,1])
	cbar.ax.tick_params(axis='y', direction='out')
	cbar.set_label('$\mathdefault{\%C_{group}}$')

	# set date format
	date_axis(ax[4],time[f],months_size=10,fraction=0.04)
	for i in [0,1,2,3,4]:
		ax[i].set_xlim(ax[4].get_xlim())
		ax[i].text(min(ax[i].get_xlim()) + (max(ax[i].get_xlim())-min(ax[i].get_xlim()))*0.01 ,240,lab_fig[i],fontweight = 'bold',ha='left')

	
#Save the figure
#	figname = '/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/fig/png/'+floats[f]+'.png'
#	plt.savefig(figname,dpi=300)
#	plt.close()

