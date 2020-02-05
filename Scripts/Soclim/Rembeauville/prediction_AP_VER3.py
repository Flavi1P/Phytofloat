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
from datetime import datetime, date, time

import warnings
warnings.simplefilter(action = 'ignore', category = RuntimeWarning)

import matplotlib
matplotlib.use('Agg')


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

#set my project wd
os.chdir('/home/flavien/Documents/these/Phytofloat')
#==============================================================================
# Perform PLS
#==============================================================================
exec(open('Scripts/Soclim/Rembeauville/PLS.py').read())

#==============================================================================
# Declare the matrices
#==============================================================================


floats = ['049b','036b','037c','107c','104c']
#floats = ['049b','037c','107c','104c']
#floats = ['036b']

#floats = ['049b']
#floats = ['104c']
#floats = ['107c']
#floats = ['036b']


#floats = ['036b']
#floats = ['019b']

#floats = ['059c']
#floats = ['064b']
#floats = ['006b']
#floats = ['019b']
#floats = ['048b']
#floats = ['050b']

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
a = []
b = []

threshold = 0.03 # Density crierion (De Boyer Montégut et al., 2004)
depth = np.arange(1000)
dep_sel = (depth<250)


#==============================================================================
# Open the data and calculate F490
#==============================================================================
#Start the loop for each float
print(floats)
for f in np.arange(len(floats)): #For each float
    print(f)
#Open the data
    datadir = 'Data/Soclim/data/'+floats[f]+'/'
    timeD = np.loadtxt(datadir+'TIME.txt',dtype=int, delimiter=',')

    for ff in range(len(timeD)):
        a=str(timeD[ff,0])+'-'+str(timeD[ff,1])+'-'+str(timeD[ff,2])+'-'+str(timeD[ff,3])+'-'+str(timeD[ff,4])
        b.append(datetime.strptime(a, "%Y-%m-%d-%H-%M"))
    time=np.transpose(b)

#	time2  = []
#	time2.append(np.load(datadir+'TIME.npy'))
    lon=np.loadtxt(datadir+'LON.txt')
    lat=np.loadtxt(datadir+'LAT.txt')
    temp=np.loadtxt(datadir+'TEMP.txt',dtype=float, delimiter=',')
    sal=np.loadtxt(datadir+'SAL.txt',dtype=float, delimiter=',')
    ox=np.loadtxt(datadir+'OX.txt',dtype=float, delimiter=',')
    sigma=np.loadtxt(datadir+'SIGMA.txt',dtype=float, delimiter=',')
    chl=np.loadtxt(datadir+'CHLA.txt',dtype=float, delimiter=',')
    bbp=np.loadtxt(datadir+'BBP.txt',dtype=float, delimiter=',')
    cp=np.loadtxt(datadir+'CP.txt',dtype=float, delimiter=',')
    par=np.loadtxt(datadir+'PAR.txt',dtype=float, delimiter=',')
    ed=np.loadtxt(datadir+'ED490.txt',dtype=float, delimiter=',')



#	bbp_spikes.append(np.load(datadir+'BBP_SPIKES.npy'))
#	chl_spikes.append(np.load(datadir+'CHL_SPIKES.npy'))
#	bbp_spikes[f][bbp_spikes[f]>0.048] = np.nan # threshold for acceptable data -> See Nathan for a justification
#	chl_spikes[f][chl_spikes[f]>0.05] = np.nan # threshold for acceptable data -> See Nathan for a justification

    mld.append(np.repeat(np.nan,sigma.shape[1])) #Empy array for MLD
    ox_sat.append(np.empty(temp.shape)) ; ox_sat[f][:] = np.nan
    aou.append(np.copy(ox_sat[f]))
    F_all.append(np.repeat(np.nan,bbp.shape[1]))

	#Start the loop for each  profile

    for p in np.arange(sigma.shape[1]): # For each profile

		#Calculate oxygen saturation
        dens = (sigma[:,p]+1000)/1000


        ox_sat[f][:,p] = (ox[:,p] * 100) / O2sol(temp[:,p],sal[:,p]*dens) *100 # Not used in the prediction !

#		aou[f][:,p] = (O2sol(temp[:,p],sal[:,p])*dens) - ox[:,p]
#		ox_sat[f][:,p] = ox[:,p] *100 # Not used in the prediction !


		#Calculate MLD
        sigma_ref = np.nanmean(sigma[19:21,p])
        dep_val = depth[np.where(sigma[:,p]> sigma_ref + threshold)]
        if any(dep_val):
            mld[f][p] = min(dep_val)
        else:
            mld[f][p] = 50 # a minmum summer value if daily superficial temperature create a fake MLD (should not be necessary)

		#Unquench chlorophyll using the Xing et al., 2012 method
        mld_sel = (depth <= mld[f][p])
		#mld_sel = (depth <= 50)
        chl_mld_values = chl[mld_sel,p]
#		dep_unquench = min(depth[np.where(chl_mld_values == np.nanmax(chl_mld_values))])
#		chl[:dep_unquench,p] = np.nanmean(chl[dep_unquench-1:dep_unquench+1,p])

		#Calculate F490 unsing Xing et al., 2011
        F_all[f][p] = F490(depth,chl[:,p],ed[:,p],100) # much deeper than Ze

	    #Set bio-optical values to 0 in the deep
        chl[:,p] = chl[:,p] - np.nanmin(chl[:,p])
        bbp[:,p] = bbp[:,p] - np.nanmin(bbp[:,p])
        cp[:,p] = cp[:,p] - np.nanmin(cp[:,p])

	#Apply a median F490 to each float
    F_all[f][F_all[f]<0.05] = np.nan # delete erroneous calculations that lead to extreme F values. Be careful with locations other than SO
    F_all[f][F_all[f]>0.8] = np.nan # delete erroneous calculations lead to extreme F values. Be careful with locations other than SO !
	#F_all[f] = F_all[f]*0.7 # At this step you might want to deliberately bias F490. I don't think it is usefull (see Roesler et al. 2017).
    F_all[f] = nonan(F_all[f])
    chl = chl * (np.nanmedian(F_all[f]))


    #	chl[f] = chl[f] *  0.24
    print(floats[f], 'F490=',np.nanmedian(F_all[f]))





        #==============================================================================
        # Calcualte ratios and predict plankton groups
        #==============================================================================
        #for f in np.arange(len(floats)): #For each float

	#Declare empty arrays for mld and chl_mld
    chl_int.append(np.repeat(np.nan,chl.shape[1]))
    chl_mld.append(np.repeat(np.nan,chl.shape[1]))
    temp_mld.append(np.repeat(np.nan,chl.shape[1]))
    bbp_mld.append(np.repeat(np.nan,bbp.shape[1]))
    par_mld.append(np.repeat(np.nan,bbp.shape[1]))

	#Calculate the ratio of bio-optical data
    chl_bbp.append(chl/bbp)
    chl_bbp[f][(chl_bbp[f]<0) | (chl_bbp[f]>2e3)] = np.nan # cleanup extreme values

    chl_cp.append(chl/cp)
    chl_cp[f][(chl_cp[f]<0) | (chl_cp[f]>1e3)] = np.nan  # cleanup extreme values

    bbp_cp.append(bbp/cp)
    bbp_cp[f][(bbp_cp[f]<0) | (bbp_cp[f]>1)] = np.nan  # cleanup extreme values



	#Declare empty vectors to store predictions
    bact.append(np.empty(chl.shape))
    bact[f][:,:] = np.nan
    pico.append(np.copy(bact[f]))
    nano.append(np.copy(bact[f]))
    diat.append(np.copy(bact[f]))



	# Predict phyto classes from the surface data--------------------------------------
    for p in np.arange(sigma.shape[1]): # For each profile

		#calculate values within mld
        chl_int[f][p]  = np.nansum(chl[mld_sel,p])
        chl_mld[f][p] = np.nanmean(chl[mld_sel,p])
        temp_mld[f][p] = np.nanmean(temp[mld_sel,p])
        bbp_mld[f][p] = np.nanmean(bbp[mld_sel,p])
        par_mld[f][p] = np.nanmean(par[mld_sel,p])

#		print lat[1]
#		print lon[1]
#		print time[1]

#		print 'TEMP'
#		print temp[:,1]
#		print 'sal'
#		print sal[:,1]
#		print 'ox'
#		print chl[:,1]
##		print 'bbp'
#		print bbp[:,1]
#		print 'cp'
#		print cp[:,1]
#		print ox_sat[f][:,1]

        X = np.vstack([depth[dep_sel],
						temp[dep_sel,p],
						sal[dep_sel,p],
						ox_sat[f][dep_sel,p],
						chl[dep_sel,p],
						bbp[dep_sel,p],
						cp[dep_sel,p],
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

    T = mdates.date2num(time)





    import numpy
    print(numpy.size(bact))
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/time'+floats[f]+'.txt', T, delimiter=",")
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/depth'+floats[f]+'.txt', depth, delimiter=",")

    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/chl'+floats[f]+'.txt', chl, delimiter=",")
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/bbp'+floats[f]+'.txt', bbp, delimiter=",")
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/cp'+floats[f]+'.txt', cp, delimiter=",")

    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/bact'+floats[f]+'.txt',bact[f], delimiter=",")

    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/pico'+floats[f]+'.txt',pico[f], delimiter=",")
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/nano'+floats[f]+'.txt',nano[f], delimiter=",")
    numpy.savetxt('Data/Soclim/DATA_TEXT_OUT/diat'+floats[f]+'.txt',diat[f], delimiter=",")
