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
	
#==============================================================================
# Perform PLS
#==============================================================================

#==============================================================================
# Declare the matrices
#==============================================================================


#floats = ['049b','036b','037c','107c','104c']
#floats = ['049b','037c','107c','104c']
floats = ['049b']
floats = ['104c']
floats = ['107c']
floats = ['037c','037c']


#floats = ['036b']


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

# Structures to store float data
time1 = []
lon1 = []
lat1 = []
temp1 = []
sal1 = []
ox1 = []
ox_sat1 = []
aou1 = []
sigma1 = []
chl1 = []
bbp1 = []
cp1 = []
par1 = []
ed1 = []

F_all = []
a = []
b = []

threshold = 0.03 # Density crierion (De Boyer Montégut et al., 2004)
depth = np.arange(1000)
depth1 = np.arange(1001)

dep_sel = (depth<250)


#Open the data

datadir = '/home/admt/PROGRAM_AP/CODE_REMBAUVILLLE/data/'+floats[1]+'/'
timeD = np.loadtxt(datadir+'TIME.txt',dtype=int, delimiter=',')
	
for ff in range(len(timeD)):
	a=str(timeD[ff,0])+'-'+str(timeD[ff,1])+'-'+str(timeD[ff,2])+'-'+str(timeD[ff,3])+'-'+str(timeD[ff,4])
	b.append(datetime.strptime(a, "%Y-%m-%d-%H-%M"))
time=np.transpose(b)			
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
	
fig, ax = plt.subplots()
	
plt.subplot(221)
plt.plot(sigma[:,4], -depth, '-')
plt.title(time[4])
	
plt.subplot(222)
plt.plot(chl[:,4], -depth, '-')
plt.title(str(lat[4])+' '+ str(lon[4]))
	
plt.subplot(223)
plt.plot(bbp[:,4], -depth, '-')

plt.subplot(224)
plt.plot(cp[:,4], -depth, '-')

figname = '/var/www/oao/bioargo/PHP/AP/REMBAU/PROF_TXT_'+floats[1]+'.png'
plt.savefig(figname,dpi=300)
plt.close()





time1.append(np.load(datadir+'TIME.npy'))
lon1.append(np.load(datadir+'LON.npy'))
lat1.append(np.load(datadir+'LAT.npy'))
temp1.append(np.load(datadir+'TEMP.npy'))
sal1.append(np.load(datadir+'SAL.npy'))
ox1.append(np.load(datadir+'OX.npy'))
sigma1.append(np.load(datadir+'SIGMA.npy'))
chl1.append(np.load(datadir+'CHLA.npy'))
bbp1.append(np.load(datadir+'BBP.npy'))
cp1.append(np.load(datadir+'CP.npy'))
par1.append(np.load(datadir+'PAR.npy'))
ed1.append(np.load(datadir+'ED490.npy'))

	
f = 0

fig, ax = plt.subplots()
	
plt.subplot(221)
plt.plot(sigma1[f][:,2], -depth1, '-')
plt.title(time1[f][2])
	
plt.subplot(222)
plt.plot(chl1[f][:,2], -depth1, '-')
plt.title(str(lat1[f][2])+' '+str(lon1[f][2]))
	
plt.subplot(223)
plt.plot(bbp1[f][:,2], -depth1, '-')
	
plt.subplot(224)
plt.plot(cp1[f][:,2], -depth1, '-')
	
figname = '/var/www/oao/bioargo/PHP/AP/REMBAU/PROF_NPY_'+floats[1]+'.png'
plt.savefig(figname,dpi=300)
plt.close()

		

fig, ax = plt.subplots()
	
plt.subplot(321)
plt.plot(temp[:,4], -depth, '-')
plt.plot(temp1[f][:,2], -depth1, '-')

plt.title(time1[f][2])
	
plt.subplot(322)
plt.plot(chl[:,4], -depth, '-')
plt.plot(chl1[f][:,2], -depth1, '-')
plt.title(str(lat1[f][2])+' '+str(lon1[f][2]))
	
plt.subplot(323)
plt.plot(bbp[:,4], -depth, '-')
plt.plot(bbp1[f][:,2], -depth1, '-')
	
plt.subplot(324)
plt.plot(cp[:,4], -depth, '-')
plt.plot(cp1[f][:,2], -depth1, '-')
	
plt.subplot(325)
plt.plot(sigma[:,4], -depth, '-')
plt.plot(sigma1[f][:,2], -depth1, '-')
	
	
figname = '/var/www/oao/bioargo/PHP/AP/REMBAU/PROF_ALL'+floats[1]+'.png'
plt.savefig(figname,dpi=300)
plt.close()







