# -*- coding: utf-8 -*-

#import matplotlib
#matplotlib.use('Agg')
import os
import numpy as np
import matplotlib.pyplot as plt
plt.interactive(True)
plt.close('all')
from scipy import stats
#Import the Garcia and Gordon oxygen solubility
from O2sol import O2sol
from plot_tools import *

os.chdir('/home/flavien/Documents/these/Phytofloat')

#Define a way to plot PLS results
def PLS_plot(pls,ax1,ax2,scalex,scaley,xlab,ylab,colors):
	plt.axhline(0,color='k',zorder='bottom')
	plt.axvline(0,color='k',zorder='bottom')
	plt.scatter(pls.x_loadings_[:,ax1]*scalex[0],pls.x_loadings_[:,ax2]*scalex[1],c='k',s=15,zorder='top')
	for i in np.arange(pls.x_loadings_.shape[0]):
		plt.text(pls.x_loadings_[i,ax1]*scalex[0]+0.05,pls.x_loadings_[i,ax2]*scalex[1],xlab[i],fontsize=10,va='center',ha='left',zorder='top')
	plt.scatter(pls.y_loadings_[:,ax1]*scaley[0],pls.y_loadings_[:,ax2]*scaley[1],c = colors,s=15,zorder='top')
	for i in np.arange(pls.y_loadings_.shape[0]):
		plt.text(pls.y_loadings_[i,ax1]*scaley[0]+0.05,pls.y_loadings_[i,ax2]*scaley[1],ylab[i],fontsize=10,va='center',ha='left',color=colors[i],zorder='top')

def nonan(data):
	return data[~np.isnan(data)]


#=========================================================================================
#Open the the microphyto phyto counts data
#=========================================================================================
data = np.genfromtxt('Data/Soclim/data/data_micro.csv',delimiter='\t',dtype='S')

fov = data[1,3:-1].astype('<f8') # Field of view, do not take into account TB1
counts = data[2:,3:-1].astype('<f8') # Raw counts, do not take into account TB1
spe = data[2:,1]# Sepcies names, avoid silicoflagellates
vol = data[2:,2].astype('<f8') # Volume for each species, avoid silicoflagellates
abund = counts/fov*100 #Calculate cells abundance (cell/L)
abund = abund.T # !!! Transpose matrix : rows = stations, columns = species
V_spe= abund*vol # Biovolume for each species
sample_micro = data[0,3:-1] # Samples labels, do not take into account TB1
cla_micro = data[2:,0]

# Calculate the total biovolume of each microphyto class
V_micro = np.empty([abund.shape[0],np.unique(cla_micro).shape[0]])
for s in np.arange(V_micro.shape[0]): # For each sample
	for c in np.arange(V_micro.shape[1]): # for each class
		V_micro[s,c] = np.sum(V_spe[s,cla_micro == np.unique(cla_micro)[c]])

#Calculate each microphyto class contribution to carbon
C_micro = np.empty(V_micro.shape)
C_micro[:,0] = 0.117 * (V_micro[:,0]**0.881) # Diatom
C_micro[:,1] = 0.760 * (V_micro[:,1]**0.819) # Dino
C_micro[:,2] = 0.216 * (V_micro[:,2]**0.939) # Alloricate Cilliate
C_micro[:,3] = 0.261 * (V_micro[:,3]**0.860) #  silicoflagellates
C_micro = C_micro/12/1e6 # Convert pgC to  �molC

# Sum Carbon of Dino + CIlliate + Silicoflagellates
C_diat = C_micro[:,0]
C_other = np.sum(C_micro[:,1:],axis=1)
C_micro = np.vstack([C_diat.T,C_other.T]).T
C_micro_tot = np.sum(C_micro,axis=1)

# Sum Volume of Dino + CIlliate + Silicoflagellates
V_diat = V_micro[:,0]
V_other = np.sum(V_micro[:,1:],axis=1)
V_micro = np.vstack([V_diat.T,V_other.T]).T
V_micro_tot = np.sum(V_micro,axis=1)

#=========================================================================================
#Open the cytometry data
#=========================================================================================
data_cyto = np.genfromtxt('Data/Soclim/data/data_cyto.csv',delimiter='\t',dtype='S')
sample_cyto = data_cyto[3:,0]
sample_sel = data_cyto[3:,1].astype('<f8')
sample_sel = (sample_sel==1)
biovol_cyto = data_cyto[1,2:].astype('<f8')
carbon_cyto = data_cyto[2,2:].astype('<f8')
abund_cyto = data_cyto[3:,2:].astype('<f8')*1e3 # Convert cell/mL to cell/L
V_cyto = abund_cyto*biovol_cyto
C_cyto = abund_cyto*carbon_cyto/12/1e6 # Convert pgC to  �molC

#sum all pico (proc+syn+'pico')
V_bact = V_cyto[:,0] ; V_pico = np.sum(V_cyto[:,1:4],axis=1) ; V_nano = V_cyto[:,-1]
C_bact = C_cyto[:,0] ; C_pico = np.sum(C_cyto[:,1:4],axis=1) ; C_nano = C_cyto[:,-1]

#=========================================================================================
#Merge phyto data and Calculate the relative contribution of each class to total biovolume and carbon
#=========================================================================================
# Define phyto classes and associated colors
plank_lab = ['Bact','Pico','Nano','Micro']
plank_colors = ['0.8','c','g','r']

V_all = np.vstack([V_bact.T,V_pico.T,V_nano.T,V_micro_tot.T]).T
C_all = np.vstack([C_bact.T,C_pico.T,C_nano.T,C_micro_tot.T]).T
C_all_tot = np.sum(C_all,axis=1)

# Calculate the relative biovolume
V_all_rel = (V_all.T/np.sum(V_all,axis=1)).T
C_all_rel = (C_all.T/np.sum(C_all,axis=1)).T  # Relative contribution of each class to total microphyto biovo

#=========================================================================================
#Open the ctd data + poc/PON data
#=========================================================================================
#data = np.genfromtxt('data/data_btl.csv',delimiter='\t',dtype='S')
data = np.genfromtxt('Data/Soclim/data/data_btl_JULIA.csv',delimiter='\t',dtype='S')


sample_btl = data[1:,0]
poc = data[1:,-3].astype('<f8')
pon = data[1:,-2].astype('<f8')
data_ctd = data[1:,4:-3].astype('<f8')
ox_sol = O2sol(data_ctd[:,1],data_ctd[:,2])*1.027 # approximate density -> better use the real one
ox_sat = (data_ctd[:,3]+11)/ox_sol*100 # Check the CTD calibration ... +11 �mol/L
data_ctd[:,3] = ox_sat # Replace oxygen concentration by oxygen saturation


#=========================================================================================
# Remove Unecessary stations
data_ctd = data_ctd[sample_sel,:]

V_all = V_all[sample_sel,:] ; V_all_rel = V_all_rel[sample_sel,:]
C_all = C_all[sample_sel,:] ; C_all_rel = C_all_rel[sample_sel,:]
C_all_tot = C_all_tot[sample_sel]
sample_btl = sample_btl[sample_sel]
sample_micro = sample_micro[sample_sel]
sample_cyto = sample_cyto[sample_sel]
poc = poc[sample_sel]
pon = pon[sample_sel]

#Optionnal : Check the corresondance between phyto samples and CTD samples
for i in np.arange(len(sample_micro)):
	print (sample_micro[i] , sample_btl[i] , sample_cyto[i])

#=========================================================================================
# Select data : do not use oxygen saturation !!!
btl_lab = np.array(['Depth','T','S','Oxsat','Chl','$\mathdefault{b_{bp}}$','$\mathdefault{c_{p}}$','Chl:$\mathdefault{b_{bp}}$','Chl:$\mathdefault{c_{p}}$','$\mathdefault{b_{bp}}$:$\mathdefault{c_{p}}$'])
btl_sel = np.array([True,True,True,False,True,True,True,True,True,True]) # this selector removes Oxy saturation as predictor
btl_lab = btl_lab[btl_sel]
X = data_ctd[:,btl_sel]
Y = C_all_rel
# Calculate standardised varaibles for correlations
X_std = (X-np.nanmean(X,axis=0))/np.nanstd(X,axis=0)
Y_std = (Y-np.nanmean(Y,axis=0))/np.nanstd(Y,axis=0)

#=========================================================================================
# Calculate correlation coeff between variables
corr = np.empty([X.shape[1],Y.shape[1]])
pval = np.empty([X.shape[1],Y.shape[1]])
for i in np.arange(corr.shape[0]):
	for j in np.arange(corr.shape[1]):
		slope, intercept, r_value, p_value, std_err = stats.linregress(X_std[:,i],C_all_rel[:,j])
		corr[i,j] = r_value
		pval[i,j] = p_value

#=========================================================================================
# Perform PLS on relative contrib to Ctot
from sklearn.cross_decomposition import PLSRegression
pls = PLSRegression(n_components=X.shape[1],scale=False).fit(X,Y) # Maximum number of components = number of potential predictors

B = pls.coef_
Y_pred = pls.predict(X)
rsqd_pls = pls.score(X,Y)
print('\n-----------------------------\nPLS score = ', rsqd_pls)

# Estimate PLS performance
print ('PLS results:')
Y_all = np.reshape(Y,Y.size)
Y_pred_all = np.reshape(Y_pred,Y_pred.size)
for i in np.arange(len(plank_lab)):
	rmse = np.mean((Y[:,i] - Y_pred[:,i])**2)**0.5
	slope, intercept, r_value, p_value, std_err = stats.linregress(Y[:,i],Y_pred[:,i])
	print(plank_lab[i],'slope=',slope,' rsqd=',r_value**2,' pval=',p_value,' rmse=',rmse)
slope, intercept, r_value_all, p_value, std_err = stats.linregress(Y_all,Y_pred_all)
print ('All pooled','slope=', slope,'rsqd=', r_value_all**2,'pval=', p_value,'rmse=', np.mean((Y_all - Y_pred_all)**2)**0.5)
print ('-----------------------------')


#=========================================================================================
# Plot the PLS figures
fig, ax = plt.subplots(1,2,figsize=(6,3))
plt.subplots_adjust(top = 0.85,bottom=0.25,wspace=0.7,left=0.12,right=0.95)
lab = ['$\mathdefault{\%C_{bact}}$','$\mathdefault{\%C_{pico}}$','$\mathdefault{\%C_{nano}}$','$\mathdefault{\%C_{micro}}$']
# 1: plot the linear correlation coeffs
plt.sca(ax[0])
cplot = ax[0].pcolor(np.flipud(corr),vmin=-1,vmax=1,cmap='RdBu_r')
cbar = plt.colorbar(cplot,shrink=0.5,aspect=10,ticks=[-1,-0.5,0,0.5,1])
cbar.ax.tick_params(direction='out')
ax[0].set_yticks(np.arange(B.shape[0])+0.5)
ax[0].set_xticks(np.arange(B.shape[1])+0.5)
ax[0].set_yticklabels(btl_lab[::-1])
ax[0].set_xticklabels(lab,rotation=45,ha='center',va='center')
ax[0].tick_params(axis='x', which='major',pad=12)
x = np.arange(Y.shape[1])
y = np.arange(X.shape[1])[::-1]
for i in np.arange(corr.shape[0]):
	for j in np.arange(corr.shape[1]):
		if pval[i,j]<0.015:
			plt.plot(x[j]+0.5,y[i]+0.5,'ok',mfc='none',markersize=6)
#Add colored points for each plankton group
pos = ax[0].get_position()
axpt = fig.add_axes([pos.x0-0.05,pos.y0-0.2,pos.width,0.2])
axpt.set_xlim(ax[0].get_xlim()) ; axpt.set_ylim([0,1])
axpt.plot(0.7,0.3,'ok',mfc='0.8')
axpt.plot(1.7,0.3,'ok',mfc='c')
axpt.plot(2.7,0.3,'ok',mfc='g')
axpt.plot(3.7,0.3,'ok',mfc='r')
axpt.axis('off')

# 2: PLot the obs vs. pred
for i in np.arange(Y.shape[1]):
	ax[1].plot(Y[:,i],Y_pred[:,i],'ok',mfc=plank_colors[i])
ax[1].plot([1e-3,1],[1e-3,1],'--k',zorder='bottom')
ax[1].set_xlabel('$\mathdefault{\%C_{group}\ observed}$')
ax[1].set_ylabel('$\mathdefault{\%C_{group}\ predicted}$')
ax[1].set_xlim([1e-2,1]) ; ax[1].set_ylim([1e-2,1])
ax[1].set_yscale('log') ; ax[1].set_xscale('log')
text = '$\mathdefault{R^{2}=\ }$'+ str(np.round(rsqd_pls,2))
ax[1].set_title(text,fontsize=10)
for i in [0,1]:
	nice_axes(ax[i])
ax[0].text(0,9.7,'a',fontsize=12)
ax[1].text(1e-2,1.5,'b',fontsize=12)
