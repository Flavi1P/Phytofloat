# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.dates as mdates

# Define default rendering parameters
from matplotlib import rcParams
rcParams['font.size'] = 10
rcParams['axes.labelsize'] = 10
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['legend.fontsize'] = 10
rcParams['axes.linewidth'] = 0.5
rcParams['lines.linewidth'] = 0.5
rcParams['lines.markersize'] = 3
rcParams['patch.linewidth'] = 0.5
rcParams['font.sans-serif'] = 'Arial'

#Define a nice way to show axes
def nice_axes(ax):
    ax.tick_params(axis='both', which='both',direction='out')
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()

#Define a nice way to suppelemntary axes
def nice_axesup(ax):
    ax.tick_params(axis='both', which='both', direction='out')
    ax.yaxis.tick_right()
    ax.xaxis.tick_bottom()

#Define a basic boxplot function
def simple_boxplot(ax,x,data,width,showdata,showmean):
	if len(np.unique(data))>1:
		median = np.median(data)
		mean = np.nanmean(data)
		q1 = np.percentile(data,25)
		q2 = np.percentile(data,75)
		interquartile = q2-q1
		m1 = min(data[data>median-(1.5*interquartile)])
		m2 = max(data[data<median+(1.5*interquartile)])
		ax.fill_between([x-width/2.,x+width/2.],[q1,q1],[q2,q2],facecolor='0.8')
		ax.plot([x-width/2.,x+width/2.],[median,median],'-k')
		ax.plot([x,x],[m1,m2],'-k')
		
		if showdata==True:
			ax.plot(np.repeat(x,len(data)),data,'ok',mfc='none',zorder=25)
		if showmean==True:
			ax.plot(x,mean,'ok')
	else:
		print 'Boxplot warning = constant data as input'
	return median,mean,q1,q2
	

def date_axis(ax,t,months_size,fraction):
	'''
	Fuction to add nice date format :months and and years
	
	inputs:
		ax: the origina axis to be formatted
		t: the datetime vector of original data
		months_size: the fontweight of the months
		fraction: the height of new axes as a fraction of the original axes
	'''
	#Be sure the date format is set to US
	import locale
	locale.setlocale(locale.LC_ALL,'en_US.UTF-8')
	
	from datetime import timedelta
	#find the first day of the month  and first day of the year from the time array
	first_day_month = np.unique(np.array([datetime(i.year, i.month, 1) for i in t]))
	last_tick = first_day_month[-1]+timedelta(31)
	first_day_month = np.hstack([first_day_month,last_tick])
	first_day_year = np.unique(np.array([datetime(i.year, 1, 1) for i in t]))
	months_names = [datetime.strftime(i,'%b')[0] for i in first_day_month]
	years_names =  [datetime.strftime(i,'%Y') for i in first_day_year]

	#Set axes limits
	pos = ax.get_position()
	ax_time = plt.gcf().add_axes([pos.x0,pos.y0-fraction,pos.width,fraction])
	ax_time.axis('off')
	ax_time.set_ylim([0,1])
	ax_time.set_xlim([first_day_month[0],first_day_month[-1]])
	ax.set_xlim(ax_time.get_xlim())

	#Set axes ticks and labels
	ax.set_xticks(first_day_month) ; ax.set_xticklabels([])
	for i in np.arange(len(months_names)-1)[1::2] :
		ax_time.text(first_day_month[i]+timedelta(14),0.7,months_names[i],va='center',ha='center',fontsize=months_size)
	for i in np.arange(len(years_names)):
		ax_time.axvline(first_day_year[i],color='k')
		if first_day_year[i]>= first_day_month[0]:
			ax_time.text(first_day_year[i]+timedelta(14),0,years_names[i],va='bottom',ha = 'left')

#define nice_boxplot:
def nice_boxplot(ax,values,positions,color):
	bp = ax.boxplot(values,positions=positions,widths = 0.6,sym='',patch_artist=True)
	plt.setp(bp['boxes'], color='black',facecolor=color,linewidth=0.5)
	plt.setp(bp['whiskers'],color='k',linestyle='-')
	plt.setp(bp['medians'], color='black')

#Define a was to add labels automatically
def add_labels(fig):
	letters = 'abcdefghijklmnopqrstuvwxyz'
	for i in np.arange(len(fig.axes)):
		curr_ax = plt.gcf().axes[i]
		xlim = curr_ax.get_xlim()
		ylim = curr_ax.get_ylim()
		curr_ax.text(min(xlim)-(max(xlim)-min(xlim))*0.2,max(ylim)+(max(ylim)-min(ylim))*0.15,letters[i],fontsize=12)

def reg_plot(x,y,alpha,confidence,ax,color):
    '''
    
    Calculates and plot a regression line + CI of regression
    
    Parameters
    ----------
    x : original x data
    y :  original y data
    alpha : probability (e.g. alpha = 0.05) for th F statistic
    confidence : if True, plot the confidence intervals as dotted lines
    ax : the axes in which the lines are plotted
    color : the color of the regressions line + CI lines (dashed by default)
    
    Returns
    -------
    slope :  the slope of the regression
    intercept : the intercept of the regression
    rsqd : Rsquared of the regression line
    pval : pval of the regression resulting from the F test of the full model
    mdl : model results
    
    '''    
    
    import numpy as np
    import statsmodels.api as sm
    from scipy.stats import t
    
    sel = (~np.isnan(x)) & (~np.isnan(y))
    x = x[sel]
    y = y[sel]
    
    #Define confidence limit and number of predicted points
    c_limit=1-alpha/2
    n_pred = 100
    
    #Perform the fit and prediciton
    X = sm.add_constant(x)
    res = sm.OLS(y,X).fit()
    y_fit = res.params[0] + res.params[1]*x
    y_err = y -y_fit
    
    #Calculate the appropriate t value
    mean_x = np.mean(x)			# mean of x
    n = len(x)				# number of samples in origional fit
    tstat = t.ppf(c_limit, n-1)         # appropriate t value
    s_err = np.sum(np.power(y_err,2))	# sum of the squares of the residuals

    # create series of new test x-values to predict for
    p_x = np.linspace(np.min(x),np.max(x),n_pred)
    p_y = res.params[0] + res.params[1]*p_x
    #Calculate confidence intervals
    confs = tstat * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
    ((np.sum(np.power(x,2)))-n*(np.power(mean_x,2))))))
   
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    ax.plot(p_x,p_y,'-',color=color,zorder='top')
    if confidence == True:
		ax.plot(p_x,lower,'--',color=color,zorder='top')
		ax.plot(p_x,upper,'--',color=color,zorder='top')
    
    rsqd = res.rsquared
    pval = res.f_pvalue
    slope = res.params[1]
    intercept = res.params[0]
    return slope,intercept,rsqd,pval,res
