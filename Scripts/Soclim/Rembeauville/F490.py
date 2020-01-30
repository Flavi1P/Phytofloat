# -*- coding: utf-8 -*-

def F490(depth,fluo,ed490, depth_max):
	'''
	Function to calculate F490 as described in Xing et al., 2011
	
	Inputs:
		- Depth (m)
		- Raw Fluorescence
		- Ed490
		- Depth max : the maximum depth to calculate the F490 value
		
	Output:
		- F490 value
	'''
	import numpy as np
	from scipy import stats
	Kw = 0.01660
	alpha = 0.0825 # Parameters of Kd490 =f(Chl) from Morel 1988 (probably not true for SO waters ... !)
	beta = 0.6529 # Parameters of Kd490 =f(Chl) from Morel 1988 (probably not true for SO waters ... !)
	
	sel1 = depth<depth_max # Select data above depth max (generally 100 m)
	sel2 = ed490<(0.5*max(ed490)) # Select data below two optical depths (avoid quenching)
	sel = sel1 & sel2
	Ed = ed490[sel]
	fluo = fluo[sel]
	if len(fluo)>10: # If there is enough datapoints ~ 10 minimum
		Cn = np.empty(Ed.shape[0])
		An = np.empty(Ed.shape[0])
		
		for j in np.arange(Ed.shape[0]): # For each depth interval
			Cn[j] = np.log(Ed[j]) + Kw # Calculate observed attenuation at 490 nm (taking into account water)
			An[j] = np.sum(alpha*(fluo[:j]**beta)) # Calculate theoritical attenuation from fluo data
		x = An
		y = Cn
		      
		slope, intercept, r_value, p_value, std_err = stats.linregress(x,y) # Perform a first regression betwen theoritical and observed attenuation
		fit = (slope*x) + intercept
		resid = np.abs(y-fit) # Calculate residuals
		
		while r_value**2 <0.98: # While R2 of regression is <0.98 (mainly outliers due to clouds)
			x = x[resid!=np.max(resid)] # discard datapoints with maximum residual value
			y = y[resid!=np.max(resid)] # discard datapoints with maximum residual value
			if (len(x)>2) & (len(y)>2):
				slope, intercept, r_value, p_value, std_err = stats.linregress(x,y) # Perform a new regression
				fit = (slope*x) + intercept 
				resid = np.abs(y-fit) # Calculate new residuals and so on ...
			else:
				slope= np.nan
				r_value=np.nan
		
		slope=-slope
		F490 = slope**(1/beta) # Store the result
		rsquared = r_value**2 # Store the R2 (optional)
	else:
		F490 = np.nan
		rsquared = np.nan
	
	return F490