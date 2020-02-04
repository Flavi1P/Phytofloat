# -*- coding: utf-8 -*-
def O2sol(T,S):
	'''
	Function to calculate oxygen solubility from Garcia and Gordon, 1992
	Inputs:
		T the temperature
		S the salinity
	Returns
		Oxygen concentration at saturation (µmol/kg)
	Check
		Solubility at T = 10°C and S = 35, O2 = 274.610
	'''
	import numpy as np
	Ts = np.log((298.15-T)*((273.15+T)**-1))
	A0 = 5.80871
	A1 = 3.20291
	A2 = 4.17887
	A3 = 5.10006
	A4 = -9.86643e-2
	A5 = 3.80369
	B0 = -7.01577e-3
	B1 = -7.70028e-3
	B2 = -1.13864e-2
	B3 = -9.51519e-3
	C0 = -2.75915e-7
	O2sat = A0 + (A1*Ts) + (A2*(Ts**2)) + (A3*(Ts**3)) + (A4*(Ts**4)) + (A5*(Ts**5)) + (S*(B0 + (B1*Ts) + (B2*(Ts**2)) + (B3*(Ts**3)))) + (C0*(S**2))
	return np.exp(O2sat)
