# PhD_Code

This file contains the codes used in my 2019 - 2022 PhD in Physics. 
The purpose of these codes is to take in an inflationary potential V and from it compute the scalar and tensor power specta using the Mukhanov equations. 

The first code (Background.py) solves the Hamilton-Jacobi equation

(H')^2 - (3/2)H^2 = -(1/2)V.

This allows us to compute the Hubble parameter and it's derivative. This allows for the computation of the slow-roll parmeters (epsilon, eta and Xi). The H derived from the background code is given in terms of phi and the results are splined with the e-fold array to use in the later codes. To compute the slow-roll parameters change the v_func(t) and vprime(t) functions to the desired inflationary models in Background.py and set the initial conditions. 

The next code aH.py solves the ODE, 

dln(aH)/dN = epsilon - 1. 

This allows us to calculate the scale factor (a) by taking the array aH and dividing by the H derived from the background code. 

The main code for the computation of the scalar power spectrum is file Muk.py. This code pulls the results from Background.py and aH.py and plots the power spectrum P(k). It does this by solving the Mukhanov equation (in terms of e-fold) for a range of wavenumbers (k) to find the value of the Mukhanov variable at horizon exit. 




The code output.py runs the full package and plots H, epsilon Vs e-fold and the two power spectra Vs K.
