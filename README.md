# PhD_Code

This file contains the codes used in my 2019 - 2022 PhD in Physics. 
The purpose of these codes is to take in an inflationary potential V and from it compute the scalar and tensor power specta using the Mukhanov equations. 

The first code (Background.py) solves the Hamilton-Jacobi equation

(H')^2 - (3/2)H^2 = -(1/2)V.

This allows us to compute the Hubble parameter and it's derivative. This allows for the computation of the slow-roll parmeters (epsilon, eta and Xi). The H derived from the background code is given in terms of phi and the results are splined with the e-fold array to use in the later codes. 

The next code aH.py solves the ODE, 

dln(aH)/dN = epsilon - 1. 

This allows us to calculate the scale factor (a) by taking the array aH and dividing by the H derived from the background code. 

The main code for the computation of the scalar power spectrum is the file Muk.py. This code pulls the results from Background.py and aH.py. 
