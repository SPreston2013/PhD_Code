#Code to compute the power spectrum from the solution of the Mukhanov equation.
import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.interpolate import UnivariateSpline
from Background import Background_Output
from aH import Ni
from aH import aH
from aH import aHt
import time
start_time = time.time()


print()

print('Running Scalar Perturbation Code...')

print()

plt.rcParams["font.size"] = 6
plt.autoscale(enable=True, axis='both', tight=None)

points = int(10.) #Plot points used in k array.
######################################################################
#H, H', N and epsilon from background code.
H_a = Background_Output[0]

Hp_a = Background_Output[1]

ep_a = Background_Output[2]

N_a = Background_Output[3]
######################################################################
#Array for k to use later.
k_a = np.geomspace(start = 0.05, stop = 1000. , num = points)
######################################################################
#Spline epsilon from the background code.
epfunc = UnivariateSpline(N_a, ep_a, s=0, k=3)

#Spline H and H' from the background code, then take the derivatives.
Hfunc =  UnivariateSpline(N_a, H_a, s=0,k=3)
Hpfunc = UnivariateSpline(N_a, Hp_a, s=0,k=3)
Hppfunc =  Hpfunc.derivative(n=1)
Hpppfunc =  Hpfunc.derivative(n=2)

#Spline aH from the Vacc code.
aHfunc = UnivariateSpline(aHt, aH, s=0, k=3)
######################################################################
def H(t):
    """Hubble Parameter from background code"""
    return Hfunc(t)[()] #Mpl
    
def H_prime(t):
    """Derivative of H wrt phi"""
    return Hpfunc(t)[()]

def aH(t):
    """aH from vacc code"""
    return aHfunc(t)[()] #Mpc^(-1)

#Define Slow roll parameters.
def ep(t):
    """Epsilon from background code"""
    return epfunc(t)[()] #No unit

def eta(t):
    return 2*(Hppfunc(t)[()]/H(t)) #No unit

def xi(t):
    return 4*((Hpfunc(t)[()]*Hpppfunc(t)[()])/(H(t))**2) #No unit

def sigma(t):
    return 2*eta(t) - 4*ep(t) #No unit

#Function used in the Mukhanov equation
def f(t):
    return (2 - 4*ep(t) - 1.5*sigma(t) - 2*ep(t)**2 + 0.25*sigma(t)**2 + xi(t)) #No unit

#Define z
def z(t):
    return 2*aH(t)*(H_prime(t)/H(t)**2) #Need unit Mpc(-1).

######################################################################
#Mukhanov equation in terms of the e-fold number.
def sol(t, Y, k):
    """Mukhanov equation in terms of e-fold"""
    return (Y[1], (1 - ep(t))*Y[1] - ((k/aH(t))**2 - f(t))*Y[0])
    
######################################################################
#Mukhanov Power spectrum
def m_spec(k, R):
    """Mukhanov Power spectrum"""
    return np.array((k**3/(2*np.pi**2)) * R)
    
######################################################################
#Create empty lists.
ur = []

ui = []

R_r = []

R_i = []

s_spec = []
######################################################################
#Create Ni function.
def Ni_find(t):
    """ODE equation from aH code."""
    return ep(t) - 1

######################################################################
#Solve Mukhanov equation for each mode k in k_a.

for k in k_a:

#Define e-fold to start and stop the code, runs from Ni to Ne
    
    area = 0.0
    
    step = 0.001
    
    left = 60.0 - np.log(k/k_a[0])
    
    right = left + step
    
    while(area>np.log(.1)):
    
        left = left + step
        
        right = right + step
    
        trap = (step/2.)*(Ni_find(left)+ Ni_find(right))
        
        area = area + trap
     
     
    Ni = (left+right)/2.
    
    Ne = 0.5*(60.0 - np.log(k/k_a[0]))
#u_k real I.C
    ur_0 = [1/(np.sqrt(2*k)), 0] #1/root Mpc^(-1)
#u_k imaginary I.C
    ui_0 = [0, -np.sqrt(k/2)*(1/aH(Ni))] #1/root Mpc^(-1)

#Use solve_ivp with the RK45 method, solr is the real part of u_k and soli is the Imaginary.
    solr = solve_ivp(sol, (Ni, Ne), ur_0, method = 'RK45', t_eval = np.geomspace(start = Ni, stop = Ne, num = 10000), vectorized = True,  dense_output=True, args=(k,))
        
    soli = solve_ivp(sol, (Ni, Ne), ui_0, method = 'RK45', t_eval = np.geomspace(start = Ni, stop = Ne, num = 10000), vectorized = True,  dense_output=True, args=(k,))

#Take the [0] solutions from the solver for u_k. Solve_ivp gives an array solution with [0] being u_k and [1] being d(u_k)/dN.
    ur_sol = solr.y[0]
    
    #Take last u_k at horizon exit
    ur_cut = ur_sol[-1]
    
    #Repeat for imaginary part.
    ui_sol = soli.y[0]
    
    ui_cut = ui_sol[-1]

#Define curvature perturbation used in the power spectrum.
    R_rs = -(ur_cut/z(Ne)) #unit Mpc(3/2)
   
    R_is = -(ui_cut/z(Ne))
    
#Slow-roll power spectrum to compare with Mukhanov solution.

    s = ((1 - epfunc(Ne)[()] + ((0.0814514 - 3)/8)*sigma(Ne))**2)*(1/(8*np.pi**2))*((H(Ne)**2)/ep(Ne))
    
    
    
#append results to lists.
    ur.append(ur_cut)
    
    ui.append(ui_cut)
    
    R_r.append(R_rs)
    
    R_i.append(R_is)
    
    s_spec.append(s)

######################################################################
#Ni to Ne array from solve_ivp. (solr.t = soli.t)
N = solr.t

#Use Numpy Arrays for easier use later if required.
ur = np.array(ur)

ui = np.array(ui)

s_spec = np.array(s_spec)

R_r = np.array(R_r)

R_i = np.array(R_i)

#|R|^2 for Mukhanov power spectrum.

R = (R_r**2 + R_i**2) #Unit Mpc(3)
    
######################################################################
#Call m_spec and input arrays into Mukhanov power spectrum.
m_spec = m_spec(k_a, R)

#Find the primordial tilt n_s by taking the gradient of the power spectrum and print to terminal.
n_s = 1 + (np.log(m_spec[1])- np.log(m_spec[0]))/(np.log(k_a[1])- np.log(k_a[0]))

print()

print('ns =', n_s)

print()

######################################################################
#Plot the power spectrum.

plt.plot(k_a, m_spec, label = 'P(k)')

plt.plot(k_a, s_spec, label = 'P(k) slow roll')

#Plot Mukhanov variable against e-fold (if required).
#plt.plot(N, R_r, '-b', label = 'Re{R}')

#plt.plot(N, R_i, '-r', label = 'Im{R}')

#plt.plot(N, np.sqrt(R), 'b', label = '|R|')

#Print run time to terminal.
print()

print ("Scalar Code took", time.time() - start_time, "seconds to run")

print()

#Plot parameters.
plt.legend()

plt.title('P(k)')

plt.xlabel('k')
plt.ylabel('P(k)')

plt.yscale('log')
plt.xscale('log')

plt.show()

#Save results to file

np.save('Muk', m_spec)




























