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
from Muk import points
from Muk import k_a
import time
start_time = time.time()

#Error check point.
print('Running Tensor Perturbation Code...')

plt.rcParams["font.size"] = 6
plt.autoscale(enable=True, axis='both', tight=None)
######################################################################
#H', H'', H'''0 plus H and epsilon from background code.
H_a = np.array(Background_Output[0])

Hp_a = np.array(Background_Output[1])

ep_a = np.array(Background_Output[2])

N_a = np.array(Background_Output[3])
######################################################################
#Spline epsilon from the background code. 
epfunc = UnivariateSpline(N_a, ep_a, s=0, k=5)

#Spline H and H' from the background code, then take the derivatives.
Hfunc =  UnivariateSpline(N_a, H_a, s=0,k=5)

#Spline aH from the Vacc code. 
aHfunc = UnivariateSpline(aHt, aH, s=0, k=5)
######################################################################
def H(t):
    """Hubble parameter"""
    return np.array(Hfunc(t)[()])

def aH(t):
    
    return np.array(aHfunc(t)[()])
    
def a(t):
    """Scale factor"""
    return np.array(aHfunc(t)[()])/(np.array(Hfunc(t)[()]))

#Define Slow roll parameters.
def ep(t):
    return np.array(epfunc(t)[()])

######################################################################
#Mukhanov equation in terms of the e-fold number.
def sol(t, Y):
    """Tensor equation"""
    return np.array((Y[1], (1 - ep(t))*Y[1] -((k/aH(t))**2 - (2 - ep(t)))*Y[0]))

######################################################################

#Tensor Power spectrum
def h_spec(k, He, h):
    """Tensor power spectrum P(t)"""
    return np.array(((2*k**3/np.pi**2) * h))
    
#Define slow-roll power spectrum.
def SlowRoll_h_Spec(H):
    """Slow-Roll power spectrum"""
    return 14*((H**2)/(8*np.pi**2))
######################################################################

#Create empty lists.
vr = []

vi = []

He = []

aHe = []

for k in k_a:
#efold to stop the code, runs from Ni to Ne
    Ni = 60.00000001 - np.log(k/k_a[0])
    
    Ne = np.array(60.0 - np.log(k/k_a[0]))
    
#u_k real I.C
    vr_0 = np.array([1/(np.sqrt(2*k)), 0])
#u_k imaginary I.C
    vi_0 = np.array([0, -np.sqrt(k/2)/aHfunc(Ni)])

#Use solve_ivp with the RK45 method, solr is the real part and soli is the Imaginary. 
    solr = solve_ivp(sol, (Ni, Ne), vr_0, method = 'RK45', t_eval = np.linspace(start = Ni, stop = Ne, num = 1000))
        
    soli = solve_ivp(sol, (Ni, Ne), vi_0, method = 'RK45', t_eval = np.linspace(start = Ni, stop = Ne, num = 1000))
    
    
    H_e = H(Ne)
    
    aH_e = aH(Ne)

#Take the solutions from the solver.
    vrs = solr.y[0]
#Take value at horrizon exit.
    vr_cut = vrs[-1]
    

#Repeat prior step for imaginary part.
    vis = soli.y[0]

    vi_cut = vis[-1]
   
    
#append lists 
    vr.append(vr_cut)
    
    vi.append(vi_cut)
    
    He.append(H_e)
    
    aHe.append(aH_e)
   


N = solr.t

#Numpy Array
vr = np.array(vr)

vi = np.array(vi)

He = np.array(He)

aHe = np.array(aHe)


#Define the tensor perturbation h
hr = 2*(vr*He/aHe)

hi = 2*(vi*He/aHe)

h = hr**2 + hi**2
######################################################################

#Input arrays
h_spec = h_spec(k_a, He, h)

Sh_spec = SlowRoll_h_Spec(He)

#Print initial value of tensor power spectrum to terminal.
print()

print('n_t =', (np.log(h_spec[1]) - np.log(h_spec[0]))/(np.log(k_a[1]) - np.log(k_a[0])))

print()

print ("Tensor Code took", time.time() - start_time, "seconds to run")

######################################################################
#Plot the power spectrum.
plt.plot(k_a, h_spec, label = 'P(t)')

plt.plot(k_a, Sh_spec, label = 'Slow-Roll P(t)')

#Plot Tensor var.
#plt.plot(N, hr, '-b', label = 'Re{R}')
#plt.plot(N, hi, '-r', label = 'Im{R}')
#plt.plot(N, np.sqrt(hr**2 + hi**2), 'b', label = '|h|')
#Plot parameters.
plt.legend()

plt.title('p(t)')

plt.xlabel('k')
plt.ylabel('p(t)')

plt.yscale('log')
plt.xscale('log')

plt.show()

np.save('Ten', h_spec)
