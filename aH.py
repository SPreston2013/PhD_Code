import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate 
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from Background import Background_Output

np.set_printoptions(threshold=np.inf)

plt.rcParams["font.size"] = 7

plt.autoscale(enable=True, axis='both', tight=None)
######################################################################
h = np.float64(1e-3)

######################################################################
#Take result array from background code.

H_a = np.array(Background_Output[0])

ep_a = np.array(Background_Output[2])

N_a = np.array(Background_Output[3])

######################################################################
#k_* to set initial conditions  (k= aH @ N = 60.0)
k = 0.05

######################################################################
#Interpolate epsilon and H with the e-fold array from background code.
epfunc = interp1d(N_a, ep_a, bounds_error = False)

Hfunc =  interp1d(N_a, H_a, bounds_error = False)

######################################################################
#Define the ODE dln(aH)/dN = epsilon - 1.
def aH(y, t):

    return (epfunc(t)[()] - 1)

######################################################################
#Initial conditions set.
yi = np.log(k)


ti = 60.0
######################################################################
#Create lists with initial conditions as the first element.

tlist = [ti]


ylist = [yi]


zlist = [yi]
######################################################################
#Start RK4 run (forward integrate).

for number in range(len(N_a) - 1):

    k1 = h*aH(yi, ti)
    
    
    k2 = h*aH(yi + 0.25 * k1, ti + 0.25 * h)
    

    k3 = h*aH(yi + (3./32.)*k1 + (9./32.)*k2, ti + (3./8.) * h)
   

    k4 = h*aH(yi + (1932./2197.)*k1 - (7200./2197.)*k2 + (7296./2197.)*k3, ti + (12./13.)*h)
    

    k5 = h*aH(yi + (439./216.)*k1 - 8.*k2 + (3680./513.)*k3 - (845./4104.)*k4, ti + h)
    

    k6 = h*aH(yi - (8./27.)*k1 + 2.*k2 - (3544./2565.)*k3 + (1859./4104.)*k4 - (11./40.)*k5, ti + 0.5*h)
    
        # READY FOR NEXT STEP


    ti = ti + h
    

    yi = yi + (25./216.)*k1 + (1408./2565.)*k3 + (2197./4101.)*k4 - (1./5.)*k5
    
    
    zi = yi + (16./135.)*k1 + (6656./12825.)*k3 + (28561./56430.)*k4 - (9./50.)*k5 + (2./55.)*k6
   

    tlist.append(ti)
    

    ylist.append(yi)
    

    zlist.append(zi)
   


kaH = np.array(np.log(k) - zlist)

j = [i for i in kaH if i <= np.log(100.)]

N_i_list = tlist[:len(j)]

y_i_list = ylist[:len(j)]

Ni = N_i_list[-1]


N_i = N_i_list[-1]


y_i = y_i_list[-1]


Nlist = [N_i]


ylist2 = [y_i]


zlist2 = [y_i]

h2 = np.float64(-15e-4)

for number in range(len(N_a) - 1):

    k_1 = h2*aH(y_i, N_i)
    
    
    k_2 = h2*aH(y_i + 0.25 * k_1, N_i + 0.25 * h2)
    

    k_3 = h2*aH(y_i + (3./32.)*k_1 + (9./32.)*k_2, N_i + (3./8.) * h2)
   

    k_4 = h2*aH(y_i + (1932./2197.)*k_1 - (7200./2197.)*k_2 + (7296./2197.)*k_3, N_i + (12./13.)*h2)
    

    k_5 = h2*aH(y_i + (439./216.)*k_1 - 8.*k2 + (3680./513.)*k_3 - (845./4104.)*k_4, N_i + h2)
    

    k_6 = h2*aH(y_i - (8./27.)*k_1 + 2.*k_2 - (3544./2565.)*k_3 + (1859./4104.)*k_4 - (11./40.)*k_5, N_i + 0.5*h2)
    
        # READY FOR NEXT STEP


    N_i = N_i + h2
    

    y_i = y_i + (25./216.)*k_1 + (1408./2565.)*k_3 + (2197./4101.)*k_4 - (1./5.)*k_5
    
    
    z_i = y_i + (16./135.)*k_1 + (6656./12825.)*k_3 + (28561./56430.)*k_4 - (9./50.)*k_5 + (2./55.)*k_6
   

    Nlist.append(N_i)
    

    ylist2.append(y_i)
    

    zlist2.append(z_i)


zlist2 = np.array(zlist2)


aH = np.exp(zlist2[::-1])


aHt = np.array(Nlist[::-1])



print()

print('N_i is:', Ni)

print()


#plt.plot(aHt, aH, '-r')
#plt.xlabel('N')
#plt.ylabel('aH')
#plt.yscale('log')
#plt.show()
