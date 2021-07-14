#Work in Planck units where c=hbar = 1.

# Masses in M_pl (reduced planck mass)

# Parametric solution X

# No backward integration...

import matplotlib.pyplot as plt


import numpy as np

from numpy import exp, log, sqrt, sinh, cosh, tanh, arccosh, arcsinh


from scipy.integrate import quad

from scipy.interpolate import interp1d

from input import m_potential

from input import n_potential

from input import v_func

from input import vprime

from input import ti


plt.rcParams["font.size"] = 6

plt.autoscale(enable=True, axis='both', tight=None)
##############################################################

print()

print('Calculating Background Variables for n = ', n_potential)

print()

#Globals

tiny = 1e-6 # when to use taylor expansion for sinh cosh tanh?

h = np.float64(-5e-4 )# Step size for ODE


##############################################################

#H-J ODE

def f(y, t):

   # y =  parametric variable X    

    """H-J equation to be solved using RK45."""

    return sqrt(1.5) - 0.5*vprime(t)/(tanh(y)*v_func(t))

##############################################################

# Choices above result in the following conditions for RK4


Hi = sqrt(v_func(ti)/2.5) # initial H, unit = Mpl, this must be between sqrt(V/3) and sqrt(V/2)



dHi = sqrt(1.5*Hi**2-0.5*v_func(ti))   # initial dH/dphi determined by HJ, unitless (in Planck units)


yi = arccosh(Hi*sqrt(3./v_func(ti)))  # initial X

##############################################################

ep = 2.*(dHi/Hi)**2

#print("Initial epsilon = ", ep)


N = 400.
# Define lists to append on forward integration

tlist_i= [ti]



Hlist_i=[Hi]

dHlist_i=[dHi]



Xlist_i = [yi]

dXlist_i = [f(yi,ti)]



eplist_i = [ep]


Nlist_i = [N]

dNlist_i=[1./sqrt(2.*ep)]


#Error
Slist = []

tol = 1e-6

while (ep<1.):



        # Do RK4

        k1 = h*f(yi, ti)
    

        k2 = h*f(yi + 0.25 * k1, ti + 0.25 * h)


        k3 = h*f(yi + (3./32.)*k1 + (9./32.)*k2, ti + (3./8.) * h)


        k4 = h*f(yi + (1932./2197.)*k1 - (7200./2197.)*k2 + (7296./2197.)*k3, ti + (12./13.)*h)


        k5 = h*f(yi + (439./216.)*k1 - 8.*k2 + (3680./513.)*k3 - (845./4104.)*k4, ti + h)


        k6 = h*f(yi - (8./27.)*k1 + 2.*k2 - (3544./2565.)*k3 + (1859./4104.)*k4 - (11./40.)*k5, ti + 0.5*h)


        # READY FOR NEXT STEP


        ti = ti+h


        yi = yi + (25./216.)*k1 + (1408./2565.)*k3 + (2197./4101.)*k4 - (1./5.)*k5


        zi = yi + (16./135.)*k1 + (6656./12825.)*k3 + (28561./56430.)*k4 - (9./50.)*k5 + (2./55.)*k6
        

        s = ((tol*-h)/(2*np.abs(zi - yi)))**0.25

        #print(s*-h)

        dy=f(yi,ti)

    
        vfun= v_func(ti)


        
        if (yi<tiny): # taylor series 

            ep = 3.*yi**2-2.*yi**4 + 17.*(yi**6)/15.

            H = sqrt(vfun/3.)*(1.+0.5*yi**2+(yi**4)/24.)

            dH = sqrt(vfun/2.)*(yi+(yi**3)/6.+(yi**5)/120.)

           

        else: 

            ep = 3.*(tanh(yi))**2

            H = sqrt(vfun/3.)*cosh(yi)

            dH = sqrt(vfun/2.)*sinh(yi)

         

                
        dN = 1./sqrt(2.*ep)
        
        N = N + h * dN
        
        
        # append


        tlist_i.append(ti)



        Xlist_i.append(yi)
        
        dXlist_i.append(dy)

        

        Hlist_i.append(H)

        dHlist_i.append(dH)



        eplist_i.append(ep)


        Nlist_i.append(N)
        
        dNlist_i.append(dN)


        Slist.append(s)
##############################################################

#End of inflation: report results from while loop.

print("Inflation ends at:  phi = ",  tlist_i[-1], '   H = ' ,  Hlist_i[-1], '   eps= ', eplist_i[-1])


Ntot = -h*sum(dNlist_i) 



if (Ntot<60.):

    print("Not enough inflation: the IC chosen only allowed N=", Ntot)
    
else: #Trim the array to find the correct initial condition

    j=-1

    dNcount =dNlist_i[j]

    dNbound=-70.0/h # total efold desired

    while (dNcount<dNbound):

        j+=-1      # keep going down the array starting from the end

        dNcount += dNlist_i[j]



    tlist_i = tlist_i[j:]  #slicing, keeping the tail of these arrays

    Xlist_i = Xlist_i[j:]

    dXlist_i = dXlist_i[j:]

    Hlist_i = Hlist_i[j:]

    dHlist_i = dHlist_i[j:]

    eplist_i = eplist_i[j:]

    dNlist_i = dNlist_i[j:]

    Nlist_i = Nlist_i[j:]

    Nlist_i = [z - Nlist_i[-1] for z in Nlist_i]

    print()

    #print("Found a good initial condition at phi_i= ", tlist_i[0], '   H = ' ,  Hlist_i[0], '   eps= ', eplist_i[0])

    #print()

    #print("Check: Ntot =", -h*sum(dNlist_i))

    #print()
##############################################################

# Tensor to scalar ratio and primordial tilt at N = 60

r = 16. * eplist_i[0]

ns = 1 - 6*eplist_i[0]

nt = -2. * eplist_i[0]


##############################################################

# Print Output.

print('The tensor to scalar ratio (slow roll approx) =' ,r)

print()

print('n_s (slow roll approx) =' ,ns)

print()

print('n_t (slow roll approx) =' ,nt)

print()

##############################################################

# Plot for test.

plt.title('$\epsilon(\phi)$ v.s N for n = ' +str(n_potential))


plt.plot(Nlist_i, eplist_i, '-r')


plt.xlabel('N')


plt.ylabel('$\epsilon$')

plt.show()

#Make N increase to use spline
N_arr = np.array(Nlist_i[::-1])
H_arr = np.array(Hlist_i[::-1])
Hp_arr = np.array(dHlist_i[::-1])
ep_arr = np.array(eplist_i[::-1])
##############################################################
#Output

Background_Output = [H_arr, Hp_arr, ep_arr, N_arr]
