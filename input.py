import numpy as np

from numpy import exp, log, sqrt, sinh, cosh, tanh, arccosh, arcsinh


from scipy.integrate import quad

from scipy.interpolate import interp1d


#Globals
alpha = 1.

m_potential = 4.44e-6

n_potential = 1.

##############################################################

# I.C. Choose initial H and initial phi


ti = 100. #0.5/m_potential #Phi_initial. [unit = Mpl].. For r<<1, need phi_i to be roughly ~1/m

#Define Potential

def v_func(t):

    """Potential for a given inflation model."""

    #return 9*(alpha**2)*m_potential*(tanh(t/sqrt(6*alpha)))**4 #Alpha attractor (ti approx 10)
    
    
    #return (((3*m_potential**2)/4)*(1 - np.exp(-np.sqrt(2/3)*t))**2) #R + R^2 (ti approx 25)
    
    
   return np.float64(0.5 * m_potential**3 * t**n_potential) # unit = Mpl^4 (power law potential) (ti approx 100)



def vprime(t):

    """Potential derivative for a given inflation model."""

    #return (6*sqrt(6)*(alpha**2)*m_potential*(tanh(t/sqrt(6*alpha)))**3 * (1./cosh(t/sqrt(6*alpha)))**2)/sqrt(alpha)
    
    
    #return (np.sqrt((3*m_potential**2)/4)*np.exp(-np.sqrt(2/3)*t)*(1 - np.exp(-np.sqrt(2/3)*t))**2) # R + R^2
    
    
    return np.float64(0.5 * m_potential**3 * n_potential*t**(n_potential-1)) # unit = Mpl^4 (power law potential)


