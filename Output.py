#Output Code run to plot the scalar and tensor power spectrum.
from Background import Background_Output
from Muk import m_spec
from Tensor import k_a, h_spec
import matplotlib.pyplot as plt
import numpy as np

######################################################################
#Inport arrays from background code.
H_a = Background_Output[0]

Hp_a = Background_Output[1]

ep_a = Background_Output[2]

N_a = Background_Output[3]

######################################################################
#Plot all the results.
fig, axs = plt.subplots(2, 2)

#Plot H Vs e-fold
axs[0, 0].plot(N_a, H_a)
axs[0, 0].set_title('H(N) Vs e-fold')
axs[0, 0].set(xlabel='N', ylabel='H')

#Plot Epsilon Vs e-fold.
axs[0, 1].plot(N_a, ep_a, 'tab:orange')
axs[0, 1].set_title('$\epsilon$ Vs e-fold')
axs[0, 1].set(xlabel='N', ylabel='$\epsilon$')

#Plot scalar power spectrum Vs wavenumber (k).
axs[1, 0].plot(k_a, m_spec, 'tab:green')
axs[1, 0].set_title('P(k) Vs Wavenumber (k)')
axs[1, 0].set(xlabel='k', ylabel='P(k)')
axs[1, 0].set_xscale('log')
axs[1, 0].set_yscale('log')

#Plot tensor power spectrum Vs wavenumber (k).
axs[1, 1].plot(k_a, h_spec, 'tab:red')
axs[1, 1].set_title('$P_h(k)$ Vs Wavenumber (k)')
axs[1, 1].set(xlabel='k', ylabel='$P_h(k)$')
axs[1, 1].set_xscale('log')
axs[1, 1].set_yscale('log')
######################################################################
#Print values of scalar and tensor power spectrum at pivot scale (k_*) to terminal.
print()

print('A_s = ', m_spec[0])

print()

print('A_t = ', h_spec[0])

print()

print('Tensor to scalar ratio is:', h_spec[0]/m_spec[0])

print()
######################################################################

plt.tight_layout()

plt.show()
