#!/usr/bin/env python
# Floquet stability analysis of Mathieu equation

###############################################################################
# configure the environment
###############################################################################

from os import environ

# show figure in macOS window
environ['matplotlib.backend'] = 'macosx'

# output figure as PDF file
#environ['matplotlib.backend'] = 'pdf'


###############################################################################
# import libraries
###############################################################################

import numpy as np
from scipy.integrate import odeint
from pltconfig import *
import matplotlib


###############################################################################
# Mathieu equation solver
###############################################################################

# Mathieu equation (for two solutions, e.g. principal matrix)
def f(y,t,a,q):
	return y[1],-(a-2.0*q*np.cos(2.0*t))*y[0],y[3],-(a-2.0*q*np.cos(2.0*t))*y[2]

# external force period
T = np.pi; t = np.linspace(0,T,64)

# calculate Floquet exponent
def mu(a, q):
	y = odeint(lambda y,t: f(y,t,a,q), [1.0,0.0,0.0,1.0], t, atol=1.0e-12)
	z = np.arccosh((y[-1,0]+y[-1,3])/2.0 + 0j)/T; return z.real

# test that it works
print(mu(0.0,0.1))


###############################################################################
# scan parameter space for stability
###############################################################################

# scan resolution
n = 128

# scan range
a0 = 5.0; q0 = 10.0; z = 1.0
#a0 = 0.3; q0 = 1.0; z = 0.1

# 2D mesh of parameters
a = np.linspace(-a0,a0,n)
q = np.linspace(0.0,q0,n)
Q,A = np.meshgrid(q, a)

# compute exponents
floquet = np.vectorize(mu)(A,Q)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot stability diagram
plt.imshow(floquet, cmap='Reds', vmin=0.0, vmax=z, extent=[0,q0,-a0,a0], origin='lower')
#plt.imshow(floquet, cmap='Blues', vmin=0.0, vmax=z, extent=[0,q0,-a0,a0], origin='upper', alpha=0.5)

# plot driving curve
#plt.plot(q,0.3*q, 'g-')

plt.xlim(0,q0); plt.ylim(-a0,a0)
plt.xlabel('q'); plt.ylabel('a')

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
