#!/usr/bin/env python
# Bessel function J_a(kr) with real parameters (via build-in ODE solver)

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
from scipy.special import gamma,jv
from pltconfig import *
import matplotlib


###############################################################################
# ODE solver
###############################################################################

# solve Bessel equation r d/dr (r d/dr J) + (kappa^2 r^2 - alpha^2) J = 0

# equation parameters
kappa = 1.0; alpha = 0.0

# split into 1st order ODE via q = J, p = rJ', pack as y=[q,p]
def f(y,r):
	return y[1]/r, -((kappa*r)**2-alpha**2)/r * y[0]

# asymptotic period of oscillations
R = 2.0*np.pi/kappa

# asymptotic at r = 0 is J_alpha(kr) ~ (kr/2)^alpha/Gamma(alpha+1)
epsilon = 1.0e-6; q0 = pow(kappa*epsilon/2,alpha)/gamma(alpha+1); p0 = alpha*q0

r = np.linspace(epsilon,10*R,10*64)
y = odeint(f, [q0,p0], r, atol=1.0e-12)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot J_alpha(kr) from scipy libraries
#plt.plot(r, jv(alpha,r), 'b-', label='', alpha=1.0)

# plot J_alpha(kr) solution in red solid line, opaque
plt.plot(r, y[:,0], 'r-', label='', alpha=1.0)

# plot residual error
#plt.plot(r, y[:,0]-jv(alpha,r), 'r-', label='', alpha=1.0)

plt.xlabel('r'); plt.ylabel(r'$J_\alpha(\kappa r)$')

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
