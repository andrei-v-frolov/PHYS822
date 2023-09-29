#!/usr/bin/env python
# build-in ODE solver demo - anharmonic oscillator

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
# ODE solver
###############################################################################

# anharmonic oscillator L = v^2/2 - x^4/4
def f(y,t):
	return y[1],-y[0]**3

# period of anharmonic oscillator
T = 7.416298709205487673735401388781040185

t = np.linspace(0,100*T,100*64)
y = odeint(f, [1.0,0.0], t, atol=1.0e-12)

# energy drift
e = y[:,0]**4/4.0 + y[:,1]**2/2.0 - 0.25


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot x(t) in red solid line, opaque
#plt.plot(t, y[:,0], 'r-', label='', alpha=1.0)

# plot energy drift in red solid line, opaque
plt.plot(t, e, 'r-', label='', alpha=1.0)

plt.xlabel('t'); plt.ylabel(r'$\delta E$')

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
