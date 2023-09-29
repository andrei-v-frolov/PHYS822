#!/usr/bin/env python
# integrate Mathieu equation [https://en.wikipedia.org/wiki/Mathieu_function]

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

a = 0.0; q = 0.1

# Mathieu equation (canonical parametrization)
def f(y,t):
	return y[1],-(a-2.0*q*np.cos(2.0*t))*y[0]

# external force period
T = np.pi

t = np.linspace(0,100*T,100*64)
y = odeint(f, [1.0,0.0], t, atol=1.0e-12)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot x(t) in red solid line, opaque
plt.plot(t, y[:,0], 'r-', label='', alpha=1.0)

plt.xlabel('t'); plt.ylabel('x')

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
