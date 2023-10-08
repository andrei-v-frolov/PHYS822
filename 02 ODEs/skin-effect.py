#!/usr/bin/env python
# skin effect current profile via build-in Bessel function J_0(kr)

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
from scipy.special import jv
from pltconfig import *
import matplotlib


###############################################################################
# skin effect current profile
###############################################################################

# penetration depth
kappa = np.sqrt(1.0j)

# asymptotic period of oscillations
R = 2.0*np.pi/kappa.real

# current profile for wire of radius R
r = np.linspace(-R,R,10*64)
J = np.vectorize(lambda x: jv(0,kappa*x))(r)

# normalize profile to unit current at the edge
J /= J[-1]


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot J(r) solution
plt.plot(r/R, np.real(J), 'c-', label='', alpha=1.0)
plt.plot(r/R, np.imag(J), 'g-', label='', alpha=1.0)
plt.plot(r/R, np.abs(J),  'r-', label='', alpha=1.0)

# configure labels
plt.xlabel(r'$r/R$'); plt.ylabel(r'$J_0(\kappa r)/J_0(\kappa R)$')
plt.xlim(-1,1)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
