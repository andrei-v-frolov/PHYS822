#!/usr/bin/env python
# numpy and matplotlib 2D polar plot demo

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
from scipy.fft import *
from pltconfig import *
import matplotlib


###############################################################################
# compute something
###############################################################################

# number of points in a grid
n = 256; m = 512

# initialize 1D rho & phi grids
r = 2.0; twopi = 2.0*np.pi
drho = r/n; dphi = twopi/m

rho = np.linspace(drho/2.0, r-drho/2.0, num=n)
phi = np.linspace(dphi/2.0, twopi-dphi/2.0, num=m)

# 2D mesh LUTs of rho and phi coordinates
Phi,Rho = np.meshgrid(phi, rho)

# compute some function
f = np.exp(Rho)*np.cos(2.0*Phi)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# set projection to polar coordinates
plt.subplot(111, polar=True)

# plot f(x) in plasma colormap, opaque
plt.contourf(phi, rho, f, cmap='jet', extent=[0,twopi,0,r], levels=50)
plt.colorbar()

# manually set plot area
plt.yticks([])

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
