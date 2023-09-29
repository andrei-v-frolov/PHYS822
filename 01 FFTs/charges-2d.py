#!/usr/bin/env python
# superposition of point charges in 2D free space

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
nx = 1024; ny = 512

# initialize 1D x & y grids
lx = 2.0; dx = 2.0*lx/nx; x0 = lx - dx/2.0
ly = 1.0; dy = 2.0*ly/ny; y0 = ly - dy/2.0

x = np.linspace(-x0, x0, num=nx)
y = np.linspace(-y0, y0, num=ny)

# 2D mesh LUTs of x and y coordinates
X,Y = np.meshgrid(x, y)

###############################################################################
# superposition of point charges
###############################################################################

q1 = 3.0; q2 = 4.0; q3 = 4.0

g = -12.0
g -= q1*np.log((X+1.0)**2+(Y+0.7)**2)
g += q2*np.log((X-1.6)**2+(Y-0.3)**2)
g += q3*np.log((X+0.2)**2+(Y-0.6)**2)

z0 = 30.0

###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
#plt.imshow(rho, cmap='seismic', vmin=-1.5, vmax=1.5, extent=[-x0,x0,-y0,y0], origin='lower')
plt.imshow(g, cmap='seismic', vmin=-z0, vmax=z0, extent=[-x0,x0,-y0,y0], origin='lower')
plt.contour(g, cmap='coolwarm', vmin=-z0, vmax=z0, levels=40, extent=[-x0,x0,-y0,y0])

# manually set plot area
plt.xlim(-lx,lx); plt.ylim(-ly,ly)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
