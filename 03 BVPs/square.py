#!/usr/bin/env python
# rectangular waveguide TM modes (analytic solution)

###############################################################################
# configure the environment
###############################################################################

from os import environ

# show figure in macOS window
#environ['matplotlib.backend'] = 'macosx'

# output figure as PDF file
environ['matplotlib.backend'] = 'pdf'


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
nx = 512; ny = 256
kx = 3; ky = 2

# initialize 1D x & y grids
lx = 2.0; dx = lx/nx
ly = 1.0; dy = ly/ny

x = np.linspace(dx/2.0, lx - dx/2.0, num=nx)
y = np.linspace(dy/2.0, ly - dy/2.0, num=ny)

# 2D mesh LUTs of x and y coordinates
X,Y = np.meshgrid(x, y)

# compute the mode function
f = np.sin(kx*np.pi*X/lx)*np.sin(ky*np.pi*Y/ly)

# symmetric color scale
z = max(abs(np.max(f)),abs(np.min(f)))
levels = np.linspace(-z,z,100)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
plt.imshow(f, cmap='jet', extent=[0,lx,0,ly], vmin=-z, vmax=z)
#plt.contourf(f, cmap='jet', extent=[0,lx,0,ly], levels=levels)

# manually set plot area
plt.title(r"$\gamma = %.6f$" % np.sqrt((np.pi*kx/lx)**2 + (np.pi*ky/ly)**2))
plt.xlim(0,lx); plt.ylim(0,ly)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
