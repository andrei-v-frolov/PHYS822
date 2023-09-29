#!/usr/bin/env python
# numpy and matplotlib 2D plots demo

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

# compute a function
sigma = 0.4; f = np.exp(-(X*X+Y*Y)/(2.0*sigma*sigma))

# exact Laplacian
laplacian = (X*X+Y*Y-2.0*sigma*sigma)/sigma**4 * f


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
plt.imshow(f, cmap='Oranges', extent=[-x0,x0,-y0,y0])
plt.contour(f, cmap='plasma', extent=[-x0,x0,-y0,y0])

# manually set plot area
plt.xlim(-lx,lx); plt.ylim(-ly,ly)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
