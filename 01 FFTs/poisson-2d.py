#!/usr/bin/env python
# solve 2D Poisson equation, with zero boundary potential

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
sigma = 0.1; f = np.exp(-(X*X+Y*Y)/(2.0*sigma*sigma))

# exact Laplacian
laplacian = (X*X+Y*Y-2.0*sigma*sigma)/sigma**4 * f

# print computed array shape (note index is in y,x order!)
print(f.shape)


###############################################################################
# Fast Fourier transform
###############################################################################

# wave numbers, note that grid SKIPS k=0
kx = (np.pi/(2.0*lx)) * np.linspace(1,nx, num=nx)
ky = (np.pi/(2.0*ly)) * np.linspace(1,ny, num=ny)

# 2D mesh LUTs of x and y wave numbers
Kx,Ky = np.meshgrid(kx, ky)

# compute Laplacian via Type II DST
# F = dst(dst(f).T).T
# L = -(Kx*Kx+Ky*Ky) * F
# g = idst(idst(L.T).T)

# solve Poisson equation with zero potential boundaries
rho  = 1.5*np.exp(-((X+1.0)**2+(Y+0.7)**2)/(2.0*sigma*sigma))
rho -= 1.0*np.exp(-((X-1.6)**2+(Y-0.3)**2)/(4.0*sigma*sigma))
rho -= 0.5*np.exp(-((X+0.2)**2+(Y-0.6)**2)/(8.0*sigma*sigma))

Rho = dst(dst(rho).T).T
Phi = Rho/(Kx*Kx+Ky*Ky)
g = idst(idst(Phi.T).T)

# peak potential value (either min or max)
z0 = max(abs(np.max(g)),abs(np.min(g)))


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
#plt.imshow(rho, cmap='seismic', vmin=-1.5, vmax=1.5, extent=[-x0,x0,-y0,y0], origin='lower')
plt.imshow(g, cmap='seismic', vmin=-z0, vmax=z0, extent=[-x0,x0,-y0,y0], origin='lower')
plt.contour(g, cmap='coolwarm', levels=20, extent=[-x0,x0,-y0,y0])

# manually set plot area
plt.xlim(-lx,lx); plt.ylim(-ly,ly)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
