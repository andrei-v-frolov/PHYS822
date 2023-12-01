#!/usr/bin/env python
# cylindrical waveguide TM modes (spectral solver for radial eigenfunctions)

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
from scipy.linalg import *
from pltconfig import *
import matplotlib


###############################################################################
# spectral cylindrical grid
###############################################################################

# number of points in a grid
n = 256; m = 512

# mode parameters and index
r = 1.0; mode = 4; l = 2

# initialize 1D rho & phi grids
twopi = 2.0*np.pi; halfpi = np.pi/2.0
dtheta = halfpi/n; dphi = twopi/m

theta = np.linspace(halfpi-dtheta/2.0,dtheta/2.0, num=n)
rho = r*np.cos(theta); phi = np.linspace(0.0, twopi, num=m)


###############################################################################
# basis functions and spectral operators
###############################################################################

def basis(k):
	return np.cos(theta)**l*np.sin(theta)**2*np.cos(2*k*theta)

def laplacian(k):
	return -4/(r*r)*(
		(k*k+l+1)*np.cos(2*k*theta) +
		k*(3.0-2*(l+2)*np.sin(theta)**2)*np.sin(2*k*theta)/np.sin(2*theta)
		)*np.cos(theta)**l


###############################################################################
# generalized eigenvalue problem
###############################################################################

L = np.zeros([n,n])
B = np.zeros([n,n])

for k in range(0,n):
	B[:,k] = basis(k)
	L[:,k] = laplacian(k)

# solve for radial eigenmodes
w,v = eig(L,B)

# index by eigenvalue
idx = np.argsort(-w.real)

# output eigenvalues and residuals
#for i in idx:
	#print(-w[i].real,w[i].imag, norm(np.matmul(L,v[:,i]) - w[i]*np.matmul(B,v[:,i])))

# 2D mesh LUTs
gamma = -w[idx[mode]].real
f = np.matmul(B,v[:,idx[mode]])
Phi,F = np.meshgrid(phi, f)

# symmetric color scale
z = max(abs(np.max(f)),abs(np.min(f)))
levels = np.linspace(-z,z,100)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# set projection to polar coordinates
plt.subplot(111, polar=True)

# plot f(x) in plasma colormap, opaque
plt.contourf(phi, rho, F*np.cos(l*Phi), cmap='jet', extent=[0,twopi,0,r], levels=levels)
#plt.colorbar()

# manually set plot annotations
plt.title(r"$\gamma = %.6f$" % np.sqrt(gamma))
plt.xticks([])
plt.yticks([])

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
