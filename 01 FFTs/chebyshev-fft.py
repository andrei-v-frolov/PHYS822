#!/usr/bin/env python
# solve 1D Poisson equation in infinite interval using rational Chebyshev basis

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
n = 2048

# characteristic scale
l = 0.4

# initialize theta grid [-pi,pi], EXCLUDING end points
dtheta = 2.0*np.pi/n; theta = np.linspace(-np.pi+dtheta/2, np.pi-dtheta/2, num=n)

# initialize Chebyshev spectral grid, and the real space one
# they cover the entire space twice, but it makes FFTs simple
y = np.cos(theta); x = l/np.tan(theta)

# define a function extracting first half of array
def half(array):
	return np.split(array,2)[0]

# compute a function
sigma = 0.1; f = np.exp(-x*x/(2.0*sigma*sigma))

# exact Laplacian
laplacian = (x*x-sigma*sigma)/sigma**4 * f


###############################################################################
# Fast Fourier transform grid
###############################################################################

# wave numbers, note that grid INCLUDES k=0
nn = n//2 + 1; k = np.linspace(0.0, np.pi/dtheta, num=nn)

###############################################################################
# compute Laplacian of a function
###############################################################################

# compute first and second derivative via FFT
#F = rfft(f); df = irfft((-1.0j*l)*k*F)/(l*l+x*x)
#G = rfft(df); g = irfft((-1.0j*l)*k*G)/(l*l+x*x)

###############################################################################
# solve 1D Poisson equation (d/dx)^2 g = -f via FFT
# note that it has linear potential asymptotic at infinity!
# strategy to deal with it: bravely ignore zero mode
###############################################################################

#k[0] = np.finfo('d').max # null zero frequency mode
#F = rfft((l*l+x*x)*f); df = irfft((1.0j/l)*F/k)
#G = rfft((l*l+x*x)*df); g = irfft((1.0j/l)*G/k)

###############################################################################
# solve spherical Poisson equation d^2u/dx^2 = -x*f, g = u/x via FFT
# this is reduction of spherical geometry, has 1/|x} asymptotic
###############################################################################

k[0] = np.finfo('d').max # null zero frequency mode
F = rfft((l*l+x*x)*f*x); df = irfft((1.0j/l)*F/k)
G = rfft((l*l+x*x)*df); g = irfft((1.0j/l)*G/k)/x

# optionally shift potential to zero at infinity
g = -g + np.max(g)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
plt.plot(half(x), half(g), 'r-', label='', alpha=1.0)

plt.xlim(-15*l,15*l)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
