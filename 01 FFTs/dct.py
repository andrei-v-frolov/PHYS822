#!/usr/bin/env python
# scipy 1D DST demo - calculate Laplacian

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
n = 1024

# initialize x grid [-l,l], EXCLUDING end points
l = 1.0; dx = 2.0*l/n; x0 = l - dx/2.0
x = np.linspace(-x0, x0, num=n)

# compute a function
sigma = 0.1; f = np.exp(-x*x/(2.0*sigma*sigma))

# exact Laplacian
laplacian = (x*x-sigma*sigma)/sigma**4 * f


###############################################################################
# Fast Fourier transform
###############################################################################

# wave numbers, note that grid INCLUDES k=0
k = (np.pi/(2.0*l)) * np.linspace(0,n-1, num=n)

# compute Type II DCT 
F = dct(f)

# print shape of the resulting array, its first and last elements
print(F.shape, F[0], F[-1])

# compute Laplacian via inverse FFT
g = idct(-k*k*F)


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
plt.plot(x, g-laplacian, 'r-', label='', alpha=1.0)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
