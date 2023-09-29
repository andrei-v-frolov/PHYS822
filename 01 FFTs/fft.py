#!/usr/bin/env python
# numpy 1D FFT demo - calculate Laplacian

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
from numpy.fft import *
from pltconfig import *
import matplotlib


###############################################################################
# compute something
###############################################################################

# number of points in a grid
n = 1024

# initialize x grid [-l,l], EXCLUDING last point
l = 1.0; dx = 2.0*l/n
x = np.linspace(-l, l-dx, num=n)

# compute a function
sigma = 0.1; f = np.exp(-x*x/(2.0*sigma*sigma))

# exact Laplacian
laplacian = (x*x-sigma*sigma)/sigma**4 * f


###############################################################################
# Fast Fourier transform
###############################################################################

# Nyquist frequency (// is integer division, % formatted output)
nn = n//2 + 1; print("Nyquist dimension is n=%i" % nn)
k = np.linspace(0.0, np.pi/dx, num=nn)

# compute real FFT
F = rfft(f)

# print shape of the resulting array, its first and last elements
print(F.shape, F[0], F[-1])

# compute Laplacian via inverse FFT
g = irfft(-k*k*F)


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
