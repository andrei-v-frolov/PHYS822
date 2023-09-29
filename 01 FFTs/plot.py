#!/usr/bin/env python
# numpy and matplotlib plots demo

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

# initialize x grid [-l,l]
l = 1.0; dx = 2.0*l/(n-1)
x = np.linspace(-l, l, num=n)

# compute a function
sigma = 0.1; f = np.exp(-x*x/(2.0*sigma*sigma))


###############################################################################
# make a figure
###############################################################################

# configure figure size
fig = plt.figure(figsize=(6, 4), frameon=False)

# plot f(x) in red solid line, opaque
plt.plot(x, f, 'r-', label='', alpha=1.0)

# tighten the whitespace (optional)
plt.tight_layout()

if environ['matplotlib.backend'] == 'macosx':
	plt.show()
else:
	plt.savefig("output.pdf", bbox_inches='tight', pad_inches=0.02, transparent=True)
