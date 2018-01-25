#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from grand_tour import Topography


# Get the topography handle
topo = Topography(latitude=42.1, longitude=86.3, path="share/topography",
                  stack_size=49)

# Load the trigger map
with open("share/events.triggers.p", "rb") as f:
    n, rate_map = pickle.load(f)
rate = rate_map[4] / n
x = numpy.arange(0., rate.shape[1] * rate_map[1], rate_map[1]) + rate_map[0]
y = numpy.arange(0., rate.shape[0] * rate_map[3], rate_map[3]) + rate_map[2]

# Compute the topography
h = numpy.zeros(rate.shape)
for i, yi in enumerate(y):
    for j, xj in enumerate(x):
        h[i, j] = topo.ground_altitude(xj, yi)

# Show the map
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

plt.figure()
# plt.contourf(x * 1E-03, y * 1E-03, h * 1E-03, 100, cmap="terrain")
# plt.pcolor(x, y, rate, norm=LogNorm())
plt.pcolor(x * 1E-03, y * 1E-03, h * 1E-03, cmap="terrain")
plt.colorbar()
plt.show()
