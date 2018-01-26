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
    n, latitude, longitude, rate = pickle.load(f)

# Compute the topography
h = numpy.zeros(rate.shape)
for i, la in enumerate(latitude):
    for j, lo in enumerate(longitude):
        h[i, j] = topo.ground_altitude(la, lo, geodetic=True)

# Show the map
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

plt.figure()
plt.contour(longitude, latitude, h, 10, colors="k")
plt.pcolor(longitude, latitude, rate, cmap="nipy_spectral", norm=LogNorm(),
           alpha=0.5)
plt.colorbar()
plt.xlabel(r"longitude (deg)")
plt.ylabel(r"latitude (deg)")
plt.title(r"trigger rate (Hz / deg$^2$)")
plt.show()
