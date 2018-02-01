#!/usr/bin/env python
import json

import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from grand_tour import Topography

# Get the topography
topo = Topography(latitude=42.1, longitude=86.3, path="share/topography",
                  stack_size=49)

# Get antennas
with open("share/setup/hotspot.json", "rb") as f:
    antennas = json.load(f)

def axisEqual3D(ax):
    extents = numpy.array([getattr(ax, 'get_{}lim'.format(dim))()
                           for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = numpy.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)

# Loop over antennas.
numpy.random.shuffle(antennas)
n_figs = 10
for ifig, antenna in enumerate(antennas):
    r0 = antenna[:3]
    theta, phi = map(numpy.deg2rad, antenna[3:])
    s = numpy.sin(theta)
    n = numpy.array((s * numpy.cos(phi), s * numpy.sin(phi), numpy.cos(theta)))
    r1 = r0 - 4.5 * n

    # Sample the ground around the antenna
    x = numpy.linspace(-100., 100., 21) + r0[0]
    y = numpy.linspace(-100., 100., 21) + r0[1]
    zg = numpy.zeros((len(y), len(x)))
    for i, yi in enumerate(y):
        for j, xj in enumerate(x):
            zg[i, j] = topo.ground_altitude(xj, yi)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = numpy.meshgrid(x, y)
    ax.plot_wireframe(X, Y, zg)
    ax.scatter(*r0, s=100, marker="o", facecolors="k", edgecolors="k")
    ax.scatter(*r1, s=100, marker="o", facecolors="k", edgecolors="k")
    ax.plot(*zip(r0, r1), color="k")
    axisEqual3D(ax)
    if ifig >= n_figs - 1:
        break

plt.show()
