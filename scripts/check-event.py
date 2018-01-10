#!/usr/bin/env python

import os
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from retro.event import EventIterator
from grand_tour import Topography

def make_rotation(direction, angle):
    """
    Create a rotation matrix corresponding to the rotation around a general
    axis by a specified angle.

    R = dd^T + cos(a) (I - dd^T) + sin(a) skew(d)

    Parameters:

        angle : float a
        direction : array d
    """
    d = numpy.array(direction, dtype=numpy.float64)
    d /= numpy.linalg.norm(d)

    eye = numpy.eye(3, dtype=numpy.float64)
    ddt = numpy.outer(d, d)
    skew = numpy.array([[    0,  d[2],  -d[1]],
                        [-d[2],     0,  d[0]],
                        [d[1], -d[0],    0]], dtype=numpy.float64)

    mtx = ddt + numpy.cos(angle) * (eye - ddt) + numpy.sin(angle) * skew
    return mtx

# Load the event
tag = "E.1e17_X.93815_Y.74409_Z.-559_T.91_P.21_D.3750214666968983.more-antennas"
topo = Topography(43, 87, "flat/10")
event = EventIterator("share/events/{:}.json".format(tag)).next()
antennas = numpy.array(event["antennas"])
energy, position, direction, altitude = map(numpy.array, event["tau_at_decay"])

gamma = numpy.deg2rad(3.)
zcmin = 14E+03
zcmax = 165E+03 * energy / 1E+09 + 55E+03
rmin = position + zcmin * direction

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(antennas[:,0], antennas[:,1], antennas[:,2], "ko")
ax.plot((position[0],), (position[1],), (position[2],), "ko")

r0 = position
r1 = position + zcmax * direction
r = numpy.vstack((r0, r1))
ax.plot(r[:,0], r[:,1], r[:,2], "r-")

u = numpy.cross(direction, (0., 0., 1.))
u /= numpy.linalg.norm(u)
R = make_rotation(u, gamma)
u = numpy.dot(R, direction)

for phi in numpy.linspace(numpy.deg2rad(-90.), numpy.deg2rad(90.), 60):
    R = make_rotation(direction, phi)
    v = numpy.dot(R, u)
    r0 = position + zcmin * v
    r1 = position + zcmax * v
    r = numpy.vstack((r0, r1))
    ax.plot(r[:,0], r[:,1], r[:,2], "g-")

plt.show()
