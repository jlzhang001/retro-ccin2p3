#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time

import numpy
from retro.event import EventIterator

from grand_tour import Topography


# Get the global topography handle
topo = Topography(latitude=42.1, longitude=86.3, path="share/topography",
                  stack_size=49)

def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def add_triggers(event, latitude, longitude, rate):
    """Extract the antenna trigger rates from an event"""
    year = 365.25 * 24. * 60. * 60.
    w = [v[0] * primary_flux(v[1]) for v in event["primaries"]]
    w = sum(w) * year / event["statistics"][1]
    n = event["statistics"][0]
    la, lo = [], []
    if w > 0.:
        for x, y, z, _, _ in event["antennas"]:
            lla = topo.local_to_lla((x, y, z))
            la.append(lla[0])
            lo.append(lla[1])
        if la:
            p, _, _ = numpy.histogram2d(lo, la, (longitude, latitude))
            rate += w * p
    return n

def antennas_density(latitude, longitude):
    """Compute the density of antennas"""
    Dx, Dy, s = 0.5 * 66.5E+03, 0.5 * 150.4E+03, 500.
    x = numpy.arange(-Dx, Dx + s, s)
    y = numpy.arange(-Dy, Dy + s, s)
    la, lo = [], []
    for i, yi in enumerate(y):
        for j, xj in enumerate(x):
            zij = topo.ground_altitude(xj, yi)
            lla = topo.local_to_lla((xj, yi, zij))
            la.append(lla[0])
            lo.append(lla[1])
    p, _, _ = numpy.histogram2d(lo, la, (longitude, latitude))
    return p

def process(path, latitude, longitude, rate):
    """Process a set of event files"""
    tstart = time.time()
    generated, data = 0, []
    for name in os.listdir(path):
        if not name.startswith("events"):
            continue
        filename = os.path.join(path, name)
        t0 = time.time()
        print "o Processing", filename
        for event in EventIterator(filename):
            generated += add_triggers(event, latitude, longitude, rate)
        print "  --> Done in {:.1f} s".format(time.time() - t0)
    if generated > 0:
        d = antennas_density(latitude, longitude)
        rate /= (generated * d)
    print rate.max()
    latitude = 0.5 * (latitude[1:] + latitude[:-1])
    longitude = 0.5 * (longitude[1:] + longitude[:-1])

    while path.endswith("/"):
        path = path[:-1]
    path += ".triggers.p"
    with open(path, "wb+") as f:
        pickle.dump((generated, latitude, longitude, rate), f, -1)
    print "o Triggers dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    latitude = numpy.linspace(41.797, 42.396, 101)
    longitude = numpy.linspace(85.387, 87.204, 101)
    rate = numpy.zeros((100, 100))
    process(sys.argv[1], latitude, longitude, rate)
