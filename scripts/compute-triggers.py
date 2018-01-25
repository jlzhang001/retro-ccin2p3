#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time

import numpy
from retro.event import EventIterator


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def add_triggers(event, rate_map):
    """Extract the antenna trigger rates from an event"""
    year = 365.25 * 24. * 60. * 60.
    w = [v[0] * primary_flux(v[1]) for v in event["primaries"]]
    w = sum(w) * year / event["statistics"][1]
    n = event["statistics"][0]
    if w > 0.:
        x0, dx, y0, dy, rate = rate_map
        rx, ry = 1. / dx, 1. / dy
        for x, y, _, _, _ in event["antennas"]:
            ix, iy = int((x - x0) * rx), int((y - y0) * ry)
            rate[iy, ix] += w
    return n


def process(path, rate_map):
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
            generated += add_triggers(event, rate_map)
        print "  --> Done in {:.1f} s".format(time.time() - t0)

    path += ".triggers.p"
    with open(path, "wb+") as f:
        pickle.dump((generated, rate_map), f, -1)
    print "o Triggers dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    rate_map = -33500., 500., -75450., 500., numpy.zeros((303, 135))
    process(sys.argv[1], rate_map)
