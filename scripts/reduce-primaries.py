#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time

from retro.event import EventIterator


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def summarise_primaries(event):
    """Extract the relevant neutrino info from an event"""
    year = 365.25 * 24. * 60. * 60.
    nrm = year / event["statistics"][1]
    data = [(v[0] * primary_flux(v[1]) * nrm, v[1]) for v in event["primaries"]]
    n = event["statistics"][0]
    return n, data


def process(path, summarise):
    """Summarise a set of event files"""
    tstart = time.time()
    generated, data = 0, []
    for name in os.listdir(path):
        if not name.startswith("events"):
            continue
        filename = os.path.join(path, name)
        t0 = time.time()
        print "o Processing", filename
        for event in EventIterator(filename):
            n, d = summarise(event)
            generated += n
            if d is not None:
                data += d
        print "  --> Done in {:.1f} s".format(time.time() - t0)

    if len(data) == 0:
        raise RuntimeError("No event found")

    while path.endswith("/"):
        path = path[:-1]
    path += ".primaries.p"
    with open(path, "wb+") as f:
        pickle.dump((generated, data), f, -1)
    print "o Events dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    process(sys.argv[1], summarise_primaries)
