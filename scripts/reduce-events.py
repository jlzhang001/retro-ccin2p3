#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time

sys.path.append("../../lib/python")
from retro.event import EventIterator


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def summarise_tau(event):
    """Extract the relevant tau info from an event"""
    year = 365.25 * 24. * 60. * 60.
    _, e, (x, y, z), _, (la, lo, h), (t, p) = event["tau_at_decay"]
    w = [v[0] * primary_flux(v[1]) for v in event["primaries"]]
    w = sum(w) * year / event["statistics"][1]
    n = event["statistics"][0]
    data = None
    if w > 0.:
        data = [w, e, x, y, z, la, lo, h, t, p]
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
                data.append(d)
        print "  --> Done in {:.1f} s".format(time.time() - t0)

    if len(data) == 0:
        raise RuntimeError("No event found")

    while path.endswith("/"):
        path = path[:-1]
    path += ".p"
    with open(path, "wb+") as f:
        pickle.dump((generated, data), f, -1)
    print "o Events dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    process(sys.argv[1], summarise_tau)
