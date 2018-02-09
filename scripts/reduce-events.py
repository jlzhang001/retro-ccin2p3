#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time
import numpy as np

from retro.event import EventIterator
from common import checkTrig

def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def summarise_tau(event,opt='sel'):
    """Extract the relevant tau info from an event"""
    year = 365.25 * 24. * 60. * 60.
    _, e, (x, y, z), _, (la, lo, h), (t, p) = event["tau_at_decay"]
    w = [v[0] * primary_flux(v[1]) for v in event["primaries"]]
    w = sum(w) * year / event["statistics"][1]
    n = event["statistics"][0]
    data = None
    if opt == 'sel':
    	[bTrig,_,_] = checkTrig(event,5,'a')
    else:
	bTrig = True
	
    if bTrig and w > 0.:
        data = [w, e, x, y, z, la, lo, h, t, p]
    return n, data


def process(path, summarise,opt='sel'):
    """Summarise a set of event files"""
    tstart = time.time()
    generated, data = 0, []
    nf = 0
    for name in os.listdir(path):
        if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        t0 = time.time()
        #print "o Processing", filename
        for event in EventIterator(filename):
            n, d = summarise(event,opt)
            generated += n
            if d is not None:
                data.append(d)
	if float(nf)/100 == np.floor(nf/100):
	    print 'Nb of json files processed:',nf
	
	nf += 1
        #print "  --> Done in {:.1f} s".format(time.time() - t0)
        if nf==1000:
	  break

    if len(data) == 0:
        raise RuntimeError("No event found")

    while path.endswith("/"):
        path = path[:-1]
    if opt=='sel':
    	path += "_sel.p"
    else:
    	path += ".p"
    
    with open(path, "wb+") as f:
        pickle.dump((generated, data), f, -1)
    print "o Events dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    print "Usage: >python reduce-events.py <path to json file> [<sel/whatever>]"
    print "sel: includes antenna response (default)"

    if len(sys.argv)==3:
        process(sys.argv[1], summarise_tau,sys.argv[2])
    else:
        process(sys.argv[1], summarise_tau,'sel')
	
