#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time
import numpy as np

from retro.event import EventIterator
from common import checkTrig
from common import checkCluster


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
    ants = np.array(event["antennas"])
    xants = ants[:,0]
    yants = ants[:,1]
    zants = ants[:,2]
    alpha = ants[:,3]
    beta = ants[:,4]
    nc = np.shape(ants)[0]  # Number of antennas in cone

    if opt == 'sel':  # Trigger analysis
    	try:
    	  volts = np.array(event["voltage"])
    	  na = np.shape(volts)[0]    # Number of simulated antennas
    	except:
    	  #print "No voltage!"
    	  data = [w, e, x, y, z, la, lo, h, t, p, nc, 1, 1]  # Unit value for na & nt allowing for log plots (see show-hotspot.py)
    	  return n, data, [nc, 0, 0]

    	_,antsIDs = checkTrig(event,'a')
    	na = len(antsIDs)   # Nb of simulated antennas
    	bTrig, bCluster = checkCluster(event,antsIDs,'a')
    	nt = sum(bCluster) # Nb of trigged antennas
    else:  # Only cone analysis
	bTrig = True
	na, nt = 1, 1
	
    if bTrig and w > 0.:  # THis shower triggered
        data = [w, e, x, y, z, la, lo, h, t, p, nc, na, nt]
    else:
        data = None
	 	
    dataAll = [nc,na,nt,alpha,beta]
    return n, data, dataAll


def process(path, summarise,opt='sel'):
    """Summarise a set of event files"""
    tstart = time.time()
    generated, data, dataAll = 0, [], []
    nf = 0
    for name in os.listdir(path):
        if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        t0 = time.time()
        #print "o Processing", filename
        for event in EventIterator(filename):
            n, d, d2 = summarise(event,opt)
            generated += n
            if d is not None:
                data.append(d)
	    dataAll.append(d2)
	if float(nf)/100 == np.floor(nf/100):
	    print 'Nb of json files processed:',nf
	
	nf += 1
        #print "  --> Done in {:.1f} s".format(time.time() - t0)
        #if nf==1000:
	#  break

    if len(data) == 0:
        raise RuntimeError("No event found")

    while path.endswith("/"):
        path = path[:-1]
    if opt=='sel':
    	path += "_sel.p"
    else:
    	path += ".p"
    
    with open(path, "wb+") as f:
        pickle.dump((generated, data, dataAll), f, -1)
    print "o Events dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    print "Usage: >python reduce-events.py <path to json file> [<sel/whatever>]"
    print "sel: includes antenna response (default)"

    if len(sys.argv)==3:
        process(sys.argv[1], summarise_tau,sys.argv[2])
    else:
        process(sys.argv[1], summarise_tau,'sel')
	
