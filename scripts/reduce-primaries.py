#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time
import numpy as np
from retro.event import EventIterator
from common import checkTrig
from common import checkCluster
from common import setStep
from common import checkCone

#sys.path.append("../retro/lib/python/retro/")


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def summarise_primaries(event,opt='sel'):
    #print """**** Extract the relevant neutrino info from an event"""
    
    year = 365.25 * 24. * 60. * 60.
    nrm = year / event["statistics"][1]
    
    n = event["statistics"][0]
    data = None
    if opt == 'sel':
    	antsIDs = checkTrig(event) 
    else:  # Cone selection
        antsIDs = checkCone(event)
    bTrig, bCluster = checkCluster(event,antsIDs)
	
    if bTrig:  # Shower was detected!
    	v2 = event["primaries"]
	data = []
	for k in range(len(v2)):
	    data.append([v2[k][0] * primary_flux(v2[k][1]) * nrm, v2[k][1]])
	        	    	  
    return n, data


def process(path,summarise,opt='sel'):
    """Summarise a set of event files"""
    tstart = time.time()
    nf = 0
    generated, data = 0, []
    for name in os.listdir(path):
        if not name.endswith("json"):
            continue
	nf += 1
        filename = os.path.join(path, name)
        t0 = time.time()
        #print "o Processing", filename
        for event in EventIterator(filename):
            n, d = summarise_primaries(event,opt)
	    generated += n
            if d is not None:
                data += d
	#print "  --> Done in {:.1f} s".format(time.time() - t0)
	if float(nf)/100 == np.floor(nf/100):
		print 'Nb of json files processed:',nf
	#if nf==19900:
	#  print "#### Warning! Stopping at",nf
	#  break

    if len(data) == 0:
        raise RuntimeError("No event found")

    while path.endswith("/"):
        path = path[:-1]
    if opt == 'sel':
    	path += ".primaries_sel.p"
    else:
    	path += ".primaries.p"
    
    with open(path, "wb+") as f:
        pickle.dump((generated, data), f, -1)
    print "o Events dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    print "Usage: >python reduce-primaries.py <path to json file> [<sel/cone>]"
    print "sel: includes antenna response (default)"
    if len(sys.argv)==3:
        process(sys.argv[1],summarise_primaries,sys.argv[2])
    else:
        process(sys.argv[1],summarise_primaries)
