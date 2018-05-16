#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time
import numpy as np

sys.path.append("../../lib/python")
from grand_tour import Topography
from retro.event import EventIterator
from common import checkTrig
from common import checkCluster


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)


def summarise_tau(event,opt='sel'):
    """Extract the relevant tau info from an event"""
    global origin, topo
    year = 365.25 * 24. * 60. * 60.
    _, e, (x, y, z), u, (la, lo, h), (t, p) = event["tau_at_decay"]   # Positions in meters
    [xx, yx, zx] = np.array([x, y, z]) + np.array(u) * 8000  # Assuming Xmax 8km after tau decay
    if event["origin"] != origin: #  Update the topography handle if required
      print "Loading new topo tile..."
      latitude, longitude = origin = event["origin"]
      topo = Topography(latitude=latitude, longitude=longitude,path="share/topography", stack_size=121)

    zg = topo.ground_altitude(x, y)   #Meters
    zgx = topo.ground_altitude(xx, yx)   #Meters
    height = z-zg
    heightx = zx-zg
    #print x,y,height, xx,yx,heightx
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

    targetw = np.loadtxt('isnotintarget.txt')
    if opt == 'sel':  # Trigger analysis
    	try:
    	  volts = np.array(event["voltage"])
    	  nv = np.shape(volts)[0]    # Number of simulated antennas
    	except:
    	  #print "No voltage for shower",event["tag"]
	  return n, None
    	antsIDs = checkTrig(event)
    else:  # Only cone analysis
	antsIDs = checkCone(event)
    nt = len(antsIDs)   # Nb of trigged antennas
    bTrig, bCluster = checkCluster(event,antsIDs)
    nl = sum(bCluster) # Nb of clustered antennas

    if bTrig and w > 0.:  # This shower triggered
        posCluster = [xants[antsIDs[bCluster]], yants[antsIDs[bCluster]],zants[antsIDs[bCluster]]]  # matrix of trigged antenna pos
	posDecayTile = np.transpose(np.tile(np.array([x, y, z]),(nl,1)))  # matrix of decay pos
	dist = np.linalg.norm(posCluster-posDecayTile,axis=0)
	dmin,dmax = np.min(dist)/1e3,np.max(dist)/1e3
	mask = np.in1d(w,targetw)
	if sum(mask)>0:
	  print event["tag"],w
        data = [w, e, x, y, z, la, lo, h, t, p, height, heightx, nc, nv, nt, nl, dmin, dmax]
    else:
        data = [0, e, x, y, z, la, lo, h, t, p, height, heightx, nc, nv, nt, 0,0,0]
    return n, data


def process(path, summarise,opt='sel'):
    """Summarise a set of event files"""
    tstart = time.time()
    global origin, topo
    origin, topo = None, None
    generated, data = 0, []
    nf = 0
    novolt = 0
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
	    else:
	      novolt += 1
	      #print 'novolt=',novolt
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
