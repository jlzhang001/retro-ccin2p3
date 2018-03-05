#!/usr/bin/env python
import cPickle as pickle
import os
import sys
import time

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors

from retro.event import EventIterator
from grand_tour import Topography
from common import checkTrig


# Get the global topography handle
topo = Topography(latitude=42.1, longitude=86.3, path="share/topography",stack_size=49)

DISPLAY = 1
step = 500
plt.style.use("../retro/plugins/display/deps/mplstyle-l3/style/l3.mplstyle")

def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-04 / (3. * e**2)

def histogram2d(*args, **kwargs):
    """Encapsulate numpy.histogram2d for matrix convention compatibility"""
    n, x, y = numpy.histogram2d(*args, **kwargs)
    n = n.T
    return n, x, y
    
def add_triggers(event, latitude, longitude, rate, opt='sel'):
    """Extract the antenna trigger rates from an event"""
    
    year = 365.25 * 24. * 60. * 60.
    w = [v[0] * primary_flux(v[1]) for v in event["primaries"]]
    w = sum(w) * year / event["statistics"][1]
    n = event["statistics"][0]

    if opt=='sel':
        [bTrig,antsIDs] = checkTrig(event) 
    else:
        bTrig = True
	x = np.array(event["antennas"])[:,0]
	antsIDs = np.array(range(len(x)))
    
    la, lo = [],[]
    xant, yant = [],[]
    if bTrig and w>0:  # Shower was detected!
      j = 0
      for x, y, z, _, _ in event["antennas"]:
          if np.any(antsIDs==j):  #This is a trigged antenna
	    lla = topo.local_to_lla((x, y, z))
            la.append(lla[0])
            lo.append(lla[1])
	    xant.append(x)
	    yant.append(y)
	  j+=1

      if la:
          p, _, _ = np.histogram2d(lo, la, (longitude, latitude))
          rate += w * p

      if DISPLAY:
        # Now do plots
        #
        xant = np.array(xant)
        yant = np.array(yant)
        z = np.array(z)

        yantr = yant[::-1] # Pointing east
        plt.figure(31)
        plt.scatter(yantr/1e3,xant/1e3, marker='o', alpha=0.75)
        plt.xlabel(r"Easting (km)")
        plt.ylabel(r"Northing (km)")
        plt.grid(True)

        la = np.array(la)
        lo = np.array(lo)
        # LatLong
        plt.figure(32)
        plt.scatter(lo,la, marker='o', alpha=0.75)
        plt.xlabel(r"longitude (deg)")
        plt.ylabel(r"latitude (deg)")
        plt.grid(True)

        fig = plt.figure(1)
        plt.scatter(yantr/1e3,xant/1e3,marker='o')
        fig = plt.figure(2)
        plt.scatter(lo,la,marker='o')
	    
    return n

def antennas_density(latitude, longitude):
    """Compute the density of antennas"""
    Dx, Dy, s = 0.5 * 66.5E+03, 0.5 * 150.4E+03, step
    x = np.arange(-Dx, Dx + s, s)
    y = np.arange(-Dy, Dy + s, s)
    z = np.zeros((len(y),len(x)))
    la, lo = [], []
    for i, yi in enumerate(y):
        for j, xj in enumerate(x):
            zij = topo.ground_altitude(xj, yi)
            lla = topo.local_to_lla((xj, yi, zij))
            la.append(lla[0])
            lo.append(lla[1])
	    z[i,j]=zij
	    
    p, _, _ = np.histogram2d(lo, la, (longitude, latitude))

    if DISPLAY:
      # Now do plots
      la = np.array(la)
      lo = np.array(lo)
      n = len(y)/2
      lay=la[n*len(x):(n+1)*len(x)]
      a = np.array(range(0,len(y)))*len(x)
      lox=lo[a]
      zt = np.transpose(z)
      yr = y[::-1]
 
      # Local referential
      plt.figure(1)
      norm = colors.Normalize(vmin=0,vmax=1800)
      plt.pcolor(yr/1e3, x/1e3,zt,cmap="terrain", alpha=0.75,norm=norm)
      plt.xlim(min(yr/1e3),max(yr/1e3) )
      plt.ylim(min(x/1e3),max(x/1e3))
      plt.xlabel(r"Easting (km)")
      plt.ylabel(r"Northing (km)")
      plt.colorbar()
      #plt.show()
 
      # LatLong
      plt.figure(2)
      plt.pcolor(lox,lay, zt, cmap="terrain", alpha=0.75,norm=norm)
      plt.xlim(min(lox),max(lox))
      plt.ylim(min(lay),max(lay))
      plt.xlabel(r"longitude (deg)")
      plt.ylabel(r"latitude (deg)")
      plt.colorbar()

      plt.show()
  
    return p

def process(path, latitude, longitude, rate,opt='sel'):
    """Process a set of event files"""
    tstart = time.time()
    generated, data = 0, []
    nf = 0
    for name in os.listdir(path):
        nf += 1
	if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        t0 = time.time()
        #print "o Processing", filename
        for event in EventIterator(filename):
            generated += add_triggers(event, latitude, longitude, rate,opt)
        #print "  --> Done in {:.1f} s".format(time.time() - t0)
	if float(nf)/100 == np.floor(nf/100):
  	  print 'Nb of json files processed:',nf
  	if nf==1000:
	  break


    if generated > 0:
        d = antennas_density(latitude, longitude)
        rate /= (generated * d)
    latitude = 0.5 * (latitude[1:] + latitude[:-1])
    longitude = 0.5 * (longitude[1:] + longitude[:-1])

    while path.endswith("/"):
        path = path[:-1]
    path += ".triggers.p"
    #ratet = np.transpose(rate)
    with open(path, "wb+") as f:
        pickle.dump((generated, latitude, longitude, rate), f, -1)
    print "o Triggers dumped to", path
    print "  --> All done in {:.1f} s".format(time.time() - tstart)


if __name__ == "__main__":
    latitude = np.linspace(41.797, 42.396, 101) # How was that computed?
    longitude = np.linspace(85.387, 87.204, 101)  # How was that computed?
    rate = np.zeros((100, 100))
    if len(sys.argv)==3:
        process(sys.argv[1],latitude,longitude,rate,sys.argv[2])
    else:
        process(sys.argv[1],latitude,longitude,rate)
