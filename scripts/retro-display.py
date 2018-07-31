#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
#  Author: Valentin NIESS (niess@in2p3.fr)
#
#  A basic event display for RETRO, based on matplotlib.
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>

import sys
import os
CC = 0
if CC==1:
  RETRODIR = "/pbs/throng/trend/soft/sim/GRANDsim/retro/"
  MAPDIR = "/sps/hep/trend/neu/maps/ASTER-GDEM2/"
else:
  RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
  MAPDIR = "/home/martineau/GRAND/GRAND/data/maps/ASTER-GDEM2/"
sys.path.append(RETRODIR)
sys.path.append(RETRODIR+"lib/python/")
from retro.event import EventIterator
from grand_tour import Topography

#from display import Display  # Not used anymore
sys.path.append("/home/martineau/GRAND/soft/neutrinos/retro/lib/python")
from display_omh import DisplayOMH # in /home/martineau/GRAND/soft/neutrinos/retro/lib/python
import pylab as pl
import numpy as np

display = DisplayOMH()

def ants(event):

  ra = np.array(event["antennas"])[:, :3] * 1E-03
  #pl.figure(1)
  #pl.plot(ra[:,0],ra[:,1],'go') 
  #pl.show()
  return ra[:,0],ra[:,1]

def dumpIni(event):

  global origin, topo
  if event["origin"] != origin:
      # Update the topography handle if required
      latitude, longitude = origin = event["origin"]
      print "Topography around (lat,long)=",latitude, longitude
      topo = Topography(latitude=latitude, longitude=longitude,  	  path="share/topography", stack_size=121)


  # Format info for txt file to be analyzed in ciompareConeIniNew.m
  _, e, (x0, y0, z0), u, (la, lo, h), (t, p) = event["tau_at_decay"]
  print "Decay pos (m) = ",x0, y0, z0
  zg = topo.ground_altitude(x0, y0)
  print "Ground altitude @ decay point (GRANDref, m)=",zg 
  print "Shower vec dir =",u
  
  # Compute the shower energy
  shower_energy = 0.
  for (pid_, momentum) in event["decay"]:
    aid = abs(pid_)
    if aid in (12, 13, 14, 16):
      continue
    shower_energy += sum(m**2 for m in momentum)**0.5
  print "shower energy (GeV) =",shower_energy
  ra = np.array(event["antennas"])[:, :3]
  nants = np.shape(ra)[0]
  print "Nants in cone = ", nants

  x = np.array([shower_energy,u[0],u[1],u[2],x0,y0,z0,zg,nants])
  np.savetxt(f,x.reshape(1, x.shape[0]),fmt="%3.4f")  
    
if __name__ == "__main__":
    i = 0
    global origin, topo
    origin, topo = None, None
    #xall = np.zeros(shape=(0,0))
    #yall = np.zeros(shape=(0,0))
    path = sys.argv[1]
    f = open('compIniGRAND.txt','ab')
    #np.savetxt(f,"#E (GeV), u(GRANDconv), x_decay(GRANDconv)",delimiter=" ",fmt="%s")
    for name in os.listdir(path):
        if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        print "o Processing", filename
        for event in EventIterator(filename):
	    i += 1
	    #xa,ya = ants(event)
	    #xall = np.append(xall,xa)
	    #yall = np.append(yall,ya)
	    #display(event)
	    #fresnel(event)
	    dumpIni(event)
	    
        #if i >1000:
        #  break
    
    f.close()
    #pl.figure()
    #pl.subplot(211)
    #pl.hist(xall,100)
    #pl.xlabel('Antenna X pos (km)')
    #pl.subplot(212)
    #pl.hist(yall,100)
    #pl.xlabel('Antenna Y pos (km)')
    #pl.show()
