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
#from display import Display  # Not used anymore
sys.path.append("/home/martineau/GRAND/soft/neutrinos/retro/lib/python")
from display_omh import DisplayOMH # in /home/martineau/GRAND/soft/neutrinos/retro/lib/python


display = DisplayOMH()

if __name__ == "__main__":
    path = sys.argv[1]
    for name in os.listdir(path):
        if not name.endswith("json"):
	    continue
	filename = os.path.join(path, name)
        print "o Processing", filename
        for event in EventIterator(filename):
            display(event)
	    #fresnel(event)
