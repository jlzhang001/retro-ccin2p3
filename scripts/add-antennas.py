#!/usr/bin/env python
import os
from retro.event import EventIterator, EventLogger
from grand_tour import Topography

# Settings
DATA_DIR = "share/events"
TAG = "E.1e17_X.93815_Y.74409_Z.-559_T.91_P.21_D.3750214666968983"

topography = Topography(latitude=43, longitude=87, path="flat/10")

filename = TAG + ".more-antennas.json"
path = os.path.join(DATA_DIR, filename)
log_event = EventLogger(path=path)

filename = TAG + ".json"
path = os.path.join(DATA_DIR, filename)
for event in EventIterator(path):
    antennas = event["antennas"]
    nx, ny = 10, 10
    for i in xrange(nx):
        x = -50E+03 + (i * 100E+03) / (nx - 1)
        for j in xrange(ny):
            y = -50E+03 + (j * 100E+03) / (ny - 1)
            z = topography.ground_altitude(x, y) + 3.
            antennas.append((x, y, z))
    log_event(**event)
