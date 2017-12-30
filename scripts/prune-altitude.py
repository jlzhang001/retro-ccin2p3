#!/usr/bin/env python

import json
import os
import time

from retro.event import EventIterator, EventLogger
from grand_tour import Topography

# Settings
DATA_DIR = "share/events/earth-skimming"

outdir = os.path.join(DATA_DIR, "pruned")
if not os.path.exists(outdir):
    os.makedirs(outdir)

topography = Topography(latitude=43, longitude=87, path="flat/10")

for filename in os.listdir(DATA_DIR):
    if not filename.startswith("events"):
        continue
    t0 = time.time()
    print "# Processing", filename
    subindex = 1
    path = os.path.join(DATA_DIR, filename)
    for i, event in enumerate(EventIterator(path)):
        if (i % 100) == 0:
            v = filename.split(".")
            v.insert(-1, str(subindex))
            newname = ".".join(v)
            subindex += 1
            path = os.path.join(DATA_DIR, "pruned", newname)
            log_event = EventLogger(path=path)
        tau = event["tau_at_decay"]
        tau.append(topography.local_to_lla(tau[1])[2])
        event.pop("previous")
        log_event(**event)
    print "  --> Done in {:.1f} s".format(time.time() - t0)
