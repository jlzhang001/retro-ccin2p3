#!/usr/bin/env python
import os
import sys
import time

from retro.event import EventIterator, EventLogger

# Format the output directory
out_prefix = sys.argv[1]
while out_prefix.endswith("/"):
    out_prefix = out_prefix[:-1]
out_prefix = os.path.basename(out_prefix)
out_dir = os.path.realpath(
    os.path.join(sys.argv[1], "..", out_prefix + "-split"))
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

# Loop over files
log_event = None
out_index, n_events = 0, 0
for filename in os.listdir(sys.argv[1]):
    if not filename.startswith("events"):
        continue
    path = os.path.join(sys.argv[1], filename)
    print "o Processing file", path
    t0 = time.time()
    # Loop over events in a file
    for event in EventIterator(path):
        event.pop("previous")
        if log_event is None:
            # Create a new output file
            out_index += 1
            path = os.path.join(out_dir, "events.{:}.json".format(out_index))
            log_event = EventLogger(path=path)
        log_event(**event)
        n_events += 1
        if n_events == 400:
            # Stop writting to the current output file
            log_event = None
            n_events = 0
    print "  --> Done in {:.1f} s".format(time.time() - t0)
