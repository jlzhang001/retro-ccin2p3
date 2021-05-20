#!/usr/bin/env python
#===============================================================================
# Grid Engine steering options
#===============================================================================
# Submit job under TREND group
#$ -P P_trend

# The job name
#$ -N retro-ww

# Merge the stdout et stderr in a single file
#$ -j y

# Files .e et .o copied to current working directory
#$ -cwd

# Notify stop and kill signals before issuing them.
#$ -notify

# Job array indices
#$ -t 1

# CPU time
#$ -l ct=12:00:00

# Memory
#$ -l h_rss=4.0G

# Disk space
#$ -l h_fsize=16.0G

# Request CentOS7
#$ -l os=cl7

# Request access to iRODS and /sps
#$ -l irods=1,sps=1
#===============================================================================
"""Generate tau decays for the hotspot array of 66x150 km2
"""

import sys
sys.path.append("lib/python")

import os
import time

import ccin2p3

# Settings
SETUP = (500., 4.5)
RETRO_HASHTAG = "cfdecde"
N_EVENTS = 10000
OUTDIR = "irods://grand/test/tau-ww"

topography = {
    "latitude" : 42.0,
    "longitude" : 86.0,
    "density" : 2.65E+03,
    "path" : "topography",
    "stack_size" : 50 }


# Install RETRO
print "# Installing RETRO ..."
t0 = time.time()
rootdir, tmpdir, tag = ccin2p3.retro_install(hashtag=RETRO_HASHTAG,
                                             topography=(topography, 5, 6))
print "  --> Done in {:.1f} s".format(time.time() - t0)


# Generate the configuration card
options = {
    "generator" : { "theta" : ["uniform", [85.0, 95.0]],
                    "energy" : [10**7.5, 10**10.5],
                    "position" : ["geodetic",
                                  [[-0.5, 0.5],
                                   [-0.5, 0.5],
                                   [0, 5E+03]]] },
    "topography" : topography,

    "selector" : {
        "vertex": { "limit": 4.0 }},

    "primary": {
        "events": 10000,
        "requested": 100,
        "longitudinal": False },

    "setup": {
        "path": "ww://{:.0f}/{:.2f}".format(*SETUP) }}


# Run RETRO
print ""
print "# Running RETRO ..."
t0 = time.time()
outfile = "events.{:}.json".format(tag)
ccin2p3.retro_run(N_EVENTS, options, path=os.path.join(
    OUTDIR, outfile))
print "  --> Done in {:.1f} s".format(time.time() - t0)
