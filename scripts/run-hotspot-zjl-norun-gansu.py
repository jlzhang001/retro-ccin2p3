#!/usr/bin/env python
#===============================================================================
# Grid Engine steering options
#===============================================================================
# Submit job under TREND group
#$ -P P_trend

# Merge the stdout et stderr in a single file
#$ -j y

# The job name
#$ -N hotspot

# Files .e et .o copied to current working directory
#$ -cwd

# Notify stop and kill signals before issuing them.
#$ -notify

# Job array indices
#$ -t 1-100
# CPU time
#$ -l ct=48:00:00

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
ARRAY_SIZE = ((-33.25E+03, 43.25E+03), (-55.52E+03, 75.2E+03))
ANTENNA_SPACING = 500.
ANTENNA_HEIGHT = 4.5
RETRO_HASHTAG = "d468302"
N_EVENTS = 200
SELECTOR_SETUP = "3deg"
#OUTDIR = "irods://grand/sim/hotspot-130x77km2/taus"
OUTDIR = "/pbs/home/j/jlzhang/sps/grandsoft/retro-ccin2p3/scripts/data"

#DEM = "SRTMGL1"
DEM = "ASTER3"
topography = {
    "latitude" : 41.36,
    "longitude" : 96.545,
    "density" : 2.65E+03,
    "path" : "topography",
    "stack_size" : 144 }


# Install RETRO
print "# Installing RETRO ..."
t0 = time.time()
rootdir, tmpdir, tag = ccin2p3.retro_install(hashtag=RETRO_HASHTAG,
                                             topography=(topography, 7, 8),
                                             dem=DEM)
print "  --> Done in {:.1f} s".format(time.time() - t0)


# Generate the configuration card
delta = 300E+03
sx, sy = map(lambda x: (x[0] - delta, x[1] + delta), ARRAY_SIZE)
options = {
    "generator" : { "theta" : ["uniform", [85.0, 95.0]],
                    "energy" : [10**7.5, 10**10.5],
                    "position" : [sx, sy, [0, 5E+03]] },
    "topography" : topography,

    "selector" : {
        "setup" : { "cone" : SELECTOR_SETUP, "xmax": False },
        "vertex": { "limit": 4.0 }},

    "primary": {
        "events": 10000,
        "requested": 100,
        "longitudinal": False}}


# Generate the antenna setup
import grand_tour
opts = topography.copy()
del opts["density"]
topo = grand_tour.Topography(**opts)

dc = ANTENNA_SPACING
nx, ny = map(lambda x: int((x[1] - x[0]) / dc) + 1, ARRAY_SIZE)
setup = []
for i in xrange(ny):
    yi = ARRAY_SIZE[1][0] + i * dc
    for j in xrange(nx):
        xj = ARRAY_SIZE[0][0] + j * dc
        uij, alpha, beta = topo.ground_normal(xj, yi, step=200.,
            angles=True)
        zij = topo.ground_altitude(xj, yi)
        setup.append((xj + uij[0] * ANTENNA_HEIGHT,
                      yi + uij[1] * ANTENNA_HEIGHT,
                      zij + uij[2] * ANTENNA_HEIGHT,
                      alpha, beta))

# Run RETRO

                   
  
