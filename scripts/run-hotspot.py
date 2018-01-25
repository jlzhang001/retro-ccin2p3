#!/usr/bin/env python
#===============================================================================
# Grid Engine steering options
#===============================================================================
# Submit job under TREND group
#$ -P P_trend

# Merge the stdout et stderr in a single file
#$ -j y

# Files .e et .o copied to current working directory
#$ -cwd

# Notify stop and kill signals before issuing them.
#$ -notify

# Job array indices
#$ -t 1-20

# CPU time
#$ -l ct=12:00:00

# Memory
#$ -l h_rss=4.0G

# Disk space
#$ -l h_fsize=16.0G

# Request CentOS7
#$ -l os=cl7

# Request access to sps
#$ -l sps=1
#===============================================================================
"""Generate tau decays for the hotspot array of 100x100 km2
"""

import sys
sys.path.append("lib/python")

import os
import time

import ccin2p3

# Settings
ARRAY_SIZE = 66.5E+03, 150.4E+03
ANTENNA_HEIGHT = 4.5
RETRO_HASHTAG = "796747f"
N_EVENTS = 1000
OUTDIR = "/sps/hep/trend/niess/retro/hotspot"

topography = {
    "latitude" : 42.1,
    "longitude" : 86.3,
    "path" : "topography",
    "stack_size" : 144 }


# Install RETRO
print "# Installing RETRO ..."
t0 = time.time()
rootdir, tmpdir, tag = ccin2p3.retro_install(hashtag=RETRO_HASHTAG,
                                             topography=(topography, 5, 6))
print "  --> Done in {:.1f} s".format(time.time() - t0)


# Generate the configuration card
sx, sy = map(lambda x: 0.5 * x + 300E+03, ARRAY_SIZE)
options = {
    "generator" :
                [[ 1.0, {
                        "theta" : ["uniform", [85.0, 95.0]],
                        "energy" : [10**7.5, 10**10.5],
                        "position" : [[-sx, sx],
                                      [-sy, sy],
                                      [0, 5E+03] ]}],
                [ 2.0, {
                        "theta" : ["linear", [90.0, 95.0]] }]],
    "topography" : topography,

    "selector" : {
        "vertex": { "limit": 3.0 }},

    "primary": {
        "events": 10000,
        "requested": 100,
        "longitudinal": False}}


# Generate the antenna setup
import grand_tour
topo = grand_tour.Topography(**options["topography"])

dc = 500.
nx, ny = map(lambda x: int(x / dc) + 1, ARRAY_SIZE)
setup = []
for i in xrange(ny):
    yi = -0.5 * ARRAY_SIZE[1] + i * dc
    for j in xrange(nx):
        xj = -0.5 * ARRAY_SIZE[0] + j * dc
        uij, alpha, beta = topo.ground_normal(xj, yi, angles=True)
        zij = topo.ground_altitude(xj, yi)        
        setup.append((xj + uij[0] * ANTENNA_HEIGHT,
                      yi + uij[1] * ANTENNA_HEIGHT,
                      zij + uij[2] * ANTENNA_HEIGHT,
                      alpha, beta))


# Run RETRO
print ""
print "# Running RETRO ..."
t0 = time.time()
outfile = "events.{:}.json".format(tag)
ccin2p3.retro_run(N_EVENTS, options, setup, outfile=os.path.join(
    OUTDIR, outfile))
print "  --> Done in {:.1f} s".format(time.time() - t0)
