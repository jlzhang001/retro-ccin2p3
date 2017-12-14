#!/usr/bin/env python
#===============================================================================
# Grid Engine steering options
#===============================================================================
## Submit job under TREND group
#$ -P P_trend

## Merge the stdout et stderr in a single file
#$ -j y

## Files .e et .o copied to current working directory
#$ -cwd

## Notify stop and kill signals before issuing them.
#$ -notify

## Job array indices
#$ -t 1-4

## CPU time
#$ -l ct=05:00:00

## Memory
#$ -l vmem=3.0G

## Disk space
#$ -l fsize=1.0G

## Request CentOS7
#$ -l os=cl7

## Request access to sps
#$ -l sps=1
#===============================================================================
"""Generate tau decays for flat array of 100x100 km2
"""

import sys
sys.path.append("lib/python")

import os
import time

import ccin2p3

# Settings
ARRAY_SIZE = 100E+03
ANTENNA_HEIGHT = 3.
RETRO_HASHTAG = "898d13c"
N_EVENTS = 5000
OUTDIR = "/sps/hep/trend/niess/retro"


# Install RETRO
print "# Installing RETRO ..."
t0 = time.time()
rootdir, tmpdir, tag = ccin2p3.retro_install(hashtag=RETRO_HASHTAG)
print "  --> Done in {:.1f} s".format(time.time() - t0)


# Generate the configuration card
s = ARRAY_SIZE / 2 + 50E+03
options = {
	"generator" : {
                "theta" : [85.0, 95.0],
                "energy" : [10**7.5, 10**11.5],
                "position" : [
                        [-s, s], [-s, s], [0, 1E+03] ]},

        "topography" : {
                "latitude" : 43,
                "longitude" : 87,
                "path" : "flat/4" },

        "primary" : {
                "events" : 10000,
                "requested" : 100,
                "longitudinal" : False },
}


# Generate the antenna setup
import grand_tour
topo = grand_tour.Topography(**options["topography"])

dc = 500.
n = int(ARRAY_SIZE / dc) + 1
positions = []
for i in xrange(n):
	yi = -0.5 * ARRAY_SIZE + i * dc 
	for j in xrange(n):
		xj = -0.5 * ARRAY_SIZE + j * dc
		zij = topo.ground_altitude(xj, yi) + ANTENNA_HEIGHT
		positions.append((xj, yi, zij))


# Run RETRO
print ""
print "# Running RETRO ..."
t0 = time.time()
outfile = "events.{:}.json".format(tag)
ccin2p3.retro_run(N_EVENTS, options, positions, outfile=os.path.join(
    OUTDIR, outfile))
print "  --> Done in {:.1f} s".format(time.time() - t0)
