"""Utilities for running RETRO at CCIN2P3
"""

import json
import math
import os
import shutil
import subprocess
import sys
import time

TOPO_PATH = { "ASTER3" : "/sps/trend/jlzhang/grandsoft/maps/e4ftl01.cr.usgs.gov/ASTT/ASTGTM.003/2000.03.01",
              "ASTER" : "/sps/trend/neu/maps/ASTER-GDEM2",
              "SRTMGL1" : "/sps/trend/niess/SRTMGL1.003" }


def system(cmd, mute=False):
    """Run a system command
    """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
                         stderr=subprocess.PIPE)
    out, err = p.communicate()
    if mute:
        return
    if err:
        raise RuntimeError(err)
    return out


def retro_install(path=None, topography=None, hashtag=None, dem="ASTER"):
    """Install RETRO locally
    """

    # Set the path
    if path is None:
        path = os.getenv("TMPDIR")
        if (path is None) or (path == "/scratch"):
            user = os.getenv("USER")
            path = os.path.join("/tmp", user)
    if not os.path.exists(path):
        os.makedirs(path)
    rootdir = os.getcwd()
    os.chdir(path)

    # Get RETRO from GitHub
    if os.path.exists("retro"):
        shutil.rmtree("retro")
    system("git clone https://github.com/grand-mother/retro", mute=True)
    if hashtag is not None:
        system("cd retro && git checkout " + hashtag, mute=True)

    # Build RETRO
    system("cd retro && ./build.sh", mute=True)
    sys.path.append(os.path.join(path, "retro", "lib", "python"))

    # Get the topography tiles from /sps
    if topography:
        sx, sy = None, None
        if not isinstance(topography, dict):
            topography, sx, sy = topography
        if os.path.exists(topography["path"]):
            shutil.rmtree(topography["path"])
        os.makedirs(topography["path"])
        lat0 = int(topography["latitude"])
       lng0 = int(topography["longitude"])
        if sx is None:
            sx = sy = (int(math.sqrt(topography["stack_size"])) - 1) / 2
        for i in xrange(-sx, sx + 1):
            lat = lat0 + i
            if lat >= 0: sn = "N"
            else: sn = "S"
            for j in xrange(-sy, sy + 1):
                lng = lng0 + j
                if lng >= 0: ew = "E"
                else: ew = "W"
                if dem == "ASTER":
                    filename = "ASTGTM2_{:}{:02d}{:}{:03d}_dem.tif".format(
                        sn, lat, ew, lng)
                elif dem == "ASTER3":
                       filename = "ASTGTMV003_{:}{:02d}{:}{:03d}_dem.tif".format(
                          sn, lat, ew, lng)
                elif dem == "SRTMGL1":
                      filename = "{:}{:02d}{:}{:03d}.SRTMGL1.hgt".format(
                         sn, lat, ew, lng)
                     filename = "{:}{:02d}{:}{:03d}.SRTMGL1.hgt".format(
                         sn, lat, ew, lng)
                else :
                     raise RuntimeError("no dem files!!!")

                path = os.path.join(TOPO_PATH[dem], filename)
                shutil.copy(path, topography["path"])

    # Build the worker tag
    tag = os.getenv("JOB_ID")
    if tag is None:
        tag = "local"
    else:
        tid = os.getenv("SGE_TASK_ID")
        if tid is not None:
            tag = ".".join((tag, tid))

    return rootdir, path, tag


def retro_run(events, options, setup=None, path=None):
    """Run RETRO with the given options and antenna positions
    """

    # Dump the configuration card
    if ((path is None) or (path.startswith("irods://")) or
        (path.startswith("hpss://"))):
        outfile = "events.json"
    else:
        outdir = os.path.dirname(path)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir)
        outfile = path
    card = options.copy()
    card["processor"] = {"requested": events}
    card["logger"] = {"path": outfile}
    if setup is not None:
        card["setup"] = { "path": "setup.json" }

    with open("card.json", "wb+") as f:
        json.dump(card, f)

    # Dump the antenna layout
    if setup is not None:
        with open("setup.json", "wb+") as f:
            json.dump(setup, f)

    # Run RETRO
    system("./retro/bin/retro card.json")

    # Copy the data if required
    if path.startswith("irods://"):
        path = path.split("irods://")[-1]
        for i in xrange(20):
            try:
                system("iput -f {:} {:}".format(outfile, path))
            except RuntimeError as err:
                print err
                sys.stdout.flush()
                time.sleep(6.)
            else:
                break
        else:
            print "error: failed to upload", outfile
    elif path.startswith("hpss://"):
        path = path.replace("hpss://", "/hpss/in2p3.fr/group/trend/")
        system("rfcp {:} {:}".format(outfile, path))
        system("rfchmod 664 {:}".format(path))
                                     
