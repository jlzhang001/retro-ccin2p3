"""Utilities for running RETRO at CCIN2P3
"""

import json
import os
import shutil
import subprocess
import sys


def system(cmd, mute=False):
	"""Run a system command
	"""
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE,
	                     stderr=subprocess.PIPE)
	out, err = p.communicate()
	if mute:
		return
	if err:
		raise RuntimeError(err)
	return out


def retro_install(path=None, hashtag=None):
	"""Install RETRO locally
	"""
	
	# Set the path
	if path is None:
		path = os.getenv("TMPDIR")
		if path == "/scratch":
			user = os.getenv("USER")
			path = os.path.join("/tmp", user)
	if not os.path.exists(path):
		os.path.makedirs(path)
	rootdir = os.getcwd()
	os.chdir(path)


	# Get RETRO from GitHub
	if os.path.exists("retro"):
		shutil.rmtree("retro")
	system("git clone https://github.com/grand-mother/retro", mute=True)
	if hashtag is not None:
		system("cd retro && git checkout " + hashtag, mute=True)


	# Build RETRO
	system("cd retro && ./install.sh", mute=True)
	sys.path.append(os.path.join(path, "retro", "lib", "python"))
	
	
	# Build the worker tag
	tag = os.getenv("JOB_ID")
	if tag is None:
		tag = "local"
	else:
		tid = os.getenv("SGE_TASK_ID")
		if tid is not None:
			tag = ".".join((tag, tid))

	return rootdir, path, tag


def retro_run(events, options, positions=None, outfile=None):
	"""Run RETRO with the given options and antenna positions
	"""

	# Dump the configuration card
	if outfile is None:
		outfile = "events.json"
	card = options.copy()
        card["processor"] = { "requested" : events }
        card["logger"] = { "path" : outfile }
        if positions is not None:
		card["antenna"] = { "position": "positions.json" }

	with open("card.json", "wb+") as f:
		json.dump(card, f)


	# Dump the antenna layout
	if positions is not None:
		with open("positions.json", "wb+") as f:
			json.dump(positions, f)
			
	# Run RETRO
	system(". retro/setup.sh && ./retro/bin/retro-run card.json")
