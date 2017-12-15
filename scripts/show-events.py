#!/usr/bin/env python

import os
import numpy
import matplotlib.pyplot as plt
from retro.event import EventIterator
from grand_tour import Topography


# Settings
DATA_DIR = "/sps/hep/trend/niess/retro"


# Load and parse the events
topo = Topography(43, 87, "flat/4")
theta, phi, altitude = [], [], []
position = []
for name in os.listdir(DATA_DIR):
	if not name.startswith("events"):
		continue
	filename = os.path.join(DATA_DIR, name)
	for event in EventIterator(filename):
		_, r, u = event["tau_at_decay"]
		t, p = topo.local_to_angular(r, u)
		_, _, a = topo.local_to_lla(r)
		theta.append(t)
		phi.append(p)
		altitude.append(a * 1E-03)
		position.append(r)
position = numpy.array(position) * 1E-03


def plot_histogram(samples):
	"""Plot a 1D histogram of sampled values
	"""
	n, x = numpy.histogram(samples, 40)
	x = 0.5 * (x[1:] + x[:-1])
	norm = 1. / numpy.trapz(n, x)
	p = n * norm
	dp = numpy.sqrt(n) * norm
	plt.figure()
	plt.errorbar(x, p, xerr=x[1] - x[0], yerr=dp, fmt="ko")


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

plot_histogram(theta)
plt.xlabel(r"$\theta$ (deg)")
plt.ylabel(r"pdf (deg$^{-1}$)")

plot_histogram(phi)
plt.axis((-200., 200., 0., 4E-03))
plt.xlabel(r"$\phi$ (deg)")
plt.ylabel(r"pdf (deg$^{-1}$)")

plot_histogram(altitude)
plt.axis((0., 1., 0., 3.))
plt.xlabel(r"decay altitude (km)")
plt.ylabel(r"pdf (km$^{-1}$)")

plt.figure()
plt.plot(position[:,0], position[:,1], "k.")
plt.xlabel(r"x (km)")
plt.ylabel(r"y (km)")

plt.show()
