#!/usr/bin/env python

import os
import numpy
import matplotlib.pyplot as plt
from retro.event import EventIterator
from grand_tour import Topography


# Settings
DATA_DIR = "share/events"
PHI0 = 2E-08
E0 = 10**7.5
E1 = 10**14.5

# Load and parse the events
topo = Topography(43, 87, "flat/4")
theta, phi, altitude = [], [], []
position, weight = [], []
generated = 0
for name in os.listdir(DATA_DIR):
    if not name.startswith("events"):
        continue
    filename = os.path.join(DATA_DIR, name)
    for event in EventIterator(filename):
        _, r, u = event["tau_at_decay"]
        t, p = topo.local_to_angular(r, u)
        _, _, a = topo.local_to_lla(r)
        w = sum(v[0] for v in event["primaries"][0]) / event["primaries"][1]
        w *= event["statistics"][0]
        generated += event["statistics"][1]
        theta.append(t)
        phi.append(p)
        altitude.append(a * 1E-03)
        position.append(r)
        weight.append(w)
position = numpy.array(position) * 1E-03
year = 365.25 * 24. * 60. * 60.
weight = numpy.array(weight) * PHI0 * (1. / E0 - 1. / E1) * year
mu = sum(weight) / generated
sigma = sum(weight**2) / generated
sigma = numpy.sqrt((sigma - mu**2) / generated)
print "Total rate = {:.3f} +-{:.3f} a^-1".format(mu, sigma)


def plot_histogram(samples):
    """Plot a 1D histogram of sampled values
    """
    n, x = numpy.histogram(samples, 40, weights=weight / generated)
    n2, _ = numpy.histogram(samples, x, weights=weight**2 / generated)
    x = 0.5 * (x[1:] + x[:-1])
    norm = 1. / (x[1] - x[0])
    p = n * norm
    dp = numpy.sqrt((n2 - n * n) / generated) * norm
    plt.figure()
    plt.errorbar(x, p, xerr=x[1] - x[0], yerr=dp, fmt="ko")


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

plot_histogram(theta)
plt.xlabel(r"$\theta_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")

plot_histogram(phi)
plt.axis((-200., 200., 0., 8E-04))
plt.xticks(numpy.linspace(-180., 180., 5))
plt.yticks(numpy.linspace(0., 8E-04, 5))
plt.xlabel(r"$\phi_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")

plot_histogram(altitude)
# plt.axis((0., 1., 0., 3.))
plt.xlabel(r"decay altitude (km)")
plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")

plt.figure()
plt.plot(position[:, 0], position[:, 1], "k.")
plt.xlabel(r"x (km)")
plt.ylabel(r"y (km)")

plt.show()
