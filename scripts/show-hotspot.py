#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import os
import numpy
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from retro.event import EventIterator


# Settings
DATA_DIR = "share/events"


def primary_flux(e):
    # Waxman-Bahcall bound with 1 / 3 of tau neutrinos
    return 2E-04 / (3. * e**2)


# Load and parse the events
year = 365.25 * 24. * 60. * 60.

def add_events(path, energy, position, theta, phi, altitude, weight):
    """Load events from a file and append them
    """
    generated = 0
    for event in EventIterator(path):
        _, e, r, u, (_, _, a), (t, p) = event["tau_at_decay"]
        w = numpy.array([v[0] * primary_flux(v[1])
                         for v in event["primaries"]])
        norm = year / event["statistics"][1]
        w = sum(w) * norm
        generated += event["statistics"][0]
        if w == 0:
            continue
        theta.append(t)
        phi.append(p)
        altitude.append(a * 1E-03)
        position.append(r)
        energy.append(e)
        weight.append(w)
    return generated

# Merge a set of decay files
theta, phi, altitude = [], [], []
position, energy, weight = [], [], []
generated = 0
for name in os.listdir(DATA_DIR):
    if not name.startswith("events"):
        continue
    filename = os.path.join(DATA_DIR, name)
    generated += add_events(filename, energy, position, theta, phi, altitude,
                            weight)

position = numpy.array(position) * 1E-03
weight = numpy.array(weight)
mu = sum(weight) / generated
sigma = sum(weight**2) / generated
sigma = numpy.sqrt((sigma - mu**2) / generated)
print "Rate = {:.3f} +-{:.3f} a^-1".format(mu, sigma)

def plot_histogram(samples, weight, generated, plot=plt.plot, clr="k",
                   new=True, factor=None):
    """Plot a 1D histogram of sampled values
    """
    if plot == plt.loglog:
        n, b = numpy.histogram(numpy.log(samples), 40,
                               weights=weight / generated)
        n2, _ = numpy.histogram(numpy.log(samples), b,
                                weights=weight**2 / generated)
        norm = 1. / (b[1] - b[0])
        x = numpy.exp(0.5 * (b[1:] + b[:-1]))
        norm /= x
        b = numpy.exp(b)
        xerr = (x - b[:-1], b[1:] - x)
    else:
        n, b = numpy.histogram(samples, 40, weights=weight / generated)
        n2, _ = numpy.histogram(samples, b, weights=weight**2 / generated)
        x = 0.5 * (b[1:] + b[:-1])
        xerr = x[1] - x[0]
        norm = 1. / xerr
    p = n * norm
    dp = numpy.sqrt((n2 - n * n) / generated) * norm
    if new:
        plt.figure()
    if factor is None:
        y = p
        yerr = dp
    else:
        factor = x**factor
        y = p * factor
        yerr = dp * factor
    if sum(y) > 0:
        plot(x, y, clr + "o")
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=clr + "o")
    return x, p


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

x, p = plot_histogram(energy, weight, generated, plot=plt.loglog)

plt.xlabel(r"energy, E$_\tau$ (GeV)")
plt.ylabel(r"E$_\tau \times$ rate (a$^{-1}$)")
plt.savefig("tau-energy.png")

plt.figure()
plt.semilogx(x, cumtrapz(p, x, initial=0.) / numpy.trapz(p, x), "k-")
plt.xlabel(r"energy limit, E$_\tau$ (GeV)")
plt.ylabel(r"ratio")
plt.savefig("tau-ratio.png")

plot_histogram(theta, weight, generated)
plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-zenith.png")

plot_histogram(phi, weight, generated)
plt.xticks(numpy.linspace(-180., 180., 5))
plt.axis((-200., 200., 0., 4E-03))
plt.yticks(numpy.linspace(0., 4E-03, 5))

plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-azimuth.png")

plot_histogram(altitude, weight, generated, plot=plt.semilogy)
plt.axis((0., 10., 1E-03, 1E+00))
plt.xlabel(r"decay altitude, h$_\tau$ (km)")
plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
plt.savefig("tau-altitude.png")

plt.figure()
plt.plot(position[:, 0], position[:, 1], "k.")
dx, dy = 0.5 * 150.4, 0.5 * 66.5
plt.plot((-dx, -dx, dx, dx, -dx), (-dy, dy, dy, -dy, -dy), "w--")
plt.xlabel(r"local x (km)")
plt.ylabel(r"local y (km)")
plt.savefig("tau-xy.png")

plt.show()
