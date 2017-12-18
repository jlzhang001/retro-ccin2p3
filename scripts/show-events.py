#!/usr/bin/env python

import os
import numpy
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from retro.event import EventIterator
from grand_tour import Topography


# Settings
DATA_DIR = "share/events"
PHI0 = 2. / 3. * 1E-04


# Load and parse the events
year = 365.25 * 24. * 60. * 60.
topo = Topography(43, 87, "flat/10")

def load(path):
    """Load and parse a set of decay files
    """
    theta, phi, altitude = [], [], []
    position, energy, weight = [], [], []
    generated = 0
    for name in os.listdir(path):
        if not name.startswith("events"):
            continue
        filename = os.path.join(path, name)
        for event in EventIterator(filename):
            e, r, u = event["tau_at_decay"]
            t, p = topo.local_to_angular(r, u)
            _, _, a = topo.local_to_lla(r)
            w = sum(v[0] for v in event["primaries"][0]) / event["primaries"][1]
            w *= event["statistics"][0] * PHI0 * year
            generated += event["statistics"][1]
            if t < 90.: continue
            theta.append(t)
            phi.append(p)
            altitude.append(a * 1E-03)
            position.append(r)
            energy.append(e)
            weight.append(w)
    position = numpy.array(position) * 1E-03
    weight = numpy.array(weight)
    mu = sum(weight) / generated
    sigma = sum(weight**2) / generated
    sigma = numpy.sqrt((sigma - mu**2) / generated)
    print "Total rate = {:.3f} +-{:.3f} a^-1".format(mu, sigma)

    return theta, phi, altitude, position, energy, weight, generated

theta, phi, altitude, position, energy, weight, generated = load(DATA_DIR)

def plot_histogram(samples, weight, generated, log=False, clr="k", new=True):
    """Plot a 1D histogram of sampled values
    """
    if log:
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
    if log:
        plt.loglog(x, p * x, clr + "o")
        plt.errorbar(x, p * x, xerr=xerr, yerr=dp * x, fmt=clr + "o")
    else:
        plt.plot(x, p, clr + "o")
        plt.errorbar(x, p, xerr=xerr, yerr=dp, fmt=clr + "o")
    return x, p


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

x, p = plot_histogram(energy, weight, generated, log=True)
plt.xlabel(r"energy, E$_\tau$ (GeV)")
plt.ylabel(r"E$_\tau \times$ rate (a$^{-1}$)")
plt.savefig("tau-energy.png")

plt.figure()
plt.semilogx(x, cumtrapz(p, x, initial=0.), "k-")
plt.xlabel(r"energy limit, E$_\tau$ (GeV)")
plt.ylabel(r"rate (a$^{-1}$)")
plt.savefig("tau-rate.png")

plot_histogram(theta, weight, generated)
plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-zenith.png")

plot_histogram(phi, weight, generated)
plt.axis((-200., 200., 0., 8E-04))
plt.xticks(numpy.linspace(-180., 180., 5))
plt.yticks(numpy.linspace(0., 8E-04, 5))
plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-azimuth.png")

plot_histogram(altitude, weight, generated)
plt.axis((0., 3., 0., 0.35))
plt.xlabel(r"decay altitude, h$_\tau$ (km)")
plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
plt.savefig("tau-altitude.png")

plt.figure()
plt.plot(position[:, 0], position[:, 1], "k.")
plt.xlabel(r"local x (km)")
plt.ylabel(r"local y (km)")
plt.savefig("tau-xy.png")

plt.show()
