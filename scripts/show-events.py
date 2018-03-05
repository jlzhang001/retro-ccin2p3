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
            K = [numpy.linalg.norm(v[3]) > 6371E+03
                 for v in event["primaries"][0]]
            J = [not k for k in K]
            K, J = map(lambda X: numpy.nonzero(X)[0], (K, J))
            wp = numpy.array([v[0] for v in event["primaries"][0]])
            norm = event["statistics"][0] * PHI0 * year / event["primaries"][1]
            wE = sum(wp[J]) * norm
            wA = sum(wp[K]) * norm
            generated += event["statistics"][1]
            if wE + wA == 0:
                continue
            theta.append(t)
            phi.append(p)
            altitude.append(a * 1E-03)
            position.append(r)
            energy.append(e)
            weight.append((wE, wA))
    position = numpy.array(position) * 1E-03
    weight = numpy.array(weight)
    def get_stats(w, n):
        mu = sum(w) / n
        sigma = sum(w**2) / n
        sigma = numpy.sqrt((sigma - mu**2) / n)
        return mu, sigma
    mE, sE = get_stats(weight[:,0], generated)
    mA, sA = get_stats(weight[:,1], generated)
    print "Earth skimming = {:.3f} +-{:.3f} a^-1".format(mE, sE)
    print "Atmospheric    = {:.3f} +-{:.3f} a^-1".format(mA, sA)
    print "Events ratio   = {:.2f}".format(
        sum(weight[:,0] > 0.) / float(len(weight)))

    if mE <= 0.:
        # Bring the atmospheric events to the front
        weight[:,[0, 1]] = weight[:,[1, 0]]

    return theta, phi, altitude, position, energy, weight, generated

theta, phi, altitude, position, energy, weight, generated = load(DATA_DIR)

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

x, p = plot_histogram(energy, weight[:, 0], generated, plot=plt.loglog)
x2, p2 = plot_histogram(energy, weight[:, 1], generated, plot=plt.loglog,
                        new=False, clr="r")
earthskimming = sum(p2) > 0.

plt.xlabel(r"energy, E$_\tau$ (GeV)")
plt.ylabel(r"E$_\tau \times$ rate (a$^{-1}$)")
plt.savefig("tau-energy.png")

plt.figure()
plt.semilogx(x, cumtrapz(p, x, initial=0.) / numpy.trapz(p, x), "k-")
if earthskimming:
    plt.semilogx(x2, cumtrapz(p2, x2, initial=0.) / numpy.trapz(p2, x2), "r-")
plt.xlabel(r"energy limit, E$_\tau$ (GeV)")
plt.ylabel(r"ratio")
plt.savefig("tau-ratio.png")

plot_histogram(theta, weight[:, 0], generated)
plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-zenith.png")

plot_histogram(phi, weight[:, 0], generated)
plt.xticks(numpy.linspace(-180., 180., 5))
if earthskimming:
    plt.axis((-200., 200., 0., 2E-02))
    plt.yticks(numpy.linspace(0., 2E-02, 5))
else:
    plt.axis((-200., 200., 0., 4E-03))
    plt.yticks(numpy.linspace(0., 4E-03, 5))

plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-azimuth.png")

plot_histogram(altitude, weight[:, 0], generated, plot=plt.semilogy)
plot_histogram(altitude, weight[:, 1], generated, plot=plt.semilogy, new=False,
               clr="r")
if earthskimming:
    plt.axis((0., 3., 1E-03, 1E+02))
else:
    plt.axis((0., 10., 1E-03, 1E+00))
plt.xlabel(r"decay altitude, h$_\tau$ (km)")
plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
plt.savefig("tau-altitude.png")

plt.figure()
plt.plot(position[:, 0], position[:, 1], "k.")
plt.plot((-50., -50., 50., 50., -50.), (-50., 50., 50., -50., -50.), "w--")
plt.xlabel(r"local x (km)")
plt.ylabel(r"local y (km)")
plt.savefig("tau-xy.png")

plt.show()
