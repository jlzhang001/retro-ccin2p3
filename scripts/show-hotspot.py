#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt

# Load the reduced events
with open("share/events.p", "rb") as f:
    n, data = pickle.load(f)
data = numpy.array(data)

# Print the total rate
mu = sum(data[:,0]) / n
sigma = sum(data[:,0]**2) / n
sigma = numpy.sqrt((sigma - mu**2) / n)
print "Rate = {:.3f} +-{:.3f} a^-1".format(mu, sigma)

def plot_histogram(samples, weight, generated, plot=plt.plot, clr="k",
                   new=True, factor=None):
    """Plot a 1D histogram of sampled values"""
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

x, p = plot_histogram(data[:, 1], data[:, 0], n, plot=plt.loglog)
plt.xlabel(r"energy, E$_\tau$ (GeV)")
plt.ylabel(r"E$_\tau \times$ rate (a$^{-1}$)")
plt.savefig("tau-energy.png")

plt.figure()
plt.semilogx(x, cumtrapz(p, x, initial=0.) / numpy.trapz(p, x), "k-")
plt.xlabel(r"energy limit, E$_\tau$ (GeV)")
plt.ylabel(r"ratio")
plt.savefig("tau-ratio.png")

plot_histogram(data[:, 8], data[:, 0], n)
plt.axis((85., 95., 0., 3.))
plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-zenith.png")

plot_histogram(data[:, 9], data[:, 0], n)
plt.xticks(numpy.linspace(-180., 180., 5))
plt.axis((-200., 200., 0., 2E-02))
plt.yticks(numpy.linspace(0., 2E-02, 5))
plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
plt.ylabel(r"rate (deg$^{-1}$ a$^{-1}$)")
plt.savefig("tau-azimuth.png")

plot_histogram(data[:, 7] * 1E-03, data[:, 0], n, plot=plt.semilogy)
plt.axis((0., 5., 1E-02, 1E+01))
plt.xlabel(r"decay altitude, h$_\tau$ (km)")
plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
plt.savefig("tau-altitude.png")

plt.figure()
plt.plot(data[:, 2] * 1E-03, data[:, 3] * 1E-03, "k.")
dx, dy = 0.5 * 150.4, 0.5 * 66.5
plt.plot((-dx, -dx, dx, dx, -dx), (-dy, dy, dy, -dy, -dy), "w--")
plt.xlabel(r"local x (km)")
plt.ylabel(r"local y (km)")
plt.savefig("tau-xy.png")

plt.show()
