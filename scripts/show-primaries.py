#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
import matplotlib.pyplot as plt


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-08 / (3. * e**2)

# Load the reduced events
with open("share/hotspot-150x67km2.primaries.p", "rb") as f:
    n, data = pickle.load(f)
data = numpy.array(data)

# Print the total rate
mu = sum(data[:,0]) / n
sigma = sum(data[:,0]**2) / n
sigma = numpy.sqrt((sigma - mu**2) / n)
print "Rate = {:.3f} +-{:.3f} a^-1".format(mu, sigma)


# Bin the primaries
def compute_spectrum(samples, weight, generated):
    """Estimate the primaries differential rate by binning"""
    lnx, w = numpy.log(samples), weight / generated
    n, b = numpy.histogram(lnx, 40, weights = w)
    n2, _ = numpy.histogram(lnx, b, weights=weight * w)
    norm = 1. / (b[1] - b[0])
    x = numpy.exp(0.5 * (b[1:] + b[:-1]))
    norm /= x
    b = numpy.exp(b)
    xerr = (x - b[:-1], b[1:] - x)
    p = n * norm
    dp = numpy.sqrt((n2 - n * n) / generated) * norm
    return x, p, xerr, dp


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")

plt.figure()
x, p, dx, dp = compute_spectrum(data[:, 1], data[:, 0], n)
plt.loglog(x, p * x, "k.")
plt.errorbar(x, p * x, xerr=dx, yerr=dp * x, fmt="k.")
plt.xlabel(r"energy, E$_\nu$ (GeV)")
plt.ylabel(r"E$_\nu \times$ rate (a$^{-1}$)")
plt.axis((1E+07, 1E+12, 1E-03, 1E+01))
plt.savefig("primaries-rate.png")

plt.figure()
norm = 1. / primary_flux(x)
y, dy = p * norm, dp * norm
plt.loglog(x, y, "k.")
plt.errorbar(x, y, xerr=dx, yerr=dy, fmt="k.")
plt.xlabel(r"energy, E$_\nu$ (GeV)")
plt.ylabel(r"exposure, 1 year (cm$^2$ s sr)")
plt.axis((1E+07, 1E+12, 1E+15, 1E+19))
plt.savefig("exposure.png")

plt.show()
