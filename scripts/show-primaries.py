#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
import matplotlib.pyplot as plt

plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")


def primary_flux(e):
    """Waxman-Bahcall bound with 1 / 3 of tau neutrinos"""
    return 2E-08 / (3. * e**2)

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

def doPlots(n,x, p, dx, dp,col="k",leg=""):
    # Show the distributions
    plt.figure(1)
    plt.loglog(x, p * x, col)
    plt.errorbar(x, p * x, xerr=dx, yerr=dp * x, fmt=col,label=leg)
    plt.xlabel(r"energy, E$_\nu$ (GeV)")
    plt.ylabel(r"E$_\nu \times$ rate (yr$^{-1}$)")
    plt.axis((1E+07, 1E+12, 1E-03, 1E+01))

    plt.figure(2)
    norm = 1. / primary_flux(x)
    y, dy = p * norm, dp * norm
    plt.loglog(x, y, col)
    plt.errorbar(x, y, xerr=dx, yerr=dy, fmt=col)
    plt.xlabel(r"energy, E$_\nu$ (GeV)")
    plt.ylabel(r"exposure, 1 year (cm$^2$ s sr)")
    plt.axis((1E+07, 1E+12, 1E+15, 1E+19))

def getRate(pfile):
    # Load the reduced events
    with open(pfile, "rb") as f:
    	n, data = pickle.load(f)
    data = numpy.array(data)

    # Print the total rate
    mu = sum(data[:,0]) / n
    sigma = sum(data[:,0]**2) / n
    sigma = numpy.sqrt((sigma - mu**2) / n)
    print "Rate = {:.3f} +-{:.3f} yr^-1".format(mu, sigma)

    return n, data

files = ["share/HS1.primaries.p","share/HS1_sel.primaries.p.5ants.4s","share/HS1_sel.primaries.p.8ants.10s",]
col = ["k","r","g","b"]
leg = ['Cone selection','Ag.: 5 ants>60$\muV$pp','Cons: 8 ants>150$\muV$pp']
for i in range(len(files)):
    print "Loading",files[i]
    n, data = getRate(files[i])
    x, p, dx, dp = compute_spectrum(data[:, 1], data[:, 0], n)
    doPlots(n,x, p, dx, dp,col[i],leg=leg[i])

eprel = 1e11*numpy.array([0.0010,0.0030,0.0100,0.0300,0.1000,0.3000,1.0000,3.0000])
exposureprel = 1e17*numpy.array([0.0260,0.2619,0.9504,1.7266,2.5552,3.4848,4.5618,5.2054])
exposureprelc = 1e+17 *numpy.array([0,0.0769,0.5324,1.2734,2.0508,2.7981,4.1219,5.0788])
exposureprel = exposureprel*10000/7500  # Scale from 7500 to 1000km2
exposureprels = exposureprelc*10000/7500  # Scale from 7500 to 1000km2

plt.figure(1)
plt.savefig("primaries-rate.png")
plt.figure(2)
plt.loglog(eprel,exposureprel,'r--',label="Prelim, ag.")
plt.loglog(eprel,exposureprelc,'g--',label="Prelim, cons.")
#plt.legend('loc=best')
plt.savefig("exposure.png")

plt.show()
