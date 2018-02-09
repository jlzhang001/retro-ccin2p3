#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
from scipy.integrate import cumtrapz
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from grand_tour import Topography


def plot_histogram(samples, weight, generated, plot=plt.plot, col="k",figID=1, factor=None):
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
    
    plt.figure(figID)
    if factor is None:
        y = p
        yerr = dp
    else:
        factor = x**factor
        y = p * factor
        yerr = dp * factor
    if sum(y) > 0:
        plot(x, y, col + "o")
        plt.errorbar(x, y, xerr=xerr, yerr=yerr, fmt=col + "o")
    return x, p


# Show the distributions
plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")
def doPlots(n, data,col="k"):
  x, p = plot_histogram(data[:, 1], data[:, 0],n,figID=1,plot=plt.loglog,col=col)
  plt.xlabel(r"energy, E$_\tau$ (GeV)")
  plt.ylabel(r"E$_\tau \times$ rate (yr$^{-1}$)")
  plt.savefig("tau-energy.png")

  plt.figure(2)
  plt.semilogx(x, cumtrapz(p, x, initial=0.) / numpy.trapz(p, x),color=col)
  plt.xlabel(r"energy limit, E$_\tau$ (GeV)")
  plt.ylabel(r"ratio")
  plt.savefig("tau-ratio.png")

  plot_histogram(data[:, 8], data[:, 0],n,col=col,figID=3)
  plt.axis((85., 95., 0., 6.))
  plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
  plt.ylabel(r"rate (deg$^{-1}$ yr$^{-1}$)")
  plt.savefig("tau-zenith.png")

  plot_histogram(data[:, 9], data[:, 0], n,col=col,figID=4)
  plt.xticks(numpy.linspace(-180., 180., 5))
  plt.axis((-200., 200., 0., 1E-01))
  plt.yticks(numpy.linspace(0., 1E-01, 5))
  plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
  plt.ylabel(r"rate (deg$^{-1}$ yr$^{-1}$)")
  #plt.savefig("tau-azimuth.png")

  plot_histogram(data[:, 7] * 1E-03, data[:, 0], n,col=col,figID=5, plot=plt.semilogy)
  plt.axis((0., 5., 1E-02, 1E+01))
  plt.xlabel(r"decay altitude, h$_\tau$ (km)")
  plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
  plt.savefig("tau-altitude.png")


def doTopoPlot(data):
  rate, xe, ye = histogram2d(data[:, 6], data[:, 5], 100,
  			     normed=True, weights=data[:, 0])
  rate *= mu
  x = 0.5 * (xe[1:] + xe[:-1])
  y = 0.5 * (ye[1:] + ye[:-1])
 
  # Compute the topography
  topo = Topography(latitude=42.1, longitude=86.3, path="share/topography",
  		    stack_size=49)
  h = numpy.zeros(rate.shape)
  for i, yi in enumerate(y):
      for j, xj in enumerate(x):
  	  h[i, j] = topo.ground_altitude(yi, xj, geodetic=True)
 
  plt.figure()
  plt.contour(x, y, h, 10, colors="k", alpha=0.75)
  plt.pcolor(x, y, rate, cmap="nipy_spectral", norm=LogNorm(), alpha=0.5)
  dx, dy = 0.5 * 66.5E+03, 0.5 * 150.4E+03
  y0, x0, _ = topo.local_to_lla((-dx, -dy, 0.))
  y1, x1, _ = topo.local_to_lla((dx, dy, 0.))
  plt.plot((x0, x0, x1, x1, x0), (y0, y1, y1, y0, y0), "w-")
  plt.colorbar()
  plt.xlabel(r"longitude (deg)")
  plt.ylabel(r"latitude (deg)")
  plt.title(r"$\tau$ rate (yr$^{-1}$ / deg$^2$)")
  plt.savefig("tau-rate-map.png")
  
# Estimate the decay density
def histogram2d(*args, **kwargs):
    """Encapsulate numpy.histogram2d for matrix convention compatibility"""
    n, x, y = numpy.histogram2d(*args, **kwargs)
    n = n.T
    return n, x, y


col = ["k","r","g"]
files = ["share/HS1.p", "share/HS1_sel.p.5ants.4s"]
for i in range(len(files)):
# Load the reduced events
  print"Opening",files[i]
  with open(files[i], "rb") as f:
    n, data = pickle.load(f)
  data = numpy.array(data)

  # Print the total rate  
  mu = sum(data[:,0]) / n 
  sigma = sum(data[:,0]**2) / n
  sigma = numpy.sqrt((sigma - mu**2) / n)
  print "Rate = {:.3f} +-{:.3f} yr^-1".format(mu, sigma)
  doPlots(n,data,col[i])
  doTopoPlot(data)

plt.show()


 

