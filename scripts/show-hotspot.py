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
from scipy import stats
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
  #plt.axis((85., 95., 0., 6.))
  plt.xlim(85., 95.)
  plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
  plt.ylabel(r"rate (deg$^{-1}$ yr$^{-1}$)")
  plt.savefig("tau-zenith.png")

  plot_histogram(data[:, 9], data[:, 0], n,col=col,figID=4)
  plt.xticks(numpy.linspace(-180., 180., 5))
  #plt.axis((-200., 200., 0., 1E-01))
  #plt.yticks(numpy.linspace(0., 1E-01, 5))
  plt.xlim(-180., 180.)
  plt.xlabel(r"azimuth, $\phi_\tau$ (deg)")
  plt.ylabel(r"rate (deg$^{-1}$ yr$^{-1}$)")
  #plt.savefig("tau-azimuth.png")

  plot_histogram(data[:, 7], data[:, 0], n,col=col,figID=5)
  #plt.axis((0., 5., 1E-02, 1E+01))
  plt.xlabel(r"decay altitude, h$_\tau$ (m)")
  plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
  plt.savefig("tau-altitude.png")

  plot_histogram(data[:, 10], data[:, 0], n,col=col,figID=6)
  plot_histogram(data[:, 11], data[:, 0], n,col='r',figID=6)
  plt.xlabel(r"decay height above ground (m)")
  plt.ylabel(r"rate (km$^{-1}$ a$^{-1}$)")
  plt.savefig("tau-height.png")
  

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
  plt.plot((x0, x0, x1, x1, x0), (y0, y1, y1, y0, y0), "w-",lw=5)
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
#files = ["share/HS1.p", "share/HS1_sel1000.p.5ants.4s"]
#files = ["share/HS1_sel500.p.5ants.4s", "share/HS1_sel1000.p.5ants.4s"]
#files = ["share/HS1_sel1000.p.5ants.3s"]
#files = ["share/HS1_sel1000.p.5ants.3s","share/HS1flat_sel1000.p.5ants.3s"]
#files = ["share/flat_sel1000.p.5ants.3s"]

files = ["/home/martineau/GRAND/GRAND/data/massProd/HS1/HS1_sel.p"]
#files = ["/home/martineau/GRAND/GRAND/data/massProd/HS1flat/HS1flat_sel.p"]
for i in range(len(files)):
# Load the reduced events
  print"Opening",files[i]
  with open(files[i], "rb") as f:
    n, data = pickle.load(f)
  data = numpy.array(data)
  theta = data[:, 8]
  ncones = data[:,12]/4 # Scale down by 4 for 500-->1000m step
  nvolts = data[:,13]/4 # Scale down by 4 for 500-->1000m step
  ntrigs = data[:,14]
  nclust = data[:,15]
  dmin = data[:,16]
  dmax = data[:,17]
  
  plt.figure(17)
  plt.subplot(221)
  plt.hist(dmin[dmin>0],100)
  plt.xlabel('Dmin (km)')
  plt.subplot(222)
  plt.plot(theta[dmin>0],dmin[dmin>0],'ob')
  plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
  plt.ylabel('Dmin (km)')
  plt.subplot(223)
  plt.hist(dmax[dmax>0],100)
  plt.xlabel('Dmax (km)')
  plt.subplot(224)
  plt.plot(theta[dmin>0],dmax[dmax>0],'ob')
  plt.xlabel(r"zenith, $\theta_\tau$ (deg)")
  plt.ylabel('Dmax (km)')
  
  print "Nb events with antennas further than 90km:",sum(dmin>90)
  
  plt.figure(23)
  plt.title('$N_{ants}$ in events (with 1+ voltage')
  plt.subplot(411)
  plt.hist(numpy.log10(ncones[ncones>0]),100)
  #plt.xlim(0.1,max(numpy.log10(ncones)))
  plt.xlabel(r"$log_{10}(N_{Cone}^*)$")
  print "\n**** All stats only for events with 1+ voltage***"
  print "Cones stats: [N,minAnts,maxAnts,<Ants>]",len(ncones),numpy.min(ncones),numpy.max(ncones),numpy.mean(ncones)
  #stats.describe(ncones)
  plt.subplot(412)
  plt.hist(numpy.log10(nvolts[nvolts>0]),100)
  #plt.xlim(0.1,max(numpy.log10(ncones)))
  plt.xlabel(r"$log_{10}(N_{RadioSim}^*)$")
  print "Radio sim stats: [N,minAnts,maxAnts,<Ants>]",len(nvolts[nvolts>0]),numpy.min(nvolts[nvolts>0]),numpy.max(nvolts[nvolts>0]),numpy.mean(nvolts[nvolts>0])
  print "Radio sim stats (events with 5+ ants): [N,minAnts,maxAnts,<Ants>]",len(nvolts[nvolts>4]),numpy.min(nvolts[nvolts>4]),numpy.max(nvolts[nvolts>4]),numpy.mean(nvolts[nvolts>4])
  
  plt.subplot(413)
  plt.hist(numpy.log10(ntrigs[ntrigs>0]),100)
  #plt.xlim(0.1,max(numpy.log10(ncones)))
  plt.xlabel(r"$log_{10}(N_{Trig})$")
  print "Trigged events stats: [N,minAnts,maxAnts,<Ants>]",len(ntrigs[ntrigs>0]),numpy.min(ntrigs[ntrigs>0]),numpy.max(ntrigs[ntrigs>0]),numpy.mean(ntrigs[ntrigs>0])
  print "Trigged events stats (events with 5+ ants): [N,minAnts,maxAnts,<Ants>]",len(ntrigs[ntrigs>4]),numpy.min(ntrigs[ntrigs>4]),numpy.max(ntrigs[ntrigs>4]),numpy.mean(ntrigs[ntrigs>4])
    
  plt.subplot(414)
  plt.hist(numpy.log10(nclust[nclust>0]),100)
  #plt.xlim(0.1,max(numpy.log10(ncones)))
  plt.xlabel(r"$log_{10}(N_{Clust})$")
  print "Clustered events stats: [N,minAnts,maxAnts,<Ants>]",len(nclust[nclust>0]),numpy.min(nclust[nclust>0]),numpy.max(nclust[nclust>0]),numpy.mean(nclust[nclust>0])
  print "Clustered events stats (events with 5+ ants): [N,minAnts,maxAnts,<Ants>]",len(nclust[nclust>4]),numpy.min(nclust[nclust>4]),numpy.max(nclust[nclust>4]),numpy.mean(nclust[nclust>4])

  # Print the total rate  
  mu = sum(data[:,0]) / n 
  sigma = sum(data[:,0]**2) / n
  sigma = numpy.sqrt((sigma - mu**2) / n)
  print "Rate = {:.3f} +-{:.3f} yr^-1".format(mu, sigma)
  doPlots(n,data,col[i])
  doTopoPlot(data)

plt.show()


 

