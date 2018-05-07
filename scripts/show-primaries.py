#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
import matplotlib.pyplot as plt

#plt.style.use("deps/mplstyle-l3/style/l3.mplstyle")


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

def doPlots(n,x, p, dx, dp,lin='-',col="k",leg=""):
    # Show the distributions
    plt.figure(1)
    plt.loglog(x, p * x, linestyle=lin,color=col,lw=4,label=leg)
    #plt.errorbar(x, p * x, xerr=dx, yerr=dp * x, fmt=col+'+')
    plt.xlabel(r"energy, E$_\nu$ (GeV)")
    plt.ylabel(r"E$_\nu \times$ rate (yr$^{-1}$)")
    plt.axis((3E+07, 1E+12, 1E-03, 1E+01))

    plt.figure(2)
    norm = 1. / primary_flux(x)
    y, dy = p * norm, dp * norm
    plt.loglog(x, y, color=col,linestyle=lin,lw=4,label=leg)
    #plt.errorbar(x, y, xerr=dx, yerr=dy, fmt=col+'+')
    plt.xlabel(r"energy, E$_\nu$ (GeV)")
    plt.ylabel(r"exposure, 1 year (cm$^2$ s sr)")
    plt.axis((3E+07, 1E+12, 1E+15, 1E+18))
    plt.grid(True)

    exi = numpy.array(range(70,120,1))/10.  # Conformed energy set
    xi = numpy.power(10,exi)
    yi = numpy.interp(xi,x,y) # Interpolated values for the comformed energy set
    
    return xi,yi
    
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

#files = ["share/old/HS1_sel1000.primaries.p.5ants.3s","share/old/HS1_sel1000.primaries.p.8ants.10s",]
files = ["/home/martineau/GRAND/GRAND/data/massProd/HS1/HS1freespace/noatt.primaries.p","share/HS1ground.primaries.5ants.5s","share/HS1freespace.primaries.5ants.5s","share/HS1freespace.primaries.5ants.5s.noAtt","share/HS1ground.primaries.5ants.3s","share/HS1freespace.primaries.5ants.3s","share/HS1freespace.primaries.5ants.3s.noAtt"]
col = ["k","r","b","g","r","b","g"]
lin = ["-","-","-","-","--","--","--"]
leg = ['Cone selection','Ground.+5$\sigma$','Free space + Att +5$\sigma$','Free space NoAtt +5$\sigma$','Ground.+3$\sigma$','Free space + Att +3$\sigma$','Free space NoAtt +3$\sigma$']
yv = numpy.zeros([len(files),50])
for i in range(len(files)):
    print "Loading",files[i]
    n, data = getRate(files[i])
    x, p, dx, dp = compute_spectrum(data[:, 1], data[:, 0], n)
    xi,yv[i,:] = doPlots(n,x, p, dx, dp,lin[i],col[i],leg[i])   
    plt.legend(loc='best')
    
eprel = 1e11*numpy.array([0.0010,0.0030,0.0100,0.0300,0.1000,0.3000,1.0000,3.0000])
expprea = 1e17*numpy.array([0.0260,0.2619,0.9504,1.7266,2.5552,3.4848,4.5618,5.2054])
expprec = 1e17 *numpy.array([0,0.0769,0.5324,1.2734,2.0508,2.7981,4.1219,5.0788])
exppreai = numpy.interp(xi,eprel,expprea)
exppreci = numpy.interp(xi,eprel,expprec)
#exppreai = exppreai*10000/7500  # Scale from 7500 to 1000km2 
#exppreci = exppreci*10000/7500  # Scale from 7500 to 1000km2

plt.figure(1)
plt.legend(loc='best')
plt.savefig("primaries-rate.png")
plt.figure(2)
plt.loglog(xi,exppreai,'r-.',lw=3,label="Agr. (prelim)")
plt.loglog(xi,exppreci,'g-.',lw=3,label="Cons. (prelim)")
plt.legend(loc='best')
plt.savefig("exposure.png")

#ra500 = yv[1,:]/yv[3,:]  # 500/1000
#ra1500 = yv[5,:]/yv[3,:]  # 1500/1000
#rc500 = yv[2,:]/yv[4,:]  # 500/1000
#rc1500 = yv[6,:]/yv[4,:]  # 1500/1000
#rac= yv[4,:]/yv[3,:]  # Conservative to Agressive
#rpa = exppreai/yv[3,:]  # Agressive prel to Agressive now
#rca = exppreci/yv[4,:]  # Conservative prel to Conservative now

#plt.figure(21)
#plt.loglog(xi, yv[3,:], color="r",lw=4,label='Agr. (1000m) ')
#plt.loglog(xi, yv[4,:], color="b",lw=4,label='Cons. (1000m)')
#plt.loglog(xi,exppreai,'r-.',lw=3,label="Prelim agr. (800m)")
#plt.loglog(xi,exppreci,'b-.',lw=3,label="Prelim cons. (800m)")
#plt.legend(loc='best')
#plt.xlabel(r"energy, E$_\nu$ (GeV)")
#plt.ylabel(r"exposure, 1 year (cm$^2$ s sr)")
#plt.axis((1E+08, 1E+12, 1E+14, 1E+18))
#plt.grid(True)
    
#plt.figure(3)
#plt.semilogx(xi,ra500,':r',lw=3,label='500m/1000m (Agr.)')
#plt.semilogx(xi,ra1500,'--r',lw=3,label='1500m/1000m (Agr.)')
#plt.semilogx(xi,rc500,':b',lw=3,label='500m/1000m (Cons.)')
#plt.semilogx(xi,rc1500,'--b',lw=3,label='1500m/1000m (Cons.)')
#plt.grid(True)
#plt.legend(loc='best')
#plt.axis((1E+08, 1E+12, 0, 4))
#plt.xlabel(r"energy, E$_\nu$ (GeV)")
#plt.ylabel(r"Exposure ratio")

#plt.figure(31)
#plt.semilogx(xi,rac,'k-',lw=3)
#plt.grid(True)
#plt.axis((1E+08, 1E+12, 0, 1))
#plt.xlabel(r"energy, E$_\nu$ (GeV)")
#plt.ylabel(r"Exposure ratio Cons./Agr.")

#plt.figure(32)
#plt.semilogx(xi,rpa,'r-',lw=3,label='Prelim agr./Agr.')
#plt.semilogx(xi,rca,'b-',lw=3,label='Prelim cons./Cons')
#plt.grid(True)
#plt.legend(loc='best')
#plt.axis((1E+08, 1E+12, 0, 4))
#plt.xlabel(r"energy, E$_\nu$ (GeV)")
#plt.ylabel(r"Exposure ratio")
    
plt.show()

