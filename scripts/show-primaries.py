#!/usr/bin/env python

# Prune obsolete site packages at CC
import sys
exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]

import cPickle as pickle

import numpy
import matplotlib.pyplot as plt
DISPLAY = 1
dur = 3 # years
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
    plt.xlabel(r"E$_\nu$ (GeV)")
    plt.ylabel(r"1-year exposure (cm$^2$ s sr)")
    plt.axis((3E+07, 1E+12, 1E+15, 5E+18))
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

if DISPLAY == 1:  # HS1 final
  files = ["share/HS1ground.primaries.5ants.2s.full.090818","share/HS1ground.primaries.5ants.2s.f50200.090818","share/HS1ground.primaries.5ants.3s.f50200.090818","share/HS1ground.primaries.5ants.5s.f50200.090818"]
  leg = ['Ground Full Aug 08 (2$\sigma$ th)','Ground 50-200MHz Aug 08 (2$\sigma$ th)','Ground 50-200MHz Aug 08 (3$\sigma$ th)','Ground 50-200MHz Aug 08 (5$\sigma$ th)']
  col = ["b","b",'magenta','cyan']
  lin = ["-","--",'-','-']
  
if DISPLAY == 2:  # Free space vs ground
  files = ["share/HS1freespace.primaries.5ants.2s.v2.50200","share/HS1freespace.primaries.5ants.2s.noAtt.50200","share/HS1ground.primaries.5ants.2s.50200","share/HS1ground.primaries.5ants.2s.full.090818","share/HS1ground.primaries.5ants.2s.f50200.090818","share/HS1ground.primaries.5ants.3s.f50200.090818","share/HS1ground.primaries.5ants.5s.f50200.090818"]
  col = ["g","g","b","k","k","r","r"]
  lin = ["-","--","-","--","-","-","--"]
  leg = ['Free space + att 50-200MHz May 21 (2$\sigma$ th)' ,'Free space + No att 50-200MHz May 21 (2$\sigma$ th)','Ground 50-200MHz May 21(2$\sigma$ th)','Ground Full Aug 08 (2$\sigma$ th)','Ground 50-200MHz Aug 08 (2$\sigma$ th)','Ground 50-200MHz Aug 08 (3$\sigma$ th)','Ground 50-200MHz Aug 08 (5$\sigma$ th)']
  
if DISPLAY == 3:  # HS1 vs flat
  files = ["share/HS1freespace.primaries.5ants.2s.v2.50200","share/HS1freespace.primaries.5ants.5s.v2.50200","share/flat_freespace.primaries.5ants.2s.v2.50200","share/HS1freespace.primaries.5ants.2s.v2.50200.500m"]
  col = ["r","g","brown",'r']
  lin = ["-","-","-","--"]
  leg = ['HS1 (2$\sigma$ th)','HS1 (5$\sigma$ th)','Flat (2$\sigma$ th)','HS1 500m (2$\sigma$ th)']
if DISPLAY == 4:  # HS1 frequency study
  files = ["share/HS1freespace.primaries.5ants.2s.v2.50200","share/HS1freespace.primaries.5ants.2s.v2.70150","share/HS1freespace.primaries.5ants.2s.v2.5090","share/HS1freespace.primaries.5ants.2s.v2.150200"]
  #files = ["share/HS1freespace.primaries.5ants.2s.v2.50200","share/HS1freespace.primaries.5ants.2s.v2.70150","share/HS1ground.primaries.5ants.2s.50200","share/jsons_att_v2f70150.primaries_sel.p"]
  col = ["g","b","r","k"]
  lin = ["-","-","-","-"]
  leg = ['HS1 50-200MHz (2$\sigma$ th)','HS1 70-150MHz (2$\sigma$ th)','HS1 Ground 50-90MHz (2$\sigma$ th)','HS1 Ground 150-200MHz (2$\sigma$ th)']
if DISPLAY == 5: # 60000km2 area
  #files = ["share/simu60kcone.primaries.5ants.agr.800m.v0","share/simu60kcone.primaries.5ants.agr.800m.v2","share/simu60kcone.primaries.5ants.noclust.agr.800m.v2","share/simu60kcone.primaries.5ants.agr.800m.v3","share/simu60kcone.primaries.5ants.noclust.agr.800m.v3","share/simu60kcone.primaries.5ants.agr.1200m.v3","share/simu60kcone.primaries.5ants.noclust.agr.1200m.v3","share/jsons.primaries_sel.p"]
  #col = ["k","b","b",'k','k','g','g',"b"]
  #lin = ["-","-","--","-","--",'-','--',":"]
  #leg = ['60k cone new 800m v0 (5 ants)','60k cone new 800m v2 (5 ants)','60k cone new 800m v2 - no cluster (5 ants)','60k cone new 800m v3 (5 ants)','60k cone new 800m v3 - no cluster (5 ants)','60k cone new 1200m v3 (5 ants)','60k cone new 1200m v3 - no cluster (5 ants)',"60k cone radio (agr)"]
  files = ["share/simu60kcone.primaries.5ants.agr.800m.v0","share/simu60kcone.primaries.5ants.agr.800m.v3","share/simu60kcone.primaries.5ants.noclust.agr.800m.v3","share/simu60.primaries.5ants.agr.v2106","share/simu60.primaries.5ants.agr.v0807","share/simu60kconepresel.primaries.5ants.agr.1000m.09082018","share/simu60.primaries.5ants.2s.full.090818","share/simu60.primaries.5ants.2s.f50200.1000m.090818","share/simu60.primaries.5ants.3s.f50200.090818","share/simu60.primaries.5ants.5s.f50200.090818"]
  col = ["k","k","k","g","r","b","b","b",'magenta','cyan']
  lin = ["--","-",":","-",'-',":","-","--",'-','-','-']
  leg = ['60k cone ana 800m June 1st (5 ants)','60k cone ana 800m June 22nd (5 ants)','60k cone ana 800m June 22d (noclust)',"60k RM June 22nd (2$\sigma$ th)","60k RM July 8 (2$\sigma$ th)",'60k cone presel 1000m Aug 08 (5 ants)','60k RM Aug 08 (2$\sigma$ th)','60k RM 50-200MHz Aug 08 (2$\sigma$ th)','60k RM 50-200MHz Aug 08 (3$\sigma$ th)','60k RM 50-200MHz Aug 08 (5$\sigma$ th)']
  
yv = numpy.zeros([len(files),50])
for i in range(len(files)):
    print "Loading",files[i]
    n, data = getRate(files[i])
    x, p, dx, dp = compute_spectrum(data[:, 1], data[:, 0], n)
    xi,yv[i,:] = doPlots(n,x, p, dx, dp,lin[i],col[i],leg[i])   
    plt.legend(loc='best')
    
    # Compute integral limit
    y = p/primary_flux(x)
    lim = 2.44/(numpy.trapz(y/x,numpy.log(x))*dur)*3
    lim200k_allflavor = lim/20
    print "3 years all flavor limit=",lim,"GeV/cm2/s/sr"
    print "GRAND200k limit projection",lim200k_allflavor,"GeV/cm2/s/sr"
    print "New/ini (agr.)",lim200k_allflavor/2e-10 

    # Differential limit
    dlim = 2.44*x/(numpy.log(10)*y*dur)*3  # All flavor
    dlim200k = dlim/20  # 
    plt.figure(17)
#    if i == 0:
      #plt.loglog(x,dlim,linestyle="--",color="k",lw=4,label="60k cone new 800m June 1st (5 ants)")
#      plt.loglog(x,dlim200k,linestyle="-",color="orange",lw=4,label="GRAND200k (2$\sigma$ th)")
#      print x,dlim200k
    if i == 1:
#      plt.loglog(x,dlim200k,linestyle="-",color="green",lw=4,label="GRAND200k (5$\sigma$ th)")
      plt.loglog(x,dlim,linestyle="-",color="k",lw=4,label="60k cone new 800m June 22nd (5 ants)")
#    if i == 3:
#      plt.loglog(x,dlim,linestyle="-",color="b",lw=4,label="60k RM July 8 (agr)")
#      plt.loglog(x,dlim200k,linestyle="--",color="brown",lw=4,label="GRAND200k - flat(2$\sigma$ th)")
    if i == 4:
      plt.loglog(x,dlim,linestyle="-",color="b",lw=4,label="60k RM July 8 50-200MHz (agr)")
            
    if i == 4:
      eHS1ini = 1.0e+11 *numpy.array([0.0030, 0.0038, 0.0048, 0.0060, 0.0075, 0.0095, 0.0119, 0.0150, 0.0189, 0.0238, 0.0300, 0.0378, 0.0475, 0.0599, 0.0754, 0.0949, 0.1194, 0.1504, 0.1893, 0.2383, 0.3000, 0.3777, 0.4755, 0.5986, 0.7536, 0.9487, 1.1943, 1.5036, 1.8929, 2.3830, 3.0000, 3.7768, 4.7547, 5.9858, 7.5357, 9.4868])
      lHS1ini = 1.0e-05 *numpy.array([0.0012, 0.0011, 0.0010, 0.0010, 0.0010, 0.0011, 0.0012, 0.0013, 0.0014, 0.0016, 0.0018, 0.0021, 0.0025, 0.0029, 0.0034, 0.0040, 0.0047, 0.0055, 0.0065, 0.0077, 0.0092, 0.0109, 0.0129, 0.0154, 0.0185, 0.0222, 0.0268, 0.0326, 0.0399, 0.0491, 0.0611, 0.0767, 0.0973, 0.1249, 0.1623, 0.2136])
      eGRANDini = 1.0e+11 *numpy.array([0.0010, 0.0013, 0.0016, 0.0020, 0.0025,0.0032, 0.0040, 0.0050, 0.0063, 0.0079, 0.0100, 0.0126, 0.0158, 0.0200, 0.0251, 0.0316, 0.0398, 0.0501, 0.0631, 0.0794, 0.1000, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1.0000, 1.2589, 1.5849, 1.9953, 2.5119, 3.1623, 3.9811, 5.0119, 6.3096, 7.9433])
      lGRANDini = 1.0e-07 *numpy.array([0.0094, 0.0059, 0.0042, 0.0032, 0.0026, 0.0023, 0.0022, 0.0021, 0.0021, 0.0022, 0.0024, 0.0026, 0.0029, 0.0032, 0.0037, 0.0042, 0.0049, 0.0056, 0.0065, 0.0076, 0.0089, 0.0104, 0.0123, 0.0145, 0.0172, 0.0205, 0.0245, 0.0295, 0.0355, 0.0429, 0.0520, 0.0630, 0.0762, 0.0918, 0.1097, 0.1297, 0.1508, 0.1716, 0.1898, 0.2026]) 
      e60kini= 1.0e+11*numpy.array([0.001000000000000, 0.001258925411794, 0.001584893192461, 0.001995262314969, 0.002511886431510, 0.003162277660168, 0.003981071705535, 0.005011872336273, 0.006309573444802, 0.007943282347243, 0.010000000000000, 0.012589254117942, 0.015848931924611, 0.019952623149689, 0.025118864315096, 0.031622776601684, 0.039810717055350, 0.050118723362727,  0.063095734448019, 0.079432823472428, 0.100000000000000, 0.125892541179417, 0.158489319246111, 0.199526231496888, 0.251188643150958, 0.316227766016838, 0.398107170553497, 0.501187233627271, 0.630957344480194, 0.794328234724282, 1.000000000000000, 1.258925411794166, 1.584893192461117, 1.995262314968883, 2.511886431509582, 3.162277660168380, 3.981071705534969, 5.011872336272735, 6.309573444801943, 7.943282347242822])
      sim60kini=1.0e-06*numpy.array([0.012353678072464,  0.007830433632699,  0.005494362808308,  0.004209641738704,  0.003478524281274,  0.003065745021102,  0.002853329522561,  0.002779796648992,  0.002812857585562,  0.002936393702096,  0.003143925737675,  0.003435219141273,  0.003814505453464,  0.004289606922992,  0.004871627763942,  0.005575060772800,  0.006418254501583,  0.007424239685458,  0.008621946117694,  0.010047863791528,  0.011748220291098,  0.013781761312723,  0.016223229492577,  0.019167628884747,  0.022735319833885,  0.027077879350227,  0.032384433309571,  0.038887739990732,  0.046868569153439,  0.056655741297852,  0.068617440512125,  0.083137065812658,  0.100564213704777,  0.121129310199709,  0.144810992546161,  0.171152219781848,  0.199039307339552,  0.226492117959622,  0.250562200288043,  0.267482798841317])
      e60kini1200 = 1.0e+11*numpy.array([0.0010, 0.0013, 0.0016, 0.0020, 0.0025, 0.0032, 0.0040, 0.0050, 0.0063, 0.0079, 0.0100, 0.0126, 0.0158, 0.0200, 0.0251, 0.0316, 0.0398, 0.0501, 0.0631, 0.0794, 0.1000, 0.1259, 0.1585, 0.1995, 0.2512, 0.3162, 0.3981, 0.5012, 0.6310, 0.7943, 1.0000, 1.2589, 1.5849, 1.9953, 2.5119, 3.1623, 3.9811, 5.0119, 6.3096, 7.9433])
      sim60kini1200 = 1.0e-06*numpy.array([0.0384, 0.0203, 0.0123, 0.0083, 0.0062, 0.0050, 0.0044, 0.0041, 0.0040, 0.0040, 0.0042, 0.0045, 0.0049, 0.0055, 0.0061, 0.0070, 0.0080, 0.0091, 0.0105, 0.0122, 0.0141, 0.0164, 0.0191, 0.0224, 0.0263, 0.0311, 0.0369, 0.0439, 0.0526, 0.0632, 0.0761, 0.0917, 0.1104, 0.1323, 0.1573, 0.1846, 0.2126, 0.2388, 0.2595, 0.2704])


      if DISPLAY==5:
        plt.loglog(e60kini,sim60kini,linestyle='--',color="brown",lw=4,label="60k ini (agr)")
      else:
        plt.loglog(eHS1ini,lHS1ini,linestyle='--',color="b",lw=4,label="HS1 ini (agr)")
      #plt.loglog(e60kini1200,sim60kini1200,linestyle='--',color="r",lw=4,label="60k ini (agr) - 1200m")
      #plt.loglog(eGRANDini,lGRANDini,linestyle='--',color="orange",lw=4,label="GRAND200k ini (agr)")
      plt.xlabel(r"E$_\nu$ (GeV)")
      plt.ylabel(r"E$^2\Phi(E)$ limit (GeV/s/cm$^{2}$/sr)")
      plt.axis((3E+07, 1E+12, 1E-10, 1E-05))
      plt.legend(loc='best')
      plt.grid(True)
    
eprel = 1e11*numpy.array([0.0010,0.0030,0.0100,0.0300,0.1000,0.3000,1.0000,3.0000])
expprea = 1e17*numpy.array([0.0260,0.2619,0.9504,1.7266,2.5552,3.4848,4.5618,5.2054])
expprec = 1e17 *numpy.array([0,0.0769,0.5324,1.2734,2.0508,2.7981,4.1219,5.0788])
exppreai = numpy.interp(xi,eprel,expprea) # Interpolate ini rate @ new energies
exppreci = numpy.interp(xi,eprel,expprec)
fact = (7500/(0.8*0.8))/10000
fact = 1
exppreai = exppreai/fact  # Scale to same nb of antennas
exppreci = exppreci/fact  

plt.figure(1)
plt.legend(loc='best')
plt.savefig("primaries-rate.png")
plt.figure(2)
if DISPLAY<5:
  plt.loglog(xi,exppreai,'-.',color='brown',lw=4,label="HS1 ini (agr)")
#if DISPLAY==1:
#  plt.loglog(xi,exppreci,'-.',color = 'magenta',lw=4,label="HS1 ini (Cons)")
plt.legend(loc='best')
plt.savefig("exposure.png")

#rpa = yv[2,:]/exppreai  # Agressive prel to Agressive now
#rca = yv[1,:]/exppreci  # Conservative prel to Conservative now
#plt.figure(3)
#plt.semilogx(xi,rpa,'red',lw=3,label='Present/ini (Agr)')
#plt.semilogx(xi,rca,'green',lw=3,label='Present/ini (Cons)')
#plt.grid(True)
#plt.legend(loc='best')
#plt.axis((1E+08, 1E+12, 0, 1.5))
#plt.xlabel(r"E$_\nu$ (GeV)")
#plt.ylabel(r"HS1 exposure ratio")


## Sensitivity limit
    
plt.show()

