import sys
import numpy as np
import pylab as pl
RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
sys.path.append(RETRODIR)
sys.path.append(RETRODIR+"lib/python/")
from retro.event import EventIterator

#noise = 7.9 #muV 150-200MHz
#noise = 7.7 #muV 50-90MHz
#noise = 12. #muV 70-150MHz
noise = 15. #muV (minimal rms noise level on HorizonAntenna in 50-200MHz)

# Simu new
step = 1000.
stepini = 500.

# Simu ini
#step = 800.
#step = 1200.
#stepini = 400.

DISPLAY = 0

def checkCluster(event,antIDs,nantsmin=4):
    #print """checkCluster: check if antennas cluster.""" 
    if len(antIDs)==0:  # Empty array
      #print "Less than min ants nb above threshold."
      if DISPLAY:
        pl.close("all")
      return False, []
      
    bTrig = False
    ants = np.array(event["antennas"])
    xants = ants[antIDs,0]
    yants = ants[antIDs,1]
    zants = ants[antIDs,2]
    posAnts = np.array([xants, yants, zants])
    nAnts = len(xants)  
    if nAnts<= nantsmin:
        return False, []
    d = 1e12*np.ones(shape=(nAnts,nAnts))
    for i in range(nAnts):
      for j in range(i,nAnts):
         d[i,j]=np.linalg.norm(posAnts[:,i]-posAnts[:,j])
         d[j,i]=d[i,j]
	 
    a = np.sort(d)
    d3rdneighbourg = a[:,nantsmin]	#  Distance of 3rd closest trigged antenna
    bCluster= d3rdneighbourg<2*step  # Antennas with 3rd neighbourgh within 2 steps
    if np.sum(bCluster)>0: # At least one
    	bTrig = True;
	
    if DISPLAY:
      print a[:,0:4]
      print d3rdneighbourg
      print bCluster
      print bTrig
      pl.figure(1)
      pl.scatter(xants[bCluster],yants[bCluster],marker='o',color='g')
      pl.show()
      raw_input()
      pl.close("all")
      
    return bTrig,bCluster
    #return True,bCluster

def checkCone(event):
    #print "checkCone"
    ants = np.array(event["antennas"])
    bAntIn = setStep(event,step,stepini)  # is it a "true" antenna pos? (ie sim with 500m step)
    antIDs = np.where(bAntIn)[0]  # Antenna ID
    xants = np.array(ants[antIDs,0])
    yants = ants[antIDs,1]
    zants = ants[antIDs,2]
    posAnts = np.array([xants, yants, zants])
    posAnts = np.transpose(posAnts)
    nAnts = len(xants)  
    
    _, e, (x0, y0, z0), u, (la, lo, h), (t, p) = event["tau_at_decay"]
    posDecay = [x0,y0,z0]
    d = np.linalg.norm(posAnts-posDecay,axis=1)
    #tAntsID = antIDs[np.where(d<90000)[0]]
    tAntsID = antIDs
    
    if DISPLAY:
      ants = np.array(event["antennas"])
      xants = np.array(ants[:,0])
      yants = ants[:,1]
      pl.figure(1)
      pl.scatter(xants[tAntsID],yants[tAntsID],marker='o',color='b')
    
    return tAntsID
   
def checkTrig(event,thresh=2):
    #print """checkTrig: determines if event trigs array. Returns list of trigged antennas.""" 
    th = thresh*noise  #
    bAntIn = setStep(event,step,stepini)  # is it a "true" antenna pos? (ie sim with 500m step)
    antIDs = np.where(bAntIn)[0]  # Antenna ID
    #print "len(antIDs):",len(antIDs)
    bTrig = False
    try:
      v = np.array(event["voltage"])
    except:
      #print "No voltage!"
      return []
    #print "len(v)",len(v)  
    antsin = []
    Ampx=[]
    Ampy=[]
    Ampz=[]
    Ampxy = [];   
    for i in range(np.shape(v)[0]):
      if np.any(antIDs==int(v[i,0])):  # This is a true antenna
        antsin.append(int(v[i,0]))  # Index of antennas with radio simulation
        Ampx.append(float(v[i,1]))  # NS arm
        Ampy.append(float(v[i,2]))  # EW arm
        Ampxy.append(float(v[i,4]))  # EW arm
        Ampz.append(float(v[i,3]))  # Vert arm

    antsin = np.array(antsin)
    #print "len(antsin)",len(antsin)
    Ampx = np.array(Ampx)
    Ampy = np.array(Ampy)
    Ampz = np.array(Ampz)
    Ampxy = np.array(Ampxy)
    trigMat = np.array([Ampx>th,Ampy>th,Ampz>th,Ampxy>th*np.sqrt(2)])
    iTrig = trigMat.sum(axis=0)>0 # Vector of antenna trigger flags (Warning: <=>index of "voltage" field)
    tAntsID =  np.array(antsin[iTrig])  # ID of trigged antennas   
    ntrigs = np.sum(iTrig)  # Nb of antennas trigged in this event
    
    if DISPLAY:
      ants = np.array(event["antennas"])
      xants = np.array(ants[:,0])
      yants = ants[:,1]
      pl.figure(1)
      print tAntsID
      if np.size(tAntsID)>0:
        pl.scatter(xants[tAntsID],yants[tAntsID],marker='o',color='b')
    
    return tAntsID

#
def setStep(event,step,stepini):
    #
    r = int(step/stepini)
    ants = np.array(event["antennas"])
    xants = np.array(ants[:,0])
    yants = np.array(ants[:,1])
    xantsr = xants-np.min(ants[:,0])  # Local coordinates, using 1st antenna as ref, otherwise rounding does not work
    yantsr = yants-np.min(ants[:,1])
    #xantsr = xants
    #yantsr = yants
    #zants = np.array(ants[:,2])
    #alpha = np.array(ants[:,3])
    
    xsteps = np.round(xantsr/stepini)
    xin = (xsteps % r == 0)
    ysteps = np.round(yantsr/stepini)
    yin = (ysteps % r == 0)
    ain = np.logical_and(xin,yin)
    
    if DISPLAY:    
      pl.figure(1)
      pl.scatter(xants,yants,marker='+',color='k')
      pl.scatter(xants[ain],yants[ain],marker='o',color='k')
      print "len(ain):",len(ain)
      #pl.figure(2)
      #pl.subplot(211)
      #pl.hist(xin,100)
      #pl.subplot(212)
      #pl.hist(yin,100)
      #sel = np.where( (xants>-115500) & (xants<-115000) )
      #pl.figure(3)
      #pl.subplot(311)
      #pl.hist(xants[sel],1000)
      #pl.subplot(312)
      #sl.hist(xsteps[sel],1000)
      #pl.subplot(313)
      #pl.hist(xin[sel],1000)
    
    return ain
    
