import sys
import numpy as np
import pylab as pl
RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
sys.path.append(RETRODIR)
sys.path.append(RETRODIR+"lib/python/")
from retro.event import EventIterator

noise = 15. #muV (minimal rms noise level on HorizonAntenna in 50-200MHz)
step = 500.
DISPLAY = 1

    
def checkCluster(event,antIDs,nantsmin=3):
    #print """checkCluster: check if antennas cluster.""" 
    if len(antIDs)==0:  # Empty array
      #print "Less than min ants nb above threshold."
      if DISPLAY:
        pl.close("all")
      return False, []
      
    bTrig = False
    ants = np.array(event["antennas"])
    xants = np.array(ants[antIDs,0])
    yants = ants[antIDs,1]
    zants = ants[antIDs,2]
    posAnts = np.array([xants, yants, zants])
    nAnts = len(xants)  
    if nAnts<= nantsmin:
        return False, []
    d = 1e12*np.ones(shape=(nAnts,nAnts))
    for i in range(nAnts):
      for j in range(nAnts):
         d[i,j]=np.linalg.norm(posAnts[:,i]-posAnts[:,j])
    
    a = np.sort(d)
    d3rdneighbourg = a[:,nantsmin]	#  Distance of 3rd closest trigged antenna
    bCluster= d3rdneighbourg<2*step  # Antennas with 3rd neighbourgh within 2 steps
    if np.sum(bCluster)>0: # At least one
    	bTrig = True;
	
    if DISPLAY:
      print a[:,0:4]
      print bCluster
      print bTrig
      pl.figure(1)
      pl.scatter(xants[bCluster],yants[bCluster],marker='o',color='g')
      pl.show()
      raw_input()
      pl.close("all")
      
    return bTrig,bCluster

def checkCone(event):
    ants = np.array(event["antennas"])
    bAntIn = setStep(event,step)  # is it a "true" antenna pos? (ie sim with 500m step)
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
    tAntsID = antIDs[np.where(d<90000)[0]]
    
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
    print "Threshold (muV) = ",th
    print "Step (m) = ",step
    bAntIn = setStep(event,step)  # is it a "true" antenna pos? (ie sim with 500m step)
    antIDs = np.where(bAntIn)[0]  # Antenna ID
    #print "len(antIDs):",len(antIDs)
    bTrig = False
    try:
      v = np.array(event["voltage"])
    except:
      #print "No voltage!"
      return False,[]
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
    tAntsID =  antsin[iTrig]  # ID of trigged antennas   
    ntrigs = np.sum(iTrig)  # Nb of antennas trigged in this event
    
    if DISPLAY:
      ants = np.array(event["antennas"])
      xants = np.array(ants[:,0])
      yants = ants[:,1]
      pl.figure(1)
      pl.scatter(xants[tAntsID],yants[tAntsID],marker='o',color='b')
      pl.show()
    return tAntsID

#
def setStep(event,step):
    r = int(step/500)
    ants = np.array(event["antennas"])
    xants = np.array(ants[:,0])
    yants = ants[:,1]
    zants = ants[:,2]
    alpha = ants[:,3]
    
    xsteps = np.floor(xants/500)
    xin = (xsteps % r == 0)
    ysteps = np.floor(yants/500)
    yin = (ysteps % r == 0)
    ain = np.logical_and(xin,yin)
    
    if DISPLAY:    
      pl.figure(1)
      pl.scatter(xants,yants,marker='+',color='k')
      pl.scatter(xants[ain],yants[ain],marker='o',color='r')
    #print "len(ain):",len(ain)
    
    return ain
    
