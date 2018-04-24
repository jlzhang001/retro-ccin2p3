import numpy as np
import pylab as pl
from retro.event import EventIterator

noise = 15. #muV (minimal rms noise level on HorizonAntenna in 50-200MHz)
step = 1000.

def checkCluster(event, antIDs,mod='a'):
    """Check if antennas cluster.""" 
    
    if len(antIDs)==0:  # Empty array
      return False, []
      
    if mod == 'a':
        ntrigthresh=5
    if mod == 'c':
        ntrigthresh=8
    bTrig = False
    ants = np.array(event["antennas"])
    xants = np.array(ants[antIDs,0])
    yants = ants[antIDs,1]
    zants = ants[antIDs,2]
    posAnts = np.array([xants, yants, zants])
    nAnts = len(xants)
    # TOBEDONE
    #for i in range(nAnts):
    #  v = np.tile(posAnts[:,i], (nAnts,1))
    #  v = np.transpose(v)
    #  d = np.linalg.norm(v-posAnts,axis=0)
    #  #ds = np.ndarray.sort(d)
    #  
    #  print d
    #  pl.figure(12)
    #  pl.plot(d)
    #  pl.show()
    
    # TO BE REPLACED (see computeCluster.m)   
    d = 1e12*np.ones(shape=(nAnts,nAnts))
    for i in range(nAnts):
      for j in range(i+1,nAnts):
         d[i,j]=np.linalg.norm(posAnts[:,i]-posAnts[:,j])
    		 
    minDist= np.ndarray.min(d,axis=1)
    bCluster = minDist<step*2  # Clustered antennas if closer than twice step size
    ntrigs = np.sum(bCluster)  # Nb of antennas trigged in this event
    if ntrigs>=ntrigthresh:
    	bTrig = True;
	
    return bTrig,bCluster
    #pl.figure(1)
    #pl.scatter(xants,yants)
    #pl.show()

def checkTrig(event,mod='a'):
    """Determines if event trigs array. Returns list of trigged antennas.""" 
    if mod  == 'a': # Agressive threshold
    	v = 3   # 6 sigmas
	ntrigthresh=5
    if mod  == 'c': # Conservative threshold
    	v = 10
	ntrigthresh=8
    th = v*noise  #2xrms is peak-peak noise level
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
    if ntrigs>=ntrigthresh:
    	bTrig = True;
    return bTrig,tAntsID

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
        
    #pl.figure(1)
    #pl.scatter(xants,yants,marker='+')
    #pl.scatter(xants[ain],yants[ain],marker='o',color='g')
    #pl.show()
    #print "len(ain):",len(ain)
    return ain
    
