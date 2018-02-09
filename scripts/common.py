import numpy as np
from retro.event import EventIterator

noise = 15. #muV (minimal rms noise level on HorizonAntenna in 50-200MHz)

def checkTrig(event,ntrigthresh=5,mod='a'):
    """Determines if event trigs array. Returns list of trigged antennas.""" 
    if mod  == 'a': # Agressive threshold
    	    v = 2
    if mod  == 'c': # Conservative threshold
    	    v = 5
    th = 2*v*noise  #2xrms is peak-peak noise level
    bTrig = False

    try:
      v = np.array(event["voltage"])
    except:
      #print "No voltage!"
      return bTrig,[],[]
      
    antsin = []
    Ampx=[]
    Ampy=[]
    Ampz=[]
    Ampxy = [];   
    for i in range(np.shape(v)[0]):
      antsin.append(int(v[i,0]))  # Index of antennas with radio simulation
      Ampx.append(float(v[i,1]))  # NS arm
      Ampy.append(float(v[i,2]))  # EW arm
      Ampxy.append(float(v[i,4]))  # EW arm
      Ampz.append(float(v[i,3]))  # Vert arm
    
    antsin = np.array(antsin)
    Ampx = np.array(Ampx)
    Ampy = np.array(Ampy)
    Ampz = np.array(Ampz)
    Ampxy = np.array(Ampxy)
    trigMat = np.array([Ampx>th,Ampy>th,Ampz>th,Ampxy>th*np.sqrt(2)])
    tAntsInd = trigMat.sum(axis=0)>0 # Vector of antenna trigger flags (<=>index of "voltage" field)
    ntrigs = np.sum(tAntsInd)  # Nb of antennas trigged in this event
    tAntsID = antsin[tAntsInd]  # trigged antennas IDs (<=>index of "antennas" field)

    if ntrigs>=ntrigthresh:
    	bTrig = True;

    return bTrig,tAntsInd,tAntsID

