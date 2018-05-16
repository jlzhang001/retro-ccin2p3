#import matplotlib
#matplotlib.use('Agg')
# Taken from Efield_2Dmap.py

import os
import sys
import time
from sys import argv
import glob
import numpy as np
import pylab as pl
import math
import scipy.interpolate as itp
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
RETRODIR = "/home/martineau/GRAND/soft/neutrinos/retro/"
sys.path.append(RETRODIR)
sys.path.append(RETRODIR+"lib/python/")
from retro.event import EventIterator, EventLogger
from common import checkTrig
sys.path.append("/home/martineau/GRAND/soft/neutrinos/simulations")
from modules import TopoToAntenna


exclude = "/usr/local/python/python-2.7/lib/python2.7/site-packages"
sys.path = [v for v in sys.path if not (exclude in v)]
pl.style.use("/home/martineau/GRAND/soft/neutrinos/retro-ccin2p3/deps/mplstyle-l3/style/l3.mplstyle")


pl.ion()
#pl.ioff()

noise = 15. #muV (minimal rms noise level on HorizonAntenna in 50-200MHz)
th = 4*noise


def amp2DPattern(zvector,tit):
    pl.figure()
    norm = colors.Normalize(vmin=zvector.min(),vmax=zvector.max())
    pl.scatter(decay_pos[0]/1e3, decay_pos[1]/1e3,  c='r', marker='h',s=100)
    pl.plot([decay_pos[0]/1e3, vend[0]/1e3],[decay_pos[1]/1e3, vend[1]/1e3],'r')
    pl.scatter(xants/1e3,yants/1e3,c=zvector,marker='o',cmap='jet',s=200,norm=norm)
    pl.axis('equal')
    pl.grid(True)
    pl.xlabel('SN (km)')
    pl.ylabel('EW (km)')
    pl.title(tit)
    cbar = pl.colorbar()
    cbar.set_label('peak-peak voltage ($\mu$V)')

def doScatterPlot(x,y,figId=-1,mark='o',xlab='X axis',ylab='Y axis'):
    if figId<0:
      fig = pl.figure()
    else:
      fig = pl.figure(figId)
    pl.scatter(x,y, marker=mark,c="k",s=50)
    pl.xlabel(xlab)
    pl.ylabel(ylab)
    pl.grid(True)

def doScatterCol(x,y,z,figId=-1,xlab='X axis',ylab='Y axis',clab='Color Axis'):
    if figId<0:
      fig = pl.figure()
    else:
      fig = pl.figure(figId)
    if figId > 80:
      norm = colors.Normalize(-100,+100)
    else:
      norm = colors.Normalize(vmin=z.min(),vmax=z.max())
    pl.scatter(x,y,c=z,marker='o',cmap='jet',s=200,norm=norm)
    pl.xlabel(xlab)
    pl.ylabel(ylab)
    cbar = pl.colorbar()
    cbar.set_label(clab)
    pl.grid(True)

if __name__ == "__main__":
 wkdir = sys.argv[1] # path where the simulation file is
 if len(sys.argv)!=2:
     print """\
     This script plot the 2D map of the Efield & voltages amplitude pattern
     Usage:  python plotShowerCaras.py [folder containing the .json file]
     """
     sys.exit(1)

 # First load json infos
 json_files =  glob.glob(wkdir+'/*'+'.voltage.json')
 for i in range(len(json_files)):  #
   json_file = json_files[i]
   print "o Processing ",json_file
   for evt in EventIterator(json_file):  # Should only be one
     # Tau decay infos
     tau =  evt["tau_at_decay"]
     decay_pos =  tau[2]
     decay_pos = np.array(decay_pos)
     tau_dir = tau[5]
     zen = tau_dir[0]
     # Compute the shower energy
     shower_energy = 0.
     for (pid_, momentum) in evt["decay"]:
    	aid = abs(pid_)
    	if aid in (12, 13, 14, 16):
    	    continue
    	shower_energy += sum(m**2 for m in momentum)**0.5
     shower_energy = shower_energy/1e9	
     print shower_energy
     #if shower_energy>1:
     #  print 'Skip E'
     #  continue	
     if zen>89:
       print 'Skip Zen'
       continue
     az = tau_dir[1]
     print  "Decay at position",decay_pos,"in direction (theta,phi)=",tau_dir
     cz = np.cos(zen*np.pi/180)
     sz = np.sin(zen*np.pi/180)
     ca = np.cos(az*np.pi/180)
     sa = np.sin(az*np.pi/180)
     k = np.array([ca*sz,sa*sz,cz])  # Direction vector
     vend = k*9e4+decay_pos  # 90km line

     # Ants infos
     ants = np.array(evt["antennas"])
     xants = ants[:,0]
     print "Nb of antennas in cone (1000m step):",float(len(xants))/4
     yants = ants[:,1]
     zants = ants[:,2]
     alpha = ants[:,3]
     beta = ants[:,4]
 
     pl.figure(222)
     pl.subplot(131)
     pl.plot(xants/1e3,yants/1e3,'k+')
     pl.plot(decay_pos[0]/1e3,decay_pos[1]/1e3,'ro')
     pl.plot([decay_pos[0]/1e3, vend[0]/1e3],[decay_pos[1]/1e3, vend[1]/1e3],'r')
     pl.grid(True)
     pl.xlabel('SN (km)')
     pl.ylabel('EW (km)')
     pl.subplot(132)
     pl.plot(xants/1e3,zants/1e3,'k+')
     pl.plot(decay_pos[0]/1e3,decay_pos[2]/1e3,'ro')
     pl.plot([decay_pos[0]/1e3, vend[0]/1e3],[decay_pos[2]/1e3, vend[2]/1e3],'r')
     pl.grid(True)
     pl.xlabel('SN (km)')
     pl.ylabel('Alt (km)')
     pl.subplot(133)
     pl.plot(yants/1e3,zants/1e3,'k+')
     pl.plot(decay_pos[1]/1e3,decay_pos[2]/1e3,'ro')
     pl.plot([decay_pos[1]/1e3, vend[1]/1e3],[decay_pos[2]/1e3, vend[2]/1e3],'r')
     pl.grid(True)
     pl.xlabel('EW (km)')
     pl.ylabel('Alt (km)')
 
     # Voltage/Trigger infos
     antsIDs = checkTrig(evt)
     if len(antsIDs)>0:
#     if 1:
 	 antsin = []
 	 Ampx=[]
 	 Ampy=[]
 	 Ampz=[]
 	 Ampxy = [];
 	 v = np.array(evt["voltage"])
 	 print "Nb of antennas with radio (1000m step):",float(np.shape(v)[0])/4
 	 for i in range(np.shape(v)[0]):
 	   antsin.append(int(v[i,0]))  # Index of antennas with radio simulation
 	   Ampx.append(float(v[i,1]))  # NS arm
 	   Ampy.append(float(v[i,2]))  # EW arm
 	   Ampz.append(float(v[i,3]))  # Vert arm
 	   Ampxy.append(float(v[i,4]))  # EW arm

 	 if len(antsIDs)>0:
	   #antsin = np.array(antsin[antsInd])
 	   bin = np.in1d(antsin,antsIDs)
	 else:
	   bin = np.arange(len(antsin))
 	 Ampx = np.array(Ampx)[bin]
 	 Ampy = np.array(Ampy)[bin]
 	 Ampz = np.array(Ampz)[bin]
 	 Ampxy = np.array(Ampxy)[bin]
	   
     else:
 	 print "Shower not triggered."
 	 continue

     print "Nb of antennas trigged:",len(antsIDs)
     print "(NS,EW,Vert,NS+EW):",np.sum(Ampx>th),np.sum(Ampy>th),np.sum(Ampz>th),np.sum(Ampxy>th*np.sqrt(2))
 
     xantsr = xants[antsin]
     yantsr = yants[antsin]
     zantsr = zants[antsin]
     if len(antsIDs)>0:
       xants = xants[antsIDs]
       yants = yants[antsIDs]
       zants = zants[antsIDs]
     else:
       xants = xantsr
       yants = yantsr
       zants = zantsr
       
     posAnts = np.array([xants,yants,zants])
 
     pl.figure(222)
     pl.subplot(131)
     pl.plot(xantsr/1e3,yantsr/1e3,'ob')
     pl.plot(xants/1e3,yants/1e3,'og')
     pl.show()
 
     # Compute shower angle
     vants = np.transpose(posAnts)
     vants = vants-decay_pos
     dDecay = np.linalg.norm(vants,axis=1)  # Distance from antennas to decay point
     angSh = np.arccos(np.dot(vants,k)/(np.linalg.norm(vants,axis=1)))*180/np.pi # Angular distae from antenna to decay point
 
     # Compute (zen, az) in antenna ref
     theta_ant = []
     phi_ant = []
     alpha = alpha[antsIDs]

     beta = beta[antsIDs]
     for j in range(len(alpha)):
       ushp = TopoToAntenna(-k,alpha[j],beta[j])  # Xmax vector in antenna frame
       zenith=np.arccos(ushp[2])*180/np.pi  # Zenith in antenna frame
       azim=math.atan2(ushp[1],ushp[0])*180/np.pi
       if azim>360:
 	   azim = azim-360
       elif azim<0:
 	   azim = azim+360
       theta_ant.append(zenith)
       phi_ant.append(azim)
 
     ################ Plots #################
 
     # Antenna array 3D-plots (to check geometry)
 
     if len(alpha)>0:
       fig = pl.figure(1)
       ax = fig.add_subplot(111, projection='3d')
       slop = alpha
       norm = colors.Normalize(vmin=slop.min(),vmax=slop.max())
       ax.scatter(xants/1e3, yants/1e3, zants,c=slop, marker='o',cmap='jet',s=100,norm=norm)
       ax.scatter(decay_pos[0]/1e3, decay_pos[1]/1e3, decay_pos[2], c='r', marker='h',s=100)
       ax.plot([decay_pos[0]/1e3, vend[0]/1e3],[decay_pos[1]/1e3, vend[1]/1e3], [decay_pos[2], vend[2]],'r')
       pl.xlabel('SN (km)')
       pl.ylabel('EW (km)')
       ax.set_zlabel('Alt (m)')
       m = cm.ScalarMappable(cmap=cm.jet,norm=norm)
       m.set_array(slop)
       pl.colorbar(m)
       #pl.zlabel('Alt (m)')
 
     # Amplitude patterns @ ground
     amp2DPattern(Ampx,'Vx (NS)')
     amp2DPattern(Ampy,'Vy (EW)')
     #amp2DPattern(Ampz,'Vz (Vert)')
     #amp2DPattern(Ampxy,'Vx+Vy')

     #doScatterPlot(dDecay/1e3,Ampx,mark='+',xlab='Decay distance (km)',ylab='AmpSN ($\mu V$)')
     #doScatterPlot(angSh,Ampx,mark='+',xlab='Shower angle (deg)',ylab='AmpSN ($\mu V$)')
     #doScatterPlot(theta_ant,Ampy,mark='+',xlab='$\Theta_{ant}$ (deg)',ylab='AmpEW ($\mu V$)')
     #doScatterPlot(theta_ant,Ampx,mark='+',xlab='$\Theta_{ant}$ (deg)',ylab='AmpSN ($\mu V$)')

     #doScatterCol(phi_ant,angSh,Ampy,xlab='$\phi_{ant}$ (deg)',ylab='Angle to axis (deg)',clab='AmpEW ($\mu V$)')
     #doScatterCol(theta_ant,angSh,Ampx,xlab='$\Theta_{ant}$ (deg)',ylab='Angle to axis (deg)',clab='AmpSN ($\mu V$)')
     #doScatterCol(theta_ant,phi_ant,Ampx,xlab='$\Theta_{ant}$ (deg)',ylab='$\phi_{ant}$ (deg)',clab='AmpSN ($\mu V$)')
     #doScatterCol(phi_ant,angSh,Ampx,xlab='$\phi_{ant}$ (deg)',ylab='Angle to axis (deg)',clab='AmpSN ($\mu V$)')
 
     #doScatterPlot(alpha,Ampx,xlab='Slope (deg)',ylab='AmpSN ($\mu V$)')
 
     #fig = pl.figure(70)
     #pl.hist(Ampx/Ampy,50)
     #pl.xlabel('AmpSN/AmpEW')
     #doScatterPlot(phi_ant,Ampx/Ampy,xlab='$\phi_{ant}$ (deg)',ylab='AmpSN/AmpEW')
 
     pl.figure()
     alphat = alpha*np.sign(np.cos(beta))
     pl.hist(alphat,100)
     pl.xlabel('Alpha (deg)')
     
     pl.show()
     raw_input()
     pl.close("all")
