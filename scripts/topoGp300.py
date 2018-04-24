import numpy
import numpy as np
import matplotlib.pyplot as plt
from grand_tour import Topography

A = numpy.array([0,0])
B = numpy.array([12000,14000])
u = B-A
traj = A+u*numpy.linalg.norm(u)


latitude, longitude = 43.5,94.0
print "Loading topography map around (lat,long)=",latitude, longitude,"..."
topo = Topography(latitude=latitude, longitude=longitude,path="share/topography", stack_size=121)
print "Done."


# Plot the topography
xmin = -50000
xmax = +50000
ymin = -50000
ymax = +50000
xt = numpy.linspace(xmin, xmax, 1001)
yt = numpy.linspace(ymin, ymax, 1001)
zt = numpy.zeros((len(yt), len(xt)))
lat1 = 43.612421
long1 = 93.820829
alt1 = 2783
ant0 = topo.lla_to_local(lat1,long1,alt1)
print "First antenna location:",ant0

for i, yi in enumerate(yt):
    for j, xj in enumerate(xt):
        zt[i, j] = topo.ground_altitude(xj, yi)        
	#lla = topo.local_to_lla([xj, yi,zt[i, j]])  # Get topo
        #latitude, longitude, altitude = lla
	#print lla
#X = SN
#Y = EW		
plt.figure(1)
plt.pcolor(xt, yt, zt, cmap="terrain", alpha=0.75)
plt.xlabel("Northing (m)")
plt.ylabel("Westing (m)")
plt.colorbar()
plt.plot(ant0[0],ant0[1],'ko')

setsA = [[-4000,-6000],[-6000,6000],[3000,20000]]
setsB = [[12000,14000],[12000,-9000],[16000,-2000]]
col = ["k-","r-","b-"]
for i in range(len(setsA)):
  A = numpy.array(setsA[i])
  B = numpy.array(setsB[i])
  u = B-A
  plt.figure(1)
  plt.plot([A[0],B[0]],[A[1],B[1]],col[i])
  normu = numpy.linalg.norm(u)
  u = u/normu
  traj = numpy.linspace(0,normu,200)
  xtraj = numpy.array(A[0]+u[0]*traj)
  ytraj = numpy.array(A[1]+u[1]*traj)
  zt = numpy.zeros(len(traj))
  for j in range(len(traj)):
    zt[j] = topo.ground_altitude(xtraj[j], ytraj[j])

  plt.figure()
  plt.plot(traj/1e3,zt,col[i])  
  plt.grid(True)
  plt.xlabel("Distance along track (km)")
  plt.ylabel("Altitude (m)")
  

plt.show()
