import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord

def ProbFuncReal(vel):
    meanV = 250 ; sigmaV = 190
    return np.sqrt(2/(np.pi))*(meanV/sigmaV**2)*np.exp(-(vel)**2/(2*sigmaV**2))

#Velocities
data = np.loadtxt("Data/VelocityOnlyPulsars.txt", comments="#", usecols=(1)) #Importing milisecond pulsar data

x = np.linspace(0, max(data), len(data))
v = np.linspace(0,1200,4000)
dx = 1200/4000
vel = ProbFuncReal(v)/np.sum(ProbFuncReal(v))

cdfvel = np.cumsum(ProbFuncReal(v) * dx)/np.sum(np.cumsum(ProbFuncReal(v)))*15

plt.title(r"Normalised Probability Distribution of Ejected Velocities $(v_{ej})$")
plt.xlabel(r"Velocity ($km s^{-1}$) ")
plt.ylabel(r"Probability")
plt.plot(v, vel, label = 'Population')
plt.plot(v, cdfvel, 'r--', label ='CDF')
plt.legend()
plt.savefig("CDFPDFVelocities.png")
plt.show()

plt.title("Histogram of Millisecond Pulsar Velocities as Observed by the ATNF")
n, bins, patches = plt.hist(data, 120, color = 'black', histtype='step')
plt.xlim(np.min(data),1300)
plt.xlabel(r"Velocity ($km s^{-1}$) ")
plt.savefig("Histogram", dpi = 300)
plt.show()

#Distances
data = np.loadtxt("Data/PulsarsDataATNF(2019-04-24).txt", comments="#", usecols=(6,7,8)) #Importing milisecond pulsar data

#data = np.stack((data[:,1], data[:,2], data[:,0]), axis=1)
#galactocentr = SkyCoord(x=data[:,0], 
#                       y=data[:,1], 
#                       z=data[:,2], 
#                       unit='kpc', 
#                       representation_type='cartesian')
#distanceMSP = galactocentr.transform_to(coord.Galactocentric)
#x_coMSP = np.array(distanceMSP .cartesian.x)
#y_coMSP = np.array(distanceMSP .cartesian.y)
#z_coMSP = np.array(distanceMSP .cartesian.z)
#MSPCoords = np.stack((z_coMSP, x_coMSP, y_coMSP), axis=1)

distancemod = [(abs(np.sqrt(np.sum(data[i]**2))-8.2)) for i in range(len(data))]

plt.title("Histogram of Millisecond Pulsar Distances \n From Galactic Center as Observed by the ATNF")
n, bins, patches = plt.hist(distancemod, 120, color = 'black', histtype='step')
plt.xlim(np.min(distancemod),20)
plt.xlabel(r"Distance (kpc) ")
plt.savefig("HistogramDistance", dpi = 300)
plt.show()

#distancemod = [np.sqrt(np.sum(MSPCoords[i]**2)) for i in range(len(MSPCoords))]
#plt.title("Histogram of Pulsar Velocities as Observed by the ATNF")
#n, bins, patches = plt.hist(distancemod, 30, color = 'black', histtype='step')
#plt.xlim(np.min(distancemod),np.max(distancemod))
#plt.xlabel(r"Velocity ($km s^{-1}$) ")
#plt.savefig("HistogramDistance", dpi = 300)
#plt.show()
