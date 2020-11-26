import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u

def ProbFuncReal(vel):
    meanV = 250 ; sigmaV = 190
    return np.sqrt(2/(np.pi))*(meanV/sigmaV**2)*np.exp(-(vel)**2/(2*sigmaV**2))

#Velocities
data = np.loadtxt("Data/VelocityOnlyPulsars.txt", comments="#", usecols=(1)) #Importing milisecond pulsar data
print(data)
x = np.linspace(0, max(data), len(data))
v = np.linspace(0,1200,4000)
dx = 1200/4000
vel = ProbFuncReal(v)/np.sum(ProbFuncReal(v)) * 0.0035/0.0012

plt.title("Histogram of Millisecond Pulsar Velocities as Observed by the ATNF")
n, bins, patches = plt.hist(data, 120, color = 'black', density=True, histtype='step')
plt.plot(v, vel, label = 'Population')
plt.xlim(np.min(data),1300)
plt.xlabel(r"Velocity ($km s^{-1}$) ")
plt.savefig("Histogram", dpi = 300)
plt.show()

#Distances
data = np.loadtxt("Data/PulsarsDataATNF(2019-04-24).txt", comments="#", usecols=(6,7,8)) #Importing milisecond pulsar data
distancemod = [(abs(np.sqrt(np.sum(data[i]**2))-8.2)) for i in range(len(data))]

plt.title("Histogram of Millisecond Pulsar Distances \n From Galactic Center as Observed by the ATNF")
n, bins, patches = plt.hist(distancemod, 120, color = 'black', histtype='step')
plt.xlim(np.min(distancemod),35)
plt.xlabel(r"Distance (kpc) ")
plt.savefig("HistogramDistance", dpi = 300)
plt.show()

beyondMW = [i for i in distancemod if i>15]

print("Total Millisecond Pulsars: ", len(data))
print("Total Millisecond Pulsars beyond the Milky Way:", len(beyondMW))