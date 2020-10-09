import matplotlib.pyplot as plt
import numpy as np
from conv_coord import conv_coordPulsars
from matplotlib import pyplot
from amuse.units import units
import random

#This script is to show how to extract data of the real pulsars
#It also shows the PDF of existing Pulsar's velocities as well as distance w.r.t
#Galactic center and the distribution of the theoretical velocities.
def VelocityProbability(vel):
    vel = vel
    meanV = 300
    sigmaV = 190
    return np.sqrt(2/np.pi)*(meanV**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

###REAL PULSARS###

#COORDINATE DATA
#Importing milisecond pulsar coordinate data
dataCords = np.loadtxt("PulsarsDataATNF(2019-04-24).txt", comments="#", usecols=(6, 7, 8))
Coordinates = conv_coordPulsars(dataCords)

#VELOCITY DATA
#Importing milisecond pulsar velocity data
dataVels = np.loadtxt("VELOCITIESPulsarsDataATNF(2019-04-24).txt", 
                      comments="#", usecols=(3, 6, -1))
                    
#Getting a histogram of current MSP velocities
#num_bins = 40
#n, bins, patches = plt.hist(dataVels[:,-1], num_bins, color='black', alpha=0.5, histtype='step')
#pyplot.show()

###HYPOTHETICAL PULSARS###
#Getting the theoretical distribution of MSP velocities
pop = 1000
vel = np.linspace(min(dataVels[:,-1]), max(dataVels[:,-1]), pop)
velocityDistribution = VelocityProbability(vel)
velocityDistr = np.stack((velocityDistribution, vel), axis=1) 

#Setting up array of initial conditions for pulsars to shoot
#I use a Gaussian distribution for the mass based on the central limit theorem
#and no reports have shown a good distribution function
meanM = 1.4
sigmaMass = 0.3 # mean and standard deviation

#Final Initial Conditions of our 1000 Millisecond-Pulsars
MassDistributionPulsars = np.random.normal(meanM, sigmaMass, 1000) | units.MSun #For every coordinate we shoot 1000 pulsars
velocitylistPulsars = random.choices(vel, weights=velocityDistribution, k =1000) | units.km/units.s
#Initialcoordinates =

#Still need to get the initial position distribution - I don't know how we should sample this
#Maybe choose existing globular positions in the LMC?