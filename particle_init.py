from amuse.lab import Particles, units
import numpy as np
import random

# initialize particle sets here. 
##### NEUTRON STAR SET #####

class neut_initialiser(object):
    def ProbFunc(vel):
        meanV = 300 ; sigmaV = 190
        return np.sqrt(2/np.pi)*(meanV**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(vrange, w, N):
        return random.choices(vrange, weights=w, 
                              k = N) | units.km/units.s   

    def massList(N):
        meanM = 1.4 ; sigmaMass = 0.3
        return np.random.normal(meanM, sigmaMass, N) | units.MSun 

    N = int(input("Give the number of millisecond pulsars you want to simulate: ", ))    
    dataVels = np.loadtxt("VELOCITIESPulsarsDataATNF(2019-04-24).txt", 
                          comments="#", usecols=(3, 6, -1))
    
    vrange = np.linspace(200, max(dataVels[:,-1]), N)
    velocityDistr = ProbFunc(vrange)
    velList = velocityList(vrange, velocityDistr, N)
    
    x = np.zeros(N) |units.kpc 
    y = np.zeros(N) | units.kpc 
    z = np.zeros(N) | units.kpc
    
    neut_InitialConditions = np.stack((x, y, z, velList, 
                                       massList(N)), axis=1)
    #def neut_stars(N_neut):
     #   return neuts

    #def neut_path_init(N_neut):
     #   for i in range(N_neut):
      #      x_neut = [] | units.kpc
       #     y_neut = [] | units.kpc
        #    z_neut = [] | units.kpc
         #   #X_neut[i] = (x_neut, y_neut, z_neut)
        #return 

#neuts = neut_path_init(10)
##### GALAXY TRACERS #####

# initialize LMC and SMC
def gal_tracer():
    gals = Particles(2)
    gal_lmc = gals[0]
    gal_lmc.mass = 2e10 | units.MSun 
    gal_lmc.position = (-0.88,-41,-27) * (1 | units.kpc)
    gal_lmc.velocity = (47,242,225) * (1 | units.kms)

    gal_smc = gals[1]
    gal_smc.mass = 2e10 | units.MSun 
    gal_smc.position = (15,-36,-42) * (1 | units.kpc)
    gal_smc.velocity = (5.35,164,136) * (1 | units.kms)
    return gals

def gal_path_init():
     x_lmc = [] | units.kpc
     y_lmc = [] | units.kpc
     z_lmc = [] | units.kpc

     x_smc = [] | units.kpc
     y_smc = [] | units.kpc
     z_smc = [] | units.kpc
     return x_lmc, y_lmc, z_lmc, x_smc, y_smc, z_smc 