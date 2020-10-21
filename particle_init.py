from amuse.lab import Particles, units
import numpy as np
import random
from conv_coord import conv_coord

# initialize particle sets here. 
##### NEUTRON STAR SET #####

class neut_initializer(object):
    
    def __init__(self):
        self.N = int(input("Give the number of millisecond pulsars you want to simulate: ", ))
        self.dataVels = np.loadtxt("VELOCITIESPulsarsDataATNF(2019-04-24).txt", 
                                   comments="#", usecols=(3, 6, -1))
    
    def N(self):
        N = self.N
        return N
        
    def ProbFunc(self, vel):
        meanV = 300 ; sigmaV = 190
        return np.sqrt(2/np.pi)*(meanV**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(self, vrange, w):
        r = [-1, 1]
        scalex = [np.random.choice(r)]
        scaley = [np.random.choice(r)]
        scalez = [np.random.choice(r)]
        vx = np.array(random.choices(vrange, weights=w, k = 1))*scalex
        vy = np.array(random.choices(vrange, weights=w, k = 1))*scaley
        vz = np.array(random.choices(vrange, weights=w, k = 1))*scalez
        return vx, vy, vz

    def massList(self):
        meanM = 1.4 ; sigmaMass = 0.3
        return np.random.normal(meanM, sigmaMass) #| units.MSun 

    def neut_path_init(self):

        x_neut = [] | units.kpc 
        y_neut = [] | units.kpc 
        z_neut = [] | units.kpc
        return x_neut, y_neut, z_neut

    def neut_stars(self):
        vrange = np.linspace(200, max(self.dataVels[:,-1]))/3
        velocityDistr = self.ProbFunc(vrange)
        neuts = Particles(self.N)
        vx, vy, vz = self.velocityList(vrange, velocityDistr)
        print(vx, vy, vz)
        for n in range(self.N):
            neut = neuts[n]
            neut.mass = self.massList() | units.MSun
            #neut.position = conv_coord(74.7875, -65.9878, 49) # NGC 1783 position
            #neut.velocity = (vx, vy, vz) | units.kms
            # Choose the following initial conditions to to check out a captured orbit
            neut.velocity = (200, 0, 0) | units.kms 
            neut.position = (-0.88, -41, -27) * (1 | units.kpc) # LMC position, slight offset from center
            neut.velocity += (47,242,225) * (1 | units.kms)
        return neuts

##### GALAXY TRACERS #####

# initialize LMC and SMC
class gal_initializer(object):
    def gal_tracer(self):
        gals = Particles(2)
        gal_lmc = gals[0]
        gal_lmc.mass = 2e10 | units.MSun 
        gal_lmc.position = conv_coord(80.8958, -69.7561, 49.97)
        gal_lmc.velocity = (47,242,225) * (1 | units.kms)
        
        gal_smc = gals[1]
        gal_smc.mass = 6e9 | units.MSun 
        gal_smc.position = conv_coord(13.1583, -72.8, 61.7)
        gal_smc.velocity = (5.35,164,136) * (1 | units.kms)
        return gals
    
    def gal_path_init(self):
        x_lmc = [] | units.kpc
        y_lmc = [] | units.kpc
        z_lmc = [] | units.kpc
        
        x_smc = [] | units.kpc
        y_smc = [] | units.kpc
        z_smc = [] | units.kpc
        return x_lmc, y_lmc, z_lmc, x_smc, y_smc, z_smc
