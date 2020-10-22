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
        return np.sqrt(1/(2*np.pi))*(1/sigmaV**2)*np.exp(-(vel-meanV)**2/(2*sigmaV**2))

    def velocityList(self, vrange, w):
        r=[-1,1]
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
        vrange = np.linspace(0, max(self.dataVels[:,-1]))/3
        velocityDistr = self.ProbFunc(vrange)
        neuts = Particles(1)
        
        #for n in range(self.N):
        vx, vy, vz = self.velocityList(vrange, velocityDistr)
        print(vx, vy, vz)
        neut = neuts[0]
        neut.mass = self.massList() | units.MSun
        neut.position = conv_coord(74.8, -66, 49) # NGC 1783 position
        neut.velocity = (vx, vy, vz) | units.kms # randomized initial  velocity
        #neut.velocity = (5.97959184, 0, 0) | units.kms # choose these init conditions for captured orbit
        #neut.position = (-0.88, -41, -27) * (1 | units.kpc) # LMC position, slight offset from center
        return neuts

##### GALAXY TRACERS #####

# initialize LMC and SMC
class gal_initializer(object):
    
    def ngc_1783_pos_init(self):
        ngc_1783_pos_init = conv_coord(74.7875, -65.9878, 49)
        return ngc_1783_pos_init[0], ngc_1783_pos_init[1], ngc_1783_pos_init[2]
    
    def orb_comp_LMC(self):
        direction = np.array([47,242,225]) 
        comp_norm = direction / np.sum(direction)
        return comp_norm
           
    def gal_tracer(self, vel_init):
        gal_MC = Particles(2)
        ngc_1783 = Particles(1)
        
        gal_lmc = gal_MC[0]
        gal_lmc.mass = 2e10 | units.MSun 
        gal_lmc.position = conv_coord(80.8958, -69.7561, 49.97)
        gal_lmc.velocity = (47,242,225) * (1 | units.kms)
        
        gal_smc = gal_MC[1]
        gal_smc.mass = 6e9 | units.MSun 
        gal_smc.position = conv_coord(13.1583, -72.8, 61.7)
        gal_smc.velocity = (5.35,164,136) * (1 | units.kms)
        
        ngc_1783.mass = 1.7e5 | units.MSun # McMillan 2017
        ngc_1783.position = self.ngc_1783_pos_init()
        ngc_1783.velocity = self.orb_comp_LMC() * vel_init
        return gal_MC, ngc_1783
    
    def gal_path_init(self):
        x_lmc = [] | units.kpc
        y_lmc = [] | units.kpc
        z_lmc = [] | units.kpc
        
        x_smc = [] | units.kpc
        y_smc = [] | units.kpc
        z_smc = [] | units.kpc
        
        x_ngc_1783 = [] | units.kpc
        y_ngc_1783 = [] | units.kpc
        z_ngc_1783 = [] | units.kpc
        return x_lmc, y_lmc, z_lmc, x_smc, y_smc, z_smc, x_ngc_1783, y_ngc_1783, z_ngc_1783
