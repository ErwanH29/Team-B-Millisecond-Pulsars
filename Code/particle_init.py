from amuse.lab import Particles, Particle, units
import numpy as np
import random
from conv_coord import conv_coord

# initialize particle sets here. 
##### NEUTRON STAR SET #####

class neut_initializer(object):
    
    def __init__(self):
        self.N = 0
        self.dataVels = np.loadtxt("VELOCITIESPulsarsDataATNF(2019-04-24).txt", 
                                   comments="#", usecols=(3, 6, -1))
    
    def N_count(self):
        N = self.N
        return int(N)
    
    def decision(self):
        return random.random() < 0.8 # from the MSP formation rate
    
    def ProbFunc(self, vel):
        meanV = 250 ; sigmaV = 190
        return np.sqrt(2/np.pi)*(meanV**2/sigmaV**3)*np.exp(-vel**2/(2*sigmaV**2))

    def velocityList(self):
        vrange = np.linspace(0, max(self.dataVels[:,-1])/3)
        r=[-1,1]
        w = self.ProbFunc(vrange)
        scalex = [np.random.choice(r)]
        scaley = [np.random.choice(r)]
        scalez = [np.random.choice(r)]
        vx = np.array(random.choices(vrange, weights=w, k = 1))*scalex
        vy = np.array(random.choices(vrange, weights=w, k = 1))*scaley
        vz = np.array(random.choices(vrange, weights=w, k = 1))*scalez
        velocity = np.concatenate((vx,vy,vz))
        return velocity

    def massList(self):
        meanM = 1.4 ; sigmaMass = 0.3
        return np.random.normal(meanM, sigmaMass) #| units.MSun 

    def neut_star_init(self, pos):
        self.N=1
        neuts = Particles(1)
        neuts.mass = self.massList() | units.MSun
        neuts.position = pos
        neuts.velocity = self.velocityList() * (1 | units.kms)
        
        return neuts

    def add_neut(self, neuts, pos, vel_sys):
        c = self.decision()
        if c==True:
            self.N += 1
            
            add_neut = Particles(1)
            add_neut.mass = self.massList() | units.MSun
            add_neut.position = pos 
            add_neut.velocity = self.velocityList() * (1 | units.kms)
            add_neut.velocity += vel_sys
            neuts.add_particle(add_neut)   
        else:
           neuts = neuts 
        return neuts

##### GALAXY TRACERS #####

# initialize LMC and SMC
class gal_initializer(object):
    
    def ngc_1783_pos_init(self):
        ngc_1783_pos_init = conv_coord(74.7875, -65.9878, 50.1) - conv_coord(80.8958, -69.7561, 49.97)
        return ngc_1783_pos_init[0], ngc_1783_pos_init[1], ngc_1783_pos_init[2]
    
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
        ngc_1783.velocity = (-0.31, 0.95, 0) * vel_init # a unit vector tangent to radius vector
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

def old_coordlist():
    LMC_coord = (-32.2877571477, -201.334477302, -194.233794365) | units.kpc
    LMC_vel = (-23.6377151397, -116.106444922, -126.743869139) | units.kms
    SMC_coord = (-35.2794369809, -140.65635819, -96.9363051162) | units.kpc
    SMC_vel = (-49.5945188999, -49.9087554726, -12.7823855815) | units.kms
    ngc_coord = (-30.3480687709, -203.481152417, -197.606147399) | units.kpc
    ngc_vel = (26.7697951421, -130.939929951, -234.017971574) | units.kms
    
    data = (LMC_coord, LMC_vel, SMC_coord, SMC_vel, ngc_coord, ngc_vel)
    return data