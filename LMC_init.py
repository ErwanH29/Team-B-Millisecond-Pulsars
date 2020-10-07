#This code is to define an LMC potential and place it at the right galactrocentric coordinates
from amuse.lab import units
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from amuse.ext.galactic_potentials import Plummer_profile, NFW_profile

cLMC = SkyCoord(ra=79.88*u.degree, dec=-69.59*u.degree, distance=49.97*u.kpc) #Data of LMC distance taken from Pieterzxynski et al 2013 in the Sun's reference frame
gcLMC_co = cLMC.transform_to(coord.Galactocentric)
x_co = np.array(gcLMC_co.cartesian.x) | units.kpc
y_co = np.array(gcLMC_co.cartesian.y) | units.kpc
z_co = np.array(gcLMC_co.cartesian.z) | units.kpc

class LMC_pot(object):
    def __init__(self):
          
        self.plum = Plummer_profile(2e10|units.MSun, 0.73|units.kpc)
        self.nfw = NFW_profile(8.18e6|units.MSun/units.kpc**3, #This value was taken from SIffert et al. 2011
                               2.6|units.kpc)
     
    def get_potential_at_point(self, eps, x, y, z):
        return self.plum.get_potential_at_point(eps, 
                                                x-x_co, 
                                                y-y_co, 
                                                z-z_co) + \
            self.nfw.get_potential_at_point(eps, 
                                                x-x_co, 
                                                y-y_co, 
                                                z-z_co)
                                            
    def get_gravity_at_point(self, eps, x, y, z):
        ax_p, ay_p, az_p = self.plum.get_gravity_at_point(eps, 
                                                x-x_co, 
                                                y-y_co, 
                                                z-z_co)
        ax_h, ay_h, az_h = self.nfw.get_gravity_at_point(eps, 
                                                x-x_co, 
                                                y-y_co, 
                                                z-z_co)
        return ax_p+ax_h, ay_p+ay_h, az_p+az_h
    

    
    
    

    
