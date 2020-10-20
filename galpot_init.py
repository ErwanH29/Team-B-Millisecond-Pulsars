# This code is to define an LMC and SMC potential and place it at the right galactrocentric coordinates
# The coordinates can be updated through d_update
from amuse.lab import units
from amuse.ext.galactic_potentials import Plummer_profile, NFW_profile

class LMC_pot(object):
    def __init__(self):
          
        self.plum = Plummer_profile(2e10|units.MSun, 0.73|units.kpc)
        self.nfw = NFW_profile(8.18e6|units.MSun/units.kpc**3, #This value was taken from SIffert et al. 2011
                               2.6|units.kpc)
        self.d = (0,0,0) | units.kpc
    def d_update(self, x_c, y_c, z_c):
        self.d[0] = x_c
        self.d[1] = y_c
        self.d[2] = z_c
        
    def get_potential_at_point(self, eps, x, y, z):
        return self.plum.get_potential_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2]) + \
            self.nfw.get_potential_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
                                            
    def get_gravity_at_point(self, eps, x, y, z):
        ax_p, ay_p, az_p = self.plum.get_gravity_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
        ax_h, ay_h, az_h = self.nfw.get_gravity_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
        return ax_p+ax_h, ay_p+ay_h, az_p+az_h
    
    def radial_force(self, x, y, z):
        r = ((x-self.d[0])**2 + (y-self.d[1])**2 + (z-self.d[2])**2).sqrt()
        return self.plum.radial_force(r) + \
            self.nfw.radial_force(r)
            
    def circular_velocity(self, x, y, z):
        r = ((x-self.d[0])**2 + (y-self.d[1])**2 + (z-self.d[2])**2).sqrt()
        fr = self.radial_force(r)
        return (-r*fr).sqrt()

class SMC_pot(object):
    def __init__(self):
          
        self.plum = Plummer_profile(4e9|units.MSun, 0.339|units.kpc) # scale radius by assuming R~m^(1/3) and comparing to LMC
        self.nfw = NFW_profile(8.18e6|units.MSun/units.kpc**3, # No values for SMC halo....
                               2.6|units.kpc)
        self.d = (0,0,0) | units.kpc
    def d_update(self, x_c, y_c, z_c):
        self.d[0] = x_c
        self.d[1] = y_c
        self.d[2] = z_c
       
    def get_potential_at_point(self, eps, x, y, z):
        return self.plum.get_potential_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2]) + \
            self.nfw.get_potential_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
                                            
    def get_gravity_at_point(self, eps, x, y, z):
        ax_p, ay_p, az_p = self.plum.get_gravity_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
        ax_h, ay_h, az_h = self.nfw.get_gravity_at_point(eps, 
                                                x-self.d[0], 
                                                y-self.d[1], 
                                                z-self.d[2])
        return ax_p+ax_h, ay_p+ay_h, az_p+az_h
    
    def radial_force(self, x, y, z):
        r = ((x-self.d[0])**2 + (y-self.d[1])**2 + (z-self.d[2])**2).sqrt()
        return self.plum.radial_force(r) + \
            self.nfw.radial_force(r)
            
    def circular_velocity(self, x, y, z):
        r = ((x-self.d[0])**2 + (y-self.d[1])**2 + (z-self.d[2])**2).sqrt()
        fr = self.radial_force(r)
        return (-r*fr).sqrt()
