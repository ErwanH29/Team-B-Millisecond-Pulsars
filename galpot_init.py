#This code is to define an LMC & SMC potential and place it at the right galactrocentric coordinates
from amuse.lab import units
from conv_coord import conv_coord, prop_motion
from amuse.ext.galactic_potentials import Plummer_profile, NFW_profile

class LMC_pot(object):
    def __init__(self, ra, dec, dis, gc = False): # gc: galactic cartesian; when gc = False: (ra, dec, dis) and when gc = True: (x, y, z)
           
        self.plum = Plummer_profile(2e10|units.MSun, 0.73|units.kpc)
        self.nfw = NFW_profile(8.18e6|units.MSun/units.kpc**3, #This value was taken from SIffert et al. 2011
                               2.6|units.kpc)
        self.ra = ra
        self.dec = dec
        self.dis = dis
        if gc == True: 
            self.d = (self.ra, self.dec, self.dis) | units.kpc
        else:
            self.d = conv_coord(self.ra, self.dec, self.dis)
        #self.d = conv_coord(79.88, -69.59, 49.97)
        #self.v = prop_motion(1.91, 0.229, 262)
        
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

class SMC_pot(object):
    def __init__(self, ra, dec, dis, gc = False):
          
        self.plum = Plummer_profile(4e9|units.MSun, 0.339|units.kpc) # scale radius by assuming R~m^(1/3) and comparing to LMC
        self.nfw = NFW_profile(8.18e6|units.MSun/units.kpc**3, # No values for SMC halo....
                               2.6|units.kpc)
        self.ra = ra
        self.dec = dec
        self.dis = dis
        if gc == True:
            self.d = (self.ra, self.dec, self.dis) | units.kpc
        else:
            self.d = conv_coord(self.ra, self.dec, self.dis)
        self.d = conv_coord(13.16, -72.8, 60.6) # taken from SIMBAD
        #self.v = prop_motion(0.772, -1.117, 145.6)
     
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
