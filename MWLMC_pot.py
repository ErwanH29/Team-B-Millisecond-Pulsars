import numpy as np
import matplotlib.pyplot as pt
from amuse.units import units, constants
from amuse.plot import plot

#Hi guys, the code for MW and LMC background potentials

#Equations and galaxy parameters from Binney & Tremaine 'Galactic Dynamics'
#Initialize attributes (masses here) for the Milkyway (self) object
class MW_galaxy(object):
    def __init__(self,
                 Mb = 1.40592e10 | units.MSun,
                 Md = 8.5608e10 | units.MSun,
                 Mh = 1.07068e11 | units.MSun):
        self.Mb = Mb
        self.Md = Md
        self.Mh = Mh
    
    #Define potential function of galaxy at point (x,y,z). eps does
    #not seem to be used here...
    def get_mwpotential_at_point(self,eps,x,y,z):
        r = (x**2+y**2+z**2)**0.5 #Define position in x,y,z-space
        R = (x**2+y**2)**0.5 #Define position in x-y plane
        #Bulge potential
        b1 = 0.3873 | units.kpc #Define characteristic length scale for pot_bulge
        pot_bulge = -constants.G*self.Mb/(r**2+b1**2)**0.5
        #Disk potential
        a2 = 5.31 | units.kpc #Disk length scale in x-y plane
        b2 = 0.25 | units.kpc #Disk length scale in z-direction
        pot_disk = -constants.G*self.Md/(R**2+(a2+(z**2+b2**2)**0.5)**2)**0.5
        #Halo potential
        a3 = 12.0 | units.kpc #Length scale for DM halo
        cut_off = 100 | units.kpc #Reach of DM halo
        d3 = r/a3
        c3 = 1 + (cut_off/a3)**1.02
        pot_halo = -constants.G*(self.Mh/a3)*d3**1.02/(1+d3**1.02) \
                   - (constants.G*self.Mh/(1.02*a3)) \
                       * (-1.02/c3 + np.log(c3) + 1.02/(1+d3**1.02) \
                          - np.log(1.0 + d3**1.02) )
        return 2 * (pot_bulge + pot_disk + pot_halo)
    #I'll find out what a rigid 
    #potential is...

#Initialize attributes for LMC - I used a plummer sphere as a model here
class LMC_galaxy(object):
    def __init__(self,
                 M0 = 10e10 | units.MSun):
        self.M0 = M0
    
    #define potential function for plummer sphere
    def get_lmcpotential_at_point(self, x, y, z, trans):
        r = (x**2 \
              + y**2 \
                  + z**2)**0.5
        #Set characteristic length scale a for plummer model
        #In this assumption, I took half the size (diameter~4.3 kpc)
        #Please, you guys think about a good/better value for this!
        a_plum = 2.0 | units.kpc
        pot_plum = -constants.G*self.M0/((r**2+a_plum**2)**0.5)
        
        return pot_plum

#Define position of LMC, in galactic coordinates. Angles and distances
#taken from SIMBAD. For ease, x will be along theta=0 deg and y will be
#along theta=270 deg
theta_LMC = 280.4652 #degrees
phi_LMC = -32.8884 #degrees
dist_LMC = 50 | units.kpc

#Define distance in x-y plane
dist_LMC_xy = dist_LMC * np.cos(np.deg2rad(phi_LMC))

#Define distances for all three coordinates
dist_LMC_y = dist_LMC_xy * np.cos(np.deg2rad(theta_LMC-270))
dist_LMC_x = dist_LMC_xy * np.sin(np.deg2rad(theta_LMC-270))
dist_LMC_z = dist_LMC * np.sin(np.deg2rad(phi_LMC))

#Translate LMC potential along axes
trans = [dist_LMC_x, dist_LMC_y, dist_LMC_z]

#Defining an x,y,z space to plot in
D = 100
x = np.linspace(-D,D,1000) | units.kpc
y = np.linspace(-D,D,1000) | units.kpc
z = np.linspace(-D,D,1000) | units.kpc
eps = 0.0 #Dont know what this eps does... really

#Checking how class functions can be called...
MW = MW_galaxy()
pot_MW = MW_galaxy.get_mwpotential_at_point(MW, eps, x, y, z)

LMC = LMC_galaxy()
pot_LMC = LMC_galaxy.get_lmcpotential_at_point(LMC, x, y, z, trans)

#Just for checking plot functionality, hehehe
plot(x, pot_MW, pot_LMC)
pt.show()       
