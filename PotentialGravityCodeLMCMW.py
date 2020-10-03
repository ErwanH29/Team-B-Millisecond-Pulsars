import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units, constants
from amuse.plot import plot
from astropy import units as u
from astropy.coordinates import SkyCoord
from amuse.ext.galactic_potentials import MWpotentialBovy2015, MiyamotoNagai_profile, NFW_profile, Plummer_profile

c = SkyCoord(ra=79.88*u.degree, dec=-69.59*u.degree, distance=49.97*u.kpc) #Data taken from Pieterzxynski et al 2013
xcordLMC = c.cartesian.x ; ycordLMC = c.cartesian.y ; zcordLMC = c.cartesian.z

#Defining an x,y,z space to plot in
D = 100
xMW = np.linspace(-D,D,1000) | units.kpc
yMW = np.linspace(-D,D,1000) | units.kpc
zMW = np.linspace(-D,D,1000) | units.kpc

eps = 0.0 #Dont know what this eps does... really

#Milky Way Statistics - Taken from Bovy 2015
disk_mass = 6.8e+10 | units.MSun
disk_size = 3 | units.kpc #Milky Way's scale radius
disk_height = 280 | units.parsec #Milky Way's scale height
MW = MWpotentialBovy2015()
pot_MW = MW.get_potential_at_point(eps, xMW, yMW, zMW) 

#Large Magellanic Cloud Statistics - Taken from Bekki 2005:
mass = 2e10 | units.MSun #LMC mass
scaleDH_radius = 2.6 | units.kpc #LMC dark matter halo scale radius
scaleB_radius = 0.73 | units.kpc #LMC bulge scale height
rhoTEST = 100 | units.MSun/units.kpc**3 #Need to find a central density value of the dark matter halo
LMCPlummerProfile = Plummer_profile(mass, scaleB_radius) #The Plummer profile for LMC
LMCNFWProfile = NFW_profile(rhoTEST, scaleDH_radius) #The NFW profile for LMC

#Calling the potentials
LMCPlummer = LMCPlummerProfile.get_potential_at_point(eps, x, y, z)
LMCNFW = LMCNFWProfile.get_potential_at_point(eps, x, y, z)
CombinedLMC = LMCNFW + LMCPlummer
pot_MW = MW.get_potential_at_point(eps, xMW, yMW, zMW)

plt.title("LMC Potentials")
plot(x, LMCPlummer, label = 'Plummer')
plot(x, LMCNFW, label = 'NFW')
plot(x, CombinedLMC, label = 'Combined')
plt.legend()
plt.show()         

plt.title("Disk Potentials")
plot(x, pot_MW)
plt.show()         