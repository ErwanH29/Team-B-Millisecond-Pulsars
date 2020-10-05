import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units, constants
from amuse.plot import plot
from astropy import units as u
from astropy.coordinates import SkyCoord
from amuse.ext.galactic_potentials import MWpotentialBovy2015, MiyamotoNagai_profile, NFW_profile, Plummer_profile

#Converting coordinate systems to Cartesian
c = SkyCoord(ra=79.88*u.degree, dec=-69.59*u.degree, distance=49.97*u.kpc) #Data of LMC distance taken from Pieterzxynski et al 2013 in the Sun's reference frame
xcordLMC = c.cartesian.x ; ycordLMC = c.cartesian.y ; zcordLMC = c.cartesian.z
#The Sun is at (8.3, 0.027, 0) kpc from the center of the Milky Way and moving at (v_R, v_theta, v_z) = (-11.1, 232.24, 7.25) km/s (From Chen et al. 2001, Gillessen et al. 2009, Bovy 2015)
#We can use astropy to convert into a galactocentric frame as well as convert the proper motion values of the HVS3

#HVS initial conditions - Taken from Gaia Dr2
HVS3RADEC = SkyCoord("04:38:12.8 -54:33:12", unit=(u.hourangle, u.deg)) #Converting R.A and Dec of HVS3 to degrees - may have to check if this is correct and if we can convert in same line to MW frame
c = SkyCoord(ra=HVS3RADEC.ra.degree*u.degree, dec=HVS3RADEC.dec.degree*u.degree, distance=61*u.kpc) #Converting degrees into cartesian coordinates in the Sun's reference frame
xcordHVS3 = c.cartesian.x ; ycordHVS3 = c.cartesian.z ; zcordHVS3 = c.cartesian.z


print("Coordinates for LMC: ", xcordLMC, ",", ycordLMC, ",", zcordLMC)
print("Coordinates for HVS3: ", xcordHVS3, ",", ycordHVS3, ",", zcordHVS3)

#Defining an x,y,z space to plot in
D = 100
xMW = np.linspace(-D,D,1000) | units.kpc ; yMW = np.linspace(-D,D,1000) | units.kpc ; zMW = np.linspace(-D,D,1000) | units.kpc
x = np.linspace(-D,D,1000) | units.kpc ; y = np.linspace(-D,D,1000) | units.kpc ; z = np.linspace(-D,D,1000) | units.kpc

eps = 0.0 #Epsilon is the softening parameter and used to remove singularities - this is 1/N^(2/3) (AMUSE book) but I can't find a reliable source for the number of stars in LMC

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
rho0 = 8.18e6 | units.MSun/units.kpc**3 #This value was taken from SIffert et al. 2011
LMCPlummerProfile = Plummer_profile(mass, scaleB_radius) #The Plummer profile for LMC
LMCNFWProfile = NFW_profile(rho0, scaleDH_radius) #The NFW profile for LMC

#Calling the potentials
LMCPlummer = LMCPlummerProfile.get_potential_at_point(eps, x, y, z)
LMCNFW = LMCNFWProfile.get_potential_at_point(eps, x, y, z)
CombinedLMC = (LMCNFW + LMCPlummer)
pot_MW = MW.get_potential_at_point(eps, xMW, yMW, zMW)

plt.title("Milky Way and LMC Potentials")
plot(x, pot_MW, label = 'Milky Way Potential')
plot(x, CombinedLMC, label = 'Combined')
plt.xlim(0,D)
plt.legend()
plt.show()   

plt.title("LMC NFW Potential")
plot(x, LMCNFW, label = 'NFW')
plt.xlim(0,D)
plt.legend()
plt.show()   

plt.title("LMC Plummer Potential")
plot(x, LMCPlummer, label = 'Plummer Profile (LMC)')
plt.xlim(0,30)
plt.legend()
plt.show()   