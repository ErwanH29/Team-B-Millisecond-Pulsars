import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units, constants
from amuse.plot import plot
from astropy import units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from amuse.ext.galactic_potentials import MWpotentialBovy2015, MiyamotoNagai_profile, NFW_profile, Plummer_profile

def potentialBH(m,x,y,z):
    return (-constants.G * mass_bh)/(x**2+y**2+z**2)**0.5

#Initialising LMC
cLMC = SkyCoord(ra=79.88*u.degree, dec=-69.59*u.degree, distance=49.97*u.kpc) #Data of LMC distance taken from Pieterzxynski et al 2013 in the Sun's reference frame
gcLMC = cLMC.transform_to(coord.Galactocentric)
print("LMC coordinates in galactocentric frame: \n x = ", gcLMC.cartesian.x, ", y = ", gcLMC.cartesian.y, ", z = ", gcLMC.cartesian.z)

#HVS initial conditions - Taken from Gaia Dr2
HVS3RADEC = SkyCoord("04:38:12.8 -54:33:12", unit=(u.hourangle, u.deg)) #Converting R.A and Dec of HVS3 to degrees
cS = coord.ICRS(ra=HVS3RADEC.ra.degree*u.degree, dec=HVS3RADEC.dec.degree*u.degree, 
                distance=61*u.kpc, pm_ra_cosdec=0.851*u.mas/u.yr, 
                pm_dec=1.936*u.mas/u.yr, radial_velocity=723*u.km/u.s)
gcS = cS.transform_to(coord.Galactocentric) #Setting HVS3 coordinates with respect to the center of the MW
print("HVS3 coordinates in galactocentric frame: \n x = ", gcS.cartesian.x, ", y = ", gcS.cartesian.y, ", z = ", gcS.cartesian.z)
print("HVS3 velocity in galactocentric frame: \n vx = ", gcS.v_x, ", vy = ", gcS.v_y, ", vz = ", gcS.v_z)

#Defining an x,y,z space to plot in
D = 100
xMW = np.linspace(-D,D,1000) | units.kpc ; yMW = np.linspace(-D,D,1000) | units.kpc ; zMW = np.linspace(-D,D,1000) | units.kpc
x = np.linspace(-D,D,1000) | units.kpc ; y = np.linspace(-D,D,1000) | units.kpc ; z = np.linspace(-D,D,1000) | units.kpc

eps = 0.0 #Epsilon is the softening parameter and used to remove singularities - this is 1/N^(2/3) (AMUSE book) but I can't find a reliable source for the number of stars in LMC

#Milky Way Statistics - Taken from Bovy 2015
disk_mass = 6.8e+10 | units.MSun
disk_size = 3 | units.kpc #Milky Way's scale radius
disk_height = 280 | units.parsec #Milky Way's scale height
mass_bh = 4*10**6 | units.MSun #Milky Way's central black hole mass
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
MWPotential2014wBH= pot_MW + potentialBH(mass_bh,xMW,yMW,zMW)

#Getting the circular velocity as a function of distance for the different potentials
r = (x**2+y**2+z**2)**0.5
vcircMW = MW.circular_velocity(r).in_(units.km/units.s)
vcircLMC = LMCPlummerProfile .circular_velocity(r).in_(units.km/units.s)+LMCNFWProfile.circular_velocity(r).in_(units.km/units.s)



#PLOTTING THE ORBITAL VELOCITY AROUND GALAXIES AND THE POTENTIALS
plt.title(r"Milky Way Circular Velocity ($v_{circ}$) vs Distance ($r$)")
plt.xlabel(r"Distance ($r$)")
plt.ylabel(r"Velocity $v_{circ}$")
plot(x, vcircMW, label = 'Milky Way Potential')
plot(x, vcircLMC, label = 'LMC Potential Potential')
plt.xlim(0,20)
plt.legend()
plt.savefig("MW+LMC Circular Velocity", dpi=300)
plt.show()  

plt.title(r"Milky Way Circular Velocity ($v_{circ}$) vs Distance ($r$)")
plt.xlabel(r"Distance ($r$)")
plt.ylabel(r"Velocity $v_{circ}$")
plot(x, vcircMW, label = 'Milky Way Potential')
plt.xlim(0,50)
plt.legend()
plt.savefig("MW Circular Velocity", dpi=300)
plt.show() 

plt.title(r"LMC Circular Velocity ($v_{circ}$) vs Distance ($r$)")
plt.xlabel(r"Distance ($r$)")
plt.ylabel(r"Velocity $v_{circ}$")
plot(x, vcircLMC, label = 'LMC Potential')
plt.xlim(0,25)
plt.legend()
plt.savefig("LMC Circular Velocity", dpi=300)
plt.show()  

plt.title("Milky Way and LMC Potentials")
plot(x, pot_MW, label = 'Milky Way Potential')
plot(x, CombinedLMC, label = 'Combined')
plt.xlim(0,D)
plt.legend()
plt.savefig("MW and LMC Potentials", dpi=300)
plt.show()   

plt.title("LMC NFW Potential")
plot(x, LMCNFW, label = 'NFW')
plt.xlim(0,D)
plt.legend()
plt.savefig("LMC NFW Potential", dpi=300)
plt.show()   

plt.title("LMC Plummer Potential")
plot(x, LMCPlummer, label = 'Plummer Profile (LMC)')
plt.xlim(0,30)
plt.legend()
plt.savefig("LMC Plummer Potential", dpi=300)
plt.show()  

plt.title("MW Black Hole Potential")
plot(x, potentialBH(mass_bh,xMW,yMW,zMW), label = 'Black Hole Potential')
plt.xlim(0,30)
plt.legend()
plt.savefig("BHPotential", dpi=300)
plt.show()  