from amuse.lab import Particles, units
from conv_coord import conv_coord

# This script is to test the updating potential fields of LMC and SMC with 1 neutron star

# initialize particle sets here. 
##### NEUTRON STAR SET #####

def neut_stars(N_neut):
    neuts = Particles(1)
    neuts.mass = 2 | units.MSun 
    neuts.position = (-0.88, -41, -27) * (1 | units.kpc)
    neuts.velocity = (-200, 0, 0) * (1 | units.kms)
    return neuts

def neut_path_init(N_neut):
    x_neut = [] | units.kpc
    y_neut = [] | units.kpc
    z_neut = [] | units.kpc
    return x_neut, y_neut, z_neut

##### GALAXY TRACERS #####

# initialize LMC and SMC
def gal_tracer():
    gals = Particles(2)
    gal_lmc = gals[0]
    gal_lmc.mass = 2e10 | units.MSun 
    gal_lmc.position = conv_coord(80.8958, -69.7561, 49.97)
    print(gal_lmc.position)
    gal_lmc.velocity = (47,242,225) * (1 | units.kms)

    gal_smc = gals[1]
    gal_smc.mass = 4e9 | units.MSun 
    gal_smc.position = conv_coord(13.1583, -72.8, 61.7)
    gal_smc.velocity = (5.35,164,136) * (1 | units.kms)
    return gals

def gal_path_init():
     x_lmc = [] | units.kpc
     y_lmc = [] | units.kpc
     z_lmc = [] | units.kpc
    
     x_smc = [] | units.kpc
     y_smc = [] | units.kpc
     z_smc = [] | units.kpc
     return x_lmc, y_lmc, z_lmc, x_smc, y_smc, z_smc
