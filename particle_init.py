from amuse.lab import Particles, units

# initialize particle sets here. 
##### NEUTRON STAR SET #####

def neut_stars(N_neut):
    return neuts

def neut_path_init(N_neut):
    for i in range(N_neut):
        x_neut = [] | units.kpc
        y_neut = [] | units.kpc
        z_neut = [] | units.kpc
        #X_neut[i] = (x_neut, y_neut, z_neut)
    return 

neuts = neut_path_init(10)
##### GALAXY TRACERS #####

# initialize LMC and SMC
def gal_tracer():
    gals = Particles(2)
    gal_lmc = gals[0]
    gal_lmc.mass = 2e10 | units.MSun 
    gal_lmc.position = (-0.88,-41,-27) * (1 | units.kpc)
    gal_lmc.velocity = (47,242,225) * (1 | units.kms)

    gal_smc = gals[1]
    gal_smc.mass = 2e10 | units.MSun 
    gal_smc.position = (15,-36,-42) * (1 | units.kpc)
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