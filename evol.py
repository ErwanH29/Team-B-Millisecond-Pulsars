import numpy as np
from particle_init_test import gal_tracer, gal_path_init, neut_stars, neut_path_init
from amuse.lab import nbody_system, units
from amuse.couple import bridge
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init import LMC_pot, SMC_pot

def neut_gal_evol():
    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()
    
    gals = gal_tracer()
    converter_1 = nbody_system.nbody_to_si(gals.mass.sum(),
                                           gals.position.length())
    
    neuts=neut_stars(1)
    converter_2 = nbody_system.nbody_to_si(neuts.mass.sum(),
                                           neuts.position.length())
    
    from amuse.community.hermite.interface import Hermite
    
    # the code for the dynamic galactic potentials
    gravity_code_1 = Hermite(converter_1)
    gravity_code_1.particles.add_particles(gals)
    ch_g2l_1 = gravity_code_1.particles.new_channel_to(gals)
    
    # the code for the neutron stars
    gravity_code_2 = Hermite(converter_2)
    gravity_code_2.particles.add_particles(neuts)
    ch_g2l_2 = gravity_code_2.particles.new_channel_to(neuts)
    
    dt = 1 #Define timestep in Myr
    
    gravity_1 = bridge.Bridge()
    gravity_1.add_system(gravity_code_1, (MWG,))
    gravity_1.timestep = dt | units.Myr
    
    gravity_2 = bridge.Bridge()
    gravity_2.add_system(gravity_code_2, (MWG, LMC, SMC))
    gravity_2.timestep = dt | units.Myr
    
    times = np.arange(0., 2000, dt) | units.Myr
    
    X_gal = gal_path_init()
    X_neut = neut_path_init(1)
    
    for time in times:
        
        gravity_1.evolve_model(time)
        gravity_2.evolve_model(time)
        
        LMC.d_update(gals[0].x, gals[0].y, gals[0].z)
        SMC.d_update(gals[1].x, gals[1].y, gals[1].z)
        
        ch_g2l_1.copy()
        
        X_gal[0].append(gals[0].x)
        X_gal[1].append(gals[0].y)
        X_gal[2].append(gals[0].z)
        
        X_gal[3].append(gals[1].x)
        X_gal[4].append(gals[1].y)
        X_gal[5].append(gals[1].z)
        
        ch_g2l_2.copy()
        
        X_neut[0].append(neuts[0].x)
        X_neut[1].append(neuts[0].y)
        X_neut[2].append(neuts[0].z)

    gravity_1.stop()   
    gravity_2.stop()
    return  X_gal, X_neut
