import numpy as np
from particle_init import neut_initializer, gal_initializer
from amuse.lab import nbody_system, units
from amuse.couple import bridge
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init import LMC_pot, SMC_pot
from amuse.plot import *
from matplotlib import pyplot

def neut_gal_evol():
    neut_code = neut_initializer()
    gal_code = gal_initializer()
    N = neut_code.N
    neuts = neut_code.neut_stars()
    L_neut = []
    for i in range(N):
        MWG = MWpotentialBovy2015()
        LMC = LMC_pot()
        SMC = SMC_pot()
        
        gals = gal_code.gal_tracer()
        converter_1 = nbody_system.nbody_to_si(gals.mass.sum(),
                                               gals.position.length())
         
        converter_2 = nbody_system.nbody_to_si(neuts[i].mass.sum(),
                                               neuts[i].position.length())
        
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
        
        times = np.arange(0., 1000, dt) | units.Myr
        
        l_gal = gal_code.gal_path_init()
        l_neut = neut_code.neut_path_init()
        
        for time in times:
            
            LMC.d_update(gals[0].x, gals[0].y, gals[0].z)
            SMC.d_update(gals[1].x, gals[1].y, gals[1].z)
            
            gravity_1.evolve_model(time)
            gravity_2.evolve_model(time)
            
            ch_g2l_1.copy()
            
            l_gal[0].append(gals[0].x)
            l_gal[1].append(gals[0].y)
            l_gal[2].append(gals[0].z)
            
            l_gal[3].append(gals[1].x)
            l_gal[4].append(gals[1].y)
            l_gal[5].append(gals[1].z)
            
            ch_g2l_2.copy()
            
            l_neut[0].append(neuts[i].x)
            l_neut[1].append(neuts[i].y)
            l_neut[2].append(neuts[i].z)
        
        print(l_neut)
        L_neut.append(l_neut)
        gravity_1.stop()   
        gravity_2.stop()
        
    return  l_gal, L_neut
