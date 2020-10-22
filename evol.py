import numpy as np
from particle_init import neut_initializer, gal_initializer
from amuse.lab import nbody_system, units, Particles
from amuse.units import quantities
from amuse.couple import bridge
from galpot_init import LMC_pot, SMC_pot, ngc_1783_pot
from amuse.ext.galactic_potentials import MWpotentialBovy2015

class drift_without_gravity(object):
    def __init__(self, convert_nbody, time= 0 |units.Myr):
        self.model_time = time
        self.convert_nbody = convert_nbody
        self.particles = Particles()
    def evolve_model(self, t_end):
        dt = t_end - self.model_time
        self.particles.position += self.particles.velocity*dt
        self.model_time = t_end
    @property
    def potential_energy(self):
        return quantities.zero
    @property 
    def kinetic_energy(self):
        return (0.5*self.particles.mass \
                   *self.particles.velocity.lengths()**2).sum()
    def stop(self):
        pass

def neut_gal_evol():
    neut_code = neut_initializer()
    gal_code = gal_initializer()
    N = neut_code.N

    L_neut = []
    for i in range(N):
        MWG = MWpotentialBovy2015()
        LMC = LMC_pot()
        SMC = SMC_pot()
        x_c, y_c, z_c = gal_code.ngc_1783_pos_init()
        vel_init = LMC.circular_velocity(x_c, y_c, z_c).in_(units.kms)
        print(vel_init)
        NGC_1783 = ngc_1783_pot()
        neuts= neut_code.neut_stars()
        gal_MC, ngc_1783 = gal_code.gal_tracer(vel_init)
        
        ngc_1783.velocity += gal_MC[0].velocity
        neuts.velocity += ngc_1783.velocity
        
        converter_1_MC = nbody_system.nbody_to_si(gal_MC.mass.sum(),
                                              gal_MC.position.sum())
        
        converter_1_ngc = nbody_system.nbody_to_si(ngc_1783.mass.sum(),
                                              ngc_1783.position.sum())
        
        converter_2 = nbody_system.nbody_to_si(neuts.mass.sum(),
                                              neuts.position.sum())


        from amuse.community.hermite.interface import Hermite
        
        # the code for the dynamic galactic potentials
        gravity_code_1_MC = Hermite(converter_1_MC)
        gravity_code_1_MC.particles.add_particles(gal_MC)
        ch_g2l_1_MC = gravity_code_1_MC.particles.new_channel_to(gal_MC)
        
        gravity_code_1_ngc = drift_without_gravity(converter_1_ngc)
        gravity_code_1_ngc.particles.add_particles(ngc_1783)
        ch_g2l_1_ngc = gravity_code_1_ngc.particles.new_channel_to(ngc_1783)

        # the code for the neutron stars
        gravity_code_2 = drift_without_gravity(converter_2)
        gravity_code_2.particles.add_particles(neuts)
        ch_g2l_2 = gravity_code_2.particles.new_channel_to(neuts)
        
        dt = 1 #Define timestep in Myr
        
        gravity_1_MC = bridge.Bridge()
        gravity_1_MC.add_system(gravity_code_1_MC, (MWG,))
        gravity_1_MC.timestep = dt | units.Myr
        
        gravity_1_ngc = bridge.Bridge()
        gravity_1_ngc.add_system(gravity_code_1_ngc, (MWG, LMC, SMC))
        gravity_1_ngc.timestep = dt | units.Myr

        gravity_2 = bridge.Bridge()
        gravity_2.add_system(gravity_code_2, (MWG, LMC, SMC, NGC_1783))
        gravity_2.timestep = dt | units.Myr
        
        times = np.arange(0., 500, dt) | units.Myr
        
        l_gal = gal_code.gal_path_init()
        l_neut = neut_code.neut_path_init()
        
        for time in times:
            
            LMC.d_update(gal_MC[0].x, gal_MC[0].y, gal_MC[0].z)
            SMC.d_update(gal_MC[1].x, gal_MC[1].y, gal_MC[1].z)
            NGC_1783.d_update(ngc_1783.x, ngc_1783.y, ngc_1783.z)
            
            gravity_1_MC.evolve_model(time)
            gravity_1_ngc.evolve_model(time)
            
            gravity_2.evolve_model(time)
            
            ch_g2l_1_MC.copy()
            
            l_gal[0].append(gal_MC[0].x)
            l_gal[1].append(gal_MC[0].y)
            l_gal[2].append(gal_MC[0].z)
            
            l_gal[3].append(gal_MC[1].x)
            l_gal[4].append(gal_MC[1].y)
            l_gal[5].append(gal_MC[1].z)
            
            ch_g2l_1_ngc.copy()
            
            l_gal[6].append(ngc_1783.x)
            l_gal[7].append(ngc_1783.y)
            l_gal[8].append(ngc_1783.z)
            
            ch_g2l_2.copy()
            
            l_neut[0].append(neuts.x)
            l_neut[1].append(neuts.y)
            l_neut[2].append(neuts.z)
        
        L_neut.append(l_neut)
        
    gravity_1_MC.stop()
    gravity_1_ngc.stop()
    gravity_2.stop()

    return  l_gal, L_neut, N
