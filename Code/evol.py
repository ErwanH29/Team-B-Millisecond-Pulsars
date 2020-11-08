import numpy as np
from coordinates import old_coordlist, capture_check
from particle_init import neut_initializer, gal_initializer
from amuse.lab import nbody_system, units, Particles
from amuse.units import quantities
from amuse.couple import bridge
from galpot_init import LMC_pot, SMC_pot, ngc_1783_pot
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from amuse.community.hermite.interface import Hermite
import pandas as pd
import pickle
import json

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

def neut_gal_evol(number_of_workers, **options):
    neut_code = neut_initializer()
    gal_code = gal_initializer()

    old_coord = old_coordlist()

    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()
    NGC_1783 = ngc_1783_pot()
    
    x_c, y_c, z_c = gal_code.ngc_1783_pos_init()
    vel_init = LMC.circular_velocity(x_c, y_c, z_c).in_(units.kms)
    gal_MC, ngc_1783 = gal_code.gal_tracer(vel_init)
    
    if options.get('old_pos') == 'current':
        
        ngc_1783.velocity += gal_MC[0].velocity
        ngc_1783.position += gal_MC[0].position
        
    if options.get('old_pos') == 'old':
        
        gal_MC[0].position = old_coord[0]
        gal_MC[0].velocity = -1*old_coord[1]
        gal_MC[1].position = old_coord[2]
        gal_MC[1].velocity = -1*old_coord[3]
        ngc_1783.position = old_coord[4]
        ngc_1783.velocity += gal_MC[0].velocity
        
    neuts = neut_code.neut_star_init(ngc_1783.position, ngc_1783.velocity)
   
    converter_1_MC = nbody_system.nbody_to_si(gal_MC.mass.sum(),
                                          gal_MC.position.sum())
    
    converter_1_ngc = nbody_system.nbody_to_si(ngc_1783.mass.sum(),
                                          ngc_1783.position.sum())
    
    converter_2 = nbody_system.nbody_to_si(neuts.mass.sum(),
                                           neuts.position.sum())

    gravity_code_1_MC = Hermite(converter_1_MC)
    gravity_code_1_MC.particles.add_particles(gal_MC)
    ch_g2l_1_MC = gravity_code_1_MC.particles.new_channel_to(gal_MC)
    
    gravity_code_1_ngc = drift_without_gravity(converter_1_ngc)
    gravity_code_1_ngc.particles.add_particles(ngc_1783)
    ch_g2l_1_ngc = gravity_code_1_ngc.particles.new_channel_to(ngc_1783)
    
    gravity_code_2 = Hermite(converter_2, number_of_workers=number_of_workers)
    gravity_code_2.particles.add_particles(neuts)
    ch_g2l_2 = gravity_code_2.particles.new_channel_to(neuts)
      
    dt = 1 #Define timestep in Myr
    
    gravity = bridge.Bridge()
    
    gravity.add_system(gravity_code_1_MC, (MWG,))
    gravity.add_system(gravity_code_1_ngc, (MWG, LMC, SMC))
    gravity.add_system(gravity_code_2, (MWG, LMC, SMC, NGC_1783))
    
    gravity.timestep = dt | units.Myr
    
    times = np.arange(0., 1000, dt) | units.Myr
    
    l_gal = gal_code.gal_path_init()
    
    neut_line = pd.DataFrame()

    for time in times:
        
        print(time)
        N = neut_code.N_count()
        
        LMC.d_update(gal_MC[0].x, gal_MC[0].y, gal_MC[0].z)
        SMC.d_update(gal_MC[1].x, gal_MC[1].y, gal_MC[1].z)
        NGC_1783.d_update(ngc_1783.x, ngc_1783.y, ngc_1783.z)
        
        gravity.evolve_model(time)
        
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

        df_conc = pd.DataFrame()
        rows = len(neut_line.index)
        ch_g2l_2.copy()
        for i in range(N):
                
            df = pd.Series({'{}'.format(time): [neuts[i].position]})
            df_conc = df_conc.append(df, ignore_index=True) 
        neut_line = neut_line.append(df_conc, ignore_index = True)
        neut_line['{}'.format(time)] = neut_line['{}'.format(time)].shift(-rows)

        #pplot(gal_MC, 0.05, 50, time, fix='y') # for plotting GIF files. WARNING: takes long!!
        if neut_code.decision(time)==True and time <= (1000 | units.Myr) :
            add_n = neut_code.add_neut(ngc_1783.position, ngc_1783.velocity)
            neuts.add_particle(add_n) 
            gravity_code_2.particles.add_particles(add_n)

    gravity.stop()
    
    neut_line = neut_line.dropna(thresh=1)
    check = capture_check(neut_line)
    f = open('working_data/check.txt', 'w')
    json.dump(check, f)
    f.close()
    
    t = open('working_data/total.txt', 'w')
    json.dump(N-1, t)
    t.close()

    neut_line.to_pickle('working_data/neut_stars_positions.pkl')
    
    with open('working_data/gal_line.pickle', 'wb') as f:
        pickle.dump(l_gal, f)
        
    return neut_line