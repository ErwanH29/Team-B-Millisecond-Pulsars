import numpy as np
from particle_init import gal_initializer
from amuse.lab import nbody_system, units, Particles
from amuse.units import quantities
from amuse.couple import bridge
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init import LMC_pot, SMC_pot

gal_code = gal_initializer()
LMC = LMC_pot()
SMC = SMC_pot()
MWG = MWpotentialBovy2015()
x_c, y_c, z_c = gal_code.ngc_1783_pos_init()
vel_init = LMC.circular_velocity(x_c, y_c, z_c).in_(units.kms)
gal_MC, ngc_1783 = gal_code.gal_tracer(vel_init)
ngc_1783.velocity += gal_MC[0].velocity
ngc_1783.position += gal_MC[0].position

ngc_1783.velocity *= -1
gal_MC[0].velocity *= -1
gal_MC[1].velocity *= -1

converter_1_MC = nbody_system.nbody_to_si(gal_MC.mass.sum(),
                                         gal_MC.position.sum())

converter_1_ngc = nbody_system.nbody_to_si(ngc_1783.mass.sum(),
                                         ngc_1783.position.sum())
from amuse.community.phigrape.interface import PhiGRAPEInterface, PhiGRAPE

gravity_code_1_MC = PhiGRAPE(converter_1_MC, MODE_GPU='gpu')
gravity_code_1_MC.particles.add_particles(gal_MC)
ch_g2l_1_MC = gravity_code_1_MC.particles.new_channel_to(gal_MC)

gravity_code_1_ngc = PhiGRAPE(converter_1_ngc, MODE_GPU='gpu')
gravity_code_1_ngc.particles.add_particles(ngc_1783)
ch_g2l_1_ngc = gravity_code_1_ngc.particles.new_channel_to(ngc_1783)

dt = 1 #Define timestep in Myr
        
gravity_1_MC = bridge.Bridge()
gravity_1_MC.add_system(gravity_code_1_MC, (MWG,))
gravity_1_MC.timestep = dt | units.Myr

gravity_1_ngc = bridge.Bridge()
gravity_1_ngc.add_system(gravity_code_1_ngc, (MWG, LMC, SMC))
gravity_1_ngc.timestep = dt | units.Myr

times = np.arange(0., 1000, dt) | units.Myr

for time in times:
    print(time)
    LMC.d_update(gal_MC[0].x, gal_MC[0].y, gal_MC[0].z)
    SMC.d_update(gal_MC[1].x, gal_MC[1].y, gal_MC[1].z)
    gravity_1_MC.evolve_model(time)
    gravity_1_ngc.evolve_model(time)
    ch_g2l_1_MC.copy()
    ch_g2l_1_ngc.copy()

gravity_1_MC.stop()
gravity_1_ngc.stop()

coord_LMC = gal_MC[0].position
vel_LMC = gal_MC[0].velocity
coord_SMC = gal_MC[1].position
vel_SMC = gal_MC[1].velocity
coord_ngc = ngc_1783.position
vel_ngc = ngc_1783.velocity
print('LMC coordinate 1000 Myr back in time is:{}'.format(coord_LMC))
print('LMC velocity 1000 Myr back in time is:{}'.format(vel_LMC))
print('SMC coordinate 1000 Myr back in time is:{}'.format(coord_SMC))
print('SMC velocity 1000 Myr back in time is:{}'.format(vel_SMC))
print('ngc coordinate 1000 Myr back in time is:{}'.format(coord_ngc))
print('ngc velocity 1000 Myr back in time is:{}'.format(vel_ngc))
