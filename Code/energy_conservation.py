import numpy as np
from matplotlib import pyplot
from amuse.plot import *
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
    gal_code = gal_initializer()
    E = []
   

    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()
    lmc_mass = 2e10 | units.MSun 
    smc_mass = 6e9 | units.MSun

    x_c, y_c, z_c = gal_code.ngc_1783_pos_init()
    vel_init = LMC.circular_velocity(x_c, y_c, z_c).in_(units.kms)
    print(vel_init)
    NGC_1783 = ngc_1783_pot()
    gal_MC, ngc_1783 = gal_code.gal_tracer(vel_init)

    ngc_1783.velocity += gal_MC[0].velocity
    #neuts.velocity += ngc_1783.velocity

    converter_1_MC = nbody_system.nbody_to_si(gal_MC.mass.sum(),
                                          gal_MC.position.sum())

    converter_1_ngc = nbody_system.nbody_to_si(ngc_1783.mass.sum(),
                                          ngc_1783.position.sum())

    En = []

    from amuse.community.hermite.interface import Hermite

    # the code for the dynamic galactic potentials
    gravity_code_1_MC = Hermite(converter_1_MC)
    gravity_code_1_MC.particles.add_particles(gal_MC)
    ch_g2l_1_MC = gravity_code_1_MC.particles.new_channel_to(gal_MC)

    gravity_code_1_ngc = drift_without_gravity(converter_1_ngc)
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

    l_gal = gal_code.gal_path_init()
 
    E_init = gravity_1_MC.kinetic_energy.value_in(units.m**2 * units.s**-2 * units.kg) +\
          lmc_mass.value_in(units.kg) *\
          (MWG.get_potential_at_point(0,gal_MC[0].x, gal_MC[0].y, gal_MC[0].z).value_in(units.m**2 * units.s**-2)+\
          SMC.get_potential_at_point(0,gal_MC[0].x, gal_MC[0].y, gal_MC[0].z).value_in(units.m**2 * units.s**-2)) +\
          smc_mass.value_in(units.kg) *\
          (MWG.get_potential_at_point(0,gal_MC[1].x, gal_MC[1].y, gal_MC[1].z).value_in(units.m**2 * units.s**-2)+\
          LMC.get_potential_at_point(0,gal_MC[1].x, gal_MC[1].y, gal_MC[1].z).value_in(units.m**2 * units.s**-2))

    for time in times:

        E_t = gravity_1_MC.kinetic_energy.value_in(units.m**2 * units.s**-2 * units.kg) +\
          lmc_mass.value_in(units.kg) *\
          (MWG.get_potential_at_point(0,gal_MC[0].x, gal_MC[0].y, gal_MC[0].z).value_in(units.m**2 * units.s**-2)+\
          SMC.get_potential_at_point(0,gal_MC[0].x, gal_MC[0].y, gal_MC[0].z).value_in(units.m**2 * units.s**-2))+\
          smc_mass.value_in(units.kg) *\
          (MWG.get_potential_at_point(0,gal_MC[1].x, gal_MC[1].y, gal_MC[1].z).value_in(units.m**2 * units.s**-2)+\
          LMC.get_potential_at_point(0,gal_MC[1].x, gal_MC[1].y, gal_MC[1].z).value_in(units.m**2 * units.s**-2))
        En.append(E_t / E_init)

        LMC.d_update(gal_MC[0].x, gal_MC[0].y, gal_MC[0].z)
        SMC.d_update(gal_MC[1].x, gal_MC[1].y, gal_MC[1].z)
        NGC_1783.d_update(ngc_1783.x, ngc_1783.y, ngc_1783.z)

        gravity_1_MC.evolve_model(time)
        gravity_1_ngc.evolve_model(time)

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

    E.append(En)
    gravity_1_MC.stop()
    gravity_1_ngc.stop()

    return  l_gal, E

X_gal, E_ratio = neut_gal_evol()

pyplot.figure(figsize=(15,5)) 
pyplot.subplot(121)
pyplot.title('Trajectory of MCs')
plot(X_gal[0][0:126], X_gal[1][0:126], lw=1.7, label='LMC_BEFORE')
plot(X_gal[0][126:1000], X_gal[1][126:1000], lw=1.7, label='LMC_AFTER')
plot(X_gal[3][0:126], X_gal[4][0:126], lw=1.7, label='SMC_BEFORE')
plot(X_gal[3][126:1000], X_gal[4][126:1000], lw=1.7, label='SMC_AFTER')

# plot the isocontour 
from amuse.ext.galactic_potentials import MWpotentialBovy2015
MW = MWpotentialBovy2015()
omega = 600 | units.km/units.s/units.km
effective_iso_potential_plot(MW, omega, center_of_rotation = [0, 0]|units.kpc,
                                  xlim = [-10, 10] | units.kpc, ylim = [-10, 10] | units.kpc,
                                  resolution = [1000, 1000], fraction_screen_filled = 0.8,
                                  quadratic_contour_levels = True, contour_kwargs = dict(linestyles = 'solid',
                                  linewidths=0.5, cmap=('Greys')), omega2 = None, center_of_rotation2 = [0, 0]|units.kpc, 
                                  fraction_screen_filled2 = 0.2, projection3D = False, number_of_contours = 15)

pyplot.xlim(-20,20)
pyplot.ylim(-40,30)

pyplot.legend()
pyplot.xlabel(r"$x$ coordinate (kpc)")
pyplot.ylabel(r"$y$ coordinate (kpc)")

pyplot.subplot(122)
pyplot.title('E_t / E_init')
pyplot.plot(E_ratio[0])
pyplot.ylim(0.8,1.2)
pyplot.xlabel("time(Myr)")
pyplot.savefig("EvolutionSystem", dpi=300)
pyplot.show()
