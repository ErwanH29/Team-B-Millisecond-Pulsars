import numpy as np
from matplotlib import pyplot
from amuse.plot import plot
from amuse.lab import units, Particles, nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init_with_parameter import LMC_pot
from galpot_init_with_parameter import SMC_pot

# Initialize Bovy potential for Milky Way
MWG = MWpotentialBovy2015()
LMC = LMC_pot(79.88, -69.59, 49.97)
SMC = SMC_pot(13.16, -72.8, 60.6)

# MCs system (as particle)
MCs = Particles(2)
lmc = MCs[0]
lmc.mass = 2e10 | units.MSun
lmc.position = (-0.88, -41, -27) * (1 | units.kpc)
lmc.velocity = (26.3984409103, 266.255095135, -515.154438166) * (1 | units.kms)
##
smc = MCs[1]
smc.mass = 2e9 | units.MSun
smc.position = (15.0553251422, -36.443450871, -42.4034312546) * (1 | units.kpc)
smc.velocity = (377.103712014, 78.374279055, 132.87362692) * (1 | units.kms)
converter = nbody_system.nbody_to_si(MCs.mass.sum(), MCs.position.length())

Lmc_gravity = Hermite(converter)
Lmc_gravity.particles.add_particles(MCs)

Smc_gravity = Hermite(converter)
Smc_gravity.particles.add_particles(MCs)

channel_lmc = Lmc_gravity.particles.new_channel_to(MCs)
channel_smc = Smc_gravity.particles.new_channel_to(MCs)

dt = 0.5  # Define timestep

# initialize star to be ejected
star = Particles(1)
star.mass = 9 | units.MSun
star.position = (-0.88, -41, -27) * (1 | units.kpc)
star.velocity = (0, -5.0, 0) * (100 | units.kms)
converter = nbody_system.nbody_to_si(star.mass.sum(),
                                     star.position.length())

gravity_code = Hermite(converter)
gravity_code.particles.add_particles(star)
ch_g2l = gravity_code.particles.new_channel_to(star)

# initial position
x = [] | units.kpc
y = [] | units.kpc

gravity_system = bridge.Bridge()
gravity_system.timestep = dt | units.Myr
gravity = bridge.Bridge()
gravity.timestep = dt | units.Myr

times = np.arange(0., 20., dt) | units.Myr
for time in times:
    gravity.evolve_model(time)
    channel_lmc.copy()
    channel_smc.copy()
    ch_g2l.copy()

    gravity_system.add_system(Smc_gravity, (MWG, LMC,))
    gravity_system.add_system(Lmc_gravity, (MWG, SMC,))

    gravity.add_system(gravity_code, (MWG, LMC))

    x_lmc = MCs[0].x.value_in(units.kpc)
    y_lmc = MCs[0].y.value_in(units.kpc)
    z_lmc = MCs[0].z.value_in(units.kpc)
    LMC = LMC_pot(x_lmc, y_lmc, z_lmc, gc=True)

    x_smc = MCs[1].x.value_in(units.kpc)
    y_smc = MCs[1].y.value_in(units.kpc)
    z_smc = MCs[1].z.value_in(units.kpc)
    SMC = SMC_pot(x_smc, y_smc, z_smc, gc=True)

    x.append(star[0].x)
    y.append(star[0].y)

plot(x, y, lw=2)
# pyplot.scatter(0,0)
pyplot.show()