import numpy as np
from matplotlib import pyplot
from amuse.plot import plot
from amuse.lab import units, Particles, nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init import LMC_pot
from galpot_init import SMC_pot

# Initialize Bovy potential for Milky Way
MWG = MWpotentialBovy2015()

# Initialize LMC & SMC potential
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
converter = nbody_system.nbody_to_si(MCs.mass.sum(),
                                     MCs.position.length())

Lmc_gravity = Hermite(converter)
Lmc_gravity.particles.add_particles(MCs)

Smc_gravity = Hermite(converter)
Smc_gravity.particles.add_particles(MCs)

channel_lmc = Lmc_gravity.particles.new_channel_to(MCs)
channel_smc = Smc_gravity.particles.new_channel_to(MCs)

dt = 0.5  # Define time step

# sl_gravity = bridge.Bridge(use_threading=False)
# sl_gravity.add_system(Smc_gravity, (Lmc_gravity,))
# sl_gravity.add_system(Lmc_gravity, (Smc_gravity,))

gravity = bridge.Bridge()
gravity.add_system(Smc_gravity, (MWG, LMC,))
gravity.add_system(Lmc_gravity, (MWG, SMC,))
gravity.timestep = dt | units.Myr

times = np.arange(0., 100., dt) | units.Myr

x_lmc = [] | units.kpc
y_lmc = [] | units.kpc

x_smc = [] | units.kpc
y_smc = [] | units.kpc

# for time in times:
#    gravity.evolve_model(time)
#    channel_lmc.copy()
#    x_lmc.append(MCs[0].x)
#    y_lmc.append(MCs[0].y)

for time in times:
    gravity.evolve_model(time)
    channel_lmc.copy()
    channel_smc.copy()
    x_lmc.append(MCs[0].x)
    y_lmc.append(MCs[0].y)
    x_smc.append(MCs[1].x)
    y_smc.append(MCs[1].y)

# pyplot.subplot(121)
plot(x_lmc, y_lmc, lw=3, label='LMC')
# pyplot.subplot(122)
pyplot.scatter(0, 0)
plot(x_smc, y_smc, lw=3, label='SMC')
pyplot.legend()
pyplot.show()
