#Code to calculate orbit in a potential. Play around with input parameters and see what happens...
import numpy as np
from matplotlib import pyplot
from amuse.plot import plot
from amuse.lab import units, Particles, nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from galpot_init import LMC_pot

#Initialize Bovy potential for Milky Way
MWG = MWpotentialBovy2015()

#Initialize LMC potential
LMC = LMC_pot()

#initialize star to be ejected
star = Particles(1)
star.mass = 9 | units.MSun
star.position = (-0.88,-41,-27) * (1 | units.kpc)
star.velocity = (0,-1.0,0) * (200 | units.kms)
converter = nbody_system.nbody_to_si(star.mass.sum(),
                                     star.position.length())

gravity_code = Hermite(converter)
gravity_code.particles.add_particles(star)
ch_g2l = gravity_code.particles.new_channel_to(star)

dt = 1 #Define timestep
gravity = bridge.Bridge(use_threading=False)
gravity.add_system(gravity_code, (MWG, LMC))
gravity.timestep = dt | units.Myr

times = np.arange(0., 100, dt) | units.Myr
x = [] | units.kpc
y = [] | units.kpc
for time in times:
    gravity.evolve_model(time)
    ch_g2l.copy()
    x.append(star[0].x)
    y.append(star[0].y)
gravity.stop()

plot(x, y, lw=1)
#pyplot.gca().set_aspect("equal", adjustable="box")
pyplot.show()
