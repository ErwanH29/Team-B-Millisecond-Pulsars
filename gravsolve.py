#Code to calculate orbit in a potential. Play around with input parameters and see what happens...
import numpy as np
from matplotlib import pyplot
from amuse.plot import plot
from amuse.lab import units, Particles, nbody_system
from amuse.couple import bridge
from amuse.community.hermite.interface import Hermite
from amuse.ext.galactic_potentials import MWpotentialBovy2015, NFW_profile, Plummer_profile

#Initialize Bovy potential for Milky Way
MWG = MWpotentialBovy2015()

#Large Magellanic Cloud Statistics - Taken from Bekki 2005:
mass = 2e10 | units.MSun #LMC mass
scaleDH_radius = 2.6 | units.kpc #LMC dark matter halo scale radius
scaleB_radius = 0.73 | units.kpc #LMC bulge scale height
rho0 = 8.18e6 | units.MSun/units.kpc**3 #This value was taken from SIffert et al. 2011
LMCPlummerProfile = Plummer_profile(mass, scaleB_radius) #The Plummer profile for LMC
LMCNFWProfile = NFW_profile(rho0, scaleDH_radius) #The NFW profile for LMC

#initialize star to be ejected
star = Particles(1)
star.mass = 9 | units.MSun
star.position = (1.0,0,0) * (8.5 | units.kpc)
star.velocity = (0,-1.0,0) * (100 | units.kms)
converter = nbody_system.nbody_to_si(star.mass.sum(),
                                     star.position.length())

gravity_code = Hermite(converter)
gravity_code.particles.add_particles(star)
ch_g2l = gravity_code.particles.new_channel_to(star)

dt = 1 #Define timestep
gravity = bridge.Bridge(use_threading=False)
gravity.add_system(gravity_code, (MWG,))
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