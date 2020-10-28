# This code needs galpot_init, particle_init_test, evol and conv_coord to run
from matplotlib import pyplot
from amuse.plot import *
from evol import neut_gal_evol
from capture_check import capture_check
from plot_potential import pplot
import pandas as pd

X_gal, X_neut, N = neut_gal_evol(old_pos='old')

plot(X_gal[0], X_gal[1], lw=0.7, label='LMC')
plot(X_gal[3], X_gal[4], lw=0.7, label='SMC')
plot(X_gal[6], X_gal[7], lw=0.7, label='NGC 1783')

check = capture_check(X_neut)

print(check)

#pyplot.gca().set_aspect("equal", adjustable="box")
# plot the isocontour 
from amuse.ext.galactic_potentials import MWpotentialBovy2015
MW = MWpotentialBovy2015()
omega = 600 | units.km/units.s/units.km
effective_iso_potential_plot(MW, omega, center_of_rotation = [0, 0] | units.kpc,
                                  xlim = [-10, 10] | units.kpc, ylim = [-10, 10] | units.kpc,
                                  resolution = [1000, 1000], fraction_screen_filled = 0.8,
                                  quadratic_contour_levels = True, contour_kwargs = dict(linestyles = 'solid',
                                  linewidths=0.5, cmap=('Greys')), omega2 = None, center_of_rotation2 = [0, 0]|units.kpc, 
                                  fraction_screen_filled2 = 0.2, projection3D = False, number_of_contours = 15)

#pyplot.colorbar()

pyplot.xlabel(r"$x$ coordinate (kpc)")
pyplot.ylabel(r"$y$ coordinate (kpc)")
pyplot.title("MSP's captured by Milky Way for N={} ejected MSP's".format(N))
#pyplot.xlim(-20,20)
#pyplot.ylim(-50,50)

pyplot.legend()
pyplot.savefig("EvolutionSystem", dpi=300)
pyplot.show()

#pplot(X_gal, X_gal, 0.05, 50, 1)