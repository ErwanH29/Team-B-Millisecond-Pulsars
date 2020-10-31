# This code needs galpot_init, particle_init_test, evol and conv_coord to run
from matplotlib import pyplot
from amuse.plot import *
from evol import neut_gal_evol
from capture_check import capture_check

X_gal, X_neut, N, E_ratio = neut_gal_evol()

pyplot.figure(figsize=(15,5)) 
pyplot.subplot(121)
plot(X_gal[0], X_gal[1], lw=0.7, label='LMC')
plot(X_gal[3], X_gal[4], lw=0.7, label='SMC')

check = capture_check(X_neut)

print(check)

#pyplot.gca().set_aspect("equal", adjustable="box")
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

#pyplot.colorbar()

pyplot.xlim(-20,20)
pyplot.ylim(-40,30)


for i in range(N):
    plot(X_neut[i][0], X_neut[i][1], lw=1, label='Pulsar {}'.format(i+1))

pyplot.legend()
pyplot.xlabel(r"$x$ coordinate (kpc)")
pyplot.ylabel(r"$y$ coordinate (kpc)")

pyplot.subplot(122)
pyplot.plot(E_ratio[0])
pyplot.savefig("EvolutionSystem", dpi=300)
pyplot.show()
