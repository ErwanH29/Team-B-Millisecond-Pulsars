import numpy as np
import matplotlib.pyplot as plt
from amuse.units import units, constants
from amuse import plot
from astropy import units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from amuse.ext.galactic_potentials import MWpotentialBovy2015

MW = MWpotentialBovy2015()
omega = 600 | units.km/units.s/units.km

plot.effective_iso_potential_plot(MW, omega, center_of_rotation = [0, 0]|units.AU,
                                  xlim = [-1.5, 1.5] | units.AU, ylim = [-1.5, 1.5] | units.AU,
                                  resolution = [1000, 1000], fraction_screen_filled = 0.33,
                                  quadratic_contour_levels = True, contour_kwargs = dict(linestyles = 'solid', 
                                  cmap=('Greys')), omega2 = None, center_of_rotation2 = [0, 0]|units.AU, 
                                  fraction_screen_filled2 = 0.2, projection3D = False)
number_of_contours = 10,
plt.colorbar()
plt.savefig("Isopotential plot of MWBovy2015()")
plt.show()