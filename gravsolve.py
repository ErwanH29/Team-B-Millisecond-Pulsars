# This code needs galpot_init, particle_init_test, evol and conv_coord to run
from matplotlib import pyplot
from amuse.plot import plot
from evol import neut_gal_evol

X_gal, X_neut = neut_gal_evol()

plot(X_gal[0], X_gal[1], lw=1, lable = 'lmc')
plot(X_gal[3], X_gal[4], lw=1, lable = 'smc')
plot(X_neut[0], X_neut[1], lw=1, lable = 'neut')
#pyplot.gca().set_aspect("equal", adjustable="box")
pyplot.legend()
pyplot.show()
