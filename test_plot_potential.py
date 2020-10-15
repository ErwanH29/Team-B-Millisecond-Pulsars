from matplotlib import pyplot
import plot_potential as pp
from evol import neut_gal_evol

X_gal, X_neut = neut_gal_evol()

lmc_position_1 = [X_gal[0][250], X_gal[1][250], X_gal[2][250]]
smc_position_1 = [X_gal[3][250], X_gal[4][250], X_gal[5][250]]
z1 = pp.pplot(lmc_position_1, smc_position_1, 0.05, 10)

lmc_position_2 = [X_gal[0][10], X_gal[1][10], X_gal[2][10]]
smc_position_2 = [X_gal[3][10], X_gal[4][10], X_gal[5][10]]
z2 = pp.pplot(lmc_position_2, smc_position_2, 0.05, 10)

pyplot.imshow(z1-z2)
pyplot.show()
