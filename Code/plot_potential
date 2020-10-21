import matplotlib.pyplot as plt
import numpy as np
from galpot_init import LMC_pot, SMC_pot
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from amuse.lab import units


def pplot(lmc_position, smc_position, grid_size_in_kpc, xy_bound):

    '''
    :param lmc_position: given the position of LMC with form as (x, y, z) in gc
    :param smc_position: given the position of SMC with form as (x, y, z) in gc
    :param grid_size_in_kpc: set grid size in units.kpc
    :param xy_bound: generate 2d grids for the x & y bounds (square image)
    '''
    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()

    # the position of LMC and SMC
    x_lmc, y_lmc, z_lmc = lmc_position[0], lmc_position[1], lmc_position[2]
    x_smc, y_smc, z_smc = smc_position[0], smc_position[1], smc_position[2]

    LMC.d_update(x_lmc, y_lmc, z_lmc)
    SMC.d_update(x_smc, y_smc, z_smc)

    # set grids
    dx, dy = grid_size_in_kpc, grid_size_in_kpc
    y, x = np.mgrid[slice(-xy_bound, xy_bound + dy, dy), slice(-xy_bound, xy_bound + dx, dx)]
    uy, ux = y | units.kpc, x | units.kpc
    z = MWG.get_potential_at_point(0, uy, ux, 0 | units.kpc).value_in(units.m**2 * units.s**-2) + \
        LMC.get_potential_at_point(0, uy, ux, 0 | units.kpc).value_in(units.m**2 * units.s**-2) + \
        SMC.get_potential_at_point(0, uy, ux, 0 | units.kpc).value_in(units.m**2 * units.s**-2)
    # remove the bound value
    z = z[:-1, :-1]

    plt.imshow(z)
    plt.colorbar()

    return z
