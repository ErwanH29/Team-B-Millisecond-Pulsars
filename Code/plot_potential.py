import matplotlib.pyplot as plt
import numpy as np
from galpot_init import LMC_pot, SMC_pot
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from amuse.lab import units


def pplot(gal_MC, grid_size_in_kpc, xy_bound, i, **options):
    
    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()
    
    if options.get('fix') == 'x':
        
        LMC.d_update((0.0 | units.kpc), gal_MC[0].y, gal_MC[0].z)
        SMC.d_update((0.0 | units.kpc), gal_MC[1].y, gal_MC[1].z)
        
        dx, dy = grid_size_in_kpc, grid_size_in_kpc
        y, x = np.mgrid[slice(-xy_bound, xy_bound + dy, dy), slice(-xy_bound, xy_bound + dx, dx)]
        uy, ux = y | units.kpc, x | units.kpcy
        p = MWG.get_potential_at_point(0, 0 | units.kpc, ux, uy).value_in(units.m**2 * units.s**-2) + \
            LMC.get_potential_at_point(0, 0 | units.kpc, ux, uy).value_in(units.m**2 * units.s**-2) + \
            SMC.get_potential_at_point(0, 0 | units.kpc, ux, uy).value_in(units.m**2 * units.s**-2)
        # remove the bound value
        p = p[:-1, :-1]
        
    if options.get('fix') == 'y':
        LMC.d_update(gal_MC[0].x, (0.0 | units.kpc), gal_MC[0].z)
        SMC.d_update(gal_MC[1].x, (0.0 | units.kpc), gal_MC[1].z)
        # set grids
        dx, dy = grid_size_in_kpc, grid_size_in_kpc
        y, x = np.mgrid[slice(-xy_bound, xy_bound + dy, dy), slice(-xy_bound, xy_bound + dx, dx)]
        uy, ux = y | units.kpc, x | units.kpc
        p = MWG.get_potential_at_point(0, ux, 0 | units.kpc, uy).value_in(units.m**2 * units.s**-2) + \
            LMC.get_potential_at_point(0, ux, 0 | units.kpc, uy).value_in(units.m**2 * units.s**-2) + \
            SMC.get_potential_at_point(0, ux, 0 | units.kpc, uy).value_in(units.m**2 * units.s**-2)
        # remove the bound value
        p = p[:-1, :-1]
        
    if options.get('fix') == 'z':
        LMC.d_update(gal_MC[0].x, gal_MC[0].y, (0.0 | units.kpc))
        SMC.d_update(gal_MC[1].x, gal_MC[1].y, (0.0 | units.kpc))
        # set grids
        dx, dy = grid_size_in_kpc, grid_size_in_kpc
        y, x = np.mgrid[slice(-xy_bound, xy_bound + dy, dy), slice(-xy_bound, xy_bound + dx, dx)]
        uy, ux = y | units.kpc, x | units.kpc
        p = MWG.get_potential_at_point(0, ux, uy, 0 | units.kpc).value_in(units.m**2 * units.s**-2) + \
            LMC.get_potential_at_point(0, ux, uy, 0 | units.kpc).value_in(units.m**2 * units.s**-2) + \
            SMC.get_potential_at_point(0, ux, uy, 0 | units.kpc).value_in(units.m**2 * units.s**-2)
        # remove the bound value
        p = p[:-1, :-1]
        
    plt.imshow(p,aspect='auto',origin='lower',extent=(x.min(),x.max(),y.min(),y.max()))
    
    if options.get('fix') == 'z':
        plt.ylabel('y in [kpc]')
        plt.xlabel('x in [kpc]')
        
    if options.get('fix') == 'x':
        plt.ylabel('z in [kpc]')
        plt.xlabel('y in [kpc]')
        
    if options.get('fix') == 'y':
        plt.ylabel('z in [kpc]')
        plt.xlabel('x in [kpc]')
    print('potential@{}'.format(i))
    plt.savefig('pot_plots/pot_evol_{}.png'.format(i),dpi=300)

    return p
