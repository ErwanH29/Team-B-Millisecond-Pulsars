import matplotlib.pyplot as plt
import numpy as np
from galpot_init import LMC_pot, SMC_pot
from amuse.ext.galactic_potentials import MWpotentialBovy2015
from amuse.lab import units
from amuse.plot import *
import ast
import pandas as pd
import pickle

def pplot(gal_MC, grid_size_in_kpc, xy_bound, **options):
    
    MWG = MWpotentialBovy2015()
    LMC = LMC_pot()
    SMC = SMC_pot()
    
    if options.get('fix') == 'x':
        
        LMC.d_update((0.0 | units.kpc), gal_MC[0].y, gal_MC[0].z)
        SMC.d_update((0.0 | units.kpc), gal_MC[1].y, gal_MC[1].z)
        
        dx, dy = grid_size_in_kpc, grid_size_in_kpc
        y, x = np.mgrid[slice(-xy_bound, xy_bound + dy, dy), slice(-xy_bound, xy_bound + dx, dx)]
        uy, ux = y | units.kpc, x | units.kpc
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
    plt.ylabel('z in [kpc]')
    plt.xlabel('x in [kpc]')
    print('potential@{}'.format(i))
    plt.savefig('pot_plots/pot_evol_{}.png'.format(i),dpi=300)

    return p

class plot_neut_pos(object):
    
    def __init__(self):
        self.dat = pd.DataFrame()

    def get_dataframe(self, dat):
        self.dat = dat
    
    def get_dataframe_from_pkl(self, dat_str):
        self.dat = pd.read_pickle(dat_str)

    def make_line(self, i):
        neut_stars = self.dat
        neut_star = neut_stars.iloc[i]
        neut_star = neut_star.replace(np.NaN, "[0, 0, 0]")
        col_len = len(neut_star)
        line_x = []
        line_y = []
        line_z = []
        for j in range(col_len):
            
            inter_t = neut_star.iloc[j]
            if inter_t == '[0, 0, 0]':
                coord = ast.literal_eval(inter_t)
                line_x.append(coord[0])
                line_y.append(coord[0])
                line_z.append(coord[0])
             
            else:
                coord = neut_star.iloc[j][0] * (1 | units.kpc**-1)
                line_x.append(coord[0])
                line_y.append(coord[1])
                line_z.append(coord[2])

        #print(line_x)
        return line_x, line_y, line_z
    
    def zero_to_nan(self,values):
        """Replace every 0 with 'nan' and return a copy."""
        return [float('nan') if x==0 else x for x in values]
    
    def plot(self, use_all, **options):
        row_len = len(self.dat.index)
        if use_all=='y':
            for i in range(row_len):
                line_x, line_y, line_z = self.make_line(i)
                print(i)

                line_x = self.zero_to_nan(line_x) | units.kpc
                line_y = self.zero_to_nan(line_y) | units.kpc
                line_z = self.zero_to_nan(line_z) | units.kpc

                if options.get('fix') == 'z':
                    plot(line_x, line_y, lw=0.5)
                if options.get('fix') == 'y':
                    plot(line_x, line_z, lw=0.5)
                if options.get('fix') == 'x':
                    plot(line_y, line_z, lw=0.5)
                    
        if use_all=='n':
            c = options.get('check')
            cl = c.read()
            cl = cl.strip()
            cl = ast.literal_eval(cl)
            for i in cl:
                print(i)
                line_x, line_y, line_z = self.make_line(i)
                
                line_x = self.zero_to_nan(line_x) | units.kpc
                line_y = self.zero_to_nan(line_y) | units.kpc
                line_z = self.zero_to_nan(line_z) | units.kpc
                
                if options.get('fix') == 'z':
                    plot(line_x, line_y, lw=0.5)
                if options.get('fix') == 'y':
                    plot(line_x, line_z, lw=0.5)
                if options.get('fix') == 'x':
                    plot(line_y, line_z, lw=0.5)
            
def plot_mw_iso():
    MW = MWpotentialBovy2015()
    omega = 600 | units.km/units.s/units.km
    effective_iso_potential_plot(MW, omega, center_of_rotation = [0, 0] | units.kpc,
                                  xlim = [-10, 10] | units.kpc, ylim = [-10, 10] | units.kpc,
                                  resolution = [1000, 1000], fraction_screen_filled = 0.8,
                                  quadratic_contour_levels = True, contour_kwargs = dict(linestyles = 'solid',
                                  linewidths=0.5, cmap=('Greys')), omega2 = None, center_of_rotation2 = [0, 0]|units.kpc, 
                                  fraction_screen_filled2 = 0.2, projection3D = False, number_of_contours = 15)

def plot_MC(**options):
    with open('gal_line.pickle', 'rb') as f:
        gal_line = pickle.load(f)
        
    if options.get('fix') == 'z':
        plot(gal_line[0], gal_line[1], lw=0.7, label='LMC', ls='--')
        plot(gal_line[3], gal_line[4], lw=0.7, label='SMC', ls='--')
        plot(gal_line[6], gal_line[7], lw=0.7, label='Globular Cluster', ls='-.')
    
    if options.get('fix') == 'y':
        plot(gal_line[0], gal_line[2], lw=0.7, label='LMC', ls='--')
        plot(gal_line[3], gal_line[5], lw=0.7, label='SMC', ls='--')
        plot(gal_line[6], gal_line[8], lw=0.7, label='Globular Cluster', ls='-.')
    
    if options.get('fix') == 'x':
        plot(gal_line[1], gal_line[2], lw=0.7, label='LMC', ls='--')
        plot(gal_line[4], gal_line[5], lw=0.7, label='SMC', ls='--')
        plot(gal_line[7], gal_line[8], lw=0.7, label='Globular Cluster', ls='-.')
