from amuse.lab import units
import numpy as np
import pandas as pd
from amuse.plot import *
from matplotlib import pyplot
import re
import math
import ast

def capture_check(L_neut):
    
    check = []
    
    x_upp = 20 | units.kpc
    x_low = -20 | units.kpc
    y_upp = 20 | units.kpc
    y_low = -20 | units.kpc
    z_upp = 20 | units.kpc
    z_low = -20 | units.kpc
    
    last_col = L_neut.iloc[:,-1]
    for i in range(len(last_col)):
        
        pos_check = last_col.iloc[i][0]

        if (x_low <= pos_check[0] <= x_upp and \
            y_low <= pos_check[1] <= y_upp and \
                z_low <= pos_check[2] <= z_upp):
            
            check.append(i)
        
    return check

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
    
    def plot(self, use_all=True, **options):
        row_len = len(self.dat.index)
        if use_all==True:
            for i in range(row_len):
                line_x, line_y, line_z = self.make_line(i)
                print(i)

                line_x = self.zero_to_nan(line_x)
                line_y = self.zero_to_nan(line_y)
                line_z = self.zero_to_nan(line_z)
                print(line_x)

                if options.get('fix') == 'z':
                    plot(line_x, line_y, lw=0.5)
                if options.get('fix') == 'y':
                    plot(line_x, line_z, lw=0.5)
                if options.get('fix') == 'x':
                    plot(line_y, line_z, lw=0.5)
                    
        if use_all==False:
            for i in options.get('check'):
                line_x, line_y, line_z = self.make_line(i)
                
                line_x = self.zero_to_nan(line_x)
                line_y = self.zero_to_nan(line_y)
                line_z = self.zero_to_nan(line_z)
                
                if options.get('fix') == 'z':
                    plot(line_x, line_y, lw=0.5)
                if options.get('fix') == 'y':
                    plot(line_x, line_z, lw=0.5)
                if options.get('fix') == 'x':
                    plot(line_y, line_z, lw=0.5)
        
        pyplot.title('Ejected MSP from LMC')
        pyplot.xlabel(r"$x$ coordinate (kpc)")
        pyplot.ylabel(r"$y$ coordinate (kpc)")
        pyplot.savefig("neutron ejection", dpi=300)
        pyplot.show()
            
        
    
    
