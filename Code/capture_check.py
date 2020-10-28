from amuse.lab import units
import numpy as np
import pandas as pd

def capture_check(L_neut):
    
    check = []
    
    x_upp = 20 | units.kpc
    x_low = -20 | units.kpc
    y_upp = 20 | units.kpc
    y_low = -20 | units.kpc
    z_upp = 10 | units.kpc
    z_low = -10 | units.kpc
    
    last_col = L_neut.iloc[:,-1]
    for i in range(len(last_col)):
        
        pos_check = last_col.iloc[i][0]

        if (x_low <= pos_check[0] <= x_upp and \
            y_low <= pos_check[1] <= y_upp and \
                z_low <= pos_check[2] <= z_upp):
            
            check.append(i)     
    return check
