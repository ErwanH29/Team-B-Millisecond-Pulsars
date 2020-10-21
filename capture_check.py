from amuse.lab import units
import numpy as np

def capture_check(L_neut):
    
    check = []
    
    x_upp = 20 | units.kpc
    x_low = -20 | units.kpc
    y_upp = 20 | units.kpc
    y_low = -20 | units.kpc
    z_upp = 10 | units.kpc
    z_low = -10 | units.kpc
    
    for j in range(len(L_neut)):
        posx_check = L_neut[j][0][-1]
        posy_check = L_neut[j][1][-1]
        posz_check = L_neut[j][2][-1]
        
        if (x_low <= posx_check <= x_upp and \
            y_low <= posy_check <= y_upp and \
                z_low <= posz_check <= z_upp):
            
            check.append(j)
            
    return check
