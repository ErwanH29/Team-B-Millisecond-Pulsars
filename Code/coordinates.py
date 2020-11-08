from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import numpy as np
from amuse.lab import units
import ast
import pandas as pd
import json

def conv_coord(right_ascension, declination, dist):

    coordinates = SkyCoord(ra=right_ascension*u.degree, 
                           dec=declination*u.degree, 
                           distance=dist*u.kpc)
    coordinates_co = coordinates.transform_to(coord.Galactocentric)
    x_co = np.array(coordinates_co.cartesian.x)
    y_co = np.array(coordinates_co.cartesian.y)
    z_co = np.array(coordinates_co.cartesian.z)
    
    d_co = (x_co, y_co, z_co) | units.kpc
    return d_co

def conv_coordPulsars(data): #I added a converting coordinate for Pulsars from heliocentric to galactocentric
    coordinates = SkyCoord(x=data[:,1], 
                           y=data[:,2], 
                           z=data[:,0], 
                           unit='kpc', 
                           representation_type='cartesian')
    coordinates_coMSP = coordinates.transform_to(coord.Galactocentric)
    x_coMSP = np.array(coordinates_coMSP.cartesian.x)
    y_coMSP = np.array(coordinates_coMSP.cartesian.y)
    z_coMSP = np.array(coordinates_coMSP.cartesian.z)
    MSPCoords = np.stack((x_coMSP, y_coMSP, z_coMSP), axis=1) | units.kpc

    return MSPCoords

def prop_motion(pm_right_ascension, pm_declination, rad_velocity):
    proper_motion = SkyCoord(pm_ra_cosdec=pm_right_ascension*u.mas/u.yr, 
                             pm_dec=pm_declination*u.mas/u.yr, 
                             radial_velocity=rad_velocity*u.km/u.s)
    proper_motion_co = proper_motion.transform_to(coord.Galactocentric)
    #vx_co = np.array(proper_motion_co.cartesian.x)
    #vy_co = np.array(proper_motion_co.cartesian.y)
    #vz_co = np.array(proper_motion_co.cartesian.z)
    #v_co = (vx_co, vy_co, vz_co) | units.km/units.s
    v_co = proper_motion_co
    return v_co

def old_coordlist():
    LMC_coord = (-32.2877571477, -201.334477302, -194.233794365) | units.kpc
    LMC_vel = (-23.6377151397, -116.106444922, -126.743869139) | units.kms
    SMC_coord = (-35.2794369809, -140.65635819, -96.9363051162) | units.kpc
    SMC_vel = (-49.5945188999, -49.9087554726, -12.7823855815) | units.kms
    ngc_coord = (-30.3480687709, -203.481152417, -197.606147399) | units.kpc
    ngc_vel = (26.7697951421, -130.939929951, -234.017971574) | units.kms
    
    data = (LMC_coord, LMC_vel, SMC_coord, SMC_vel, ngc_coord, ngc_vel)
    return data

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

def get_final_coord(check):
    
    data = pd.read_pickle('working_data/neut_stars_positions.pkl')

    mag_tot = []
    check = check.read()
    check = check.strip()
    check = ast.literal_eval(check)
    
    last_col = data.iloc[:,-1]
    
    for c in check:
        last_coord = last_col.iloc[c][0]
        mag = (last_coord[0]**2 + last_coord[1]**2 + last_coord[2]**2).sqrt()
        mag = mag * (1 | units.kpc**-1)
        mag_tot.append(mag)
    f = open('working_data/distances.txt', 'w')
    json.dump(mag_tot, f)
    f.close()
    print(' ------------------------------------------------------------\
------------------------------------------------------------ \n \
distances.txt has been generated, containing the galactocentric distances of the close neutron stars in kpc. \n \
------------------------------------------------------------\
------------------------------------------------------------')
