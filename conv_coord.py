from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import numpy as np
from amuse.lab import units

# IMPORTANT: use heliocentric input units (degree, degree, kiloparsec)!!
# Returns: galactocentric x, y and z distance in kiloparsec

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
