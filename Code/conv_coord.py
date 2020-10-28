# IMPORTANT: use heliocentric input units (degree, degree, kiloparsec)!!
# Returns: galactocentric x, y and z distance in kiloparsec
from astropy import units as u
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
import numpy as np
from amuse.lab import units

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

