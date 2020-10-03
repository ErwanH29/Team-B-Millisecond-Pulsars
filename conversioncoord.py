from astropy import units as u
from astropy.coordinates import SkyCoord

c = SkyCoord(ra=79.88*u.degree, dec=-69.59*u.degree, distance=49.97*u.kpc)

xcord = c.cartesian.x
ycord = c.cartesian.y
zcord = c.cartesian.z

print("x = ", xcord)
print("y = ", ycord)
print("z = ", zcord)

dist = (xcord**2 + ycord**2 + zcord**2)**0.5

print(type(xcord))
print("Distance to LMC: ", dist)