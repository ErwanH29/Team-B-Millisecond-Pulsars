from particle_init import neut_initialiser
import numpy as np
import itertools

p = neut_initialiser.velocityList
pindiv = p[-1]
print("Extracting all the x, y, z velocities of the stars: \n", p)
print("Extracting the x, y, z velocities of the last star: \n", pindiv)

q = neut_initialiser.neutCords
qindiv = q[-1]
print("Extracting all the x, y, z coords of the stars: \n", q)
print("Extracting all the x, y, z coords of the stars: \n", qindiv)