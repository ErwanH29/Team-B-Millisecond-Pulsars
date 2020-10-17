from particle_init import neut_initialiser
import numpy as np
import itertools

p = neut_initialiser.neut_InitialConditions
print(p, np.shape(p))

print("Extracting MSP 1:\n", neut_initialiser.neutlist_0)
print("Extracting MSP 4:\n", neut_initialiser.neutlist_3)
print("Extracting MSP 10:\n", neut_initialiser.neutlist_9)