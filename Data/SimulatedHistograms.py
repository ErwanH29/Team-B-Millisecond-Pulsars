import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import re
from amuse.units import units

inputUser = int(input("Which simulation do you want to check? "))

file_in = open('PickleFiles/neut_star_1Gyr_frac=0.55_' + str(inputUser) + '/check.txt', 'r')
for y in file_in:
    trimed_line = y.strip()
indices = re.findall(r"[-+]?\d*\.\d+|\d+", trimed_line)
indices = [int(i) for i in indices]

data = pd.read_pickle("PickleFiles/neut_star_1Gyr_frac=0.55_" + str(inputUser) + "/neut_stars_positions.pkl")

distancemod = [] 
for i in range(len(indices)):
    distance = abs(np.sqrt(np.sum(data.iloc[indices[i]][-1][0]**2) ) )
    distancemod.append((distance))

plt.title("Histogram of Millisecond Pulsar Distances \n From Galactic Center as Observed by the ATNF")
n, bins, patches = plt.hist(distancemod, 120, color = 'black', histtype='step')
plt.xlim(np.min(distancemod),20)
plt.xlabel(r"Distance (kpc) ")
plt.savefig("HistogramDistance", dpi = 300)
plt.show()