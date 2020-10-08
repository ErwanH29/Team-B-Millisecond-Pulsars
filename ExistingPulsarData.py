import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Desktop\PulsarsDataATNF(2019-04-24).txt", comments="#", usecols=(6,7, 8)) #Importing milisecond pulsar data
for i in range(0,len(data)):
    plt.scatter(data[i,1], data[i,2])
plt.show()