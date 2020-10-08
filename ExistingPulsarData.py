import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("Desktop\PulsarsDataATNF(2019-04-24).txt", comments="#") #Importing milisecond pulsar data
print(data[2,3])
for i in range(0,len(data)):
    plt.scatter(data[i,3], data[i,4])
plt.show()