The code contained in this repository simulates a system consisting of the Large and Small Magellanic Clouds, a globular cluster and the Milky Way.
The positions of the Large and Small Magellanic cloud are first calculated 1 Gyr back in time using current position and velocity coordinates.
Upon running, pulsars will be generated along the trajectory of the globular cluster, who is simultaneously affected by the Large Magellanic Cloud potential, which in itself is bridged with the Milky Way potential giving for a more realistic model of the galactic neighborhood.
After running, data can be generated from the simulation.

A simulation can be started by running interface.py from your favourite python IDE. the script will ask what you want to do (ranging from running a new simulation to plotting older data). For proper execution, the number_of_workers variable in interface.py should not exceed the amount of CPU cores present on your system.

![](https://i.imgur.com/LuFyQi8.png)

For visualization of the code, these three GIF files illustrate the dynamics of the Large and Small Magellanic Cloud potentials in the Milky Way potential over a time of 270 million years, starting at present time, as calculated by our code. The imaged planes are (x,y), (x,z) and (y,z) respectively.

![](https://github.com/ErwanH29/Team-B-Millisecond-Pulsars/blob/master/xy_0_300%25_3f.gif) 
![](https://github.com/ErwanH29/Team-B-Millisecond-Pulsars/blob/master/xz_0_300%25_3f.gif)
![](https://github.com/ErwanH29/Team-B-Millisecond-Pulsars/blob/master/yz_0_300%25_3f.gif)
