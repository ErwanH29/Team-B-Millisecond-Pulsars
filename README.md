## This branch is for the scripts the report and simulation are based off of.

**conv_coord.py**
- This script takes the positional (deg,deg,kpc) and velocity coordinates of any object and rotates it into the galactocentric frame into an (x,y,z) coordinate in -- kpc

**evol.py**
- This script traces the evolution of a particular MSP over time as it moves across the galactic systems and gets shot out at some arbitrary velocity. Also the 
- LMC and SMC center of mass trajectories are calculated here

**galpot_init.py**
- This script defines the LMC and SMC potentials. For each timestep defined in evol.py, the potentials of both galaxies are re-calculated at the new positions.

**gravsolve.py**
- This script runs the codes and plots the output. Currently, the output comprises the LMC and SMC center of mass trajectories, the MSP trajectory and the milky
- way isopotential. The isopotential is only plottable in (x,y). The rest in any combination of two coordinates. The LMC and SMC center of mass particles interact
- and both feel the milky way background potential. The MSP feels the static milky way and the updating LMC and SMC

**particle_init.py**
- This script initialises the N amount of milisecond pulsars we aim to shoot. It has the initial mass, velocity and coordinates of the stars based on various - --models.
- It also traces the path of the LMC and SMC center of mass such that we can track them as their positions evolve over time. This is done by representing the
- center of masses as particles
