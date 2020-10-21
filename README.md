This GitHub incorporates the code utilising AMUSE which analyses the possibility of millisecond pulsars shot out from the Large Magellanic Cloud being captured inside the Milky Way.

The Data directory shows the main references used in the project, observations of existing pulsars as well as small test scripts used to manipulate such data and the other repository is the code the research is based off of.

The project proposal gives a good introduction on the aims of our research and how we intend to achieve it and lies with the final report in the directory Reports.

![alt text](https://imgur.com/a/caK97wQ)

## This folder is for the scripts the report and simulation are based off of.

**conv_coord.py**
- This script takes the positional (deg,deg,kpc) and velocity coordinates of any object and rotates it into the galactocentric frame into an (x,y,z) coordinate in -- kpc

**evol.py**
- This script traces the evolution of a particular MSP over time as it moves across the galactic systems that are bridged with one another and gets shot out at some arbitrary velocity.
- The LMC and SMC center of mass trajectories are also calculated here to provide a more accurate representation of the system

**galpot_init.py**
- This script defines the LMC and SMC potentials. For each timestep defined in evol.py, the potentials of both galaxies are re-calculated at the new positions.

**gravsolve.py**
- This script runs the codes and plots the output. Currently, the output comprises the LMC and SMC center of mass trajectories, the MSP trajectory and the Milky Way isopotential. 

- The isopotential is only plottable in (x,y). The rest in any combination of two coordinates. The LMC and SMC center of mass particles interact and both feel the Milky Way background potential. The motion of the MSP is tracked using the reference frame of the Milky Way and depends on the updating LMC and SMC trajectories

**particle_init.py**
- This script initialises the N amount of milisecond pulsars we aim to shoot. It has the initial mass, velocity and coordinates of the stars based on various - --models.
- It also traces the path of the LMC and SMC center of mass such that we can track them as their positions evolve over time. This is done by representing the
- center of masses as particles
