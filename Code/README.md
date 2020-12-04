## This directory incorporates the code that the simulation is based off


**bias.py**:

* This script returns the minimum variance unbiased estimator of the distribution of MSP distances to the galactic centre and the value of minimum variance. The variance can describe the bias of the distribution.


**coordinates.py**: 

* This script converts the coordinates of the galactic systems (the LMC, SMC and Milky Way) from equatorial coordinate systems to cartesian coordinates which is done under the function conv_coord. This is done with the help of astropy which also allows a conversion of the frames from a helicoentric frame to a galactocentric one providing not only the right spatial coordinates, but velocity coordinates as well.
* The function conv_coordPulsars changes the reference frames of the millisecond pulsars from a heliocentric one into a galactocentric one. 
* The function get_final_coord calculates the final distances of the ejected pulsars to the center of the Milky Way which is used for data manipulation.

**energy_conservation.py**:

* This script describes the total energy (potential energy + kinetic energy) of the MCs-MW system in the evolution process.

**evol.py**:

* This script traces the evolution of a particular MSP over time as it moves across the galactic systems that are bridged with one another and gets shot out at some arbitrary velocity.
* The LMC and SMC center of mass trajectories are also calculated here to provide a more accurate representation of the system

**file_logistics.py**:

* This script handles files in order to use the full script. It takes care of saving simulation results and transfers them to appropriate folders.

**galpot_init.py**:

* This script defines the LMC and SMC potentials. For each timestep defined in evol.py, the potentials of both galaxies are re-calculated at the new positions.

**interface.py**:

* Upon running this script, a simulation can be started.


**obtain_old_coord.py**:

* This script evolves the LMC and SMC back in time, inside their mutual interaction and the Milky Way background potential, to obtain the coordinates 1 Gyr back in time.

**particle_init.py**:

* This script initialises the N amount of milisecond pulsars we aim to shoot. It has the initial mass, velocity and coordinates of the stars based on various - --models.
* It also traces the path of the LMC and SMC center of mass such that we can track them as their positions evolve over time. This is done by representing the
center of masses as particles

**plotters.py**:

* This script contains all the functions that plot something with the simulation data.
