## This branch is for the scripts the report and simulation are based off of.

**conv_coord.py**
- This script takes the positional and velocity coordinates of any object and rotates it into the galactocentric frame
- We still need to figure out how to manipulate the velocity coordinates into the galactocentric frame - issues with unit manipulations

**galpot_init.py**
- This script defines the SMC and LMC class potentials. Potentials based on different papers and following various models to make the gravitational field felt by pulsars
  all the more intricate.
- We still need to find the values of the Plummer profile and NFW profile for the SMC

**gravsolve.py**
- This script takes a pre-defined star and shoots it based on initial velocities and position to track it's motion when succumbed to the MW and LMC potential
- We still need to add the SMC potential as well as the relative motion of the SMC and LMC w.r.t the MW

**motion_of_MCs.py**
- This script fixes the Milky Way and treats LMC and SMC as particles with their own masses. Then track their trajectory under the influence of the MW potential.
- We still need to consider the interaction of LMC and SMC.

**test_in_dynamic_potential.py**
- This script describes the possible trajectory of a test star given initial conditions moving in the dynamic MW-LMC-SMC system potential. 
- We still need to simplify the code and reduce the occurrence of MPI_ERR_SPAWN.
