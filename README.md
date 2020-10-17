## This branch is for the scripts the report and simulation are based off of.

**conv_coord.py**
- This script takes the positional and velocity coordinates of any object and rotates it into the galactocentric frame

**evol.py**
- This script traces the evolution of a particular MSP over time as it moves across the galactic systems and gets shot out at some arbitrary velocity.

**galpot_init.py**
- This script defines the SMC and LMC class potentials. Potentials based on different papers and following various models to make the gravitational field felt by pulsars all the more intricate.

**gravsolve.py**
- This script takes a pre-defined star and shoots it based on initial velocities and position to track it's motion when succumbed to the MW and LMC potential
- We still need to add the SMC potential as well as the relative motion of the SMC and LMC w.r.t the MW

**particle_init.py**
- This script initialises the N amount of milisecond pulsars we aim to shoot. It has the initial mass, velocity and coordinates of the stars based on various models.
- It also traces the path of the LMC and SMC such that we can track them as their positions evolve over time.

**test_in_dynamic_potential.py**
- This script describes the possible trajectory of a test star given initial conditions moving in the dynamic MW-LMC-SMC system potential.
