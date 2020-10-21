# The Possibility of Exotic Millisecond Pulsars Living in the Milky Way

Team members: Shiyang Zhang, Arend Moerman, Erwan Hochart

### Summary

This report aims to take a statistical approach in analysing the possibility whether some millisecond pulsars living in the Milky Way may have originated from the Large Magellanic Cloud.

This project begins by assuming that the dynamical evolution of the Large Magellanic Cloud (LMC), Small Magellanic Cloud (SMC) and Milky Way to be negligible as millisecond pulsars get continuously shot out existing globular clusters existing within the LMC in a time span of 1Gyr.
Many more assumptions are going to be made during the report and this will be highlighted under the assumptions section of the report, including the velocity distribution and velocity mechanism of the millisecond pulsars.

### Description of the proposed work

**The Research Question:** 

Given the theoretical population distribution of millisecond pulsars living in the Milky Way, what is the fraction that some may have originated from the Large Magellanic Cloud?

For setting up the system, one main formula is used to give a velocity distribution of the milisecond pulsars. This distribution follows the works by Hansen (1997) who described the velocity of milisecond pulsars which took into account observational issues giving bias towards the lower velocities. 
Due to the velocity needing an excess of ~200 km/s to be unbound from the LMC. The distribution was taken with a lower limit of 200 km/s and an upper limit based on the faster milisecond pulsar observed to date (based on the ATNF database, last updated 2019). This heavily filters out the milisecond pulsar population as  most have a velocity of around 80 km/s (Lorimer 2008). 

Furthermore a large assumption was taken in the potentials of the Milky Way, LMC and SMC respectively as these systems will remain rigid, even as they move under the influence of one another. Nevertheless, their potential is used by the help of the AMUSE package galactics_potential which stems it's code based on galpy (Bovy 2015) for the Milky Way and a simple Plummer and NFW potential for the LMC and SMC.
By bridging the SMC and LMC to the Milky Way as well as the simulated milisecond pulsars to the globular cluster it originates from which is bridged with the LMC, we allow ourselves to obtain a nice representation of the system as they all evolve under each others influence.

To convert the information of position and velocity of all three systems (HVS, Milky Way and LMC) the python package astropy will be used. We will convert the system into the Milky Way frame of reference using this package as well. 
The initial conditions for HV3 are taken from the Gaia Dr2 data set and are;

**Why Bother?**

After looking extensively through NASA ADS' catalogue, there doesn't seem to be research that follows a familiar goal and so this can give new light to the pulsar population of the Milky Way. A celestial object integral to physics with it's many uses. In fact it seems to be poorly understood why a significant fraction of Milky Way millisecond pulsars are not in a binary system and such a report may help shed light towards this (Lorimer 2001).
Furthermore pulsars are an integral part of physics as there are many uses to it, ranging from utilising their periodic pulses they to detect the first extra-solar system planet detected, to testing general relativity with their strong gravitational fields.

**Literature references:**

Brad M. S. Hansen and E. Sterl Phinney. The pulsar kick velocity distribution. , 291(3):569–577, November 1997. doi: 10.1093/mnras/291.3.569.

Duncan R. Lorimer. Binary and Millisecond Pulsars at the New Millennium.Living Reviewsin Relativity, 4(1):5, June 2001. doi: 10.12942/lrr-2001-5.

R. N. Manchester, G. B. Hobbs, A. Teoh, and M. Hobbs. The Australia Telescope NationalFacility Pulsar Catalogue. , 129(4):1993–2006, April 2005. doi: 10.1086/428488.

Simon Portegies Zwart and Steve McMillan.Astrophysical Recipes;  The art of AMUSE.2018. doi: 10.1088/978-0-7503-1320-9.

