## This branch looks at manipulating both theoretical Millisecond Pulsars (MSP) as well as observational ones.

**AllPulsarsData.py:**

- This script extracts the position and velocities of the observed MSP data taken from the ATNF Pulsar Catalogue (last updated April 2019)
- It can be used for a comparison of our simulated MSP with real life objects - i.e comparison of distance distribution w.r.t MW, off plane pulsar population ...
- The script also defines a probability for a given velocity stemmed from a asymmetrical kick modelled after Hansen 1997. With this it provides the velocities for 1000
  different pulsars based on the weighted probability for a given velocity (where the minimum and maximum are based on the minimum and maximum velocity values observed).
- The script also initiates 1000 different mass values using a Gaussian distribution following the CLT with the standard deviation being based on Antoniadis 2016.
- The script still needs to initialise the initial positions for our simulated pulsars.

**ExistingPulsarData.py:**
- This script will be imported onto another script at a later time, but for now it plots all observed pulsars onto the galactic plane to help with a visual comparison of simulated
  pulsars with the observed ones.
  
**PulsarsDATAATNF(2019-04-24).txt:**
 - This script includes the data for all observed milisecond Pulsars following Lorimer 2008's definition where period T < 30ms.
 - It only includes the distance
 
 **VELOCITIESPulsarsDataATNF(2019-04-24).txt:**
 - This script includes the data for all observed MSP with a defined velocity (PMRA, PMDEC and VTRANS) and could possibly be used to look at the origin of some pulsars.
