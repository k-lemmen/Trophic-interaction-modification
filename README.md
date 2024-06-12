# Data and code for "Food web context modifies predator foraging and weakens trophic interaction strength"

This repository contains code and data needed to reproduce: https://doi.org/10.1101/2024.03.04.583297

## Instructions
All analyses were done in R. 

All data needed to reproduce the analyses are available in the '01 - Data' folder.

Functions used for fitting the ODEs are in the '02 - Functions' folder.

R scripts for fitting the ODEs are in the '03 - ODE Fitting' folder.

R scripts for community dynamics are in the '04 - Community Dynamics Simulation' folder. 

Outputs can be found in the “05 – Output” Folder, this includes (i) Search grids of initial starting points and the corresponding final parameter estimates used to find a global minimum (Appendix S1d) and (ii) fits of the 13 different functional response models within the ODE system. Rdata files are provided for convenience, but all outputs can be independently derived from the provided code and CSV files.

R scripts and outputs for analysis found in the supplementary materials are in the '06 - Appendix Analysis" folder. 

Metadata can be found in the “10 – Metadata” folder and provides a description and units for all variables.
