# Simulations
The output folder contains the R objects created during the simulations study.
These are read and analysed in Simulation_results.R.
One RDS object is provided for each of the 6 settings: linear, logistic and Cox regression with a Gaussian or lognormal covariate.
These objects can be inspected, and provide finer detail about the results than those possible to report in the paper.

The files required to replicate the simulation study are also included.
These files took days to run on a SLURM server using 12 or more nodes, therefore considerable computing power will be required to recreate them.
