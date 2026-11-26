# COVID-19 vaccine efficacy example
The COV002 trial data is available via request on the Vivli platform https://doi.org/10.25934/PR00010299.
We do not provide the data as part of this repository, therefore the code will produce an error.
We are yet to verify whether this dataset matches exactly the data we used, so the code may still not run even with this dataset.
Since intermediate R objects are too large to save to GitHub, and for privacy reasons, we do not provide these.

## Saving the data
The required datasets should be saved in the Data folder.

## Running the code
The code should be run in the following order:
0_Data_cleaning.R
1_Longitudinal_antibody_model.R
2_Cox_infection_model.R
3_COVID_example.R

## Computation time
0_Data_cleaning.R and 3_COVID_example.R take seconds to run, or a few minutes at most.
1_Longitudinal_antibody_model.R and 2_Cox_infection_model.R take much longer and require considerable computing power:
1_Longitudinal_antibody_model.R
The longitudinal model was run in Stan with 4 chains, in parallel on a SLURM compute cluster. This
took 2h 23m (9h 33m total walltime), and required 52GB of RAM to generate 15000 iterations per
chain, after 5000 iterations of warmup. 
2_Cox_infection_model.R
The Cox model was run in parallel using the snowfall package in R, on the imputed datasets 
from each iteration from the longitudinal model. This was run on a SLURM compute cluster 
using 24 parallel nodes, taking 2h 37m to run (2 days 12h, 6m total walltime). 
The Cox analysis required 41GB of RAM. 
