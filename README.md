# Code for simulating TB recurrence and screening interventions, in a cohort of recently treated patients in India. 

"aftermath analyses for manuscript.R" sets and samples parameter value ranges, runs cohort and intervention simulations, and generates manuscript results. 
Sources the aftermath cohort sim functions.R file. 

"aftermath cohort sim functions.R" contains the functions to generate a cohort with a given set of parameter values, 
check consistency of the cohort with asymptomatic-TB-prevalence targets, run simulated interventions, and calculate impact and cost outcomes. 

"aftermath disease course data analysis.R" was used to estimate certain model parameters from primary data; can be ignored. 

"gamma fitting.py" (old) and "weibull fitting.py" fit a distribution of time to TB onset, 
to summary statistics of the observed truncated distribution within the first 18 months after treatment completion. 

"Aftermath for 2.0 simulations.R" is a stripped-down version of the overall model with addition of death vs routine notification endpoints; 
can be used for power simulation for future trials. 

