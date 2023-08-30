# Mexican Wolf Bayesian hierarchical population model

This document provides guidance for the implementation of the Bayesian hierarchical population model associated with "Factors affecting recovery of Mexican wolves in southwestern USA" by Stewart W. Breck, Amy Davis, John K. Oakleaf, David L. Bergman, Jim deVos, J. Paul Greer, and Kim Pepin 2023 with Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.14483. Two R scrips are provided and one csv datafile.

## Files included:

1. Wolf_to_run.R
2. wolf_pop_MCMC.R
3. MexWolfDataFile.csv

## Data needs

The data needed to fit this model are annual wolf population data including:

1.  year = Year of data
2.  pop = Minimum population size by year (count of wolf population)
3.  Reproduction = Number of reproduction events by year
4.  PupRecruitment = Minimum pup recruitment number
5.  Removals = Number of individuals removed by management by year
6.  Releases = Number of management released individuals
7.  Translocations = Number of management translocations by year
8.  mort = Number of deaths from known or legal causes
9.  poaching = Number of deaths from illegal causes


## Steps to run the code
The Wolf_to_run.R code runs with the data file included in this respository. To use your own data, you may need to make some changes (including the ones listed above in the data needs section).  The following steps take you through how to run this code with the data provideds and how modify Wolf_to_run.R to work with your data. This code was created under R version 4.1.2. 

1.	Save both R scrips and the data file in the same working directory (or modify the location of the wolf_pop_MCMC.R on Line 12).
2.	Ensure you have all of the packages loaded on your machine (from both scripts).
3.  If using your own data, Change the name of your data set on Line 15.
4.  Set the design matrices for mortality, reproduction, and proportion mortality on Lines 28-30. These can be modified to compare other covariate structures than the one shown here.
5.  Change the number of MCMC iterations you would like to use (Line 22).
6.  The tuning parameters may need to be changed to fit your data better (Lines 19-21).
7.  If the different mortality types have a different probability of going unobserved then change the 'mortmis' vector on Line 24.
8.  The code to run the MCMC is on Lines 34-36.

## Model output

The model will output the full posterior data in a list for the following parameters:

1. npredsave = the abundance estimates per year
2. Mtsave = the annual mortality rate estimates
3. Rtsave = the annual reproduction rate estimates
4. betarsave = the covariate estimates for reproduction
5. pmsave = the annual posterior proportion of mortality from illegal causes
6. poacht = the annual estimated number of moralities from illegal causes
7. poachratet = the annual estimated illegal mortality rate
8. betamsave = the covariate estimates for mortality
9. betapsave = the covariate estimates for proportion of mortality from different causes


The model produces diagnostic plots for: betas for reproduction rates, betas for mortality rates, betas for mortality proportions, and abundance.  These include trace plots (using all of the data) and posterior distribution plots assuming a burn-in of 50%. 

Lines 39-66 calculate the posterior means and 95% credible intervals of the output data and compile the results into a data frame called "combdat". 
 
 




