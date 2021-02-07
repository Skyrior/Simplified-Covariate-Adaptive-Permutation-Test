# Simplified-Covariate-Adaptive-Permutation-Test
This code is part of an assessment for an RA application.
It is a simplified version of the Covariate-Adaptive Permutation Test as outlined in (Bugni, Canay, Shaikh 2017).

## -- Provides a parallel-processed Monte-Carlo simulation for applying a
##    two-sample permutation test with arbitrary test statistics on an
##    arbitrarily generated data set with no stratification, and calculates
##    the rejection probability.
##
##    The provided code only implements the two test statistics described
##    in the assessment and uses rnorm(), but with slight modification
##    can be extended to any arbitrary distribution and arbitrary test
##    statistics by supplying the functions as arguments to the Monte Carlo
##    main function. The parallel processing allows the simulation to scale
##    with logical cores very well by Amdahl's Law, as most of the computations
## -- are parallelized.

## -----------------
## QUICK START
## -----------------

## Simply select all the code and run. The main function is the mc (Monte Carlo)
## function that will create a dataset of randomly generated observations
## of input distributions, and use the given test statistics to estimate the
## probability of rejection based on the two-sample permutation test.
