## Charles Shi
## Question 2 for the Research Professional Assessment for Shaikh
## & Tabord-Meehan.
## github repository:
## https://github.com/Skyrior/Simplified-Covariate-Adaptive-Permutation-Test
## R version: 4.0.3

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

## -- WARNING: The R interpreter is not thread-safe. Care has been taken to
##    ensure thread safety but if issues arise, a quick replacement of
##    %%dopar%% to %%do%% should do the trick.

## -----------------
## QUICK START
## -----------------

## Simply select all the code and run. The main function is the mc (Monte Carlo)
## function that will create a dataset of randomly generated observations
## of input distributions, and use the given test statistics to estimate the
## probability of rejection based on the two-sample permutation test.

library(purrr) ## for lfold (reduce) use install.packages("purrr") if necessary
library(parallel)
library(doParallel)
library(foreach)
plan(multicore)
cores <- detectCores(logical=TRUE)
registerDoParallel(ifelse(cores>=4, cores-1, 1)) ## keep 1 of the cores idle.
## The program is significantly faster when utilizing basic multithreading.
## Change %%dopar%% to %%do%% if parallel processing causes issues.

## -- seed --
## Unfortunately, the parallel process does not respect the RNG seed.
## Uncomment if using the single-threaded loops for replication.
## set.seed(NULL)

## -- debug --
## Set to true to turn on most messages in functions to trace
## what is happening to the data. Warning: only turn on debug
## and run small tests. Can crash if used on the Monte Carlo function
## with >5 replications directly. Using it on a single replication
## to trace inputs and errors is ideal.

debug = FALSE

## -- gm, gn & greps --
## Default values for generated datasets.

gm=10 ## the default, global value for m.
gn=10 ## ... for n.
greps=500 ## The number of random permutations to generate (B in step 3).

### -- deprecated
### Set suppress = TRUE if the loop messages is not needed.
### suppress <- FALSE ### does not work for parallel processes.

## -- MAIN --

## The Monte-Carlo simulation
## -- mc (int m, int n, function xfunc,
##        function yfunc, function test,
##        List xdist, List ydist, int reps) --
## Performs a Monte Carlo simulation of test statistics on
## randomly generated data according to a given distribution and
## parameters.
## 
## -- Arguments --
## test: A vector of test statistics to use. For example, to use
##       functions t1 and t2, test=c(t1, t2)
##
## -- Optional Arguments
## m: Length of the vector X = (X_i)
## n: Length of the vector Y = (Y_i)
## xfunc: Distribution for X. Specifically, the Monte Carlo simulation
##        will call xfunc(m, xdist). If m=20, xdist=list(0, 1), then that means
##        we will call xfunc(20, 0, 1).
## yfunc: Distribution for Y. Calls yfunc(n, ydist).
## alpha: rejection threshold
## xdist: Parameters passed to xfunc's second argument and on.
## ydist: Parameters passed to yfunc's second argument and on.
## reps: Passed to the first argument of xfunc and yfunc,
##       For most R standard functions this would be number of reps.
## seed: Run the simulation with the given seed for replication
##       purposes.
##
## -- Returns --
## r - Probability of rejection.
mc <- function(m=gm, n=gn, xfunc=rnorm, yfunc=rnorm, test, alpha=0.05,
               xdist=list(0, 0), ydist=list(0, 0), reps=500){
  
  print(system.time({
  
  ## input checks
  if((!is.function(xfunc)) || (!is.function(yfunc))) stop(
    "Supplied xfunc or yfunc is not a function."
  )
  if((!is.list(xdist)) || (!is.list(ydist))) stop(
    "Supplied xdist or ydist is not a list of parameters. Note:",
    " numerical arrays are not lists."
  )
  
  ## Create the matrix (reps by m+n) of observations of the random variables
  ## X and Y.
  
  X <- foreach(i=1:reps, .combine=rbind, .packages="purrr") %dopar% {
    do.call(xfunc, c(m, xdist))
  }
  
  Y <- foreach(i=1:reps, .combine=rbind, .packages="purrr") %dopar% {
    do.call(yfunc, c(n, ydist))
  }
  
  Z <<- as.data.frame(cbind(X, Y))
  
  message("Randomly generated data created!")
  
  ## Calculate the permutation test statistic for each observation:
  ## Start with creating an empty list
  
  L <<- list(1)
  
  ## Might cause issue passing functions and arguments when running
  ## in parallel. Use the commented block if there are issues.
  
  for(j in 1:length(test)){
    L[[j]] <<- foreach(i=1:reps, .combine=cbind, .packages="purrr") %dopar% {
      ## .inorder is TRUE by default, so we don't have to worry about
      ## out of order binding with parallel foreach.
      apply.tpermute(unlist(unname(Z[i,])), test[[j]], m=m, n=n)
    }
  }
  
  ## Uncomment if the for loop is causing problems, and comment the block above.
  ## Or simply change %dopar% to %do%!
  ##for(j in 1:length(test)){
  ##  L[[j]] <- foreach(i=1:reps, .combine=cbind, .packages="purrr") %do% {
  ##    ## .inorder is TRUE by default, so we don't have to worry about
  ##    ## out of order binding with parallel foreach.
  ##    apply.tpermute(unlist(unname(Z[i,])), test[[j]], m=m, n=n)
  ##  }
  ##}
  
  message("Test statistics values calculated for data!")
  
  ## rejections
  
  RJ <<- list(1)
  
  for(j in 1:length(test)){
    RJ[[j]] <<- foreach(i=1:reps, .combine=cbind, .packages="purrr") %dopar% {
      unname(reject(apply.teststatistic(Z[i,], test[[j]], m=m, n=n),
                    ## the test stat for the original
                    testvalues = L[[j]][,i],
                    threshold = alpha, method = quantile))
    }  
  }
  
  ## Uncomment if the for loop is causing problems. Or simply change %dopar% to
  ## %do%!
  ##for(j in 1:length(test)){
  ##  RJ[[j]] <- foreach(i=1:reps, .combine=cbind, .packages="purrr") %do% {
  ##    unname(reject(apply.teststatistic(Z[i,], test[[j]], m=m, n=n),
  ##                  ## the test stat for the original
  ##                  testvalues = L[[j]][,i],
  ##                  threshold = alpha, method = quantile))
  ##  }  
  ##}
  
  ## finally, return the probability of rejections for each test statistic used
  
  message("Rejections processed!")
  
  }))
  
  return(
    lapply(RJ, mean)
  )
  
}

## -- TEST STATISTICS IMPLEMENTATION --
## 
## 

## -- t1 (List data, int m, int n) --
## T_{m, n} = | sqrt(N) (X_m - Y_n) |, where X_m, Y_n are the sample means.
## 
## -- Arguments --
## data: Z = (X_1, X_2, ..., X_m, Y_1, Y_2, ..., Y_n)
##
## -- Optional Arguments
## m: Demarcates where to cutoff the data to get (X_1, ..., X_m). Default
##    is the gm variable.
## n: Sanity check to ensure m+n adds up. Default is the gn variable.
##
## -- Returns --
## r - The computed test statistic.

t1 <- function(data, m=gm, n=gn){
  
  ## since we will be accessing specific indices, let's transofrm
  ## the data into a vector to speed up.
  dat <- unlist(unname(data))
  
  if(debug) message("Calling test statistic 1 with ", dat, ", ", m,
                    " and ", n, " with data length ", length(dat))
  
  ## input checks
  if(!reduce(dat,
             function(a, b){a && is.numeric(b)},
             .init=TRUE)) stop("The data supplied to test statistic 1",
                               " is not numeric.")
  if(!is.numeric(m) || !is.numeric(n)) stop("m or n supplied ",
                                            "to test statistic 1 is not",
                                            " numeric.")
  if(!((m+n)==length(dat))) message("Warning: m+n =/= N in test",
                                     " statistic 1.")
  
  b <- splitAt(dat, m+1)
  x <- b[1]
  y <- b[2]
  
  x.mean <- mean(unlist(x))
  y.mean <- mean(unlist(y))
  
  r <- abs((sqrt(length(dat)))*(x.mean-y.mean))
  
  if(debug) message("Computed t1 statistic: ", r)
  
  return(r)
  
}

## -- t2 (List data, int m, int n) --
## T_{m, n} = t1_{m, n} / sqrt(N/m V(X_1, ..., X_m) + N/m V(Y_1, ..., Y_n)),
## where V is the variance operator.
## 
## -- Arguments --
## data: Z = (X_1, X_2, ..., X_m, Y_1, Y_2, ..., Y_n)
##
## -- Optional Arguments
## m: Demarcates where to cutoff the data to get (X_1, ..., X_m). Default
##    is the gm variable.
## n: Sanity check to ensure m+n adds up. Default is the gn variable.
##
## -- Returns --
## r - The computed test statistic.

t2 <- function(data, m=gm, n=gn){
  
  ## since we will be accessing specific indices, let's transofrm
  ## the data into a vector.
  dat <- unlist(unname(data))
  
  if(debug) message("Calling test statistic 2 with ", dat, ", ", m,
                    " and ", n, " with data length ", length(dat))
  
  ## input checks
  if(!reduce(dat,
             function(a, b){a && is.numeric(b)},
             .init=TRUE)) stop("The data supplied to test statistic 2",
                               " is not numeric.")
  if(!is.numeric(m) || !is.numeric(n)) stop("m or n supplied ",
                                            "to test statistic 2 is not",
                                            " numeric.")
  if(!((m+n)==length(dat))) message("Warning: m+n =/= N in test",
                                     " statistic 2.")
  
  b <- splitAt(dat, m+1)
  x <- b[1]
  y <- b[2]
  N <- length(dat)
  
  t1 <- t1(dat, m, n)
  x.var <- var(unlist(x))
  y.var <- var(unlist(y))
  
  r <- t1/(sqrt((N/m)*x.var + (N/n)*y.var))
  
  if(debug) message("Computed t2 statistic: ", r)
  
  return(r)
  
}


## -- AUXILLARY FUNCTIONS BELOW --
## 
## 

## -- get.permute (int N,) --
## A permutation is a mapping from {1, ..., N} to {1, ..., N}.
## The get.permute function constructs a random permutation.
## 
## -- Arguments --
## N: integer in [1, \infty).
##
## -- Returns --
## vector[] - a list of the second components of the permutation, in order.
## For example if a permutation is {(1, 3), (2, 1), (3, 2)}
## then this function outputs {3, 1, 2}.

get.permute <- function(N){
  
  if(debug) message("Calling get.permute() with arguments ", N)
  
  ## input checks
  tryCatch({tempvar <<- ((N%%1 == 0) && (N>=1))},
           error = function(a){
             stop("Supplied 'N' in get.permute() is not a positive integer.")
           })
  if(!tempvar) stop("Supplied 'N' in get.permute() is not a positive integer.")
  if(debug) message("Input in get.permute validated.")
  n <- as.integer(N)
  
  ## create a vector from 1 to N
  li <- 1:N
  
  ## produces a random, unbiased permutation with Fisher-Yates shuffle
  r <- sample(li)
  
  if(debug) message("Permutation generated in get.permute: ", r)
  
  ## output the permutation
  return(r)
}

## -- apply.permute (Vector list, Vector permutation) --
## Applies a permutation on the given list (vector), effectively shuffling it.
## 
## -- Arguments --
## list: Any generic list or vector.
##
## -- Optional Arguments --
## permutation: If supplied, applies the given permutation on the
## list. If not supplied, the function generates a random permutation.
##
## -- Returns --
## Vector[] - The shuffled list.

apply.permute <- function(list, permutation=NULL){
  
  ## since we will be accessing specific indices, let's transform
  ## the data into a vector.
  li <- unlist(unname(list))
  
  if(debug) message("Calling apply.permute() with arguments ",
                    li, ", ", permutation)
  
  ## input checks
  if(!(is.vector(li))) stop("Supplied a non-vector in apply.permute()")
  
  perm <- if(!is.null(permutation)){
    permutation
  } else {
    get.permute(length(li))
  }
  tryCatch({if(!isTRUE(all(perm == floor(perm)))){
    stop("Supplied permutation in apply.permute() is not a list of",
         "integers.")}},
    error = function(a){
      stop("Supplied permutation in apply.permute() is not a",
           "list of numbers.")
    })
  if(debug) message("Permutation used in apply.permute: ", perm)
  
  ## apply the permutation on the given list.
  
  s <- rep(NA, length(li))
  
  for(var in 1:length(li)){
    s[var] <- li[perm[var]]
  }
  
  return(s)
  
}

## -- apply.teststatistic (List list, function teststatistic) --
## Applies the given test statistic on the given list. The function
## will throw an error if there is a type mismatch: for example, if the
## test statistic is a function that processes strings but the given list
## is a list of booleans. In essence this wrapper is an error handler.
## 
## -- Arguments --
## list: Any generic list.
## teststatistic: Any appropriate function as long as teststatistic takes
## one list input.
## m: Demarcates where to cutoff the data to get (X_1, ..., X_m). Default
##    is the gm variable.
## n: Sanity check to ensure m+n adds up. Default is the gn variable.
##
## -- Returns --
## x - The output of the test statistic as applied to the given list. This
##     can be anything, depending on the supplied test statistic. If the
##     output is NOT a number, also generates a warning (but does not abort).

apply.teststatistic <- function(list, teststatistic,
                                m=gm, n=gn){
  
  d <- unlist(list)
  if(!is.function(teststatistic)) stop("Supplied test statistic",
                                       "is not a function!")
  r <- tryCatch({teststatistic(data=d, m, n)},
                error = function(a){
                  stop("Failed to apply the test statistic on the data!",
                       "\n", "Supplied list: ", d, "\n",
                       "Error: ", a)
                })
  
  if(!is.numeric(r)) message("Warning: The test statistic generated, ", r,
                             " is not a number.")
  
  if(length(r)>1) message("Warning: The test statistic generated, ", r,
                          " is a list!")
  return(r[1]) ##does not make sense to have multiple test statistic, will warn.
}

## -- apply.tpermute (List data, function teststatistic, int rep) --
## Applies the given test statistic on rep numbers of random permutations of
## the data.
## 
## -- Arguments --
## data: Any generic list.
## teststatistic: Any appropriate function as long as teststatistic takes
## one vector input.
##
## -- Optional Arguments --
## rep: Number of replications. (Default: 500)
## according to the given seed.
## m: Demarcates where to cutoff the data to get (X_1, ..., X_m). Default
##    is the gm variable.
## n: Sanity check to ensure m+n adds up. Default is the gn variable.
##
## -- Returns --
## r - The list of outputs of the test statistic as applied to the data
##     rep number of times.

apply.tpermute <- function(data, teststatistic, rep=500,
                           m=gm, n=gn){
  if(debug) message("Calling apply.tpermute() with arguments ",
                    data, ", ", "with test statistic ",
                    deparse(substitute(teststatistic)), ", ",
                    rep)
  
  permutations <- replicate(500, vector(length=40), simplify=FALSE)
  for(var in 1:rep){
    permutations[[var]] <- apply.permute(data) ##the data will be flattened
                                              ## into a vector. No worries.
  }
  
  r <- rep(NA, rep)
  for(var in 1:rep){
    r[var] <- apply.teststatistic(permutations[var],
                                  teststatistic, m=m, n=n)
  }
  
  return(r)
}

## -- reject (int testorig, List testvalues, int threshold, func method) --
## Determine whether to reject the original data (with test statistics value
## testorig) based on the quantile of test statistics applied on
## random permutations on an observation. Specifically, we reject if
## T_{m, n}(Z_1, ..., Z_N) > c_{1-alpha}, where T_{m, n} is the test statistic,
## c_{1-alpha} is the (1-alpha)-quantile of the test statistic matrix.
## To generate the test statistics values, use apply.tpermute
## 
## -- Arguments --
## testorig: The test statistic for the original observation.
## testvalues: Test values generated by apply.tpermute, or self-supplied.
## threshold: the alpha in (1-alpha)
## method: The method applied on testvalues using threshold.
##          Make sure to curry the method if necessary!
##          This function will call method(testvalues, 1-threshold) to get
##          c_{1-alpha}.
##
## -- Returns --
## r - 1 if rejected, 0 if not.

reject <- function(testorig, testvalues, threshold, method){
  if(debug) message("Calling reject() with arguments ",
                    data, ", ", "with test statistic ",
                    deparse(substitute(teststatistic)), ", ",
                    rep, " ", alpha)
  
  if(!is.numeric(testorig)) stop("Supplied non-numerical value testorig.",
                                 " in: reject()")
  if(!is.numeric(threshold) && !(0<=threshold && threshold<=1))
    stop("threshold supplied in reject() is not valid!")
  
  ## Apply the method on relevant parameters.
  c <- do.call(method, list(testvalues, 1-threshold))
  
  if(debug) message("Computed (1-", alpha, ")-method: ", c)
  
  if(debug) message("Original Test Statistic: ", testorig)
  
  ## decide whether to reject
  rej <- ifelse(testorig>c, 1, 0)
  return(rej)
}

## -- splitAt(vector, pos) --
## Splits vector at position.
## Ex: (1, 3, 5, 7, 9) with pos=2 => (1), (3, 5, 7, 9)

splitAt <- function(vector, pos){
  unname(split(vector, cumsum(seq_along(vector) %in% pos)))
}

## -- BASIC RUNS & DEMONSTRATIONS --

## -- part (a) (r1, r2, r3) & (b) (s1, s2, s3)

## reps = 500 demonstrates working code, but reps at around 10k
## demonstrates convergence and is decently fast! (around 1 min for r1,
## and 10~20 minutes roughly for r2, r3 on 24 logical cores (AMD 5900x))

r1 <- mc(m=20, n=20, test=c(t1, t2), xdist=list(0, 1), ydist=list(0, 1),
          reps = 500)

r2 <- mc(m=200, n=200, test=c(t1, t2), xdist=list(0, 1), ydist=list(0, 1),
         reps = 500)

r3 <- mc(m=500, n=100, test=c(t1, t2), xdist=list(0, 1), ydist=list(0, 1),
         reps = 500)

s1 <- mc(m=20, n=20, test=c(t1, t2), xdist=list(0, 5), ydist=list(0, 1),
         reps = 500)

s2 <- mc(m=200, n=200, test=c(t1, t2), xdist=list(0, 5), ydist=list(0, 1),
         reps = 500)

s3 <- mc(m=500, n=100, test=c(t1, t2), xdist=list(0, 5), ydist=list(0, 1),
         reps = 500)

## To view the data for r1, use View(r1). To reference a specific result,
## use r1[[x]], where x is the x-th test statistic employed.