## Charles Shi
## Question 2 for the Research Professional Assessment for Shaikh
## & Tabord-Meehan.
## github repository:
## https://github.com/Skyrior/Simplified-Covariate-Adaptive-Permutation-Test
## R version: 4.0.3

library(purrr) ## for lfold (reduce) use install.packages("purrr") if necessary
library(foreach)
library(doParallel) ## speed up the loop if we have multicore.
## for reference, on a 12-core system, the speedup
## is around 8x.

cores <- detectCores(logical = FALSE)
cl <- makeCluster(cores)
registerDoParallel(cl, cores=cores)

## -- debug --
## Set to true to turn on most messages in functions to trace
## what is happening to the data. Warning: only turn on debug
## and run small tests.

debug = FALSE

## -- gm & gn --
## Default values for generated datasets.

gm=10
gn=10



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
  if(debug) message("Calling test statistic 1 with ", data, ", ", m,
                    " and ", n)
  
  ## input checks
  if(!reduce(data,
             function(a, b){a && is.numeric(b)},
             .init=TRUE)) stop("The data supplied to test statistic 1",
                               " is not numeric.")
  if(!is.numeric(m) || !is.numeric(n)) stop("m or n supplied ",
                                            "to test statistic 1 is not",
                                            " numeric.")
  if(!((m+n)==length(data))) message("Warning: m+n =/= N in test",
                                     " statistic 1.")
  
  b <- splitAt(data, m+1)
  x <- b[[1]]
  y <- b[[2]]
  
  x.mean <- mean(x)
  y.mean <- mean(y)
  
  r <- abs((sqrt(length(data)))*(x.mean-y.mean))
  
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
  if(debug) message("Calling test statistic 2 with ", data, ", ", m,
                    " and ", n)
  
  ## input checks
  if(!reduce(data,
             function(a, b){a && is.numeric(b)},
             .init=TRUE)) stop("The data supplied to test statistic 2",
                               " is not numeric.")
  if(!is.numeric(m) || !is.numeric(n)) stop("m or n supplied ",
                                            "to test statistic 2 is not",
                                            " numeric.")
  if(!((m+n)==length(data))) message("Warning: m+n =/= N in test",
                                     " statistic 2.")
  
  b <- splitAt(data, m+1)
  x <- b[[1]]
  y <- b[[2]]
  N <- length(data)
  
  t1 <- t1(data, m, n)
  x.var <- var(x)
  y.var <- var(y)
  
  r <- t1/(sqrt((N/m)*x.var + (N/m)*y.var))
  
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
## List[] - a list of the second components of the permutation, in order.
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
  
  ## create a list from 1 to N
  li <- 1:N
  
  ## produces a random, unbiased permutation with Fisher-Yates shuffle
  r <- sample(li)
  
  if(debug) message("Permutation generated in get.permute: ", r)
  
  ## output the permutation
  return(r)
}

## -- apply.permute (List list, List permutation) --
## Applies a permutation on the given list, effectively shuffling it.
## 
## -- Arguments --
## list: Any generic list.
##
## -- Optional Arguments --
## permutation: If supplied, applies the given permutation on the
## list. If not supplied, the function generates a random permutation.
##
## -- Returns --
## List[] - The shuffled list.

apply.permute <- function(list, permutation=NULL){
  
  if(debug) message("Calling apply.permute() with arguments ",
                    list, ", ", permutation)
  
  ## input checks
  if(!(is.vector(list))) stop("Supplied a non-vector in apply.permute()")
  
  perm <- if(!is.null(permutation)){
    permutation
  } else {
    get.permute(length(list))
  }
  tryCatch({if(!isTRUE(all(perm == floor(perm)))){
              stop("Supplied permutation in apply.permute() is not a list of",
                   "integers.")}},
           error = function(a){
             stop("Supplied permutation in apply.permute() is not a",
                  "list of numbers.")
           })
  if(debug) message("Permutation used in apply.permute: ", perm)
  
  r <- rep(NA, length(list))
  
  ## permutes the input list
  for(v in 1:length(list)){
    r[v] <- list[perm[v]]
  }
  
  return(r)
  
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
  if(debug) message("Calling apply.teststatistic() with arguments ",
                    list, ", ", "with test statistic ",
                    deparse(substitute(teststatistic)))
  if(!is.function(teststatistic)) stop("Supplied test statistic",
                                       "is not a function!")
  if(debug) message("Test statistic in apply.teststatistic takes",
                    " the following arguments: ", formalArgs(teststatistic))
  r <- tryCatch({teststatistic(list, m, n)},
                error = function(a){
                  stop("Failed to apply the test statistic on the data!",
                       "\n", "Supplied list: ", list, "\n",
                       "Supplied test statistic: ", teststatistic, "\n",
                       "Error: ", a)
                })
  
  if(!is.numeric(r)) message("Warning: The test statistic generated, ", r,
                             " is not a number.")
  
  return(r)
}

## -- apply.tpermute (List data, function teststatistic, int rep) --
## Applies the given test statistic on rep numbers of random permutations of
## the data.
## 
## -- Arguments --
## data: Any generic list.
## teststatistic: Any appropriate function as long as teststatistic takes
## one list input.
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
  
  ## Number of loops each core should handle at most:
  
  chunk.size <- ceiling(rep/cores)
  
  ## Parallel processing gives noticeable improvement
  ## when m, n gets large.
  
  r <- foreach(i=1:cores, .combine='c') %dopar%
    {
      res <- rep(NA, chunk.size)
      for(x in ((i-1)*chunk.size+1):(i*chunk.size)){
        if(x <= rep){
          res[x-(i-1)*chunk.size] <- 
            apply.teststatistic(apply.permute(data),
                                teststatistic, m=m, n=n)
        }
        
      }
      res
    }
  
  return(r)
}

## -- reject (List data, function teststatistic, int rep, 
##             int alpha) --
## Determine whether to reject the original data based on the
## two-sample permutation test. Specifically, we reject if
## T_{m, n}(Z_1, ..., Z_N) > c_{1-alpha}, where T_{m, n} is the test statistic,
## c_{1-alpha} is the (1-alpha)-quantile of the genereated list by applying
## T_{m, n} on random permutations of the original data, (Z_1, ..., Z_n).
## 
## -- Arguments --
## data: Any generic list.
## teststatistic: Any appropriate function as long as teststatistic takes
## one list input.
##
## -- Optional Arguments --
## rep: Number of replications. (Default: 500)
## alpha: The threshold of rejection.
## according to the given seed.
## m: Demarcates where to cutoff the data to get (X_1, ..., X_m). Default
##    is the gm variable.
## n: Sanity check to ensure m+n adds up. Default is the gn variable.
##
## -- Returns --
## r - The list of outputs of the test statistic as applied to the data
##     rep number of times.

reject <- function(data, teststatistic, rep=500, alpha=0.05,
                   m=gm, n=gn){
  if(debug) message("Calling reject() with arguments ",
                    data, ", ", "with test statistic ",
                    deparse(substitute(teststatistic)), ", ",
                    rep, " ", alpha)
  
  ## get the test statistics on the permuted data
  tperms <- apply.tpermute(data, teststatistic, rep=rep,
                           m=m, n=n)
  
  ## check if every entry is a numeric value.
  if(!reduce(tperms,
            function(a, b){a && is.numeric(b)},
            .init=TRUE)) stop("The test statistic generated",
                                  " a non-numeric output.")
  
  ## get the (1-alpha)-quantile. Default type 7 (continuous quantile).
  c <- quantile(tperms, probs=(1-alpha), na.rm=TRUE)
  
  if(debug) message("Computed (1-", alpha, ")-quantile: ", c)
  
  ## get the test statistic for the original data.
  origt <- apply.teststatistic(data, teststatistic, m=m, n=n)
  
  if(debug) message("Original Test Statistic: ", origt)
  
  ## decide whether to reject
  rej <- ifelse(origt>c, 1, 0)
  return(rej)
}

## -- splitAt(vector, pos) --
## Splits vector at position.
## Ex: (1, 3, 5, 7, 9) with pos=2 => (1), (3, 5, 7, 9)

splitAt <- function(vector, pos){
  unname(split(vector, cumsum(seq_along(vector) %in% pos)))
}
