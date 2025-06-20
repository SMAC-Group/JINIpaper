## Here we present how to perform the simulation studies in @zhang2022flexible for the beta regression, with and without rounding. We warn the readers that the code presented here may take time to run on a personal computer. Indeed, for the paper, the computations were performed on the "Baobab" and "Yggdrasil" HPC clusters provided by the University of Geneva. 

## Setting

## We use the following setting for the beta regression, with and without rounding.

## packages
library(IBpaper)
library(JINIpaper)
library(ib)
library(betareg)

## simulation specifics
MC <- 1000 # number of simulations
n <- 200 # sample size
p <- 30 # number of coefficients
# Non-zero coefficients:
# intercept is -0.5
beta_nz <- c(-.5,rep(1,5),rep(-1.5,5),rep(2,5))
# Coefficients:
beta <- c(beta_nz,rep(c(-.1,.1),p-15)[1:(p-15)])
gamma <- log(5) # precision parameter
theta <- c(beta,gamma)

# seeds for random number generator
set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC) # for the response
seed$design <- sample.int(1e7,MC) # for the design
seed$ib <- sample.int(1e7,MC) # for the iterative bootstrap

# IB specifics
H <- 200 # number of simulated estimators
bbc_control <- ibControl(H=H, maxit=1L, constraint=FALSE) # control for bootstrap bias correction (BBC)
ib_control <- ibControl(H=H, maxit=50L, constraint=FALSE) # control for iterative bootstrap (IB)

## Beta regression without rounding

## Here is the code for the simulation of the beta regression computed on the actual responses without rounding. 

## For saving the results
res <- list(mle = matrix(ncol=p+2,nrow=MC),
            bbc = matrix(ncol=p+2,nrow=MC),
            jini = matrix(ncol=p+2,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p, sd = 4/sqrt(n)), nrow = n) # simulate the design
  betareg_object <- make_betareg(x, theta) # see ?make_betareg
  # simulate beta responses
  y <- simulation(betareg_object, 
                  control = list(seed=seed$process[m])) 
  
  ##------ MLE estimation ----------------
  # Note: here the MLE is the initial estimator
  fit_mle <- NULL
  try(fit_mle <- betareg(y ~ x), silent=T)
  if(is.null(fit_mle)) next
  res$mle[m,] <- coef(fit_mle)
  
  ##------ Bootstrap bias correction -----
  bbc_control$seed <- seed$ib[m] # update the seed
  fit_bbc <- ib(fit_mle, control = bbc_control) # compute BBC , see ? ib::ib for more details
  res$bbc[m,] <- getEst(fit_bbc) # retrieve estimator
  
  ##------ IB bias correction ------------
  ib_control$seed <- seed$ib[m] # update the seed
  fit_jini <- ib(fit_mle, control = ib_control) # compute IB
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}

## Beta regression with rounded responses

We first need to specify the rounding mechanism, i.e. how the responses are rounded. 

## Function to simulate responses with rounding
## see ?ib::simulation and ?ib::control for more details 
rounding <- function(object, control, extra=NULL){
  simulate_betareg <- getFromNamespace("simulate_betareg", ns = "ib")
  y <- simulate_betareg(object, control$H)
  n <- length(y) / control$H
  y <- round(y,1)
  y <- (y*(n-1) + 0.5)/n
  matrix(y, ncol=control$H)
}
ib_control <- ibControl(H = H, maxit = 50L, sim = rounding) # control for iterative bootstrap (IB)

## Here is the code for the simulation of the beta regression computed on the rounded responses. 

## For saving the results
res <- list(mle = matrix(ncol=p+2,nrow=MC),
            initial = matrix(ncol=p+2,nrow=MC),
            jini = matrix(ncol=p+2,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p, sd = p^(-.5)), nrow = n) # simulate the design
  betareg_object <- make_betareg(x, theta) # see ?make_betareg
  # simulate beta regression with rounded responses
  y <- simulation(betareg_object, 
                  control = list(seed = seed$process[m], 
                                 sim = rounding)) 
  
  ##------ Initial estimator (inconsistent) ----------------
  fit_initial <- NULL
  try(fit_initial <- betareg(y ~ x), silent=T)
  if(is.null(fit_initial)) next
  res$mle[m,] <- coef(fit_initial)
  
  ##------ MLE estimation ----------------
  fit_mle <- NULL
  sv <- fit_initial # starting values
  sv[p+2] <- log(sv[p+2]) # log transform
  try(fit_mle <- em(y,cbind(1,x),sv), silent=T) # see ?em
  if(!is.null(fit_em)){
    res$mle[m,] <- fit_mle$par
    res$mle[m,p+2] <- exp(fit_mle$par[p+2])
  }
  
  ##------ IB bias correction ------------
  ib_control$seed <- seed$ib[m] # update the seed
  fit_jini <- ib(fit_initial, control = ib_control) # compute IB
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}