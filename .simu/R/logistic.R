## "Logistic Regression"

## Here we present how to perform the simulation studies in @zhang2022flexible for the logistic regression, with and without misclassification. We warn the readers that the code presented here may take time to run on a personal computer. Indeed, for the paper, the computations were performed on the "Baobab" and "Yggdrasil" HPC clusters provided by the University of Geneva. 

## Setting

## We use the following setting for the logistic regression, with and without rounding.

# packages
library(JINIpaper)
library(stocapp)

## simulation specifics
MC <- 1000 # number of simulations
n <- 2000 # sample size
p <- 200 # number of coefficients
## Non-zero coefficients:
## intercept is 1
beta_nz <- c(1,rep(3,5),rep(-5,5),rep(7,5))
## Coefficients:
beta <- c(beta_nz,rep(0,p-15))

## seeds for random number generator
set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC) # for the response
seed$design <- sample.int(1e7,MC) # for the design
seed$sc <- sample.int(1e7,MC) # for the stochastic approximation

## IB specifics
H <- 200 # number of simulated estimators
M <- 50 * H # maximum number of iteration for stochastic approximation
bbc_control <- stocappControl(H=H, maxit=1L) # control for bootstrap bias correction (BBC)
jini_control <- stocappControl(maxit=M, verbose = TRUE) # control for iterative bootstrap (IB)

## Logistic regression without misclassification 

## Here is the code for the simulation of the logistic regression computed on the actual responses without misclassification. 

## For saving the results
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            bbc = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p,sd=p^(-.5)), nrow = n) # simulate the design
  logistic_object <- make_logistic(x, beta) # see ?make_logistic
  # simulate logistic responses
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m])) 
  
  ##------ MLE estimation ----------------
  # Note: here the MLE is the initial estimator
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)
  
  ##------ Bootstrap bias correction -----
  bbc_control$seed <- seed$sc[m] # update the seed
  fit_bbc <- stocapp(fit_mle, control = bbc_control) # compute BBC
  res$bbc[m,] <- getEst(fit_bbc) # retrieve estimator
  
  ##------ IB bias correction ------------
  jini_control$seed <- seed$sc[m] # update the seed
  fit_jini <- stocapp(fit_mle, control = jini_control) # compute IB
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}

## Logistic regression with misclassified responses

## We first need to specify the misclassification mechanism, i.e. how the responses are misclassified. 

## Random misclassification:
## False Positive (FP) follows rbeta(1.5,50)
## False Negative (FN) follows rbeta(1.1,10)
fp <- 1.5 / 51.5 # mean FP
fn <- 1.1 / 21.1 # mean FN

## Function to simulate responses with random misclassification
## see ?ib::simulation and ?ib::control for more details 
random_misclassification <- function(object, control, extra){
  ftd <- fitted(object)
  n <- length(ftd)
  N <- n * control$H
  u1 <- rbeta(N, shape1 = 1.5, shape2 = 50)
  u2 <- rbeta(N, shape1 = 1.1, shape2 = 20)
  mu_star <- u1 * (1 - rep(ftd,control$H)) + (1 - u2) * rep(ftd,control$H)
  matrix(ifelse(runif(N) > mu_star,0,1), ncol=control$H)
}

jini_control <- stocappControl(maxit = M, sim = random_misclassification) # control for stochastic approximation

## Here is the code for the simulation of the logistic regression computed on the misclassified responses. 

## For saving the results
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p,sd=p^(-.5)), nr=n) # simulate the design
  logistic_object <- make_logistic(x, beta) # see ?make_logistic
  # simulate logistic with misclassified responses
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m], sim=random_misclassification)) 
  
  ##------ MLE estimation ----------------
  fit_mle <- logistic_misclassification_mle(x,y,fp,fn) # see ?logistic_misclassification_mle
  res$mle[m,] <- fit_mle
  
  ##------ Initial estimator (inconsistent) ----------------
  fit_initial <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$initial[m,] <- coef(fit_mle)
  
  ##------ IB bias correction ------------
  jini_control$seed <- seed$sc[m] # update the seed
  fit_jini <- stocapp(fit_initial, control = jini_control) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}