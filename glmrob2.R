## "Robust Logistic Regression"

## Here we present how to perform the simulation studies in @zhang2022flexible for the robust logistic regression, with and without outliers We warn the readers that the code presented here may take time to run on a personal computer. Indeed, for the paper, the computations were performed on the "Baobab" and "Yggdrasil" HPC clusters provided by the University of Geneva. 

## Setting

## We use the following setting for the robust logistic regression, with and without outliers.

# packages
library(JINIpaper)
library(stocapp)
library(glmrob2)

## simulation specifics
MC <- 1000 # number of simulations
n <- 200 # sample size
p <- 20 # number of coefficients
## Non-zero coefficients:
## intercept is 1
beta_nz <- c(1,rep(3,1),rep(-5,1),rep(7,1))
## Coefficients:
beta <- c(beta_nz,rep(0,p-3))

## seeds for random number generator
set.seed(895724)
seed <- vector("list",3)
seed$process <- sample.int(1e7,MC) # for the response
seed$design <- sample.int(1e7,MC) # for the design
seed$sc <- sample.int(1e7,MC) # for the stochastic approximation

## IB specifics
H <- 200 # number of simulated estimators
M <- 100 * H # maximum number of iteration for stochastic approximation
jini_control <- stocappControl(maxit=M) # control for iterative bootstrap (IB)

## Logistic regression without outliers

## Here is the code for the simulation of the robust logistic regression computed without outliers. 

## We consider Tukey's bisquare loss function (without correction for Fisher consistency).
## For comparison, we consider robust logistic with Huber loss (with correction for Fisher consistency) and Branco-Yohai specific robust logistic.

## For saving the results
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC),
            robHub = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p,sd=p^(-.5)), nrow = n) # simulate the design
  logistic_object <- make_logistic(x, beta, robust=TRUE) # see ?make_logistic
  # simulate logistic responses
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m])) 
  
  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)
  
  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  fit_robHub <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle")
  res$robHub[m,] <- coef(fit_robHub)
  
  ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  res$robBY[m,] <- coef(fit_robBY)
  
  ##------ Initial estimator (inconsistent, Tukey bisquare) ----------------
  fit_initial <- glmrob2(y ~ x, family=binomial(link="logit"), control = glmrob2.control(maxit=200))
  res$initial[m,] <- coef(fit_initial)
  
  ##------ Stochastic approximation bias correction ------------
  jini_control$seed <- seed$sc[m] # update the seed
  fit_jini <- stocapp(fit_initial, control = jini_control) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}

## Logistic regression with outliers on responses

## We first need to specify the outliers mechanism. 

## Function to simulate responses with outliers
## see ?ib::simulation and ?ib::control for more details 
outlying_mechanism <- function(object, control, extra){
  ftd <- fitted(object)
  n <- length(ftd)
  N <- n * control$H
  n_out <- round(n * control$eps)
  y <- matrix(ifelse(runif(N) > rep(ftd,control$H), FALSE, TRUE), ncol=control$H)
  select_outliers <- order(abs(fitted(object)-0.5), decreasing = T)[1:n_out]
  y[select_outliers] <- !y[select_outliers]
  as.integer(y)
}

jini_control_out <- stocappControl(maxit = M, sim = outlying_mechanism, eps = 0.02) # control for stochastic approximation
jini_control <- stocappControl(maxit = M) # control for stochastic approximation

## Here is the code for the simulation of the logistic regression computed on the misclassified responses. 

## For saving the results
res <- list(mle = matrix(ncol=p+1,nrow=MC),
            initial = matrix(ncol=p+1,nrow=MC),
            jini = matrix(ncol=p+1,nrow=MC),
            jini_oracle = matrix(ncol=p+1,nrow=MC),
            robHub = matrix(ncol=p+1,nrow=MC),
            robBY = matrix(ncol=p+1,nrow=MC))

for(m in 1:MC){
  ##------- simulate the process ---------
  set.seed(seed$process[m]) # set the seed
  x <- matrix(rnorm(n*p,sd=p^(-.5)), nr=n) # simulate the design
  logistic_object <- make_logistic(x, beta, robust = TRUE) # see ?make_logistic
  # simulate logistic with misclassified responses
  y <- simulation(logistic_object, 
                  control = list(seed=seed$process[m], sim=outlying_mechanism, eps=0.02)) 
  
  ##------ MLE estimation ----------------
  fit_mle <- glm(y ~ x, family=binomial(link="logit"), control = glm.control(maxit=200))
  res$mle[m,] <- coef(fit_mle)
  
  ##------ Robust estimator (Cantoni-Ronchetti) ----------------
  fit_robHub <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobMqle.control(maxit=200), method = "Mqle")
  res$robHub[m,] <- coef(fit_robHub)
  
  ##------ Robust estimator (Branco-Yohai) ----------------
  fit_robBY <- glmrob(y ~ x, family=binomial(link="logit"), control = glmrobBY.control(maxit=200), method = "BY")
  res$robBY[m,] <- coef(fit_robBY)
  
  ##------ Initial estimator (inconsistent, Tukey bisquare) ----------------
  fit_initial <- glmrob2(y ~ x, family=binomial(link="logit"), control = glmrob2.control(maxit=200))
  res$initial[m,] <- coef(fit_initial)
  
  ##------ Stochastic approximation bias correction ------------
  jini_control$seed <- seed$sc[m] # update the seed
  fit_jini <- stocapp(fit_initial, control = jini_control) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  res$jini[m,] <- getEst(fit_jini) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  ##------ Stochastic approximation bias correction (oracle) ------------
  jini_control_out$seed <- seed$sc[m] # update the seed
  fit_jini2 <- stocapp(fit_initial, control = jini_control_out) # compute JINI using stochastic approximation, see ? stocapp::stocapp for more details
  res$jini_oracle[m,] <- getEst(fit_jini2) # retrieve estimator
  # getIteration(fit_jini), if one wants to retrieve the number of iterations
  
  # print the iteration
  cat(m,"\n")
}