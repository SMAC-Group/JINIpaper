#' @title Find tuning constant for robust logistic (Mqle) with Tukey's weight
#' @description
#' Find tuning constant in the robust logistic with Tukey's weight
#' for a given desired level of efficiency depending on the design and
#' the true parameter beta.
#' @param x design matrix
#' @param beta parameter
#' @param c_min minimum tuning parameter value
#' @param c_max maximum tuning parameter value
#' @param eff level of efficiency (95% by default)
#' @param inconsistent efficiency at the inconsistent estimator?
#' @export
find_tuning_constantMqle <- function(x, beta, c_min, c_max, eff=0.95, inconsistent=FALSE){
  if(inconsistent) {
    of <- function(c, x, beta, eff) abs(roblogisticMqleVar(x, beta, c)$efficiency1 - eff) 
  } else {
    of <- function(c, x, beta, eff) abs(roblogisticMqleVar(x, beta, c)$efficiency - eff) 
  }
  opt <- optimise(of, interval=c(c_min, c_max), x=x, beta=beta, eff=eff)
  opt$minimum
}

#' @title Find tuning constant for robust logistic (Wmle) with Tukey's weight
#' @description
#' Find tuning constant in the robust logistic with Tukey's weight
#' for a given desired level of efficiency depending on the design and
#' the true parameter beta.
#' @param x design matrix
#' @param beta parameter
#' @param c_min minimum tuning parameter value
#' @param c_max maximum tuning parameter value
#' @param eff level of efficiency (95% by default)
#' @param inconsistent efficiency at the inconsistent estimator?
#' @export
find_tuning_constantWmle <- function(x, beta, c_min, c_max, eff=0.95, inconsistent=FALSE){
  if(inconsistent) {
    of <- function(c, x, beta, eff) abs(roblogisticWmleVar(x, beta, c)$efficiency1 - eff) }
  else {
    of <- function(c, x, beta, eff) abs(roblogisticWmleVar(x, beta, c)$efficiency - eff) 
  }
  opt <- optimise(of, interval=c(c_min, c_max), x=x, beta=beta, eff=eff)
  opt$minimum
}


#' @title Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
roblogisticWmle1_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  differences <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    tmp_pi <- matrix(NA_real_, nrow=p, ncol=H)
    for(h in seq_len(H)){
      seed1 <- seed + h
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- roblogisticWmle1(sim, x, t0, c)
      iter <- 1L
      while(fit_tmp$conv == 1 && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_logistic(t0, x, seed1)
        fit_tmp <- roblogisticWmle1(sim, x, t0, c)
        iter <- iter + 1L
      }
      if(fit_tmp$conv == 1) next
      tmp_pi[,h] <- fit_tmp$coefficients
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) differences[k] <- test_theta
    if(!is.finite(test_theta)) {
      warning("non-finite difference")
      break
    }
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # Alternative stopping criteria, early stop :
    # if(control$early_stop){
    #   if(tt_old <= test_theta){
    #     warning("Algorithm stopped because the objective function does not reduce")
    #     break
    #   }
    # }
    
    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- differences[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta,
       boot = tmp_pi)
}

#' @title Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
logistic_wmle_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  diff <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    tmp_pi <- matrix(NA_real_, nrow=p, ncol=H)
    for(h in seq_len(H)){
      seed1 <- seed + h
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- logistic_wmle(sim, x, c)
      iter <- 1L
      while(fit_tmp$status != 0 && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_logistic(t0, x, seed1)
        fit_tmp <- logistic_wmle(sim, x, c)
        iter <- iter + 1L
      }
      if(fit_tmp$status != 0) next
      tmp_pi[,h] <- fit_tmp$coefficients
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) diff[k] <- test_theta
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # Alternative stopping criteria, early stop :
    # if(control$early_stop){
    #   if(tt_old <= test_theta){
    #     warning("Algorithm stopped because the objective function does not reduce")
    #     break
    #   }
    # }
    
    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- diff[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta,
       boot = tmp_pi)
}

#' @title Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
logistic_mqle_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  diff <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    tmp_pi <- matrix(NA_real_, nrow=p, ncol=H)
    for(h in seq_len(H)){
      seed1 <- seed + h
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- logistic_mqle(sim, x, c)
      iter <- 1L
      while(fit_tmp$status != 0 && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_logistic(t0, x, seed1)
        fit_tmp <- logistic_mqle(sim, x, c)
        iter <- iter + 1L
      }
      if(fit_tmp$status != 0) next
      tmp_pi[,h] <- fit_tmp$coefficients
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) diff[k] <- test_theta
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # Alternative stopping criteria, early stop :
    # if(control$early_stop){
    #   if(tt_old <= test_theta){
    #     warning("Algorithm stopped because the objective function does not reduce")
    #     break
    #   }
    # }
    
    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- diff[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta,
       boot = tmp_pi)
}

#' @title Stochastic approximation for robust logistic regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
roblogisticWmle1_stocapp <- function(x, thetastart, c = 4.685061, maxit=1e4, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  diff <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    seed1 <- seed + k
    sim <- r_logistic(t0, x, seed1)
    fit_tmp <- roblogisticWmle1(sim, x, t0, c)
    iter <- 1L
    while(fit_tmp$conv == 1 && iter < 10L){
      seed1 <- seed + H * k + iter
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- roblogisticWmle1(sim, x, t0, c)
      iter <- iter + 1L
    }
    if(fit_tmp$conv == 1) next
    pi_star <- fit_tmp$coefficients

    # update value
    delta <- pi0 - pi_star
    alpha <- 1.0 / (k + 1.0)
    t1 <- t0 + alpha * delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) diff[k] <- test_theta
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta)
}

#' @title Iterative bootstrap for robust logistic regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
roblogisticMqle1_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
  p <- length(thetastart)
  pi0 <- t0 <- thetastart
  
  # test diff between thetas
  test_theta <- tol + 1
  
  # iterator
  k <- 0L
  differences <- rep(NA_real_, maxit)
  
  # Iterative bootstrap algorithm:
  while(test_theta > tol && k < maxit){
    # update object for simulation
    
    # approximate
    tmp_pi <- matrix(NA_real_, nrow=p, ncol=H)
    for(h in seq_len(H)){
      seed1 <- seed + h
      sim <- r_logistic(t0, x, seed1)
      fit_tmp <- roblogisticMqle1(sim, x, t0, c)
      iter <- 1L
      while(fit_tmp$conv == 1 && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_logistic(t0, x, seed1)
        fit_tmp <- roblogisticMqle1(sim, x, t0, c)
        iter <- iter + 1L
      }
      if(fit_tmp$conv == 1) next
      tmp_pi[,h] <- fit_tmp$coefficients
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) differences[k] <- test_theta
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
    # Alternative stopping criteria, early stop :
    # if(control$early_stop){
    #   if(tt_old <= test_theta){
    #     warning("Algorithm stopped because the objective function does not reduce")
    #     break
    #   }
    # }
    
    # Alternative stopping criteria, "statistically flat progress curve" :
    if(k > 10L){
      try1 <- differences[k:(k-10)]
      try2 <- k:(k-10)
      if(var(try1)<=1e-3) break
      mod <- lm(try1 ~ try2)
      if(summary(mod)$coefficients[2,4] > 0.2) break
    }
    
    # update increment
    k <- k + 1L
    
    # Print info
    if(verbose){
      cat("Iteration:",k,"Norm between theta_k and theta_(k-1):",test_theta,"\n")
    }
    
    # update theta
    t0 <- t1
    
    # update test
    tt_old <- test_theta
  }
  # warning for reaching max number of iterations
  if(k>=maxit) warning("maximum number of iteration reached")
  
  list(iteration = k,
       of = sqrt(drop(crossprod(delta))),
       estimate = t0,
       test_theta = test_theta,
       boot = tmp_pi)
}