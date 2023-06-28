#' @title Iterative bootstrap for robust Pareto regression with inconsistent initial estimator with Tukey's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Tukey's weight (default value is 4.685061)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
paretoWmle1_ib <- function(x, thetastart, c = 4.685061, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
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
      sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
      fit_tmp <- paretoWmle1(sim, x, t0, c)
      iter <- 1L
      fit_ok <- fit_tmp$conv == 0 && all(abs(fit_tmp$coefficients)<1e2)
      while(!fit_ok && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
        fit_tmp <- paretoWmle1(sim, x, t0, c)
        iter <- iter + 1L
        fit_ok <- fit_tmp$conv == 0 && all(abs(fit_tmp$coefficients)<1e2)
      }
      if(!fit_ok) next
      tmp_pi[1:(p-1),h] <- fit_tmp$coefficients
      tmp_pi[p,h] <- fit_tmp$scale
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # Make sure scale remains positive
    if(t1[p]<.001) t1[p] <- exp(log(t0[p]) + log(delta[p]))
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) differences[k] <- test_theta
    if(!is.finite(test_theta)) {
      warning("non-finite difference")
      break
    }
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
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

#' @title Iterative bootstrap for robust Pareto regression with inconsistent initial estimator with Huber's weights
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Huber's weight (default value is 1.345)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
paretoWmle1H_ib <- function(x, thetastart, c = 1.345, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
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
      sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
      fit_tmp <- paretoWmle1H(sim, x, t0, c)
      iter <- 1L
      fit_ok <- fit_tmp$conv == 0 && all(abs(fit_tmp$coefficients)<1e2)
      while(!fit_ok && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
        fit_tmp <- paretoWmle1H(sim, x, t0, c)
        iter <- iter + 1L
        fit_ok <- fit_tmp$conv == 0 && all(abs(fit_tmp$coefficients)<1e2)
      }
      if(!fit_ok) next
      tmp_pi[1:(p-1),h] <- fit_tmp$coefficients
      tmp_pi[p,h] <- fit_tmp$scale
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # Make sure scale remains positive
    if(t1[p]<.001) t1[p] <- exp(log(t0[p]) + log(delta[p]))
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) differences[k] <- test_theta
    if(!is.finite(test_theta)) {
      warning("non-finite difference")
      break
    }
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
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

#' @title Iterative bootstrap for robust Pareto regression
#' @param x a n x p matrix of design
#' @param thetastart an inconsistent estimator (also used as starting values)
#' @param c tuning parameter for Huber's weight (default value is 1.345)
#' @param H number of estimators for Monte Carlo approximation
#' @param maxit max number of iteration for IRWLS
#' @param tol tolerance for stopping criterion
#' @param verbose print info
#' @param seed for random number generator
#' @export
#' @importFrom VGAM vglm Coef
pareto_vglm_ib <- function(x, thetastart, H = 200, maxit=200, tol=1e-7, verbose=FALSE, seed=321){
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
      sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
      fit_tmp <- tryCatch(error = function(cnd) NULL, {vglm(sim ~ x, paretoff)})
      iter <- 1L
      while(is.null(fit_tmp) && iter < 10L){
        seed1 <- seed + H * h + iter
        sim <- r_pareto(t0[1:(p-1)], t0[p], x, seed1)
        fit_tmp <- tryCatch(error = function(cnd) NULL, {vglm(sim ~ x, paretoff)})
        iter <- iter + 1L
      }
      if(is.null(fit_tmp)) next
      tmp_pi[1:(p-1),h] <- Coef(fit_tmp)
      tmp_pi[p,h] <- fit_tmp@extra$scale
    }
    pi_star <- rowMeans(tmp_pi, na.rm=TRUE)
    
    # update value
    delta <- pi0 - pi_star
    t1 <- t0 + delta
    
    # Make sure scale remains positive
    if(t1[p]<.001) t1[p] <- exp(log(t0[p]) + log(delta[p]))
    
    # test diff between thetas
    test_theta <- sum(delta^2)
    if(k>0) differences[k] <- test_theta
    if(!is.finite(test_theta)) {
      warning("non-finite difference")
      break
    }
    
    # initialize test
    if(!k) tt_old <- test_theta+1
    
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