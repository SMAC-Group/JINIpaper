r_unif_trunc <- function(n, theta, kappa) {
  x <- runif(n,max=theta)
  t <- runif(n,max=theta*kappa)
  y <- x
  y[t<=x] <- t[t<=x]
  y
}

# ll <- function(theta, y, alpha) {
#   sum(log(1/theta + 1/alpha - 2 * y / theta / alpha))
# }

Likelihood <- function(theta, y, kappa) {
  if(theta < max(y) / kappa){ 
    return(0) } else {
    1/theta^length(y) * prod(1+1/kappa-2*y/kappa/theta)
  }
}


kappa <- .9
theta0 <- 2
y <- r_unif_trunc(10, theta0, kappa)
B <- 1e3
theta <- seq(0,10,length.out=B)
f <- rep(NA_real_, B)
for(i in seq_len(B)) f[i] <- Likelihood(theta[i], y, kappa)
plot(theta, f, type="l")
abline(v=theta0)
abline(v=max(y)/kappa)

fit_mle <- optimise(f = Likelihood, interval = c(max(y) / kappa - 1e2, 10), y = y, kappa = kappa, maximum = TRUE)
abline(v=fit_mle$maximum, col="red")
