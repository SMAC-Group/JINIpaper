# test roblogistic implementation
eps <- 1e-7
mu <- 1/4
cc <- 5
set.seed(1)
y <- rbinom(1e6, 1, mu)
r <- (y - mu) / sqrt(mu * (1 - mu))


# check derivative of psi
(psi(mu + eps, cc) - psi(mu - eps, cc)) / 0.2e1 / eps
psi1(mu, cc)

# check Expectation of psi
psiV <- Vectorize(psi, vectorize.args = "x")
mean(psiV(r, cc))
Epsi(mu, cc)

# check First derivative of expectation of psi(r)
(Epsi(mu + eps, cc) - Epsi(mu - eps, cc)) / 0.2e1 / eps
dEpsi(mu, cc)

# check Expectation of first derivative of psi(r)
psi1V <- Vectorize(psi1, vectorize.args = "x")
mean(psi1V(r, cc))
Edpsi(mu, cc)

# check Expectation of first derivative of psi(r) times r
psi1r <- function(x, c) psi1(x, c) * x
psi1rV <- Vectorize(psi1r, vectorize.args = "x")
mean(psi1rV(r, cc))
Edpsir(mu, cc)

# check Expectation of psi(r) squared
mean(psiV(r, cc) * psiV(r, cc))
Epsi2(mu, cc) 


eps <- 1e-7
mu <- 1/4
cc <- 3
set.seed(1)
y <- rbinom(1e6, 1, mu)
xp <- 1
s <- abs(y - mu) * xp

(wc(mu + eps, cc) - wc(mu - eps, cc)) / 0.2e1 / eps
wc1(mu, cc)

# check Expectation of Tukey's weight
wcV <- Vectorize("wc", "x")
mean(wcV(s, cc))
Ew(mu, cc, xp)

# check Expectation of first derivative of Tukey's weight times -1^y(y-mu)
wc1V <- Vectorize("wc1", "x")
mean(wc1V(s, cc) * (y - mu) * (-1)^y)
Edwymu(mu, cc, xp)

# check Expectation of Tukey's weight times (y-mu)
mean(wcV(s, cc) * (y - mu))
Ewymu(mu, cc, xp)

# check First derivative of expectation of Tukey's weight times (y-mu)
(Ewymu(mu + eps, cc, xp) - Ewymu(mu - eps, cc ,xp)) / 0.2e1 / eps
dEwymu(mu, cc, xp)

# check  Expectation of Tukey's weight times (y-mu) all squared
mean((wcV(s, cc) * (y-mu))^2)
Ewymu2(mu, cc, xp)


library(JINIpaper)
library(stocapp)

# test implemenetation
beta_nz <- c(1,-3)
n <- 100
p <- 2
beta <- c(beta_nz,rep(0,p-1))
mat_design <- "diag"
if(mat_design=="toeplitz") Sigma <- t(chol(toeplitz(0.8^(0:(p-1))))) else Sigma <- diag(p)
set.seed(87432984)
# simulate the design
# x <- matrix(rnorm(n*p,sd=p^(-.5)), nr=n)
x <- t(Sigma %*% matrix(rnorm(n*p,sd= 4 / sqrt(n)), nr=p))
logistic_object <- make_logistic(x, beta, robust=TRUE)
# simulate logistic
y <- simulation(logistic_object, control = list(seed=321))


check(y, x, beta, 10)

(check(y,x,beta+c(eps,0,0),10)$of - check(y,x,beta-c(eps,0,0),10)$of) / 0.2e1 / eps
(check(y,x,beta+c(0,eps,0),10)$of - check(y,x,beta-c(0,eps,0),10)$of) / 0.2e1 / eps
(check(y,x,beta+c(0,0,eps),10)$of - check(y,x,beta-c(0,0,eps),10)$of) / 0.2e1 / eps

fit_mle <- glm(y ~ x, family = binomial())
# fit_rob <- 
roblogisticMqle1(y, x, coef(fit_mle), verbose = T, c= 14)
roblogisticMqle(y, x, coef(fit_mle), verbose = T, c= 14)
roblogisticMqleVar(x, beta, c = 14) 

roblogisticWmle1(y, x, coef(fit_mle), verbose = T, c= 3)
roblogisticWmle(y, x, coef(fit_mle), verbose = T, c= 3)
roblogisticWmleVar(x, beta, c = 10) 

fit_initial <- logistic_wmle(y, x, 3)
logistic_wmle_ib(x,fit_initial$coefficients,3,H=30,maxit=30,verbose=T)
logistic_wmle_stocapp(x, beta_hat, 3, maxit = 10000, verbose = T)

cc <- seq(1,4,by=.1)
eff <- eff1 <- rep(NA_real_, length(cc))
for(i in seq_along(cc)){
  tmp <- roblogisticWmleVar(x, beta, cc[i])
  eff[i] <- tmp$efficiency
  eff1[i] <- tmp$efficiency1
}
plot(cc,eff,type="l")
lines(cc,eff1,lty=2)


find_tuning_constant <- function(x, beta, c_min, c_max, eff=0.95){
  of <- function(c, x, beta, eff) abs(roblogisticMqleVar(x, beta, c)$efficiency - eff)
  opt <- optimise(of, interval=c(c_min, c_max), x=x, beta=beta, eff=eff)
  opt$minimum
}

find_tuning_constant(x, beta, 4, 10)
