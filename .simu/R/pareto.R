library(VGAM)
library(JINIpaper)
# library(RcppDE)
# library(MASS)

starting_values <- function(y, p){
  x <- log(y[-which.min(y)])
  stat <- boxplot.stats(x)
  c(exp(1/mean(x[!x%in%stat$out])), rep(0,p))
}
# 
starting_values <- function(y, x){
  y0 <- log(log(y[-which.min(y)]))
  x0 <- x[-which.min(y),]
  stat <- boxplot.stats(y0)
  ind <- !y0%in%stat$out
  y0 <- y0[ind]
  x0 <- x0[ind,]
  fit_vglm <- vglm(y0 ~ x0, paretoff)
  coef(fit_vglm)
}


# Pareto
n <- 150
p0 <- 12
p <- p0 - 2
x0 <- 5
set.seed(9237423)
x <- matrix(rnorm(n * p, sd= 1 / sqrt(n)),nr=n)
beta0 <- c(1.5,-2,-2,2,2,rep(c(-.1,.1),(p-4)/2)[1:(p-4)])
theta0 <- c(beta0,x0)
means <- exp(cbind(1,x) %*% beta0)
# summary(means)
B <- 200
fit_estimators0 <- fit_estimators1 <- fit_estimators2 <- fit_estimators3 <- fit_estimators4 <- matrix(nr=B,nc=p0-1)
scale_est <- matrix(nr=B,nc=3)
cc <- 13

for(b in 1:B) {
  set.seed(48374+b)
  y <- rpareto(n, scale = x0, shape = means)
  # y[1] <- 1e6
  # y[sample.int(n,5)] <- rpareto(5, scale = 5 * x0, shape = 3)
  sv <- starting_values(y, p0-2)
  # fit_glm <- rlm(log(y) ~ x)
  
  
  fit1 <- paretoMle(y, x, start = sv, verbose = F)
  fit1 <- vglm(y ~ x, paretoff, coefstart = sv)
  # sv <- DEoptim(fn = paretoWmle_of, lower = rep(-3,p0-2), upper = rep(3,p0-2), y=y, x=x, c=cc, control = DEoptim.control(trace = FALSE, NP = 100 * length(beta0)))
  # paretoWmle1(y, x, start = c(1.7,rep(0,p0-3)), c = cc, verbose=T)
  # paretoWmle1H(y, x, start = c(1.7,rep(0,p0-3)), c = cc, verbose=T)
  # paretoWmle1(y, x, start = fit1$coefficients, c = cc, verbose=T)
  # paretoWmle1(y, x, start = beta0, c = cc, verbose=T)
  # paretoWmle1(y, x, start = sv$optim$bestmem, c = cc, verbose=T)
  
  fit2 <- paretoWmle1(y, x, start = sv, c = cc, verbose=F)
  if(fit2$conv !=0 ) next
  if(any(abs(fit2$coefficients)>1e3)) next
  fit_estimators2[b,] <- fit2$coefficients
  
  
  # fit0 <- vglm(y ~ x, paretoff)
  # fit_estimators0[b,] <- coef(fit0)

  fit_estimators1[b,] <- fit1$coefficients
  
  
  fit3 <- paretoWmle1_ib(x, c(fit2$coefficients, fit2$scale), c = cc, H=100, seed = 98473+b, verbose=T)
  # fit_estimators3[b,] <- fit3$estimate[1:(p0-1)]
  # scale_est[b,1] <- min(y)
  # scale_est[b,3] <- fit3$estimate[p0]
  # 
  # fit4 <- paretoWmle1_ib(x, c(fit1$coefficients, fit1$scale), c = Inf, H=100, seed = 98473+b)
  # fit_estimators4[b,] <- fit4$estimate[1:(p0-1)]
  # scale_est[b,2] <- fit4$estimate[p0]
  
  fit4 <- pareto_vglm_ib(x, thetastart = c(Coef(fit1), fit1@extra$scale), H = 100, verbose = T, seed = 98473+b)
  
  cat(b,"\n")
}

y_lim <- c(-1,1)
par(mfrow = c(1,2))
# boxplot(fit_estimators0-matrix(beta0,nc=p0-1,nr=B,byrow=T),ylim=y_lim);abline(h=0)
boxplot(fit_estimators1-matrix(beta0,nc=p0-1,nr=B,byrow=T),ylim=y_lim,main="MLE");abline(h=0)
boxplot(fit_estimators2-matrix(beta0,nc=p0-1,nr=B,byrow=T),ylim=y_lim);abline(h=0)
# boxplot(fit_estimators3-matrix(beta0,nc=p0-1,nr=B,byrow=T),ylim=y_lim,main="IB on WMLE");abline(h=0)
par(mfrow = c(1,1))

colSums(is.na(fit_estimators2))
abs(colMeans(fit_estimators1, na.rm=T) - beta0)
abs(colMeans(fit_estimators3, na.rm=T) - beta0)
abs(colMeans(fit_estimators4, na.rm=T) - beta0)

boxplot(scale_est)
sum(apply(fit_estimators4,2,sd,na.rm=T)) / sum(apply(fit_estimators3,2,sd,na.rm=T))


# Truncated Pareto
n <- 200
p0 <- 6
p <- p0 - 2
x0 <- 0.5
set.seed(9237423)
x <- matrix(rnorm(n * p, sd= 1 / sqrt(n)),nr=n)
beta0 <- c(1,-2,-2,2,2)
theta0 <- c(beta0,x0)
means <- exp(cbind(1,x) %*% beta0)
B <- 200
fit_estimators <- fit_estimators2 <- matrix(nr=B,nc=p0-1)

for(b in seq_len(B)) {
  set.seed(48374+b)
  y <- rtruncpareto(n, lower = x0, upper= 3, shape = means)
  fit <- vglm(y ~ x, truncpareto(lower=x0, upper = 3))
  fit_estimators[b,1:(p+1)] <- coef(fit)

  y[sample.int(n,10)] <- rtruncpareto(n = 10, lower = x0, upper= 3, shape = 3)
  fit <- vglm(y ~ x, truncpareto(lower=x0, upper = 3))
  fit_estimators2[b,1:(p+1)] <- coef(fit)
  cat(b,"\n")
}

par(mfrow = c(2,1))
boxplot(fit_estimators-matrix(beta0,nc=p0-1,nr=B,byrow=T));abline(h=0)
boxplot(fit_estimators2-matrix(beta0,nc=p0-1,nr=B,byrow=T));abline(h=0)
par(mfrow = c(1,1))



n <- 200
x0 <- 0.5
beta0 <- 3
theta0 <- c(beta0,alpha)
B <- 200
fit_estimators <- fit_estimators2 <- matrix(nr=B,nc=2)

for(b in seq_len(B)) {
  set.seed(48374+b)
  y <- rpareto(n, scale = x0, shape = 3)
  fit <- vglm(y ~ 1, paretoff)
  fit_estimators[b,1] <- Coef(fit)
  fit_estimators[b,2] <- fit@extra$scale
  
  y[sample.int(n,10)] <- rpareto(10, scale = 10 * x0, shape = 3)
  fit <- vglm(y ~ 1, paretoff)
  fit_estimators2[b,1] <- Coef(fit)
  fit_estimators2[b,2] <- fit@extra$scale
  cat(b,"\n")
}

par(mfrow = c(2,1))
boxplot(fit_estimators-matrix(theta0,nc=2,nr=B,byrow=T));abline(h=0)
boxplot(fit_estimators2-matrix(theta0,nc=2,nr=B,byrow=T));abline(h=0)
par(mfrow = c(1,1))

fit@extra # The estimate of alpha is here
head(fitted(fit))
coef(fit, matrix = TRUE)

x <- rnorm(1000, sd =.1)
alpha <- 2; kay <- exp(3 + x)
pdata <- data.frame(y = rpareto(n = 1000, scale = alpha, shape = kay))
fit <- vglm(y ~ x, paretoff, data = pdata, trace = TRUE)
fit@extra # The estimate of alpha is here
head(fitted(fit))
with(pdata, mean(y))

summary(fit) # Standard errors are incorrect!!

# Here, alpha is assumed known
fit2 <- vglm(y ~ 1, paretoff(scale = alpha), data = pdata, trace = TRUE)
fit2@extra # alpha stored here
head(fitted(fit2))
coef(fit2, matrix = TRUE)
summary(fit2) # Standard errors are okay


