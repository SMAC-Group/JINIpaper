mu <- seq(0,1,by=.01)
FP <- 0.0
FN <- .1
mu_star <- FP * (1.0 - mu) + (1.0 - FN) * mu
eps <- .1
select_outliers <- order(abs(mu-.5), decreasing = TRUE)[1:10]
df <- rep("normal", length(mu))
df[select_outliers] <- "outlier"
cbind(mu, mu_star, df)
