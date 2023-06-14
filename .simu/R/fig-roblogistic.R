# +---------+------+-----+-------+
# | Setting |   n  |  p  | p/n   |
# +---------+------+-----+-------+
# |       1 | 2000 | 200 | 0.100 |
# |       2 | 4000 | 300 | 0.075 |
# |       3 | 8000 | 400 | 0.050 |
# +---------+------+-----+-------+

# Non-zero coefficients:
# intercept is 1
# beta = c(rep(3,5),rep(-5,5),rep(7,5))

source(".simu/R/figures.R")

# -----------------
# Inconsistent case
# -----------------
n_setting <- 3
n_methods <- 6
p <- c(10,15,20) +1
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
setting2 = array(dim=c(MC,p[2],n_methods)),
setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
# i=3
load(paste0(".simu/data/roblogistic_1_setting_",i,".rds"))
# super_res[[i]][,,1] <- res$mle
# i=1
super_res[[i]][,,1] <- res$jimi
# super_res[[i]][,,2] <- res$consistent
super_res[[i]][,,2] <- res$robCR
super_res[[i]][,,3] <- res$robBY
super_res[[i]][,,4] <- res$mle
super_res[[i]][,,5] <- res$br
}
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","CR","BY","MLE","BR")
lab_main <- c("$n=250$","$n=500$","$n=1000$")
# theta0 <- c(1, 2, -3, 4, 0)
theta0 <- c(.3, -2, -4, 0)
# ind <- list(setting1 = list(1,2,3,4:p[1]))
# theta0 <- c(-1,2,-3,4,0)
ind <- list(setting1 = list(1,2,3,4:p[1]),
setting2 = list(1,2,3,4:p[2]),
setting3 = list(1,2,3,4:p[3]))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic_outliers", p = p, y_lim=y_lim)

# -----------------
# Consistent case
# -----------------

n_setting <- 3
n_methods <- 6
p <- c(10,15,20) +1
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
setting2 = array(dim=c(MC,p[2],n_methods)),
setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  # i=3
  load(paste0(".simu/data/roblogistic_2_setting_",i,".rds"))
  # super_res[[i]][,,1] <- res$mle
  # i=1
  super_res[[i]][,,1] <- res$jimi
  # super_res[[i]][,,2] <- res$consistent
  super_res[[i]][,,2] <- res$robCR
  super_res[[i]][,,3] <- res$robBY
  super_res[[i]][,,4] <- res$mle
  super_res[[i]][,,5] <- res$br
}
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","CR","BY","MLE","BR")
lab_main <- c("$n=250$","$n=500$","$n=1000$")
# theta0 <- c(1, 2, -3, 4, 0)
theta0 <- c(.3, -2, -4, 0)
# ind <- list(setting1 = list(1,2,3,4:p[1]))
# theta0 <- c(-1,2,-3,4,0)
ind <- list(setting1 = list(1,2,3,4:p[1]),
            setting2 = list(1,2,3,4:p[2]),
            setting3 = list(1,2,3,4:p[3]))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic", p = p, y_lim=y_lim)
