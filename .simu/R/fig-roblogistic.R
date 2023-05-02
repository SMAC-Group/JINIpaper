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
n_methods <- 4
p <- 2:4 * 100+1
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/roblogistic_1_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$initial
  super_res[[i]][,,2] <- res$jini
  super_res[[i]][,,3] <- res$consistent
  super_res[[i]][,,4] <- res$robBY
}
lab_beta <- c("$\\beta_0$","$\\beta_{1:5}$","$\\beta_{6:10}$",
              "$\\beta_{11:15}$","$\\beta_{16:p}$")
lab_methods <- c("Initial","JINI","CR","BY")
lab_main <- c("$n=2000$", "$n=4000$", "$n=8000$")
theta0 <- c(1, 3, -5, 7, 0)
ind <- list(setting1 = list(1,2:6,7:11,12:16,17:p[1]),
            setting2 = list(1,2:6,7:11,12:16,17:p[2]),
            setting3 = list(1,2:6,7:11,12:16,17:p[3]))
cols <- gg_color_hue(4, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic_outliers", p = p, y_lim=y_lim)

# -----------------
# Consistent case
# -----------------

n_setting <- 3
n_methods <- 4
p <- 2:4 * 100+1
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/roblogistic_2_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$initial
  super_res[[i]][,,2] <- res$jini
  super_res[[i]][,,3] <- res$consistent
  super_res[[i]][,,4] <- res$robBY
}
lab_beta <- c("$\\beta_0$","$\\beta_{1:5}$","$\\beta_{6:10}$",
              "$\\beta_{11:15}$","$\\beta_{16:p}$")
lab_methods <- c("Initial","JINI","CR","BY")
lab_main <- c("$n=2000$", "$n=4000$", "$n=8000$")
theta0 <- c(1, 3, -5, 7, 0)
ind <- list(setting1 = list(1,2:6,7:11,12:16,17:p[1]),
            setting2 = list(1,2:6,7:11,12:16,17:p[2]),
            setting3 = list(1,2:6,7:11,12:16,17:p[3]))
cols <- gg_color_hue(4, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic", p = p, y_lim=y_lim)
