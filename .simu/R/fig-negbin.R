# +---------+-----+----+-------+
# | Setting |  n  |  p | p/n   |
# +---------+-----+----+-------+
# |       1 | 200 | 40 | 0.200 |
# |       2 | 400 | 50 | 0.125 |
# |       3 | 800 | 60 | 0.075 |
# +---------+-----+----+-------+

source(".simu/R/figures.R")

# -----------------
# Inconsistent case
# -----------------
n_setting <- 3
n_methods <- 2
p <- 4:6 * 10+2
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/negbin_1_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$consistent
  super_res[[i]][,,2] <- res$jini
}
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$", "$\\alpha$")
lab_methods <- c("MLE","JINI")
lab_main <- c("$n=200$", "$n=400$", "$n=800$")
theta0 <- c(2, 1, -1, 0, .7)
ind <- list(setting1 = list(1,2,3,4:(p[1]-1),p[1]),
            setting2 = list(1,2,3,4:(p[2]-1),p[2]),
            setting3 = list(1,2,3,4:(p[3]-1),p[3]))
cols <- gg_color_hue(3, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-.6,.6)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "negbin_inconsistent", p = p, y_lim=y_lim)

# -----------------
# Consistent case
# -----------------

n_setting <- 3
n_methods <- 2
p <- 4:6 * 10+2
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/negbin_2_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$mle
  super_res[[i]][,,2] <- res$jini
}
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$", "$\\alpha$")
lab_methods <- c("MLE","JINI")
lab_main <- c("$n=200$", "$n=400$", "$n=800$")
theta0 <- c(2, 1, -1, 0, .7)
ind <- list(setting1 = list(1,2,3,4:(p[1]-1),p[1]),
            setting2 = list(1,2,3,4:(p[2]-1),p[2]),
            setting3 = list(1,2,3,4:(p[3]-1),p[3]))
cols <- gg_color_hue(3, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-.6,.6)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "negbin_consistent", p = p, y_lim=y_lim)
