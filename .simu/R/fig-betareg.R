# +---------+------+-----+-------+
# | Setting |   n  |  p  | p/n   |
# +---------+------+-----+-------+
# |       1 | 2000 | 300 | 0.150 |
# |       2 | 4000 | 400 | 0.100 |
# |       3 | 8000 | 500 | 0.063 |
# +---------+------+-----+-------+

source(".simu/R/figures.R")

# -----------------
# Inconsistent case
# -----------------
n_setting <- 3
n_methods <- 2
p <- 3:5 * 100+2
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/betareg_1_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$consistent
  super_res[[i]][,,2] <- res$jini
}
lab_beta <- c("$\\beta_{\\scriptscriptstyle{0}}$",
              "$\\beta_{\\scriptscriptstyle 1:5}$",
              "$\\beta_{\\scriptscriptstyle 6:10}$",
              "$\\beta_{\\scriptscriptstyle 11:15}$",
              "$\\beta_{\\scriptscriptstyle 16:p}$",
              "$\\phi$")
lab_methods <- c("MLE","JINI")
lab_main <- c("$n=2000$", "$n=4000$", "$n=8000$")
theta0 <- c(-.5, 1, -1.5, 2, 0, 5)
ind <- list(setting1 = list(1,2:6,7:11,12:16,17:(p[1]-1),p[1]),
            setting2 = list(1,2:6,7:11,12:16,17:(p[2]-1),p[2]),
            setting3 = list(1,2:6,7:11,12:16,17:(p[3]-1),p[3]))
cols <- gg_color_hue(3, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-1,1)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "betareg_inconsistent", p = p, y_lim=y_lim)

# -----------------
# Consistent case
# -----------------

n_setting <- 3
n_methods <- 2
p <- 3:5 * 100+2
MC <- 1e3
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)),
                  setting2 = array(dim=c(MC,p[2],n_methods)),
                  setting3 = array(dim=c(MC,p[3],n_methods)))
for(i in 1:n_setting){
  load(paste0(".simu/data/betareg_2_setting_",i,".rds"))
  super_res[[i]][,,1] <- res$mle
  super_res[[i]][,,2] <- res$jini
}
lab_beta <- c("$\\beta_{\\scriptscriptstyle{0}}$",
              "$\\beta_{\\scriptscriptstyle 1:5}$",
              "$\\beta_{\\scriptscriptstyle 6:10}$",
              "$\\beta_{\\scriptscriptstyle 11:15}$",
              "$\\beta_{\\scriptscriptstyle 16:p}$",
              "$\\phi$")
lab_methods <- c("MLE","JINI")
lab_main <- c("$n=2000$", "$n=4000$", "$n=8000$")
theta0 <- c(-.5, 1, -1.5, 2, 0, 5)
ind <- list(setting1 = list(1,2:6,7:11,12:16,17:(p[1]-1),p[1]),
            setting2 = list(1,2:6,7:11,12:16,17:(p[2]-1),p[2]),
            setting3 = list(1,2:6,7:11,12:16,17:(p[3]-1),p[3]))
cols <- gg_color_hue(3, alpha = 1)
pchs <- c(15:16,15:16)
y_lim <- c(-1,1)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "betareg_consistent", p = p, y_lim=y_lim)
