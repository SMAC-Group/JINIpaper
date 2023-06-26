
source(".simu/R/figures.R")

# -----------------
# Consistent case
# -----------------

# load(".simu/data/roblogistic_2_n_200_p_20.rds")
load(".simu/data/roblogistic_2_setting_2.rds")
sum(apply(res$br,2,var)) / sum(apply(res$jimi,2,var))
n_setting <- 1
n_methods <- length(res$setting$method)-1
p <- as.integer(res$setting$p) +1
MC <- as.integer(res$setting$MC)
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)))
super_res[[1]][,,1] <- res$jimi
super_res[[1]][,,2] <- res$robCR
super_res[[1]][,,3] <- res$robBY
super_res[[1]][,,4] <- res$mle
super_res[[1]][,,5] <- res$br
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","CR","BY","MLE","BR")
lab_main <- paste0("$n=",res$setting$n,"$, $p=",res$setting$p,"$")
theta0 <- c(.3, -2, -4, 0)
ind <- list(setting1 = list(1,2,3,4:p))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic", p = p, y_lim=y_lim)

# -----------------
# Consistent case
# -----------------

# load(".simu/data/roblogistic_2_n_200_p_20.rds")
load(".simu/data/roblogistic_1_setting_2.rds")
sum(apply(res$br,2,var)) / sum(apply(res$jimi,2,var))
n_setting <- 1
n_methods <- length(res$setting$method)-1
p <- as.integer(res$setting$p) +1
MC <- as.integer(res$setting$MC)
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)))
super_res[[1]][,,1] <- res$jimi
super_res[[1]][,,2] <- res$robCR
super_res[[1]][,,3] <- res$robBY
super_res[[1]][,,4] <- res$mle
super_res[[1]][,,5] <- res$br
lab_beta <- c("$\\beta_0$","$\\beta_{1}$","$\\beta_{2}$",
              "$\\beta_{4:p}$")
lab_methods <- c("JIMI","CR","BY","MLE","BR")
lab_main <- paste0("$n=",res$setting$n,"$, $p=",res$setting$p,"$")
theta0 <- c(.3, -2, -4, 0)
ind <- list(setting1 = list(1,2,3,4:p))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "roblogistic_outliers", p = p, y_lim=y_lim)
