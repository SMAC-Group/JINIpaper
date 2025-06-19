source(".simu/R/figures.R")

# -----------------
# Consistent case
# -----------------
load(".simu/data/pareto_2_setting_1.rds")
sum(apply(res$jimi_mle,2,sd,na.rm=T)) / sum(apply(res$jimi_rob,2,sd,na.rm=T))
n_setting <- 1
n_methods <- 2#length(res$setting$method)
method_names <- gsub(" ==.*","",res$setting$method)[c(1,4)]
p <- as.integer(res$setting$p) +2
MC <- as.integer(res$setting$MC)
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)))
for(i in seq_len(n_methods)) super_res[[1]][,,i] <- res[[method_names[i]]]
lab_beta <- c(paste0("$\\beta_",0:4,"$"))
lab_methods <- gsub("_"," ",toupper(method_names))
lab_main <- paste0("$n=",res$setting$n,"$, $p=",res$setting$p,"$")
beta <- str2expression(res$setting$beta)
eval(beta)
theta0 <- c(beta[1:5])
ind <- list(setting1 = list(1,2,3,4,5))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "pareto", p = p, y_lim=y_lim)

 # -----------------
# Consistent case
# -----------------

load(".simu/data/pareto_1_setting_1.rds")
n_setting <- 1
n_methods <- 2#length(res$setting$method)
method_names <- gsub(" ==.*","",res$setting$method)[c(1,4)]
p <- as.integer(res$setting$p) +2
MC <- as.integer(res$setting$MC)
super_res <- list(setting1 = array(dim=c(MC,p[1],n_methods)))
for(i in seq_len(n_methods)) super_res[[1]][,,i] <- res[[method_names[i]]]
lab_beta <- c(paste0("$\\beta_",0:4,"$"))
lab_methods <- gsub("_"," ",toupper(method_names))
lab_main <- paste0("$n=",res$setting$n,"$, $p=",res$setting$p,"$")
beta <- str2expression(res$setting$beta)
eval(beta)
theta0 <- c(beta[1:5])
ind <- list(setting1 = list(1,2,3,4,5))
cols <- gg_color_hue(n_methods, alpha = 1)
pchs <- c(15:(15+n_methods-1))
y_lim <- c(-4,4)
allplot(results = super_res, index = ind, theta = theta0, lab_par = lab_beta, 
        lab_mth = lab_methods, lab_set = lab_main, cols = cols, pchs = pchs, 
        path = paste0(getwd(),"/.simu/"), extension = "pareto_outliers", p = p, y_lim=y_lim)
