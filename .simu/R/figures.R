library(tikzDevice)
my_rmse <- function(x,theta){
  stat <- boxplot.stats(x)
  sqrt(var(x[!x%in%stat$out],na.rm=T) + mean(x[!x%in%stat$out]-theta,na.rm=T)^2)
  # sqrt(var(x,na.rm=T) + mean(x-theta,na.rm=T)^2)
}
my_bias <- function(x,theta){
  stat <- boxplot.stats(x)
  abs(mean(x[!x%in%stat$out],na.rm=T)-theta)
  # abs(mean(x,na.rm=T)-theta)
  # abs(median(x,na.rm=T)-theta)
}
my_average <- function(x){
  # mean(x,na.rm=T)
  stat <- boxplot.stats(x)
  mean(x[!x%in%stat$out],na.rm=T)
} 

gg_color_hue <- function(n, alpha) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100, alpha = alpha)[1:n]
}

summaries <- function(results,index,theta){
  dims <- dim(results)
  res <- list(bias=numeric(),rmse=numeric(),average=numeric())
  bias <- rmse <- average <- matrix(nrow=length(theta),ncol=dims[3])
  
  for(i in seq_along(theta)){
    for(j in seq_len(dims[3])){
      x <- as.vector(results[,index[[i]],j])
      bias[i,j] <- my_bias(x,theta[i])
      average[i,j] <- my_average(x)
      rmse[i,j] <- my_rmse(x,theta[i])
    }
  }
  
  res$bias <- bias
  res$rmse <- rmse
  res$average <- average
  
  return(res)
}

allplot <- function(results,index,theta,lab_par,lab_mth,lab_set,cols,pchs,path=getwd(),extension=NULL,p,y_lim){
  
  n_settings <- length(results)
  n_theta <- length(theta)
  n_methods <- dim(results[[1]])[3]
  
  if(length(lab_par) != n_theta) stop("length of 'theta' and 'lab_par' should match")
  
  pkgs <- c("\\usepackage{tikz}",
            "\\usepackage[active,tightpage,psfixbb]{preview}",
            "\\PreviewEnvironment{pgfpicture}",
            "\\setlength\\PreviewBorder{0pt}",
            "\\usepackage{amssymb}",
            "\\usepackage{bm}",
            "\\usepackage{amsthm}",
            "\\usepackage{amsbsy}",
            "\\usepackage{amsfonts}")
  
  # range_bias <- extendrange(res$bias, f = .05)
  # range_rmse <- extendrange(res$rmse, f = .05)
  
  dir.create(file.path(path,"figures"), showWarnings = FALSE)
  
  if(n_settings==1){
    tikz(paste0(path,"/figures/plot.tex"), 
         width = 6.5, 
         height = 3, 
         standAlone = TRUE, 
         packages = pkgs)
    
    par(mfrow = c(1,3), 
        mai = c(.05,.2,.05,.2), 
        oma=c(4.5,4,2,1))
  }
  
  if(n_settings==2){
    tikz(paste0(path,"/figures/plot.tex"), 
         width = 6.5, 
         height = 5.5, 
         standAlone = TRUE, 
         packages = pkgs)
    
    par(mfrow = c(3,2), 
        mai = c(.05,.2,.05,.2), 
        oma=c(3,4,2,1))
  }
  
  if(n_settings==3){
    tikz(paste0(path,"/figures/plot.tex"), 
         width = 8, 
         height = 6.67, 
         standAlone = TRUE, 
         packages = pkgs)
    
    par(mfrow = c(3,3),
        mai = c(0.1, 0.1, 0.1, 0.3)/1.6, 
        oma = c(4, 4.6, 2, 0.2))
  }
  
  if(n_settings==4){
    tikz(paste0(path,"/figures/plot.tex"), 
         width = 10, 
         height = 6.67, 
         standAlone = TRUE, 
         packages = pkgs)
    
    par(mfrow = c(3,4),
        mai = c(0.1, 0.1, 0.1, 0.3)/1.6, 
        oma = c(4, 4.6, 2, 0.2))
  }
  
  if(n_settings==5){
    tikz(paste0(path,"/figures/plot.tex"), 
         width = 10, 
         height = 6.67, 
         standAlone = TRUE, 
         packages = pkgs)
    
    par(mfrow = c(3,5),
        mai = c(0.1, 0.1, 0.1, 0.3)/1.6, 
        oma = c(4, 4.6, 2, 0.2))
  }
  
  if(n_settings>5){
    stop("ask Yuming")
  }
  
  # Boxplot
  n_bxp <- n_methods * n_theta + n_theta - 1
  for(k in seq_len(n_settings)){
    # compute summary statistics
    # summ <- summaries(results[[k]],index[[k]],theta)
    
    plot(NA, xlim = c(0.7,n_bxp+.3), ylim = y_lim, type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
    if(k==1){
      mtext("Centered Estimates", side = 2, line = 3.2, cex = 1.35)
      axis(2, cex.axis=1.6)
    }
    
    grid(nx=NA,ny=NULL)
    for(j in seq_len(n_methods)){
      for(i in seq_len(n_theta)){
        x <- c(results[[k]][,index[[k]][[i]],j]) - theta[i]
        tmp <- boxplot(x, plot = F, outline = F)
        bxp(tmp, add = TRUE, at = seq(j,n_bxp,by=(n_methods+1))[i], show.names = F, outline = F, boxfill = cols[j], xaxt = "n", yaxt = "n", axes = F, lty=1)
        points(x=seq(j,n_bxp,by=(n_methods+1))[i],y=mean(x,na.rm=T),pch = 16, cex = 1.5)
      }
    }
    abline(h = 0, lwd = 4)
    mtext(lab_set[k], side = 3, line = 0.5, cex = 1.6)
  }
  
  # Absolute bias
  for(k in seq_len(n_settings)){
    # compute summary statistics
    summ <- summaries(results[[k]],index[[k]],theta)
    if(k==1) {
      range_bias <- extendrange(summ$bias, f = .05)
      range_bias[1] <- 0.0
      # range_bias[2] <- 0.3
    }
    
    plot(NA, xlim = c(0.7,n_theta+.3), ylim = range_bias, type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
    if(k==1){
      mtext("Absolute Bias", side = 2, line = 3.2, cex = 1.35)
      # axis(2, at = seq(0,1,0.25), cex.axis=1.6)
      axis(2, las = 2, cex.axis = 1.6)
    }
    # grid(nx=NA,ny=NULL)
    grid(col = "lightgray", lty = "dotted")
    for(j in seq_len(n_methods)){
      lines(1:n_theta,summ$bias[,j],col=cols[j],pch=pchs[j],type="b",cex=2.8,lwd=2)
    }
  }
  
  # RMSE
  for(k in seq_len(n_settings)){
    # compute summary statistics
    summ <- summaries(results[[k]],index[[k]],theta)
    
    if(k==1){
      range_rmse <- extendrange(summ$rmse, f = .05)
      range_rmse[1] <- 0.0
      # range_rmse[2] <- 0.3
    }
    
    plot(NA, xlim = c(0.7,n_theta+.3), ylim = range_rmse, type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
    if(k==1){
      mtext("RMSE", side = 2, line = 3.2, cex = 1.35)
      # axis(2, at = seq(0,4,1), cex.axis=1.6)
      axis(2, las = 2, cex.axis = 1.6)
    }
    # grid(nx=NA,ny=NULL)
    grid(col = "lightgray", lty = "dotted")
    for(j in seq_len(n_methods)){
      lines(1:n_theta,summ$rmse[,j],col=cols[j],pch=pchs[j],type="b",cex=2.8,lwd=2)
    }
    axis(1, at=seq(1, n_theta, by=1), labels = FALSE)
    labs <- gsub("p}",paste0((p[k]-1),"}"),lab_par)
    mtext(text = labs, side = 1, at = seq(1, n_theta, by=1), line = 1, cex = 1.2)
    
    if(k==n_settings){
      legend("topleft",
             legend = lab_mth,
             bty = "n", 
             pt.cex = 2, 
             col = cols, 
             pch = pchs,
             pt.bg = cols, cex = 1.5, lwd = 1, lty = 1)
    }
    
  }
  
  
  dev.off()
  
  current_path <- getwd()
  setwd(paste0(path,"/figures/"))
  chain_cmd <- c("latex plot.tex",
                 "dvipdf plot.dvi",
                 "rm *.aux",
                 "rm *.log",
                 "rm *.dvi",
                 "rm *.tex")
  mapply(FUN=function(x)system(x),chain_cmd)
  if(!is.null(extension)) system(paste0("mv plot.pdf ",extension,".pdf"))
  setwd(current_path)
}


betaregplot <- function(results,index,theta,lab_par,lab_mth,lab_set,cols,pchs,path=getwd(),extension=NULL,p){
  
  n_settings <- length(results)
  n_theta <- lapply(theta, length)
  n_methods <- dim(results[[1]])[3]
  
  if(length(lab_par) != length(n_theta)) stop("length of 'theta' and 'lab_par' should match")
  
  pkgs <- c("\\usepackage{tikz}",
            "\\usepackage[active,tightpage,psfixbb]{preview}",
            "\\PreviewEnvironment{pgfpicture}",
            "\\setlength\\PreviewBorder{0pt}",
            "\\usepackage{amssymb}",
            "\\usepackage{bm}",
            "\\usepackage{amsthm}",
            "\\usepackage{amsbsy}",
            "\\usepackage{amsfonts}")
  
  # range_bias <- extendrange(res$bias, f = .05)
  # range_rmse <- extendrange(res$rmse, f = .05)
  
  dir.create(file.path(path,"figures"), showWarnings = FALSE)
  
  for(l in 1:length(n_theta)){
    if(n_settings==1){
      tikz(paste0(path,"/figures/plot.tex"), 
           width = 6.5, 
           height = 3, 
           standAlone = TRUE, 
           packages = pkgs)
      
      par(mfrow = c(1,3), 
          mai = c(.05,.2,.05,.2), 
          oma=c(4.5,4,2,1))
    }
    
    if(n_settings==2){
      tikz(paste0(path,"/figures/plot.tex"), 
           width = 6.5, 
           height = 5.5, 
           standAlone = TRUE, 
           packages = pkgs)
      
      par(mfrow = c(3,2), 
          mai = c(.05,.2,.05,.2), 
          oma=c(3,4,2,1))
    }
    
    if(n_settings==3){
      tikz(paste0(path,"/figures/plot.tex"), 
           width = 8, 
           height = 6.67, 
           standAlone = TRUE, 
           packages = pkgs)
      
      par(mfrow = c(3,3),
          mai = c(0.1, 0.1, 0.1, 0.3)/1.6, 
          oma = c(4, 4.6, 2, 0.2))
    }
    
    if(n_settings==4){
      tikz(paste0(path,"/figures/plot.tex"), 
           width = 10, 
           height = 6.67, 
           standAlone = TRUE, 
           packages = pkgs)
      
      par(mfrow = c(3,4),
          mai = c(0.1, 0.1, 0.1, 0.3)/1.6, 
          oma = c(4, 4.6, 2, 0.2))
    }
    
    if(n_settings>4){
      stop("ask Yuming")
    }
    
    # Boxplot
    n_bxp <- n_methods * n_theta[[l]] + n_theta[[l]] - 1
    for(k in seq_len(n_settings)){
      # compute summary statistics
      # summ <- summaries(results[[k]],index[[k]],theta)
      
      plot(NA, xlim = c(0.7,n_bxp+.3), ylim = c(-1, 1), type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
      if(k==1){
        mtext("Centered Estimates", side = 2, line = 3.2, cex = 1.35)
        axis(2, cex.axis=1.6)
      }
      
      grid(nx=NA,ny=NULL)
      for(j in seq_len(n_methods)){
        for(i in seq_len(n_theta[[l]])){
          x <- c(results[[k]][,index[[l]][[k]][[i]],j]) - theta[[l]][i]
          tmp <- boxplot(x, plot = F, outline = F)
          bxp(tmp, add = TRUE, at = seq(j,n_bxp,by=(n_methods+1))[i], show.names = F, outline = F, boxfill = cols[j], xaxt = "n", yaxt = "n", axes = F, lty=1)
          points(x=seq(j,n_bxp,by=(n_methods+1))[i],y=mean(x,na.rm=T),pch = 16, cex = 1.5)
        }
      }
      abline(h = 0, lwd = 4)
      mtext(lab_set[k], side = 3, line = 0.5, cex = 1.6)
    }
    
    if(l==1){
      summ <- summaries(results[[1]],index[[l]][[1]],theta[[l]])
      range_bias <- extendrange(summ$bias, f = .05)
      range_rmse <- extendrange(summ$rmse, f = .05)
    }
    if(l==2){
      summ <- summaries(results[[n_settings]],index[[l]][[n_settings]],theta[[l]])
      range_bias <- extendrange(summ$bias, f = .05)
      range_rmse <- extendrange(summ$rmse, f = .05)
    }
    # Absolute bias
    for(k in seq_len(n_settings)){
      # compute summary statistics
      # summ <- summaries(results[[k]],index[[l]][[k]],theta[[l]])
      if(k==1) range_bias <- extendrange(summ$bias, f = .05)
      
      plot(NA, xlim = c(0.7,n_theta[[l]]+.3), ylim = range_bias, type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
      if(k==1){
        mtext("Absolute Bias", side = 2, line = 3.2, cex = 1.35)
        # axis(2, at = seq(0,1,0.25), cex.axis=1.6)
        axis(2, las = 2, cex.axis = 1.6)
      }
      # grid(nx=NA,ny=NULL)
      grid(col = "lightgray", lty = "dotted")
      for(j in seq_len(n_methods)){
        lines(1:n_theta[[l]],summ$bias[,j],col=cols[j],pch=pchs[j],type="b",cex=2.8,lwd=2)
      }
    }
    
    # RMSE
    for(k in seq_len(n_settings)){
      # compute summary statistics
      summ <- summaries(results[[k]],index[[l]][[k]],theta[[l]])
      
      # if(k==1) range_rmse <- extendrange(summ$rmse, f = .05)
      
      plot(NA, xlim = c(0.7,n_theta[[l]]+.3), ylim = range_rmse, type="n", xaxt = "n", yaxt = "n", axes = T, xlab = " ", ylab = " ")
      if(k==1){
        mtext("RMSE", side = 2, line = 3.2, cex = 1.35)
        # axis(2, at = seq(0,4,1), cex.axis=1.6)
        axis(2, las = 2, cex.axis = 1.6)
      }
      # grid(nx=NA,ny=NULL)
      grid(col = "lightgray", lty = "dotted")
      for(j in seq_len(n_methods)){
        lines(1:n_theta[[l]],summ$rmse[,j],col=cols[j],pch=pchs[j],type="b",cex=2.8,lwd=2)
      }
      axis(1, at=seq(1, n_theta[[l]], by=1), labels = FALSE)
      labs <- gsub("p}",paste0((p[[l]][k]-1),"}"),lab_par[[l]])
      mtext(text = labs, side = 1, at = seq(1, n_theta[[l]], by=1), line = 1, cex = 1.2)
      
      if(k==n_settings){
        legend("topleft",
               legend = lab_mth,
               bty = "n", 
               pt.cex = 2, 
               col = cols, 
               pch = pchs,
               pt.bg = cols, cex = 1.5, lwd = 1, lty = 1)
      }
      
    }
    
    
    dev.off()
    
    current_path <- getwd()
    setwd(paste0(path,"/figures/"))
    chain_cmd <- c("latex plot.tex",
                   "dvipdf plot.dvi",
                   "rm *.aux",
                   "rm *.log",
                   "rm *.dvi",
                   "rm *.tex")
    mapply(FUN=function(x)system(x),chain_cmd)
    if(!is.null(extension)) system(paste0("mv plot.pdf ",extension,"_",l,".pdf"))
    setwd(current_path)
  }
}