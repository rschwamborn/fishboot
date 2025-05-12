#' @title Standard bootstrapped length-converted catch curve
#' 
#' @param mids A vector of mid points for each size class
#' @param catch A vector of catches for each size class
#' @param K von Bertalanffy growth coefficient
#' @param Linf von Bertalanffy asymptotic length
#' @param t0 von Bertalanffy theoretical age at length zero
#' @param binsize Size class width
#' @param ex.points Number of points after highest point to exclude from the regression
#' @param plot Logical, if TRUE, a plot of the length-converted catch curve is produced
#' @param plot.yup.lim Upper limit for y-axis
#' @param plot.ylow.lim Lower limit for y-axis
#' @param plot.xlow.lim Lower limit for x-axis
#' @param plot.xup.lim Upper limit for x-axis
#' @param N_runs Number of bootstrap iterations
#'
#' @examples 
#' md = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200)
#' ct = c(100, 200, 300, 400, 350, 320, 260, 200, 160, 100, 80, 60, 50, 40, 30, 20, 10, 5, 2, 1)
#'
#'
#'standard_iLCCC(mids = md,
#'               catch = ct,
#'               K = 0.05,
#'               Linf = 200,
#'               t0 = 0,
#'               binsize = 10,
#'               ex.points = 1,
#'               plot = TRUE,
#'               plot.yup.lim = NULL,
#'               plot.ylow.lim = NULL,
#'               plot.xlow.lim = NULL,
#'               plot.xup.lim = NULL,
#'               N_runs = 10000)
#'@export
standard_iLCCC <- function(mids = NULL,
                  catch = NULL,
                  K = NULL,
                  Linf = NULL,
                  t0 = NULL,
                  binsize = NULL,
                  ex.points = 1,
                  plot = FALSE,
                  plot.yup.lim = NULL,
                  plot.ylow.lim = NULL,
                  plot.xlow.lim = NULL,
                  plot.xup.lim = NULL,
                  N_runs = 1000) {
  
  if (is.matrix(catch)) {
    stop(noquote("Catch must be arranged as a vector"))
  }
  
  if (length(mids) != length(catch)) {
    stop(noquote("Catch and midlengths must be vectors of same length"))
  }
  
  mids <- as.vector(mids)
  catch <- as.vector(catch)
  #upper and lower limits of each size class
  class.min <- mids - (binsize/2)
  class.max <- mids + (binsize/2)
  
  
  if (t0 <- FALSE) {
    t0 <- exp(-0.3922 - 0.2752*log(Linf) - 1.038*log(K))
  } else {
    t0 <- t0
  }
  
  rel.age <- ifelse(mids <= Linf, t0 - (1/K)*log(1 - mids/Linf), NA) # age at any given size class
  dt <- ifelse(mids <= Linf, (t0 - (1/K)*log(1 - class.max/Linf)) - (t0 - (1/K)*log(1 - class.min/Linf)), NA) # amount of time from one size class to another
  # if a midlength is larger than the asympotic size, it is simply excluded from subsequent analysis
  # as calculating dt with midlenths > Linf leads to NaN values
  
  lnNdt <- ifelse(is.na(dt), NA, log((catch+1)/dt)) # only calculate logged catch if the divisor is != NA
  
  df <- data.frame(lnNdt, rel.age)
  df <- na.exclude(df)
  df[df == 0] <- NA #exclude zero values
  
  yvar <- as.numeric(df$lnNdt)
  xvar <- df$rel.age
  
  selection <- c(which(yvar == max(yvar)) + 1,
                 which(xvar == max(xvar)) - ex.points) # select the data rows after the mode (highest point) to fit the line
  # if you suspect of undersampling, modify to c(which(yvar == max(yvar)) + 1,
  # which(xvar == max(xvar))-1) or -2 or more, this should exclude the largest size classes in the sample from the regression
  
  df.cc <- as.data.frame(cbind(xvar,yvar))
  df.selec.cc <- df.cc[selection[1]:selection[2],]# creates a data frame with only the selected rows
  
  df.selec.cc[is.finite(rowSums(df.selec.cc)),] # select only finite values
  df.selec.cc[is.numeric(rowSums(df.selec.cc)),] # exclude NaN measurements
  
  # some Infs and NaNs can pop up during the process when logarithms end up as imaginary or complex numbers
  
  df.selec.cc <- subset(df.selec.cc, df.selec.cc$yvar > 1)
  
  lm.cc <- lm(yvar ~ xvar,
              data = df.selec.cc)
  
  reg.output <- summary(lm.cc)
  intercept.reg <- reg.output$coefficients[1]
  slope.reg <- -reg.output$coefficients[2]
  se.slope.reg <- reg.output$coefficients[4]
  
  # bootstrap to obtain C.I.s for the total mortality
  
  boot.Z <- NULL # empty object to store the Z replicates
  
  for (i in 1:N_runs) {
    
    resamp.df = df.selec.cc[sample(1:nrow(df.selec.cc),
                                   nrow(df.selec.cc),
                                   replace = TRUE),]
    
    resamp.df <- resamp.df[apply(resamp.df, 1,
                                 function(row) all(row != 0)),]
    resamp.df[is.na(resamp.df) | resamp.df=="Inf"] <- NA
    resamp.df[is.na(resamp.df) | resamp.df=="NaN"] <- NA
    
    lm.boot <- lm(yvar ~ xvar,
                  data = resamp.df,
                  na.action = na.omit)
    
    boot.Z <- c(boot.Z, lm.boot$coefficients[2])
    
  }
  
  Z <- median(-boot.Z, na.rm = TRUE)
  Z.range <- range(-boot.Z)
  Z.CI <- quantile(-boot.Z, probs = c(0.05, 0.95), na.rm = TRUE)
  boot.Z <- as.vector(boot.Z)
  
  output <- list(intercept.reg,
                 Z, se.slope.reg,
                 Z.CI, Z.range,
                 as.vector(-boot.Z),
                 df.selec.cc)
  names(output) <- c("Intercept (a)", "Z",
                     "Standard Error for Z",
                     "Z bootstrapped 90% C.I.",
                     "Range",
                     "Z posterior distribution",
                     "Regression input")
  
  preds <- predict(lm.cc, pred.df = data.frame(x = df.selec.cc$xvar),
                   interval = "confidence")
  
  par(mfrow = c(1,2))
  
  if (plot == TRUE) {
    
    plot(yvar ~ xvar, xlab = "Relative age (years)",
         ylab = "ln(Catch)/dt",
         main = paste("Length-converted catch curve"),
         cex.main = 0.8, ylim = c(plot.ylow.lim, plot.yup.lim),
         xlim = c(plot.xlow.lim, plot.xup.lim))
    mtext(paste("Median Z =", round(Z,3),
                ", 95% C.I. =", round(Z.CI[1],3), "to", round(Z.CI[2],3)),
          side = 3, cex = 0.7)
    points(df.selec.cc$yvar ~ df.selec.cc$xvar, col = "black",
           pch = 20)
    abline(lm.cc)
    lines(df.selec.cc$xvar, preds[,3], lty = "dashed")
    lines(df.selec.cc$xvar, preds[,2], lty = "dashed")
    
    
    hist(-boot.Z, main = "Posterior distribution of Z",
         xlab = "Z", breaks = 40, col = 0, cex.main = 0.8)
    
  }
  
  return(output)
  
  }