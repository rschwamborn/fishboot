#' The interquantile range test
#'
#' @param 'data.A' a dataset (vector) of the class `numeric`.
#' @param 'data.B' a dataset (vector) of the class `numeric`.
#' @param n.boot numeric; the number of bootstrap runs for the Qanova test (Default: `n.boot = 500`)

#' @description The function tests for differences in 95% Inter-Quantile Range between two datasets     
#' Performs a non-parametric Harrell-Davis quantile test on standardized data. 
#' This test is intended to compare the precision of statistical procedures, 
#' ie. to test for diffences in 95% confidence intervals. 
#' The 95% quantile ranges (i.e., the 95% confidence interval) of two posteriors, obtained as an output from bootstrapping (e.g., obtained from from  ELEFAN_GA_boot.R) are compared. 
#' A simple non-parametric statistical test (the Harrell-Davis quantile test, see Wilcox, 2012) is conducted to verify whether there are significant differences in 95% quantile ranges of posteriors, i.e., to verify whether there are significant differences in 95% confidence intervals of the original parameter estimates. 
#' the function Qanova within the R package WRS2 (Mair & Wilcox, 2017) is used, and applied to standardized posteriors. Standardized posteriors are calculated as: standardized_data = raw_data - 0.025 quantile.
#' Small "p" values (p < 0.05) mean that 95% interquantile ranges are signifficantly different between two datasets.
#' 
#'@return a single result ("p-value") of class `character` 
#'
#'
#' @examples
#'
#'
#' # Test for differences in 95% interquantile range between pairs of data 
#' # (e.g., pairs of posterior distributions)
#' 
#' library(WRS2)
#' 
#' A <- rnorm (2000,mean = 1, sd = 5 ) # standard deviation = 5
#' B <- rnorm (2000,mean = 3000, sd = 5.02 ) # standard deviation =  basically the same as A
#' C <- rnorm (2000,mean = 1000, sd = 50 ) # much larger standard deviation
#' 
#' interquant_r.test(A,B, 500) # performs the interquant_r.test between "A" and "B", with 500 runs 
#'  # the large "p" value (p> 0.05) means that 95% interquantile ranges are not signifficantly different between datasets "A' and "B".
#' 
#' interquant_r.test(A,C, 500) # performs the interquant_r.test between "A" and "C", with 500 runs 
#' # the small "p" value (p < 0.05) means that 95% interquantile ranges are signifficantly different between datasets "A' and "C".
#' 
#' 
#' @export
#'
interquant_r.test <- function(dat.A, dat.B, n.boot = 500) {
  
  
  dat.A.stand <-   dat.A - quantile(dat.A, 0.025 )                             
  
  dat.B.stand <-   dat.B - quantile(dat.B, 0.025 )                             
  
  
  dat.frameAB <- data.frame(dat.A.stand, dat.B.stand)
  dat.frameAB.stack <-  stack(dat.frameAB)
  
  
  res <- WRS2::Qanova(values ~ ind, dat.frameAB.stack, q =  0.975, nboot = n.boot)
  
  res.p.value <-   res$p.value[1,2]
  
  
  if (res.p.value == 0)
    
    paste("p-value < 0.0001")
  
  else
    
    paste ("p-value = ", res.p.value )
  
}



