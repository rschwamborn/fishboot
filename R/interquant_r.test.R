#' @title The interquantile range test
#'
#' @description Pending...
#'
#' @param dat.A,dat.B \code{numeric} vectors to be compared.
#' @param n.boot \code{numeric}; the number of bootstrap runs for the Qanova
#' test (default: \code{n.boot = 500}).
#'
#' @description This function tests for differences in 95% Inter-Quantile Range
#' between two datasets. It performs a non-parametric Harrell-Davis quantile
#' test on standardized data. This test is intended to compare the precision of
#' statistical procedures, ie. to test for differences in 95% confidence
#' intervals. The 95% quantile ranges (i.e., the 95% confidence interval) of two
#' posteriors, obtained as an output from bootstrapping (e.g., obtained from
#' ELEFAN_GA_boot.R) are compared. A simple non-parametric statistical test (the
#' Harrell-Davis quantile test, see Wilcox, 2012) is conducted to verify whether
#' there are significant differences in 95% quantile ranges of posteriors, i.e.,
#' to verify whether there are significant differences in 95% confidence
#' intervals of the original parameter estimates. The function
#' \link[WRS2]{Qanova} within the R package WRS2 (Mair & Wilcox, 2017) is used,
#' and applied to standardized posteriors. Standardized posteriors are
#' calculated as: \deqn{standardized_data = raw_data - Q_{0.025}}.
#' Small "p" values (p < 0.05) means that 95% interquantile ranges are
#' significantly different between two datasets.
#'
#' @return A single \code{character} indicating the p-value result.
#'
#' @export
#'
#'
#' @examples
#' # Test for differences in 95% interquantile range between pairs of data
#' # (e.g., pairs of posterior distributions)
#'
#' A <- rnorm(2000,mean = 1, sd = 5) # standard deviation = 5
#' B <- rnorm(2000,mean = 3000, sd = 5.02) # a sd similar to A
#' C <- rnorm(2000,mean = 1000, sd = 50) # a far different sd
#'
#' # performs the interquant_r.test between "A" and "B", with 500 runs
#' interquant_r.test(dat.A = A, dat.B = B, n.boot = 500)
#'
#' # The large p-value (p> 0.05) means that 95% interquantile ranges are not
#' # significantly different between datasets "A' and "B".
#'
#' # performs the interquant_r.test between "A" and "C", with 500 runs
#' interquant_r.test(dat.A = A, dat.B = C, n.boot = 500)
#'
#' # The small p-value (p < 0.05) means that 95% interquantile ranges are
#' # significantly different between datasets "A' and "C".
interquant_r.test <- function(dat.A, dat.B, n.boot = 500){

  dat.A.stand <- dat.A - quantile(dat.A, 0.025 )
  dat.B.stand <- dat.B - quantile(dat.B, 0.025 )

  dat.frameAB <- data.frame(dat.A.stand, dat.B.stand)
  dat.frameAB.stack <- stack(dat.frameAB)

  res <- Qanova(formula = values ~ ind, data = dat.frameAB.stack,
                q = 0.975, nboot = n.boot)

  res.p.value <- res$p.value[1,2]

  if(res.p.value < 1e-6) paste("p-value < 1e-6") else sprintf("p-value = %.6f", res.p.value)
}
