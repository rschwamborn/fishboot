#' @title Bootstrap-based methods for the study of fish stocks and aquatic populations
#' @description The fishboot package contains a suite of new bootstrap-based models
#' and software tools for the study of fish stocks and aquatic populations.
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach "%dopar%"
#' @importFrom grDevices colorRampPalette adjustcolor blues9 rgb
#' @importFrom graphics par layout image points box hist rect contour abline legend lines text mtext polygon segments
#' @importFrom stats runif cov quantile complete.cases
#' @importFrom utils modifyList stack
#' @importFrom ks kde Hpi
#' @importFrom TropFishR ELEFAN_SA ELEFAN_GA VBGF lfqRestructure growth_length_age
#' @importFrom fishmethods grotag
## usethis namespace: end
NULL
