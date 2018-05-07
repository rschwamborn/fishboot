#' @name alba_boot
#'
#' @title Bootstrapped VBGF estimates for the alba lfq dataset
#'
#' @description Bootstrapped VBGF estimates for the alba lfq dataset as estimated by
#'   \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#'
#' @docType data
#'
#' @format A list of class `lfqBoot`
#'
#'
#' @usage data(alba_boot)
#' @keywords data dataset length-frequency bootstrap
#'
#' @examples
#'
#' \donttest{
#' #### For documentation only - How data was produced ####
#' library(parallel)
#' library(TropFishR)
#' library(ks)
#'
#' data(alba)
#'
#' # ELEFAN_GA_boot settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' seasonalised <- FALSE
#' popSize <- 100
#' maxiter <- 50
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 200
#'
#' # Bootstrapped estimates
#' alba_boot <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = seasonalised,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
#'   seed = 20
#' )
#' }
#'
#'
#' data(alba_boot)
#'
#' head(alba_boot$bootRaw)
#' head(alba_boot$seed)
#'
#' # plot
#' univariate_density(alba_boot)
#' LinfK_scatterhist(alba_boot, phi.contour = TRUE)
#'
#'
NULL
