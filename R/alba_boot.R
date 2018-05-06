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
