#' @title Bootstrapped VBGF estimates for the alba lfq dataset
#'
#' @description Bootstrapped VBGF estimates for the alba length frequency data
#' set as estimated by \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#' @format ## `alba_boot`
#' A \code{lfqBoot} object with two levels: \code{$bootRaw} and \code{$seed}.
#' \describe{
#'   \item{\code{$bootRaw}}{A \code{data.frame} of fitted VBGF parameters
#'   (columns) by resampling (rows).}
#'   \item{\code{$seed}}{A \code{numeric} vector of seed values set prior to each
#'   resampling call to \link[fishboot]{lfqResample}.}
#' }
"alba_boot"
