#' @title Bootstrapped VBGF estimates for the \link[TropFishR]{alba} data set
#' @md
#'
#' @description Bootstrapped VBGF estimates for the alba length frequency data
#' set as estimated by \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#' @format ## `alba_boot`
#' A \code{lfqBoot} object with two levels:
#' \describe{
#'   \item{\code{$bootRaw}}{A \code{data.frame} of fitted VBGF parameters
#'   (columns) by resampling (rows).}
#'   \item{\code{$seed}}{A \code{numeric} vector of seed values set prior to each
#'   resampling call to \link[fishboot]{lfqResample}.}
#' }
"alba_boot"

#' @title Bootstrapped growth estimates for the bonito data set
#' @md
#'
#' @description Bootstrapped growth estimates for the \link[fishmethods]{bonito}
#' age-growth data calculated using the \link{grotag_boot} as this:
#'
#' ```{r eval=FALSE}
#' data("bonito", package = "fishmethods")
#'
#' grotag_boot(input.data = bonito,
#'             alpha = 35, beta = 55,
#'             design  = list(nu = 1, m = 1,p = 1, sea = 1),
#'             stvalue = list(sigma = 0.9, nu = 0.4, m = -1,
#'                            p = 0.2, u = 0.4, w = 0.4),
#'             upper   = list(sigma = 5, nu = 1, m = 2,
#'                            p = 0.5, u = 1, w = 1),
#'             lower   = list(sigma = 0, nu = 0, m = -2,
#'                            p = 0.0, u = 0, w = 0),
#'             control = list(maxit = 1e4),
#'             seed = 1234, nresamp = 100, cc.only = TRUE)
#' ```
#'
#' @format ## `bonito_boot`
#' A \code{grotagBoot} object with four columns: \code{Linf}, \code{K},
#' \code{PhiL}, \code{u}, \code{w} and \code{seed}.
"bonito_boot"
