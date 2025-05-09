#' @title Bootstraped ELEFAN_SA
#'
#' @description This function performs a bootstrapped fitting of a von
#' Bertalanffy growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_SA}}
#' function. Most of the arguments are simply passed to the function within many
#' permutations (resampling) of the original \code{lfq} data. Partial
#' (repeated fitting on original data) and full bootstrap (with resampling)
#' routines are possible, depending on \code{resample}.
#'
#'
#' @param lfq a length frequency object of the class \code{lfq} (see
#' \link[TropFishR]{lfqCreate}).
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#' growth function should be applied (default: FALSE).
#' @param init_par a list providing the Initial values for the components to be
#' optimized. When set to NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is the maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.5),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corresponds to peak spawning month (range: 0 to 1,
#'   default: 0.5),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param low_par a list providing the lower bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corresponds to peak spawning month (range: 0 to 1,
#'   default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the upper bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corresponds to peak spawning month (range: 0 to 1,
#'   default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param SA_time numeric; Maximum running time in seconds (default : 60 * 1).
#' @param maxit Integer. Maximum number of iterations of the algorithm. Default
#' is NULL.
#' @param nb.stop.improvement Integer. The program will stop when there is no
#' any improvement in 'nb.stop.improvement' steps. Default is NULL
#' @param SA_temp numeric; Initial value for temperature (default : 1e5).
#' @param MA number indicating over how many length classes the moving average
#' should be performed (defalut: 5, for more information see
#' \link[TropFishR]{lfqRestructure}).
#' @param addl.sqrt Passed to \link[TropFishR]{lfqRestructure}. Applied an
#' additional square-root transformation of positive values according to Brey
#' et al. (1988). (default: FALSE, for more information see
#' \link[TropFishR]{lfqRestructure}).
#' @param agemax maximum age of species; default NULL, then estimated from
#' \eqn{L_{inf}}.
#' @param seed seed value for random number reproducibility.
#' @param nresamp \code{numeric}, the number of permutations to run (by default
#' \code{nresamp = 10}).
#' @param resample \code{logical}. Do you want that \code{lfq} object be
#' resampled (\code{TRUE} by default).
#' @param parallel Whether a \code{logical} or \code{integer} argument
#' specifying the configuration of parallel computing. If \code{parallel = TRUE},
#' the number of threads will be defined as \code{parallel::detectCores() - 2}.
#' @param outfile \code{character}; path of the file which will register the
#' progress of the permutation completions. If it is set as \code{false},
#' \code{NA} or \code{NULL}, no file will be created.
#'
#' @details
#' If \code{resample = TRUE}, a \strong{full non-parametric bootstrap} is
#' performed with resampling from the original length-frequencies by using the
#' function \link{lfqResample}. Otherwise, if \code{resample = FALSE}, a partial
#' bootstrap is performed, reflecting solution variation due only to the search
#' algorithm, with repeated fitting to the original data (no resampling is
#' performed).
#'
#' @references \itemize{
#'  \item Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency
#' analysis: a revised and expanded user's guide to ELEFAN 0, 1 and 2.
#'  \item Efron, B., & Tibshirani, R. (1986). Bootstrap methods for standard
#' errors, confidence intervals, and other measures of statistical accuracy.
#' Statistical Science, 54-75.
#'  \item Mildenberger, T., Taylor, M. H., & Wolff, A. M., 2017. TropFishR: an R
#' package for fisheries analysis with length-frequency data. Methods in Ecology
#' and Evolution, 8(11), 1520-1527.
#'  \item Pauly, D. 1981. The relationship between gill surface area and growth
#' performance in fish: a generalization of von Bertalanffy's theory of growth.
#' Meeresforsch. 28:205-211.
#'  \item Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the
#' objective extraction of growth parameters from length-frequency data.
#' Meeresforschung, 28(4):205-211.
#'  \item Schwamborn, R., Mildenberger, T. K., & Taylor, M. H., 2019. Assessing
#' sources of uncertainty in length-based estimates of body growth in
#' populations of fishes and macroinvertebrates with bootstrapped ELEFAN.
#' Ecological Modelling, 393, 37-51.
#'  \item Schwamborn, R., Freitas, M. O., Moura, R. L., & Aschenbrenner, A. 2023.
#' Comparing the accuracy and precision of novel bootstrapped length-frequency
#' and length-at-age (otolith) analyses, with a case study of lane snapper
#' (\emph{Lutjanus synagris}) in the SW Atlantic. Fisheries Research, 264, 106735.
#'  \item Scrucca, L., 2013. GA: a package for genetic algorithms in R. Journal
#' of Statistical Software, 53(4), 1-37.
#'  \item von Bertalanffy, L., 1938. A quantitative theory of organic growth.
#' Human Biology 10, 181-213.
#' }
#'
#' @return An object of class \code{lfqBoot} containing 2 levels:
#' \describe{
#'   \item{\code{$bootRaw}}{A \code{data.frame} of fitted VBGF parameters
#'   (columns) by resampling (rows).}
#'   \item{\code{$seed}}{A \code{numeric} vector of seed values set prior to each
#'   resampling call to \link[fishboot]{lfqResample}.}
#' }
#'
#' @export
#'
#' @examples
#' # load data
#' data("alba", package = "TropFishR")
#'
#' # Define settings (for demo only, fast settings)
#' MA       <- 7
#' low_par   <- list(Linf = 9, K = 0.4, t_anchor = 0.5, C = 0, ts = 0)
#' up_par    <- list(Linf = 11, K = 0.6, t_anchor = 0.8, C = 1, ts = 1)
#' init_par <- list(Linf = 10, K = 0.5, t_anchor = 0.6, C = 0.5, ts = 0.5)
#' SA_time  <- 1.5
#' SA_temp  <- 1e5
#' nresamp  <- 2
#'
#' # Non-parallel bootstrapped curve fiting
#' res <- ELEFAN_SA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       init_par = init_par, up_par = up_par, low_par = low_par,
#'                       SA_time = SA_time, SA_temp = SA_temp,
#'                       nresamp = nresamp)
#'
#' res
#'
#' \donttest{
#' # Define settings (for demo only)
#' MA       <- 7
#' low_par  <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' init_par <- list(Linf = 12, K = 0.5, t_anchor = 0.5, C = 0.5, ts = 0.5)
#' up_par   <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' SA_time  <- 10
#' SA_temp  <- 1e5
#' nresamp  <- 12
#'
#' # parallel version
#' res <- ELEFAN_SA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       init_par = init_par, up_par = up_par, low_par = low_par,
#'                       SA_time = SA_time, SA_temp = SA_temp,
#'                       nresamp = nresamp,
#'                       parallel = TRUE)
#'
#' res
#' }
ELEFAN_SA_boot <- function(lfq,
                           seasonalised = FALSE,
                           init_par = list(Linf = 50,
                                           K = 0.5,
                                           t_anchor = 0.5,
                                           C = 0, ts = 0),
                           low_par = NULL, up_par = NULL,
                           SA_time = 60 * 1, maxit = NULL,
                           nb.stop.improvement = NULL,
                           SA_temp = 1e+05, MA = 5, addl.sqrt = FALSE,
                           agemax = NULL,
                           seed = NULL, nresamp = 10, resample = TRUE,
                           parallel = FALSE, outfile = NA){

  if(is.null(seed)) seed <- as.numeric(Sys.time())

  if(!is.logical(parallel) && !is.integer(parallel)) stop("'parallel' must be logical or integer.")

  no_cores <- if(isTRUE(parallel)) parallel::detectCores() - 2 else pmax(as.integer(parallel), 1)

  if(no_cores > 1){

    # Registering cluster
    cl <- makeCluster(spec = no_cores)
    registerDoParallel(cl = cl)

    # Run multithread process
    res <- foreach(x = seq(nresamp), .inorder = FALSE, .packages = "TropFishR") %dopar% {
      lfq_ELEFAN_SA(lfq = lfq, x = x, resample = resample, seed = seed,
                    seasonalised = seasonalised,
                    init_par = init_par,
                    low_par = low_par, up_par = up_par,
                    SA_time = SA_time, maxit = maxit,
                    nb.stop.improvement = nb.stop.improvement,
                    SA_temp = SA_temp, MA = MA, addl.sqrt = addl.sqrt,
                    agemax = agemax)
    }

    # Finish cluster
    stopCluster(cl)
  }else{

    # Empty results list
    res <- vector(mode = "list", length = nresamp)

    for(x in seq(nresamp)){
      res[[x]] <- lfq_ELEFAN_SA(lfq = lfq, x = x,
                                resample = resample, seed = seed,
                                seasonalised = seasonalised,
                                init_par = init_par,
                                low_par = low_par, up_par = up_par,
                                SA_time = SA_time, maxit = maxit,
                                nb.stop.improvement = nb.stop.improvement,
                                SA_temp = SA_temp, MA = MA,
                                addl.sqrt = addl.sqrt, agemax = agemax)
    }
  }

  # Combine results into a data.frame
  res <- do.call(what = rbind, args = res) |> as.data.frame()

  # Build output structure
  res <- list(bootRaw = res[,-ncol(res)],
              seed = res$seed)

  # Set class of output object
  class(res) <- "lfqBoot"

  res
}


# Engine function
lfq_ELEFAN_SA <- function(lfq, x, resample, seed,
                          seasonalised, init_par, low_par, up_par, SA_time,
                          maxit, nb.stop.improvement, SA_temp, MA, addl.sqrt,
                          agemax){

  seed <- seed + x
  set.seed(seed)

  # resample data
  if(resample){
    lfqb <- lfqResample(lfq = lfq)
  }else{
    lfqb <- lfq
  }

  fitboot <- ELEFAN_SA(lfq = lfqb,
                       seasonalised = seasonalised,
                       init_par = init_par,
                       low_par = low_par, up_par = up_par,
                       SA_time = SA_time, maxit = maxit,
                       nb.stop.improvement = nb.stop.improvement,
                       SA_temp = SA_temp, MA = MA, addl.sqrt = addl.sqrt,
                       agemax = agemax,
                       flagging.out = FALSE, plot = FALSE, plot.score = FALSE,
                       verbose = FALSE)

  # return result
  c(unlist(fitboot$par), seed = seed)
}
