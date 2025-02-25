#' @title Bootstraped ELEFAN_GA
#'
#' @description This function performs a bootstrapped fitting of a von
#' Bertalanffy growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_GA}}
#' function. Most of the arguments are simply passed to the function within many
#' permutations (resampling) of the original \code{lfq} data. As the original
#' function, \code{ELEFAN_GA} also conducts Electronic LEngth Frequency ANalysis
#' using a genetic algorithm (GA) to estimate growth parameters. Partial
#' (repeated fitting on original data) and full bootstrap (with resampling)
#' routines are possible, depending on \code{resample}.
#'
#'
#' @param lfq a length frequency object of the class \code{lfq} (see
#' \link[TropFishR]{lfqCreate}).
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#' growth function should be applied (default: FALSE).
#' @param low_par a list providing the minimum of the search space in case
#' of real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#' \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corresponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the maximum of the search space in case of
#' real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#' \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 1),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corresponds to peak spawning month (range: 0 to 1, default: 1),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 1),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 1);
#' }
#' @param popSize the population size. Default: 50
#' @param maxiter the maximum number of iterations to run before the
#' GA search is halted. default:100
#' @param run the number of consecutive generations without any improvement
#' in the best fitness value before the GA is stopped. Default: equals
#' \code{maxiter}.
#' @param parallel Whether a \code{logical} or \code{integer} argument
#' specifying the configuration of parallel computing. See \code{\link[GA]{ga}}
#' for details.
#' @param pmutation the probability of mutation in a parent chromosome.
#' Usually mutation occurs with a small probability, and by default is set to 0.1.
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value and by default is set to 0.8.
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for
#'    more information see \link[TropFishR]{lfqRestructure})
#' @param addl.sqrt additional square root transformation of positive values
#' according to Brey et al. (1988) (default: FALSE, for
#'    more information see \link[TropFishR]{lfqRestructure})
#' @param agemax maximum age of species; default \code{NULL}, then estimated
#' from \eqn{L_{inf}}.
#' @param ... additional parameters to pass to \code{\link[GA]{ga}}.
#' @param seed seed value for random number reproducibility.
#' @param nresamp \code{numeric}, the number of permutations to run (by default
#' \code{nresamp = 10}).
#' @param resample \code{logical}. Do you want that \code{lfq} object be
#' resampled (\code{TRUE} by default).
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
#'
#' @examples
#' # load data
#' data("alba", package = "TropFishR")
#'
#' # Define settings (for demo only, fast settings)
#' MA        <- 7
#' low_par   <- list(Linf = 9, K = 0.4, t_anchor = 0.5, C = 0, ts = 0)
#' up_par    <- list(Linf = 11, K = 0.6, t_anchor = 0.8, C = 1, ts = 1)
#' popSize   <- 12
#' maxiter   <- 5
#' run       <- 4
#' pmutation <- 0.1
#' nresamp   <-  2
#'
#' # Non-parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       up_par = up_par, low_par = low_par,
#'                       parallel = FALSE,
#'                       popSize = popSize, maxiter = maxiter,
#'                       run = run, pmutation = pmutation,
#'                       nresamp = nresamp)
#'
#' res
#'
#'
#' \donttest{
#' # Define settings (for demo only)
#' MA        <- 7
#' low_par   <- list(Linf = 8, K = 0.2, t_anchor = 0, C = 0, ts = 0)
#' up_par    <- list(Linf = 12, K = 0.9, t_anchor = 1, C = 1, ts = 1)
#' popSize   <- 40
#' maxiter   <- 30
#' run       <- 10
#' pmutation <- 0.2
#' nresamp   <-  3
#'
#' # Parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       up_par = up_par, low_par = low_par,
#'                       parallel = TRUE,
#'                       popSize = popSize, maxiter = maxiter,
#'                       run = run, pmutation = pmutation,
#'                       nresamp = nresamp)
#'
#' res
#' }
ELEFAN_GA_boot <- function(lfq,
                           seasonalised = FALSE,
                           low_par = NULL, up_par = NULL,
                           popSize = 50, maxiter = 100, run = maxiter,
                           parallel = FALSE,
                           pmutation = 0.1, pcrossover = 0.8,
                           elitism = base::max(1, round(popSize * 0.05)),
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           ...,
                           seed = NULL, nresamp = 10, resample = TRUE,
                           outfile = NA){

  if(is.null(seed)) seed <- as.numeric(Sys.time())

  # Empty results list
  res <- vector(mode = "list", length = nresamp)

  for(x in seq(nresamp)){
    res[[x]] <- lfq_ELEFAN_GA(lfq = lfq, x = x, resample = resample, seed = seed,
                              seasonalised = seasonalised,
                              low_par = low_par, up_par = up_par,
                              popSize = popSize, maxiter = maxiter, run = run,
                              pmutation = pmutation, pcrossover = pcrossover,
                              elitism = elitism, parallel = parallel,
                              MA = MA, addl.sqrt = addl.sqrt, agemax = agemax,
                              ...)
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
lfq_ELEFAN_GA <- function(lfq, x, resample, seed, seasonalised, low_par, up_par,
                          popSize, maxiter, run, pmutation, pcrossover,
                          elitism, MA, addl.sqrt, agemax, parallel, ...){

  set.seed(x + seed)

  # resample data
  if(resample){
    lfqb <- lfqResample(lfq = lfq)
  }else{
    lfqb <- lfq
  }

  fitboot <- ELEFAN_GA(lfqb,
                       seasonalised = seasonalised,
                       low_par = low_par, up_par = up_par,
                       popSize = popSize, maxiter = maxiter, run = run,
                       pmutation = pmutation, pcrossover = pcrossover,
                       elitism = elitism, parallel = parallel,
                       MA = MA, addl.sqrt = addl.sqrt, agemax = agemax,
                       flagging.out = FALSE, monitor = FALSE,
                       plot = FALSE, plot.score = FALSE, seed = x + seed, ...)

  # return result
  c(unlist(fitboot$par), seed = seed + x)
}
