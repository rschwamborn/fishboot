#' @title Bootstraped ELEFAN_GA
#'
#' @description This function performs a bootstrapped fitting of von Bertalanffy
#' growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_GA}} function.
#' Most of the arguments are simply passed to the function within many
#' permutations (resampling) of the original \code{lfq} data.
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
#' in the best fitness value before the GA is stopped. Default: equals maxiter
#' @param pmutation the probability of mutation in a parent chromosome.
#' Usually mutation occurs with a small probability, and by default is set to 0.1.
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value and by default is set to 0.8.
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for
#'    more information see \link[TropFishR]{lfqRestructure})
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE, for
#'    more information see \link[TropFishR]{lfqRestructure})
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param seed an integer value containing the random number generator state. This
#' argument can be used to replicate the results of a GA search. Note that
#' if parallel computing is required, the doRNG package must be installed.
#' (Default: 'seed = NULL').
#' @param ... additional parameters to pass to \code{\link[GA]{ga}}.
#' @param seed seed value for random number reproducibility.
#' @param nresamp \code{numeric}, the number of permutations to run (by default
#' \code{nresamp = 200}).
#' @param resample \code{logical}. Do you want that \code{lfq} object be
#' resampled (\code{TRUE} by default).
#' @param no_cores positive integer. If \code{no_cores} > 1, a 'parallel' package
#' cluster with that many cores is created.
#' @param outfile \code{character}; path of the file which will register the
#' progress of the permutation completions. If it is set as \code{false},
#' \code{NA} or \code{NULL}, no file will be created.
#'
#' @details
#' If \code{resample = FALSE}, a \strong{partial bootstrap} is performed,
#' reflecting solution variation due only to the search algorithm.
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
#' \dontrun{
#' # load data
#' data("alba", package = "TropFishR")
#'
#' # Define settings (for demo only)
#' MA        <- 7
#' low_par   <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par    <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' popSize   <- 60
#' maxiter   <- 50
#' run       <- 10
#' pmutation <- 0.2
#' nresamp   <- 12
#'
#'
#' # Parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       up_par = up_par, low_par = low_par,
#'                       popSize = popSize, maxiter = maxiter,
#'                       run = run, pmutation = pmutation,
#'                       nresamp = nresamp, seed = 1,
#'                       parallel = TRUE, no_cores = parallel::detectCores() - 2)
#'
#' res
#'
#'
#' # Non-parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       up_par = up_par, low_par = low_par,
#'                       popSize = popSize, maxiter = maxiter,
#'                       run = run, pmutation = pmutation,
#'                       nresamp = nresamp, seed = 1,
#'                       parallel = FALSE)
#'
#' res
#'
#' # Plot scatterhist of Linf and K
#' LinfK_scatterhist(res = res)
#' }
ELEFAN_GA_boot <- function(lfq,
                           seasonalised = FALSE,
                           low_par = NULL, up_par = NULL,
                           popSize = 50, maxiter = 100, run = maxiter,
                           pmutation = 0.1, pcrossover = 0.8,
                           elitism = base::max(1, round(popSize * 0.05)),
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           ...,
                           seed = NULL, nresamp = 200, resample = TRUE,
                           no_cores = 1, outfile = NA){

  if(is.null(seed)) seed <- as.integer(x = runif(n = 1, min = 0, max = 1e6))

  if(no_cores > 1){

    # Registering cluster
    cl <- makeCluster(spec = no_cores)
    registerDoParallel(cl = cl)

    # Run multithread process
    res <- foreach(x = seq(nresamp), .inorder = FALSE, .packages = "TropFishR") %dopar% {
      lfq_ELEFAN_GA(lfq = lfq, x = x, resample = resample, seed = seed,
                    seasonalised = seasonalised,
                    low_par = low_par, up_par = up_par,
                    popSize = popSize, maxiter = maxiter, run = run,
                    pmutation = pmutation, pcrossover = pcrossover,
                    elitism = elitism,
                    MA = MA, addl.sqrt = addl.sqrt, agemax = agemax, ...)
    }

    # Finish cluster
    stopCluster(cl)
  }else{

    # Empty results list
    res <- vector(mode = "list", length = nresamp)

    for(x in seq(nresamp)){
      res[[x]] <- lfq_ELEFAN_GA(lfq = lfq, x = x, resample = resample, seed = seed,
                                seasonalised = seasonalised,
                                low_par = low_par, up_par = up_par,
                                popSize = popSize, maxiter = maxiter, run = run,
                                pmutation = pmutation, pcrossover = pcrossover,
                                elitism = elitism,
                                MA = MA, addl.sqrt = addl.sqrt, agemax = agemax,
                                ...)
    }
  }

  res <- do.call(what = rbind, args = res) |> as.data.frame()

  res <- list(bootRaw = res[,-ncol(res)],
              seed = res$seed)

  class(res) <- "lfqBoot"

  res
}


# Engine function
lfq_ELEFAN_GA <- function(lfq, x, resample, seed, seasonalised, low_par, up_par,
                          popSize, maxiter, run, pmutation, pcrossover,
                          elitism, MA, addl.sqrt, agemax, ...){

  set.seed(x + seed)

  # resample data
  if(resample){
    lfqb <- lfqResample(lfq = lfq)
  }else{
    lfqb <- lfq
  }

  fitboot <- ELEFAN_GA(lfq,
                       seasonalised = seasonalised,
                       low_par = low_par, up_par = up_par,
                       popSize = popSize, maxiter = maxiter, run = run,
                       pmutation = pmutation, pcrossover = pcrossover,
                       elitism = elitism,
                       MA = MA, addl.sqrt = addl.sqrt, agemax = agemax,
                       flagging.out = FALSE, monitor = FALSE,
                       plot = FALSE, plot.score = TRUE, ...)

  # return result
  c(unlist(fitboot$par), seed = seed + x)
}
