#' Bootstrapped ELEFAN_GA
#'
#' @param lfq a length frequency object of the class `lfq`
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param low_par a list providing the minimum of the search space in case
#' of real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the maximum of the search space in case of
#' real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 1),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 1),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 1),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 1);
#' }
#' @param parallel logical; should parallelized computing be used. This differs from the
#'    `parallel` argument in `ELEFAN_GA` in that it is not used within the `ga` function for
#'    calculation at the population level, but rather for permutations. Depending on platform
#'    operating system, the argument `clusterType` can be adjusted (see argument description for
#'    details). (Default: `parallel = TRUE`)
#' @param nresamp numeric; the number of permutations to run (Default: `nresamp = 10`)
#' @param no_cores numeric (Default: `no_cores = detectCores() - 1`)
#' @param clusterType (Default: `clusterType = "PSOCK"`)
#' @param outfile character; text file name (Default: `outfile = "output.txt"`) which will
#'    records the progress of the permutation completions.
#' @param popSize the population size. (Default: `popSize = 60`)
#' @param maxiter the maximum number of iterations to run before the GA search is halted.
#'   (Default: `maxiter = 50`)
#' @param run the number of consecutive generations without any improvement
#'   in the best fitness value before the GA is stopped. (Default: `run = 10`)
#' @param pmutation numeric. A small fraction of 1.0. The probability of mutation in a
#'   parent chromosome. Usually mutation occurs with a small probability.
#'   (Default: `pmutation = 0.2`)
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value.  (Default: `pcrossover = 0.8`)
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for more information see \link{lfqRestructure})
#' @param addl.sqrt logical. Should counts be square root transformed prior to restructuring.
#' @param agemax numeric. maximum age
#' @param flagging.out logical Should flagging out be done
#' @param resample logical. Should `lfq` object be resampled (Default: `resample = TRUE`).
#'   When `resample = FALSE`, a `partial bootstrap` is performed, reflecting solution
#'   variation due only to the search algorithm.
#' @param seed seed value for random number reproducibility (Default: NULL)
#'
#' @description `ELEFAN_GA_boot` performs a bootstrapped fitting of
#'   von Bertalanffy growth function (VBGF) via the \link{ELEFAN_GA} function. Most of the arguments
#'   are simply passed to the function within many permutations (resampling) of the original
#'   lfq data.
#'
#' @return a data.frame of fitted VBGF parameters (columns) by permutation (rows).
#' @export
#'
#' @examples
#' \donttest{
#' # load data
#' library(TropFishR)
#' data(alba)
#'
#' # settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' popSize <- 60
#' maxiter <- 50
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 12
#'
#'
#' # parallel version
#' library(parallel)
#' t1 <- Sys.time()
#' res <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
#'   seed = 1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#'
#' # non-parallel version
#' t1 <- Sys.time()
#' res <- ELEFAN_GA_boot(lfq=alba, seasonalised = FALSE,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run,
#'   pmutation = pmutation, nresamp = nresamp, MA = MA, parallel = FALSE,
#'   seed = 1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#' # plot resulting distributions
#' univariate_density(res, use_hist = TRUE)
#' }
#'
ELEFAN_GA_boot <- function(lfq, seasonalised = FALSE, low_par = NULL, up_par = NULL,
  parallel = TRUE, nresamp = 10, no_cores = detectCores() - 1, clusterType = "PSOCK",
  outfile = "output.txt",
  popSize = 60, maxiter = 50, run = 10,
  pmutation = 0.2, pcrossover = 0.8,
  elitism = base::max(1, round(popSize * 0.05)),
  MA = 5, addl.sqrt = FALSE, agemax = NULL, flagging.out = TRUE,
  resample = TRUE, seed = NULL
){

  if(!is.null(outfile)){unlink(outfile)} # delete old outfile

  if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6))}

  if(parallel){ # Parallel version
    ARGS <- list(
      "lfqResample",
      "lfq", "seasonalised", "low_par", "up_par",
      "parallel", "nresamp", "no_cores",
      "popSize", "maxiter", "run",
      "pmutation", "pcrossover", "elitism",
      "MA", "addl.sqrt", "agemax", "flagging.out",
      "resample", "seed"
    )

    parFun <- function(x){
      set.seed(seed + x)

      # resample data
      if(resample){
        lfqb <- lfqResample(lfq)
      }else{
        lfqb <- lfq
      }

      # call ELEFAN_GA
      fitboot <- TropFishR::ELEFAN_GA(
        lfqb, seasonalised = seasonalised,
        low_par = low_par,
        up_par = up_par,
        popSize = popSize, maxiter = maxiter, run = run,
        pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
        MA = MA, parallel = FALSE, addl.sqrt = addl.sqrt,
        agemax = agemax, flagging.out = flagging.out,
        plot.score = FALSE, seed = NULL
      )

      # print output (for checking progress in output.txt)
      if(!is.null(outfile)){
        sink(file=outfile, append = TRUE)
        print(paste("resamp", x, "completed @", Sys.time()))
        sink()
      }

      # return result
      return( c(unlist(fitboot$par), seed + x) )
    }

    cl <- parallel::makeCluster(no_cores, type=clusterType)
    nn <- split(1:nresamp, 1:nresamp)
    parallel::clusterExport(cl, varlist = ARGS, envir=environment())
    res <- parLapply(cl, nn, parFun)
    stopCluster(cl)
  }

  if(!parallel){ # Non-parallel version
    res <- vector("list", nresamp) # empty results list
    for(x in seq(res)){
      set.seed(x + seed)

      # resample data
      if(resample){
        lfqb <- lfqResample(lfq)
      }else{
        lfqb <- lfq
      }

      # call ELEFAN_GA
      fitboot <- TropFishR::ELEFAN_GA(
        lfqb, seasonalised = seasonalised,
        low_par = low_par,
        up_par = up_par,
        popSize = popSize, maxiter = maxiter, run = run,
        pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
        MA = MA, parallel = FALSE, addl.sqrt = addl.sqrt,
        agemax = agemax, flagging.out = flagging.out,
        plot.score = FALSE, seed = NULL
      )

      if(!is.null(outfile)){
        sink(file=outfile, append = TRUE)
        print(paste("resamp", x, "completed @", Sys.time()))
        sink()
      }

      # return result
      res[[x]] <- c(unlist(fitboot$par), seed + x)
    }

  }

  tmp0 <- as.data.frame(do.call("rbind", res))
  tmp <- tmp0[,-ncol(tmp0)]

  ret <- list()
  ret$bootRaw <- tmp
  ret$seed <- as.numeric(tmp0[,ncol(tmp0)])
  class(ret) <- "lfqBoot"

  return(ret)
}

