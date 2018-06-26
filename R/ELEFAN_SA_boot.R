#' Bootstraped ELEFAN_SA
#'
#' @param lfq a length frequency object of the class `lfq`
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param init_par a list providing the Initial values for the components to be
#' optimized. When set to NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is the maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.5),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0.5),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
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
#' @param parallel logical; should parallelized computing be used. Depending on platform
#'    operating system, the argument `clusterType` can be adjusted (see argument description for
#'    details). (Default: `parallel = TRUE`)
#' @param nresamp numeric; the number of permutations to run (Default: `nresamp = 200`)
#' @param no_cores numeric (Default: `no_cores = detectCores() - 1`)
#' @param clusterType (Default: `clusterType = "PSOCK"`)
#' @param outfile character; text file name (Default: `outfile = "output.txt"`) which will
#'    records the progress of the permutation completions.
#' @param SA_time numeric; Maximum running time in seconds (default : 60 * 1).
#' @param SA_temp numeric; Initial value for temperature (default : 1e5).
#' @param maxit numeric; default: maxit = NULL.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for more information see \code{\link[TropFishR]{lfqRestructure}}.
#' @param addl.sqrt logical. Should counts be square root transformed prior to restructuring.
#' @param agemax numeric. maximum age
#' @param flagging.out logical Should flagging out be done
#' @param resample logical. Should `lfq` object be resampled (Default: `resample = TRUE`).
#'   When `resample = FALSE`, a `partial bootstrap` is performed, reflecting solution
#'   variation due only to the search algorithm.
#' @param seed seed value for random number reproducibility (Default: NULL)
#'
#' @description `ELEFAN_SA_boot` performs a bootstrapped fitting of
#'   von Bertalanffy growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_SA}} function.
#'   Most of the arguments are simply passed to the function within many
#'   permutations (resampling) of the original lfq data.
#'
#' @return a list of class `lfqBoot` containing 2 levels: `$bootRaw` - a data.frame of fitted VBGF parameters
#' (columns) by resampling (rows), `$seed` - a vector of seed values set prior to each resampling
#' call to `lfqResample`.
#'
#'
#' @examples
#' \donttest{
#' # load data
#' library(TropFishR)
#' data(alba)
#'
#' # settings (these settings may not be optimal - for demo only)
#' MA <- 7
#' init_par <- NULL
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1)
#' SA_time <- 3
#' SA_temp <- 1e5
#' nresamp <- 12
#'
#'
#' # parallel version
#' library(parallel)
#' t1 <- Sys.time()
#' res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
#'   init_par = init_par, up_par = up_par, low_par = low_par,
#'   SA_time = SA_time, SA_temp = SA_temp,
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
#' set.seed(1)
#' res <- ELEFAN_SA_boot(lfq=alba, seasonalised = FALSE,
#'   init_par = init_par, up_par = up_par, low_par = low_par,
#'   SA_time = SA_time, SA_temp = SA_temp,
#'   nresamp = nresamp, MA=MA, parallel = FALSE,
#'   seed = 1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#' # plot resulting distributions
#' univariate_density(res, use_hist = TRUE)
#'
#' # plot scatterhist of Linf and K
#' LinfK_scatterhist(res)
#'
#' }
#'
#' @importFrom stats quantile runif
#' @importFrom parallel detectCores parLapply stopCluster detectCores
#'
#' @export
#'
ELEFAN_SA_boot <- function(lfq, seasonalised = FALSE,
  init_par = NULL, low_par = NULL, up_par = NULL,
  parallel = TRUE, nresamp = 200, no_cores = detectCores() - 1, clusterType = "PSOCK",
  outfile = "output.txt",
  SA_time = 60, SA_temp = 1e5, maxit = NULL,
  MA = 5, addl.sqrt = FALSE, agemax = NULL, flagging.out = TRUE,
  resample = TRUE, seed = NULL
){

  if(!is.null(outfile)){unlink(outfile)} # delete old outfile

  if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6))}

  if(parallel){ # Parallel version
    ARGS <- list(
      "lfqResample",
      "lfq", "seasonalised",
      "init_par", "low_par", "up_par",
      "parallel", "nresamp", "no_cores",
      "SA_temp", "SA_time", "maxit",
      "MA", "addl.sqrt", "agemax", "flagging.out",
      "resample", "seed", "outfile"
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
      fitboot <- TropFishR::ELEFAN_SA(
        lfqb, seasonalised = seasonalised,
        init_par = init_par,
        low_par = low_par,
        up_par = up_par,
        SA_time = SA_time, SA_temp = SA_temp, maxit = maxit,
        MA = MA, addl.sqrt = addl.sqrt,
        agemax = agemax, flagging.out = flagging.out,
        plot.score = FALSE, verbose = FALSE
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
    res <- parallel::parLapply(cl, nn, parFun)
    parallel::stopCluster(cl)
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
      fitboot <- TropFishR::ELEFAN_SA(
        lfqb, seasonalised = seasonalised,
        init_par = init_par,
        low_par = low_par,
        up_par = up_par,
        SA_time = SA_time, SA_temp = SA_temp, maxit = maxit,
        MA = MA, addl.sqrt = addl.sqrt,
        agemax = agemax, flagging.out = flagging.out,
        plot.score = FALSE, verbose = FALSE
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

