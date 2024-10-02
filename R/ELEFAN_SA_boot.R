#' @title Bootstraped ELEFAN_SA
#'
#' @description This function performs a bootstrapped fitting of von Bertalanffy
#' growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_SA}} function.
#' Most of the arguments are simply passed to the function within many
#' permutations (resampling) of the original \code{lfq} data.
#'
#'
#' @param lfq a length frequency object of the class \code{lfq} (see
#' \link[TropFishR]{lfqCreate}).
#' @param ... Extra arguments passed to \link[TropFishR]{ELEFAN_SA}:
#' \code{seasonalised}, \code{init_par}, \code{low_par}, \code{up_par},
#' \code{SA_time}, \code{maxit}, \code{nb.stop.improvement}, \code{SA_temp},
#' \code{verbose}, \code{MA}, \code{addl.sqrt}, \code{agemax} and
#' \code{flagging.out}. Default values remain the same and just \code{plot} and
#' \code{plot.score} will be set as \code{FALSE}.
#' @param nresamp \code{numeric}, the number of permutations to run (by default
#' \code{nresamp = 200}).
#' @param resample \code{logical}. Do you want that \code{lfq} object be
#' resampled (\code{TRUE} by default).
#' @param seed seed value for random number reproducibility.
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
#' @return An object of class `lfqBoot` containing 2 levels:
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
#' \dontrun{
#' # load data
#' data("alba", package = "TropFishR")
#'
#' # Define settings (for demo only)
#' MA        <- 7
#' low_par   <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par    <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' SA_time   <- 3
#' SA_temp   <- 1e5
#' nresamp   <- 12
#'
#'
#' # Parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       init_par = NULL, up_par = up_par, low_par = low_par,
#'                       SA_time = SA_time, SA_temp = SA_temp,
#'                       nresamp = nresamp, seed = 1,
#'                       parallel = TRUE, no_cores = parallel::detectCores() - 2)
#'
#' res
#'
#'
#' # Non-parallel version
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'                       init_par = NULL, up_par = up_par, low_par = low_par,
#'                       SA_time = SA_time, SA_temp = SA_temp,
#'                       nresamp = nresamp, seed = 1,
#'                       parallel = FALSE)
#'
#' res
#'
#' # Plot scatterhist of Linf and K
#' LinfK_scatterhist(res = res)
#' }
ELEFAN_SA_boot <- function(lfq, ..., nresamp = 200, resample = TRUE, seed = NULL,
                           no_cores = 1, outfile = NA){

  if(is.null(seed)) seed <- as.integer(x = runif(n = 1, min = 0, max = 1e6))

  if(no_cores > 1){

    # Registering cluster
    cl <- makeCluster(spec = no_cores)
    registerDoParallel(cl = cl)

    # Run multithread process
    res <- foreach(x = seq(nresamp), .inorder = FALSE, .packages = "TropFishR") %dopar% {
      lfq_ELEFAN_SA(lfq = lfq, x = x, resample = resample, seed = seed, ...)
    }

    # Finish cluster
    stopCluster(cl)
  }else{

    # Empty results list
    res <- vector("list", nresamp)

    for(x in seq(res)){
      res[[x]] <- lfq_ELEFAN_SA(lfq = lfq, x = x, resample = resample, seed = seed, ...)
    }
  }

  res <- do.call(what = rbind, args = res) |> as.data.frame()

  res <- list(bootRaw = res[,-ncol(res)],
              seed = res$seed)

  class(res) <- "lfqBoot"

  res
}

lfq_ELEFAN_SA <- function(lfq, x, resample, seed, ...){
  set.seed(x + seed)

  # resample data
  if(resample){
    lfqb <- lfqResample(lfq = lfq)
  }else{
    lfqb <- lfq
  }

  fitboot <- formals(ELEFAN_SA) |>

    # Combine arguments
    modifyList(val = append(list(lfq = lfqb), list(...))) |>

    as.list() |>

    # Set plotting arguments
    modifyList(val = list(verbose = FALSE,
                          plot = FALSE,
                          plot.score = FALSE)) |>

    as.list() |>

    # Call ELEFAN_GA
    do.call(what = ELEFAN_SA)

  # return result
  c(unlist(fitboot$par), seed = seed + x)
}
