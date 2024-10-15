#' @title Bootstrapped length-at-age analysis
#'
#' @description
#' This function estimates growth parameters from length-at-age data.
#' Since it internally uses the function \link[TropFishR]{growth_length_age},
#' this function allows to perform different methods: Gulland and Holt, Ford
#' Walford, Chapman, Bertalanffy, or non linear least squares method (LSM).
#'
#'
#' @param param a \code{list} (or \code{data.frame}) consisting of following
#' parameters (levels/columns):
#' \itemize{
#'   \item \strong{age}: age measurements,
#'   \item \strong{length}: corresponding lengths in cm.
#' }
#' @param method indicating which of following methods should be applied:
#' \code{"GullandHolt"}, \code{"FordWalford"}, \code{"Chapman"},
#' \code{"BertalanffyPlot"}, or \code{"LSM"} (default).
#' @param Linf_est BertalanffyPlot requires an estimate for \eqn{L_{inf}} to
#' derive \eqn{K} and \eqn{t_0} (for more information see Details).
#' @param Linf_init initital parameter of \eqn{L_{inf}} for non-linear sqaures
#' fitting (default 10).
#' @param K_init initital parameter of \eqn{K} for non-linear sqaures fitting
#' (default 0.1).
#' @param t0_init initital parameter of \eqn{L_0} for non-linear squares fitting
#' (default 0).
#' @param seed seed value for random number reproducibility.
#' @param nresamp \code{numeric}; the number of permutations to run (Default:
#' \code{nresamp = 200}).
#' @param no_cores positive integer. If \code{no_cores} > 1, a 'parallel' package
#' cluster with that many cores is created.
#' @param cc.only \code{logical} Return complete cases only:  Do you want the
#' output object to omit results with failed estimates (with NaNs)? Remember
#' that if you define \code{cc.only = TRUE} the number of elements returned may
#' be less than specified in \code{nresamp}.
#'
#' @details
#' CI and plotting arguments are set as \code{FALSE} for each bootstrap call
#' here.
#'
#' By default, \link[TropFishR]{growth_length_age} generates a plot when is
#' called, so internally \code{grolenage_boot} executes a
#' \link[graphics]{dev.off} call in order to prevent it.
#'
#' It is important to take into account the particular considerations of each
#' method regarding the required parameters, so it is recommended to read the
#' Details of the documentation of \link[TropFishR]{growth_length_age}.
#'
#' Some words about \eqn{t_{anchor} = t_0}...
#'
#'
#' @return An object of class \code{lfqBoot} containing 2 levels:
#' \describe{
#'   \item{\code{$bootRaw}}{A \code{data.frame} of fitted VBGF parameters
#'   (columns) by resampling (rows).}
#'   \item{\code{$seed}}{A \code{numeric} vector of seed values set prior to each
#'   resampling call.}
#' }
#'
#' \code{t_anchor} (or \code{t0} in \link[TropFishR]{growth_length_age}) will be
#' only available for Bertalanffy plot and LSM method, otherwise a vector of
#' zeroes will be returned instead.
#'
#' @export
#'
#' @examples
#' # Synthetical length at age data
#' dat <- list(age = rep(x = 1:7,each = 5),
#'             length = c(rnorm(n = 5, mean = 25.7, sd = 0.9),
#'                        rnorm(n = 5, mean = 36.0, sd = 1.2),
#'                        rnorm(n = 5, mean = 42.9, sd = 1.5),
#'                        rnorm(n = 5, mean = 47.5, sd = 2.0),
#'                        rnorm(n = 5, mean = 50.7, sd = 0.4),
#'                        rnorm(n = 5, mean = 52.8, sd = 0.5),
#'                        rnorm(n = 5, mean = 54.2, sd = 0.7)))
#'
#' growth_length_age(param = dat, method = "GullandHolt")
#'
#' res <- grolenage_boot(param = dat, method = "GullandHolt")
#'
#' # Plot scatterhist of Linf and K
#' LinfK_scatterhist(res = res)
#'
#' # Bertalaffy plot
#' res <- grolenage_boot(dat, method = "BertalanffyPlot", Linf_est = 50,
#'                       nresamp = 100)
#'
#' vbgfCI_time(res = res)
grolenage_boot <- function(param, method = "LSM",
                           Linf_est = NA, Linf_init = 10,
                           K_init = 0.1, t0_init = 0,
                           seed = as.numeric(Sys.time()),
                           nresamp = 200, no_cores = 1, cc.only = FALSE){

  # Coerce param to a list
  param <- as.list(param)

  # Standardize names to lowercase
  names(param) <- tolower(names(param))

  # Select only length and age variables
  param <- param[c("length", "age")]

  if(no_cores > 1){
    # Multicore way
    # Registering cluster
    cl <- makeCluster(spec = no_cores)
    registerDoParallel(cl = cl)

    # Run multithread process
    res <- foreach(x = seq(nresamp), .inorder = FALSE, .packages = "TropFishR") %dopar% {
      tryCatch({
        grolenage_internal(param = param, method = method, seed = seed, x = x,
                           Linf_est = Linf_est, Linf_init = Linf_init,
                           K_init = K_init, t0_init = t0_init)
      }, error = \(e){NULL})
    }

    # Finish cluster
    stopCluster(cl)
  }else{
    # Single core way
    res <- vector(mode = "list", length = nresamp)
    for(x in seq(nresamp)){
      tryCatch({
        res[[x]] <- grolenage_internal(param = param, method = method,
                                       seed = seed, x = x,
                                       Linf_est = Linf_est, Linf_init = Linf_init,
                                       K_init = K_init, t0_init = t0_init)
      }, error = \(e){NULL}) |> suppressWarnings()
    }
  }


  # Build the output object
  res <- list(bootRaw = lapply(X = res, FUN = \(x) x$out) |>

                do.call(what = rbind) |> as.data.frame(),
              seed = lapply(X = res, FUN = \(x) x$seed) |> do.call(what = c))

  # If cc.only is TRUE, indexing rows with non-NaN values
  if(isTRUE(cc.only)){

    index <- complete.cases(res$bootRaw)

    if(sum(is.na(index)) == nresamp){
      warning("All the resamples contain NaN values, try changing 'seed' of setting a larger value for 'nresamp'.")
    }

    res$bootRaw <- res$bootRaw[index,]
    res$seed <- res$seed[index]
  }

  message("t_anchor = t0")

  # Set the class
  class(res) <- "lfqBoot"

  res
}

grolenage_internal <- function(param, method, seed, x,
                               Linf_est, Linf_init, K_init, t0_init){

  # Set seed
  set.seed(seed + x)

  # Resample values of length and age
  out <- lapply(X = param, FUN = \(x, index) x[index],
                index = sample(x = length(param$length), replace = TRUE)) |>

    # Execute growth_length_age
    growth_length_age(method = method,
                      Linf_est = Linf_est, Linf_init = Linf_init,
                      K_init = K_init, t0_init = t0_init,
                      CI = FALSE, age_plot = FALSE, do.sim = FALSE)

  # Flush the default plot
  dev.off()

  # Prepare output
  list(out = c(Linf = out$Linf,
               K = out$K,
               t_anchor = if(is.null(out$t0)) NA else out$t0,
               phiL = log10(out$K) + 2*log10(out$Linf)),
       seed = seed + x)
}
