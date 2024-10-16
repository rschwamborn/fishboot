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
#' @param nan_action \code{character} that defines the action that the function
#' will execute if there is a row with NaN:
#' \itemize{
#'  \item \code{nothing}: the function will return the results including the NaNs
#'  (default).
#'  \item \code{nanrm} or \code{narm}: after having the results, it will only
#'  returns the rows without NaNs. For an intercompatibility issue, for this
#'  case \code{narm} and \code{nanrm} are equivalent, but it should be noted
#'  that the function will look for and omit the NaNs (and not the NAs). See
#'  Details.
#'  \item \code{force}: The function will start an iterative process changing
#'  the internal \code{seed} values until it fulfills the \code{nresamp}. It
#'  works just together \code{time_lim} argument. See Details.
#' }
#' @param time_lim If \code{nan_action = "force"}, it defines the maximum time
#' (in seconds) that the function will last resampling until it achieves a
#' result output with no-NaN rows.
#'
#' @details
#' CI and plotting arguments are set as \code{FALSE} for each bootstrap call
#' here.
#'
#' By default, \link[TropFishR]{growth_length_age} generates a plot when is
#' called, so internally \code{grolenage_boot} executes a
#' \link[grDevices]{dev.off} call in order to prevent it.
#'
#' It is important to take into account the particular considerations of each
#' method regarding the required parameters, so it is recommended to read the
#' Details of the documentation of \link[TropFishR]{growth_length_age}.
#'
#' \code{nan_action = "force"} should be used carefully, as it is not always due
#' to bootstrap data selection factors, but also to an inadequate selection of
#' the estimation parameters that the \code{NaN} values are obtained. Also, the
#' search time may depend on the size of the input set, if you have many
#' thousands of individuals or if (in addition) the value of \code{nresamp} is
#' high, it is possible that the function will take a long time before obtaining
#' complete results. \code{time_lim} avoids falling into an infinite loop by
#' limiting the time used by this process to 5 minutes, but this value is
#' referential and may be insufficient due to the factors mentioned above.
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
#' TropFishR::growth_length_age(param = dat, method = "GullandHolt")
#'
#' res <- grolenage_boot(param = dat, method = "GullandHolt")
#'
#' # Plot scatterhist of Linf and K
#' LinfK_scatterhist(res = res)
grolenage_boot <- function(param, method = "LSM",
                           Linf_est = NA, Linf_init = 10,
                           K_init = 0.1, t0_init = 0,
                           seed = as.numeric(Sys.time()), nresamp = 200,
                           nan_action = "nothing", time_lim = 5*60){

  # Coerce param to a list
  param <- as.list(param)

  # Standardize names to lowercase
  names(param) <- tolower(names(param))

  # Select only length and age variables
  param <- param[c("length", "age")]

  res <- list()
  nan_action <- tolower(nan_action)[1]
  time_0 <- Sys.time()
  x <- 0
  while(length(res) < nresamp){
    x <- x + 1

    out <- grolenage_internal(param = param, method = method,
                              seed = seed, x = x,
                              Linf_est = Linf_est, Linf_init = Linf_init,
                              K_init = K_init, t0_init = t0_init) |>

      suppressWarnings()


    outnan <- any(is.nan(out$out))
    if(!outnan || nan_action == "nothing"){
      # Adding output as is
      res <- c(res, list(out))
    }else{
      if(nan_action == "force"){
        # Check time diff
        t_diff <- difftime(time1 = Sys.time(), time2 = time_0, units = "secs")
        if(t_diff >= time_lim) break
      }else{
        if(x == nresamp) break
      }
    }
  }


  # Build the output object
  res <- list(bootRaw = lapply(X = res, FUN = \(x) x$out) |>

                do.call(what = rbind) |> as.data.frame(),
              seed = lapply(X = res, FUN = \(x) x$seed) |> do.call(what = c))

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
