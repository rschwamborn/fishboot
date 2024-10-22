#' @title Bootstrapped length-at-age analysis
#'
#' @description
#' This function obtains growth parameter estimates from length-at-age data.
#' Since it internally uses the function \link[TropFishR]{growth_length_age},
#' this function allows to perform different methods: Gulland and Holt, Ford
#' Walford, Chapman, Bertalanffy, or non linear least squares method (LSM).
#'
#' This function performs bootstrapped fitting of the von Bertalanffy growth
#' function with estimated growth parameters (\eqn{L_{inf}}, \eqn{K} and
#' \eqn{t_0}) from length-at-age data, based on the function
#' \link[TropFishR]{growth_length_age}. The output is an object containing the
#' parameters \eqn{L_{inf}}, \eqn{K} and \eqn{t_0} (named here \code{t_anchor})
#' as well as the growth performance index \eqn{Phi’} (named \code{PhiL}).
#'
#'
#' @param param a \code{list} (or \code{data.frame}) consisting of following
#' parameters (levels/columns):
#' \itemize{
#'   \item \strong{age}: age measurements (e.g. from otoliths),
#'   \item \strong{length}: corresponding length measurements.
#' }
#' @param method indicating which of following methods should be applied:
#' \code{"GullandHolt"}, \code{"FordWalford"}, \code{"Chapman"},
#' \code{"BertalanffyPlot"}, or \code{"LSM"} (it corresponds to the non-linear
#' least squares fitting method, and is the default, which is recommended for
#' bootstrapping growth).
#' @param Linf_est BertalanffyPlot requires an estimate for \eqn{L_{inf}} to
#' derive \eqn{K} and \eqn{t_0} (for more information see Details).
#' @param Linf_init initital parameter of \eqn{L_{inf}} for non-linear sqaures
#' fitting (default 100).
#' @param K_init initital parameter of \eqn{K} for non-linear sqaures fitting
#' (default 0.1).
#' @param t0_init initital parameter of \eqn{L_0} for non-linear squares fitting
#' (default 0).
#' @param seed seed value for random number reproducibility (if it \code{NULL}
#' by default, it will set internally as \code{seed = as.numeric(Sys.time())}).
#' @param nresamp \code{numeric}; the number of permutations to run (Default:
#' \code{nresamp = 200}).
#' @param nan_action \code{character} that defines the action that the function
#' will execute if there is a row with \code{NaN} (growth rate parameters
#' inestimable for that resample):
#' \itemize{
#'  \item \code{nothing}: the function will return the results including the
#'  \code{NaN}s (default).
#'  \item \code{nanrm} or \code{narm}: after having the results, it will only
#'  returns the rows without \code{NaN}s. For this case \code{narm} and
#'  \code{nanrm} are equivalent, but it should be noted that the function will
#'  look for and omit the \code{NaN}s (and not the \code{NA}s). See Details.
#'  \item \code{force}: The function will start an iterative process changing
#'  the internal \code{seed} values until it fulfills the \code{nresamp}. It
#'  only works together with the \code{time_lim} argument. See Details.
#' }
#' @param time_lim If \code{nan_action = "force"}, it defines the maximum time
#' (in seconds) that the function will last resampling until it achieves a
#' result output with no-\code{NaN} rows.
#'
#' @details
#' It is important to take into account the particular considerations of each
#' method regarding the required parameters, so it is recommended to read the
#' Details of the documentation of \link[TropFishR]{growth_length_age}.
#'
#' CI and plotting arguments (of \link[TropFishR]{growth_length_age}) are set as
#' \code{FALSE} for each bootstrap call here. By default,
#' \link[TropFishR]{growth_length_age} generates a plot when it is called, so
#' internally \code{grolenage_boot} executes a \code{dev.off} call in order to
#' prevent it.
#'
#' \code{nan_action = "force"} should be used carefully, as estimated \code{NaN}
#' VBGF parameter values are not always a result of bootstrap data selection
#' factors. Few resamples should first be tested with different \code{Linf_init},
#' \code{K_init} and \code{t0_init} values. No selection of the realistic initial
#' parameters may also result in \code{NaN} values being obtained. The search
#' time may depend on the size of the input set. For example, if you have many
#' thousands of individuals or if (in addition) the value of \code{nresamp} is
#' too high, it is possible that the function will take a long time before
#' obtaining complete results. Even though \code{time_lim} avoids falling into
#' an infinite loop by limiting the time used by this process to 5 minutes, this
#' value is referential and might be insufficient due to the factors mentioned
#' above.
#'
#' \code{t_anchor} is the true \eqn{t_0} estimate in the case of true
#' length-at-age data, but it will only be available from "BertalanffyPlot" or
#' "LSM" methods. For the other methods, a vector of \code{NA}s will be returned
#' instead.
#'
#' @references \itemize{
#'  \item Efron, B., & Tibshirani, R., 1986. Bootstrap methods for standard
#'  errors, confidence intervals, and other measures of statistical accuracy.
#'  Statistical Science, 54-75.
#'  \item Pauly, D. 1981. The relationship between gill surface area and growth
#'  performance in fish: a generalization of von Bertalanffy's theory of growth.
#'  Meeresforsch. 28:205-211.
#'  \item Schwamborn, R., Mildenberger, T. K., & Taylor, M. H., 2019. Assessing
#'  sources of uncertainty in length-based estimates of body growth in
#'  populations of fishes and macroinvertebrates with bootstrapped ELEFAN.
#'  Ecological Modelling, 393, 37-51.
#'  \item Schwamborn, R., Freitas, M. O., Moura, R. L., & Aschenbrenner, A. 2023.
#'  Comparing the accuracy and precision of novel bootstrapped length-frequency
#'  and length-at-age (otolith) analyses, with a case study of lane snapper
#'  (\emph{Lutjanus synagris}) in the SW Atlantic. Fisheries Research, 264,
#'  106735.
#'  \item von Bertalanffy, L., 1938. A quantitative theory of organic growth.
#'  Human Biology 10, 181-213.
#' }
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
#' @export
#'
#' @examples
#' # Synthetical length at age data
#' dat <- list(age = rep(x = 1:7,each = 15),
#'             length = c(rnorm(n = 15, mean = 4.6, sd = 4),
#'                        rnorm(n = 15, mean = 22.8, sd = 7),
#'                        rnorm(n = 15, mean = 35, sd = 7),
#'                        rnorm(n = 15, mean = 43, sd = 7),
#'                        rnorm(n = 15, mean = 49, sd = 5),
#'                        rnorm(n = 15, mean = 53, sd = 5),
#'                        rnorm(n = 15, mean = 57, sd = 3)))
#'
#' # Perform bootstrapped curve fitting with grolenage_boot
#' res <- grolenage_boot(param = dat, nresamp = 70)
#'
#' # Plot scatter histograms of Linf and K
#' LinfK_scatterhist(res = res)
#'
#' # Plot univariate density plots of all parameters
#' univariate_density(res = res)
#'
#' # Plot swarm plots of all n bootstraps
#' vbgfCI_time(res = res)
#'
#' # Extract data.frame with all parameter estimates
#' # for comparisons of posterior distributions
#' print(res$bootRaw)
grolenage_boot <- function(param, method = "LSM",
                           Linf_est = NA, Linf_init = 100,
                           K_init = 0.1, t0_init = 0,
                           seed = NULL, nresamp = 200,
                           nan_action = "nothing", time_lim = 5*60){

  # Coerce param to a list
  param <- as.list(param)

  # Standardize names to lowercase
  names(param) <- tolower(names(param))

  # Select only length and age variables
  param <- param[c("length", "age")]

  # Set seed
  if(is.null(seed)) seed <- as.numeric(Sys.time())

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
