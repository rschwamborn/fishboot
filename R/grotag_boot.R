#' @title Bootstrapped tag-and-recapture growth analysis
#'
#' @description
#' This function performs bootstrapped fitting of the von Bertalanffy growth
#' function (VBGF) with estimated growth parameters (\eqn{L_{inf}}, \eqn{K} and
#' \eqn{t_0}) from tag-and-recapture data, based on the function
#' \link[fishmethods]{grotag}, that estimates VBGF parameters according to
#' Francis (1988). The output is an object containing the parameters
#' \eqn{L_{inf}} and \eqn{K}, as well as the growth performance index
#' \eqn{Phi’} (named \code{PhiL}).
#'
#'
#' This function resamples the \code{input.data} data by rows (i.e.,
#' by recapture date) several times (\code{nresamp} times, default:
#' \code{nresamp = 200}). Then, a VBGF curve is fitted to each resampled data set.
#' The output (a \code{list} of class \code{lfqBoot}) will store results (e.g.,
#' VGBGF function parameters K and Linf) in a \code{data.frame} accessible
#' through \code{$bootRaw}. The \code{$bootRaw} table also includes the growth
#' performance index \strong{Phi'}, seasonal parameters \strong{u} and \strong{w}
#' (sensu Francis, 1988), which are equal to \strong{C} (sensu Pauly and Gaschütz,
#' 1979) and and \strong{ts} (sensu Mildenberger et al., 2017). The
#' \code{$bootRaw} table also includes seed values and system time.
#'
#' @param L1,L2,T1,T2 Name of the columns to be extracted from \code{input.data}
#' and used by the \link[fishmethods]{grotag} function for the arguments
#' \code{L1}, \code{L2}, \code{T1} and \code{T2}, respectively. See Details.
#' @param alpha \code{numeric} value giving an arbitrary length alpha.
#' @param beta \code{numeric} value giving an arbitrary length beta
#' (\code{beta} > \code{alpha}).
#' @param design \code{list} specifying the design of the model to estimate. Use
#' 1 to designate whether a parameter(s) should be estimated. Type of parameters
#' are:
#' \itemize{
#'  \item \code{nu}: growth variability (1 parameter).
#'  \item \code{m}: bias parameter of measurement error (1 parameter).
#'  \item \code{p}: outlier probability (1 parameter).
#'  \item \code{sea}: seasonal variation (2 parameters: u and w).
#' }
#' Model 1 of Francis is the default settings of 0 for \code{nu}, \code{m},
#' \code{p} and \code{sea}.
#' @param stvalue Starting values of sigma(s) and depending on the \code{design}
#' argument, \code{nu}, \code{m}, \code{p}, \code{u}, and \code{w} used as input
#' in the nonlinear estimation (function \link[stats]{optim}) routine.
#' @param upper,lower Upper and lower limits of the model parameters' (\code{nu},
#' \code{m}, \code{p}, \code{u}, and \code{w}) region to be investigated.
#' @param gestimate \code{logical} specifying whether starting values of
#' \strong{ga} and \strong{gb} (growth increments of \code{alpha} and \code{beta})
#' should be estimated automatically. \code{TRUE} by default.
#' @param st.ga,st.gb If \code{gestimate=FALSE}, user-specified starting value
#' for ga and gb respectively.
#' @param st.galow,st.gaup If \code{gestimate=FALSE}, user-specified lower and
#' upper limits for \code{st.ga} used in optimization.
#' @param st.gblow,st.gbup If \code{gestimate=FALSE}, user-specified lower and
#' upper limits for \code{st.gb} used in optimization.
#' @param control Additional controls passed to the optimization function
#' \link[stats]{optim}.
#' @param input.data A growth increment object of the class \code{data.frame}.
#' @param seed seed value for random number reproducibility (if it \code{NULL}
#' by default, it will set internally as \code{seed = as.numeric(Sys.time())}).
#' @param nresamp \code{numeric}; the number of permutations to run (Default:
#' \code{nresamp = 200}).
#' @param na_action \code{character} that defines the action that the function
#' will execute if there is a row with NA:
#' \itemize{
#'  \item \code{nothing}: the function will return the results including the
#'  \code{NA}s (default).
#'  \item \code{narm}: after having the results, it will only returns the rows
#'  without \code{NA}s. See Details.
#'  \item \code{force}: The function will start an iterative process changing
#'  the internal \code{seed} values until it fulfills the \code{nresamp}. It
#'  works just together \code{time_lim} argument. See Details.
#' }
#' @param time_lim If \code{na_action = "force"}, it defines the maximum time
#' (in seconds) that the function will last resampling until it achieves a
#' result output with no-\code{NaN} rows.
#'
#' @details
#' There are 2 ways to specify the main input arguments (related to the size and
#' timing of the mark-recapture): (1) in the classical way, i.e. by defining
#' \code{L1}, \code{L2}, \code{T1} and \code{T2} as \code{numeric} vectors as
#' indicated in the \link[fishmethods]{grotag} documentation or (2) through a
#' \code{data.frame} indicated in the \code{input.data} argument. In the latter
#' case, the arguments \code{L1}, \code{L2}, \code{T1} and \code{T2} must be
#' 1-length \code{character} vectors and they will serve to indicate the column
#' names of the corresponding variables. If only one value is specified for
#' \code{input.data} and any of the other arguments is NULL, a default name
#' equal to the variable name will be assigned (e.g. \code{L1 <- “L1”}).
#'
#' \code{na_action = "force"} should be used carefully, as it is not always due
#' to bootstrap data selection factors, but also to an inadequate selection of
#' the estimation parameters that the \code{NA} values are obtained. Also, the
#' search time may depend on the size of the input set, if you have many
#' thousands of individuals or if (in addition) the value of \code{nresamp} is
#' high, it is possible that the function will take a long time before obtaining
#' complete results. \code{time_lim} avoids falling into an infinite loop by
#' limiting the time used by this process to 5 minutes, but this value is
#' referential and may be insufficient due to the factors mentioned above.
#'
#' @references \itemize{
#'  \item Efron, B., & Tibshirani, R., 1986. Bootstrap methods for standard
#'  errors, confidence intervals, and other measures of statistical accuracy.
#'  Statistical Science, 54-75.
#'  \item Francis, R.I.C.C., 1988. Maximum likelihood estimation of growth and
#'  growth variability from tagging data. New Zealand Journal of Marine and
#'  Freshwater Research, 22, p.42-51.
#'  \item Pauly, D., 1981. The relationship between gill surface area and growth
#'  performance in fish: a generalization of von Bertalanffy's theory of growth.
#'  Meeresforsch. 28:205-211.
#'  \item Schwamborn, R., Mildenberger, T. K., & Taylor, M. H., 2019. Assessing
#'  sources of uncertainty in length-based estimates of body growth in
#'  populations of fishes and macroinvertebrates with bootstrapped ELEFAN.
#'  Ecological Modelling, 393, 37-51.
#'  \item Schwamborn, R. & Schwamborn, D. F. M. C. Growth and mortality of the
#'  endangered land crab \emph{Cardisoma guanhumi} assessed through tagging with
#'  PITs and novel bootstrapped methods. Pan-American Journal of Aquatic
#'  Sciences, 16(1): 57-78.
#'  \item von Bertalanffy, L., 1938. A quantitative theory of organic growth.
#'  Human Biology 10, 181-213.
#' }
#'
#'
#' @return A \code{data.frame} of fitted VBGF parameters (columns) by resampling
#' (rows). It includes a column (\code{seed}) with seed values set prior to each
#' resampling call.
#'
#' @export
#'
#' @examples
#' # Load example DB from fishmethods package
#' data(bonito, package = "fishmethods")
#'
#' # Run the example cited on ?grotag
#' fishmethods::grotag(L1 = bonito$L1,
#'                     L2 = bonito$L2,
#'                     T1 = bonito$T1,
#'                     T2 = bonito$T2,
#'                     alpha   = 35, beta = 55,
#'                     design  = list(nu = 1, m = 1,p = 1, sea = 1),
#'                     stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.2, u = 0.4, w = 0.4),
#'                     upper   = list(sigma = 5, nu = 1, m = 2, p = 0.5, u = 1, w = 1),
#'                     lower   = list(sigma = 0, nu = 0, m = -2, p = 0.0, u = 0, w = 0),
#'                     control = list(maxit = 1e4))
#'
#' # Run the example using grotag_boot
#' res <- grotag_boot(L1 = bonito$L1,
#'                    L2 = bonito$L2,
#'                    T1 = bonito$T1,
#'                    T2 = bonito$T2,
#'                    alpha   = 35, beta = 55,
#'                    design  = list(nu = 1, m = 1,p = 1, sea = 1),
#'                    stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.2, u = 0.4, w = 0.4),
#'                    upper   = list(sigma = 5, nu = 1, m = 2, p = 0.5, u = 1, w = 1),
#'                    lower   = list(sigma = 0, nu = 0, m = -2, p = 0.0, u = 0, w = 0),
#'                    control = list(maxit = 1e4),
#'                    nresamp = 3, na_action = "narm")
#'
#' res
grotag_boot <- function(L1 = NULL, L2 = NULL, T1 = NULL, T2 = NULL,
                        alpha = NULL, beta = NULL,
                        design = list(nu = 0, m = 0, p = 0, sea = 0),
                        stvalue = list(sigma = 0.9, nu = 0.4, m = -1,
                                       p = 0.1, u = 0.4, w = 0.4),
                        upper = list(sigma = 5, nu = 1, m = 2,
                                     p = 1, u = 1, w = 1),
                        lower = list(sigma = 0, nu = 0, m = -2,
                                     p = 0, u = 0, w = 0),
                        gestimate = TRUE, st.ga = NULL,
                        st.gb = NULL, st.galow = NULL,
                        st.gaup = NULL, st.gblow = NULL,
                        st.gbup = NULL, control = list(maxit = 10000),
                        input.data = NULL,
                        seed = NULL, nresamp = 200,
                        na_action = c("nothing", "narm", "force"),
                        time_lim = 5*60){

  # Tolowerize na_action value and take just the 1st element
  na_action <- tolower(na_action)[1] |>

    # ...removing any posible punctuation character
    gsub(pattern = "[[:punct:]]", replacement = "")

  # Check na_action value and send an error message if correspond
  if(!na_action %in% c("nothing", "narm", "force")){
    stop("'na_action' value must be whether 'nothing', 'narm' or 'force'. See ?grotag_boot")
  }

  # If input.data is NULL
  if(is.null(input.data)){

    # Combine the vectors L1, L2, T1 and T2 on a list
    input.data <- list(L1 = L1, L2 = L2, T1 = T1, T2 = T2) |>

      # ...and then into a data.frame
      do.call(what = cbind.data.frame)
  }else{

    # If some of the arguments are NULL, set it with a default name
    if(is.null(L1)) L1 <- "L1"
    if(is.null(L2)) L2 <- "L2"
    if(is.null(T1)) T1 <- "T1"
    if(is.null(T2)) T2 <- "T2"

    # Then, extract that names from the input.data
    input.data <- input.data[,c(L1, L2, T1, T2)]
  }

  # Set seed
  if(is.null(seed)) seed <- as.numeric(Sys.time())

  # Empty results list
  res <- list()

  # Initialize x = 0
  x <- 0

  # Catch the time just before to start the loop
  time_0 <- Sys.time()

  # Start a while loop that will only run as long as the length of res is less
  # than the value specified in nresamp
  while(length(res) < nresamp){
    # Increasing x
    x <- x + 1

    # Run engine function
    out <- tryCatch({
      grotag_internal(input.data = input.data, x = x, seed = seed,
                      alpha = alpha, beta = beta,
                      design = design, stvalue = stvalue,
                      upper = upper, lower = lower,
                      gestimate = gestimate, st.ga = st.ga,
                      st.gb = st.gb, st.galow = st.galow,
                      st.gaup = st.gaup, st.gblow = st.gblow,
                      st.gbup = st.gbup, control = control)
    }, error = \(e){NULL}) |> suppressWarnings()

    # If result is not NULL or if the na_action is set as 'nothing'
    outnull <- is.null(out)
    if(!outnull || na_action == "nothing"){
      # ...adding output as is
      res <- c(res, list(out))
    }else{
      # ...otherwise, if na_action is 'force'
      if(na_action == "force"){
        # Check time diff
        t_diff <- difftime(time1 = Sys.time(), time2 = time_0, units = "secs")

        # If the current time is higher than the limit, force to break out the
        # loop
        if(t_diff >= time_lim) break
      }else{
        # ...otherwise (if na_action if 'narm'), just break out if 'x' reaches
        # the 'nresamp' value
        if(x == nresamp) break
      }
    }
  }


  # Compile results in a data.frame
  res <- lapply(X = res, FUN = \(x) if(is.null(x$res)) rep(NA, 6) else x$res) |>

    do.call(what = rbind) |>

    as.data.frame() |>

    cbind.data.frame(sapply(X = res, FUN = \(x) if(is.null(x$time)) NA else x$time))

  # Define column names
  colnames(res) <- c( "Linf", "K", "PhiL", "u", "w", "seed", "time")

  # Set class of output object
  class(res) <- c("grotagBoot", class(res))

  res
}

grotag_internal <- function(input.data, x, seed,
                            alpha, beta, design, stvalue,
                            upper, lower, gestimate,
                            st.ga, st.gb,
                            st.galow, st.gaup,
                            st.gblow, st.gbup,
                            control){

  set.seed(seed + x)

  # Resampling rows
  samplerows <- sample(x = seq(nrow(input.data)),
                       size = nrow(input.data),
                       replace = TRUE)

  input.data_p <- input.data[samplerows,]

  # Apply grotag function
  fitboot <- grotag(L1 = input.data_p$L1,
                    L2 = input.data_p$L2,
                    T1 = input.data_p$T1,
                    T2 = input.data_p$T2,
                    alpha = alpha, beta = beta,
                    design = design, stvalue = stvalue,
                    upper = upper, lower = lower,
                    gestimate = gestimate, st.ga = st.ga,
                    st.gb = st.gb, st.galow = st.galow,
                    st.gaup = st.gaup, st.gblow = st.gblow,
                    st.gbup = st.gbup, control = control)

  Linf <- as.numeric(fitboot$VBparms$Estimate[2])
  K <- as.numeric(fitboot$VBparms$Estimate[3])

  list(res = c(Linf = Linf,
               K    = K,
               PhiL = log10(K) + 2 * log10(Linf),
               u    = fitboot$table$Estimate[4],
               w    = fitboot$table$Estimate[5],
               seed = seed + x),
       time = as.character(Sys.time()))
}
