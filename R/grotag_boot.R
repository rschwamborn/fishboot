#' @title Resampling of growth increment data (\eqn{\frac{dL}{dt}}) from
#' mark-recapture (tagging) studies
#'
#' @description This function resamples the \code{input.data} data by rows (i.e.,
#' by recapture date) several times (\code{nresamp} times, default:
#' \code{nresamp = 200}). Then, a VBGF curve is fitted to each resapled data set.
#' The output (a \code{list} of class \code{lfqBoot}) will store results (e.g.,
#' VGBGF function parameters K and Linf) in a \code{data.frame} accessible
#' through \code{$bootRaw}. The \code{$bootRaw} table also includes the growth
#' performance index \strong{Phi'}, seaonal parameters \strong{u} and \strong{w}
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
#' @param seed seed value for random number reproducibility.
#' @param nresamp \code{numeric}; the number of permutations to run (Default:
#' \code{nresamp = 200}).
#' @param na_action \code{character} that defines the action that the function
#' will execute if there is a row with NA:
#' \itemize{
#'  \item \code{nothing}: the function will return the results including the NAs
#'  (default).
#'  \item \code{narm}: after having the results, it will only returns the rows
#'  without NAs. See Details.
#'  \item \code{force}: The function will start an iterative process changing
#'  the internal \code{seed} values until it fulfills the \code{nresamp}. It
#'  works just together \code{time_lim} argument. See Details.
#' }
#' @param time_lim If \code{na_action = "force"}, it defines the maximum time
#' (in seconds) that the function will last resampling until it achieves a
#' result output with no-NaN rows.
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
#'                     alpha = 35, beta = 55,
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
#'                    alpha = 35, beta = 55,
#'                    design  = list(nu = 1, m = 1,p = 1, sea = 1),
#'                    stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.2, u = 0.4, w = 0.4),
#'                    upper   = list(sigma = 5, nu = 1, m = 2, p = 0.5, u = 1, w = 1),
#'                    lower   = list(sigma = 0, nu = 0, m = -2, p = 0.0, u = 0, w = 0),
#'                    control = list(maxit = 1e4),
#'                    nresamp = 10, na_action = "narm")
#'
#' LinfK_scatterhist(res = res)
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
                        seed = as.numeric(Sys.time()), nresamp = 200,
                        na_action = "nothing", time_lim = 5*60){

  if(is.null(input.data)){
    input.data <- list(L1 = L1, L2 = L2, T1 = T1, T2 = T2) |>

      do.call(what = cbind.data.frame)
  }else{
    if(is.null(L1)) L1 <- "L1"
    if(is.null(L2)) L2 <- "L2"
    if(is.null(T1)) T1 <- "T1"
    if(is.null(T2)) T2 <- "T2"

    input.data <- input.data[,c(L1, L2, T1, T2)]
  }

  # Empty results list
  res <- list()
  na_action <- tolower(na_action)[1]
  time_0 <- Sys.time()
  x <- 0
  while(length(res) < nresamp){
    x <- x + 1

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


    outnull <- is.null(out)
    if(!outnull || na_action == "nothing"){
      # Adding output as is
      res <- c(res, list(out))
    }else{
      if(na_action == "force"){
        # Check time diff
        t_diff <- difftime(time1 = Sys.time(), time2 = time_0, units = "secs")
        if(t_diff >= time_lim) break
      }else{
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

  # Define class
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
