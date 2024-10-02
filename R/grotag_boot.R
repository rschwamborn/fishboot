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
#' @param ... Extra arguments passed to \link[fishmethods]{grotag}:
#' \code{alpha}, \code{beta}, \code{design}, \code{stvalue}, \code{upper},
#' \code{lower}, \code{st.ga}, \code{st.gb}, \code{st.galow}, \code{st.gaup},
#' \code{st.gblow}, \code{st.gbup} and \code{control}. It is important for the
#' user to indicate at least values of \code{alpha} and \code{beta}.
#' @param input.data A growth increment object of the class \code{data.frame}.
#' @param nresamp \code{numeric}; the number of permutations to run (Default:
#' \code{nresamp = 200}).
#' @param seed seed value for random number reproducibility.
#' @param no_cores positive integer. If \code{no_cores} > 1, a 'parallel' package
#' cluster with that many cores is created.
#' @param cc.ony \code{logical} Do you want the output object to omit results
#' with failed estimates (with NAs)? Remember that if you define
#' \code{cc.only = TRUE} the number of elements returned may be less than
#' specified in \code{nresamp}.
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
#' grotag_boot(L1 = bonito$L1,
#'             L2 = bonito$L2,
#'             T1 = bonito$T1,
#'             T2 = bonito$T2,
#'             alpha = 35, beta = 55,
#'             design  = list(nu = 1, m = 1,p = 1, sea = 1),
#'             stvalue = list(sigma = 0.9, nu = 0.4, m = -1, p = 0.2, u = 0.4, w = 0.4),
#'             upper   = list(sigma = 5, nu = 1, m = 2, p = 0.5, u = 1, w = 1),
#'             lower   = list(sigma = 0, nu = 0, m = -2, p = 0.0, u = 0, w = 0),
#'             control = list(maxit = 1e4),
#'             nresamp = 10)
grotag_boot <- function(L1 = NULL, L2 = NULL, T1 = NULL, T2 = NULL, ...,
                        input.data = NULL, nresamp = 200, seed = 123,
                        no_cores = 1, cc.ony = FALSE){

  if(is.null(input.data)){
    input.data <- do.call(what = cbind.data.frame,
                          args = list(L1, L2, T1, T2))

    colnames(input.data) <- c("L1", "L2", "T1", "T2")
  }else{
    if(is.null(L1)) L1 <- "L1"
    if(is.null(L2)) L2 <- "L2"
    if(is.null(T1)) T1 <- "T1"
    if(is.null(T2)) T2 <- "T2"

    input.data <- input.data[,c(L1, L2, T1, T2)]
  }

  if(no_cores > 1){

    # Registering cluster
    cl <- makeCluster(spec = no_cores)
    registerDoParallel(cl = cl)

    # Run multithread process
    res <- foreach(x = seq(nresamp), .inorder = FALSE, .packages = "fishmethods") %dopar% {
      tryCatch({
        grotag_internal(input.data = input.data, x = x, seed = seed, ...)
      }, error = \(e){NULL})
    }

    # Finish cluster
    stopCluster(cl)
  }else{

    # Empty results list
    res <- vector("list", nresamp)

    for(x in seq(res)){
      tryCatch({
        res[[x]] <- grotag_internal(input.data = input.data, x = x, seed = seed, ...)
      }, error = \(e){NULL}) |> suppressWarnings()
    }
  }

  if(all(sapply(X = res, FUN = is.null))){
    warning("All the resamples got fitting problems with grotag. See ?grotag_boot.")

    return(NULL)
  }

  # Compile results in a data.frame
  res <- lapply(X = res, FUN = \(x) if(is.null(x$res)) rep(NA, 6) else x$res) |>

    do.call(what = rbind) |>

    as.data.frame() |>

    cbind.data.frame(sapply(X = res, FUN = \(x) if(is.null(x$time)) NA else x$time))

  # Define column names
  colnames(res) <- c( "Linf", "K", "PhiL", "u", "w", "seed", "time")

  # Remove rows with NAs
  if(isTRUE(cc.ony)) res <- res[complete.cases(res),]

  # Define class
  class(res) <- c("grotagBoot", class(res))

  res
}

grotag_internal <- function(input.data, x, seed, ...){

  seed <- seed + x
  set.seed(seed)

  # Resampling rows
  samplerows <- sample(x = seq(nrow(input.data)),
                       size = nrow(input.data),
                       replace = TRUE)

  input.data_p <- input.data[samplerows,]

  # Apply grotag function
  fitboot <- formals(grotag) |>

    modifyList(val = append(as.list(input.data_p), list(...))) |>

    as.list() |>

    do.call(what = grotag)

  Linf <- as.numeric(fitboot$VBparms$Estimate[2])
  K <- as.numeric(fitboot$VBparms$Estimate[3])

  list(res = c(Linf = Linf,
               K    = K,
               PhiL = log10(K) + 2 * log10(Linf),
               u    = fitboot$table$Estimate[4],
               w    = fitboot$table$Estimate[5],
               seed = seed),
       time = as.character(Sys.time()))
}
