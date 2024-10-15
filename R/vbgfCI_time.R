#' @title VBGF plot and CI
#'
#' @description Pending...
#'
#' @param res Object with \eqn{L_{inf}}, \eqn{K} and \eqn{t_0}, it could be a
#' \code{data.frame}, a \code{tbl_df}, a \code{list} or a \code{lfqBoot} object.
#' @param CI \code{numeric}. Confidence interval in \% (default: 95).
#' @param agemax \code{numeric} values indicating the maximum number of years to
#' project.
#' @param plot \code{logical}. If \code{TRUE} (default), a plot is returned,
#' otherwise just a list with levels \code{limCI}, \code{inCI}, \code{density}
#' and \code{max_dens}. See Value for a detailed description of each one.
#' @param add_legend \code{logical}. Should CI and max. density legend be added
#' (Default: `add_legend = TRUE`).
#' @param add_max_dens_legend logical. Should maximum density line be added
#' (Default: `add_max_dens_legend = TRUE`).
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param perm.col,perm.lwd Color and width for each resample estimate line.
#' @param ci.col,ci.lty,ci.lwd Color, type and width for CI line.
#' @param maxd.col,maxd.lty,maxd.lwd Color, type and width for maximum density line..
#' @param ... Extra arguments passed to the main plot function.
#'
#' @details
#' Additional details...
#'
#'
#'
#' @return A \code{list} containing:
#' \describe{
#'  \item{\code{$limCI}}{A \code{data.frame} with CI limits by time.}
#'  \item{\code{$inCI}}{A \code{data.frame} with logical values defining whether
#'  bootstrapping samples are within each of the defined CIs.}
#'  \item{\code{$density}}{The multivariate kernel density estimates for each
#'  sample.}
#'  \item{\code{$max_dens}}{A \code{list} with the VBGF parameter combination
#'  having the maximum density estimate.}
#' }
#'
#' @export
#'
#' @examples
#' data(alba_boot) # lfqBoot object
#' vbgfCI_time(res = alba_boot)
#'
#' vbgfCI_time(res = alba_boot, CI = c(50, 95),
#'             ci.col = c("red", "orange"))
#'
#' data(bonito_boot) # grotagBoot object
#' LinfK_scatterhist(res = bonito_boot)
vbgfCI_time <- function(res, CI = 95, agemax = NULL, plot = TRUE,
                        add_legend = TRUE, add_max_dens_legend = TRUE,
                        xlab = "Relative time", ylab = "Length",
                        perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
                        ci.col = "black", ci.lty = 2, ci.lwd = 1,
                        maxd.col = "black", maxd.lty = 1, maxd.lwd = 2, ...){

  # Extract values of Linf, K and t_anchor from res
  x <- get_LinfKtanchor(x = res)

  if(is.null(agemax)){
    agemax <- max(ceiling((1/-x$K)*log(1-((x$Linf*0.95)/x$Linf))))
  }

  # remove phi' if included in res
  # phiLcol <- which(names(res) == "phiL")
  # x <- if(length(phiLcol > 0)) res[,-phiLcol] else res

  # First a fitting round to look for shifts in t_anchor
  age <- seq(from = 0, to = agemax, by = 0.01)
  Lt0 <- matrix(data = NA, nrow = length(age), ncol = nrow(x))
  Lt_minus <- matrix(data = NA, nrow = length(age), ncol = nrow(x))
  Lt_plus <- matrix(data = NA, nrow = length(age), ncol = nrow(x))

  for(i in seq(ncol(Lt0))){
    par0 <- par_minus <- par_plus <- as.list(x[i,])
    par_minus$t_anchor <- par0$t_anchor - 1
    par_plus$t_anchor <- par0$t_anchor + 1

    Lt0[,i] <- VBGF(param = par0, t = age)
    Lt_minus[,i] <- VBGF(param = par_minus, t = age)
    Lt_plus[,i] <- VBGF(param = par_plus, t = age)
  }

  # replace negative lengths with NA
  Lt0 <- replace(x = Lt0, list = Lt0 < 0, values = NA)
  Lt_minus <- replace(x = Lt_minus, list = Lt_minus < 0, values = NA)
  Lt_plus <- replace(x = Lt_plus, list = Lt_plus < 0, values = NA)

  # determine if a positive or negative shift improves overall covariance for each permutation
  cov0 <- cov(x = Lt0, use = "pair")
  shift <- 0 * seq(ncol(Lt0))
  for(i in seq(ncol(Lt0))){
    cov_minus <- cov(x = Lt_minus[,i], y = Lt0, use = "pair")
    cov_plus <- cov(x = Lt_plus[,i], y = Lt0, use = "pair")
    shift[i] <- c(0, -1, 1)[which.max(c(sum(cov0[i,]), sum(cov_minus), sum(cov_plus)))]
  }

  # fixed predictions
  agenew <- seq(from = 0 + min(shift),
                to = agemax + max(shift),
                by = 0.01)
  Lt <- matrix(data = NA, nrow = length(agenew), ncol = nrow(x))

  for(i in seq(ncol(Lt))){
    par0 <- as.list(x[i,])
    par0$t_anchor <- par0$t_anchor + shift[i]

    Lt[,i] <- VBGF(param = par0, t = agenew)
  }

  # multivariate kernel density estimate for max. density
  H <- Hpi(x = x, nstage = 1)
  fhat <- kde(x = x, H = H, eval.points = x)

  # predict density estimate of original data
  x$estimate <- fhat$estimate

  # maximum density
  max_dens <- fhat$eval.points[which.max(fhat$estimate),]

  if(isTRUE(plot)){
    # Empty canvas
    plot(x = 1, y = 1, type = "n",
         xlim = c(min(x$t_anchor + shift), max(agenew)) + c(-0.1, 0),
         ylim = c(0, max(Lt)*1.05),
         xlab = xlab, ylab = ylab,
         xaxs = "i", yaxs = "i", ...)

    # Drawing resampling lines
    apply(X = Lt, MARGIN = 2, FUN = lines, x = agenew,
          col = perm.col, lwd = perm.lwd)

    box()
  }

  # Determine which resamples are in the CI
  limCI <- vector(mode = "list", length(CI))
  inCI <- vector(mode = "list", length(CI))
  names(inCI) <- names(limCI) <- paste0("CI", CI)

  for(j in seq(limCI)){
    inCI[[j]] <- x$estimate > quantile(x = x$estimate, probs = (100 - CI[j])/100)
    limCI[[j]] <- data.frame(t = agenew, min = NA, max = NA)

    for(i in seq(agenew)){
      limCI[[j]]$min[i] <- min(Lt[i, which(inCI[[j]]) ], na.rm = TRUE)
      limCI[[j]]$max[i] <- max(Lt[i, which(inCI[[j]]) ], na.rm = TRUE)
    }

    if(isTRUE(plot)){
      # Expand CI line attributes to length of CI (if needed, values are recycled)
      ci.col <- rep_len(x = ci.col, length.out = length(CI))
      ci.lty <- rep_len(x = ci.lty, length.out = length(CI))
      ci.lwd <- rep_len(x = ci.lwd, length.out = length(CI))

      # Drawing CI lines
      apply(X = limCI[[j]][,c("min", "max")], MARGIN = 2, FUN = lines,
            x = limCI[[j]]$t, col = ci.col[j], lwd = ci.lwd[j], lty = ci.lty[j])
    }
  }

  inCI <- as.data.frame(inCI)

  if(isTRUE(plot)){
    lines(x = agenew, y = VBGF(param = as.list(max_dens), t = agenew),
          col = maxd.col, lwd = maxd.lwd, lty = maxd.lty)

    if(add_legend){
      legend("bottomright", legend = c(paste0("CI = ", CI, "%"), "Max. Dens."),
             bty = "n", col = c(ci.col, maxd.col), lty = c(ci.lty, maxd.lty),
             lwd = c(ci.lwd, maxd.lwd))
    }

    if(add_max_dens_legend){
      legend("topleft",
             legend = sprintf(fmt = "%s = %.2f", names(max_dens), max_dens),
             bty = "n", title = "Max. Dens.\nparameters:",
             inset = c(0, 0.1))
    }
  }

  out <- list(limCI = limCI,
              inCI = inCI,
              density = x$estimate,
              max_dens = as.list(max_dens))

  # Output
  if(isTRUE(plot)) invisible(out) else out
}

get_LinfKtanchor <- function(x){

  # Define functions for extracting Linf, K and t_anchor values from several
  # objects (classes)
  getData <- list("tbl_df" = \(obj) obj,
                  "tbl" = \(obj) obj,
                  "data.frame" = \(obj) obj,
                  # "grotagBoot" = \(obj) data.frame(obj, t_anchor = 0),
                  "grotagBoot" = \(obj){
                    data.frame(obj, t_anchor = rnorm(n = nrow(obj),
                                                     mean = 0,
                                                     sd = 1e-4))
                  },
                  "lfqBoot" = \(obj) obj$bootRaw)

  # Searching the class of 'x' within the getData definitions
  index <- match(x = class(x), table = names(getData))[1]

  # If there is not a defined way (function) to extract Linf, K and t_anchor
  # variables, return an error msg
  if(is.na(index)){
    sprintf(fmt = "Internal funtion 'get_LinfKtanchor' do not know how to extract Linf or K from a '%s' object.",
            class(x)) |> stop()
  }

  # Extracting data an coerce to data.frame
  out <- getData[[index]](x) |> as.data.frame()

  # Set the column names as lowercase
  colnames(out) <- tolower(colnames(out))

  # Indexing Linf and K only and set standard column names
  out[,c("linf", "k", "t_anchor")] |> setNames(c("Linf", "K", "t_anchor"))
}
