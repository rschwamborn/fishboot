#' Univariate kernal density estimate plot of VBGF parameter from bootstrapping results
#'
#' @param res list of class `lfqBoot`
#' @param CI numeric. Confidence interval level (Default: 95)
#' @param use_hist logical. Plot histogram in addition to smoothed kernel density.
#' @param nbreaks numeric. Number of breaks in the histogram.
#' @param mar vector. Inner margins settings. See `?par`.
#' @param oma vector. Outer margins settings.See `?par`.
#' @param mgp vector. See `?par`.
#' @param tcl vector. See `?par`.
#' @param cex numeric. See `?par`.
#' @param ... Additional arguments passed to `par`.
#'
#' @return plot
#' @export
#'
#' @examples
#'
#' data(alba_boot)
#' univariate_density(alba_boot)
#'
#'
univariate_density <- function(res, CI=95, use_hist = FALSE, nbreaks = 10,
  mar = c(1.5,2,2,0), oma = c(1.5,0,0,0.5),
  mgp = c(2,0.5,0), tcl = -0.25, cex = 1,
  ...
){

  res <- res$bootRaw

  op <- par(no.readonly = TRUE)
  par(
    # mfcol = c(floor(sqrt(ncol(res))), ceiling(sqrt(ncol(res)))),
    mfcol = c(1, ncol(res)),
    mar = mar, oma = oma,
    mgp = mgp, tcl = tcl, cex = cex, ...
  )

  VARS <- list(
    Linf = expression(italic(L)[infinity]),
    K = expression(italic(K)),
    t_anchor = expression(italic(t)[anchor]),
    C = expression(italic(C)),
    ts = expression(italic(t)[s]),
    phiL = expression(paste(phi,"'"))
  )

  # univariate plots
  for(i in seq(ncol(res))){
    x <- ks::kde(res[,i])

    h = hist(res[,i], plot=FALSE, breaks = nbreaks)

    xlim <- c(0, max(x$estimate))
    if(use_hist){
      xlim <- range(c(xlim, max(h$density)))
    }
    xlim <- xlim * c(0,1.1)

    plot(x$estimate, x$eval.points, t="n",
      xlim = xlim,
      xaxs = "i",
      ylab="", xlab="", col=1, lty=1
    )
    usr <- par()$usr

    CItxt <- paste0(round(100-CI), "%")
    inCI <- rle( x$estimate > x$cont[CItxt] )
    start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
    end.idx <- cumsum(inCI$lengths)
    limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])

    in1 <- which(x$estimate > x$cont["99%"])
    mean1 <- mean(x$eval.points[in1])

    if(use_hist){
      rect(
        xleft = 0, ybottom = h$breaks[-length(h$breaks)],
        xright = h$density, ytop = h$breaks[-1],
        col = "grey90", border = "grey50"
      )
    }else{
      for(j in seq(inCI$lengths)){
        if(inCI$values[j]){
          polygon(
            y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
            x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
            col = "grey90", #col = rgb(0,1,0,0.2),
            border = NA, lty = 3, lwd = 1
          )
        }
      }
    }

    # abline(v = x$cont[CItxt], lty=2, col="grey50")
    lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")

    # rug
    segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)

    # range of CI
    abline(h = limCI, lty = 1, lwd=1, col = 1)
    text(y =c(limCI), x = mean(usr[1:2]),
      labels = paste(sprintf("%.2f", round(c(limCI),2))),
      pos = c(1,3), offset = 0.25, col = 1
    )
    abline(h = mean1, lty = 1, lwd=1, col = 1)
    text(y =  mean1, x = mean(usr[1:2]),
      labels = sprintf("%.2f", round(mean1,2)),
      pos = 3,
      offset = 0.25, col = 1
    )

    varlab <- VARS[[match(names(res)[i], names(VARS))]]
    mtext(varlab, line=0.25, side=3)

    box()
  }
  mtext("Density", side = 1, line = 0, outer = TRUE)
  par(op)
}
