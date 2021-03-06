#' Link/K scatterplot of bootstrapping results
#'
#' @param res list of class `lfqBoot`
#' @param Linf.breaks vector. Breaks for Linf histogram.
#' @param K.breaks vector. Breaks for K histogram.
#' @param gridsize vector. 2 values for defining the resolution of the grid
#' @param H object from \code{\link[ks]{Hpi}} (Default: `H = ks::Hpi(res[,c("Linf", "K")])`)
#' @param shading logical. Should 2d field of density estimates be colored with colors
#'   specified by `shading.cols` argument (Default: `shading = TRUE`).
#' @param shading.cols vector. Colors for background shading of 2d field of density estimates
#'   (Default: `shading.cols = colorRampPalette(c("white", blues9))(50)`).
#' @param dens.contour logical. Should contour lines be added (Default: `dens.contour = TRUE`)
#' @param probs vector. Density probability cutoffs to be plotted by contours when
#'   `dens.contour = TRUE` (Default: `probs = c(25,50,75,95)`).
#' @param phi.contour logical. Should phi prime isolines be displayed
#'   (Default: `phi.contour = FALSE`)
#' @param phi.levels vector. Phi prime values to display when `phi.contour = TRUE`
#'   (Default: `phi.levels = NULL`). When not defined (`phi.levels = NULL`), the default
#'   levels are chosen automatically by the \code{\link[graphics]{contour}} function.
#' @param phi.contour.col vector. Color to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.lty vector. Line type to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.lwd vector. Line width to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.labcex vector. Labels for the contour lines.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param pt.pch pch value to use for resampling points
#' @param pt.col color to use for resampling points
#' @param pt.cex size value to use for resampling points
#' @param pt.bg background color to use for resampling points
#' @param xlab label for x-axis
#' @param ylab lavel for y axis
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#'
#' @return plot
#' @export
#'
#' @examples
#' data(alba_boot)
#' LinfK_scatterhist(alba_boot)
#'
LinfK_scatterhist = function(
  res, Linf.breaks = "Sturges", K.breaks = "Sturges",
  gridsize = rep(151, 2), H = NULL,
  shading = TRUE, shading.cols = colorRampPalette(c("white", blues9))(50),
  dens.contour = TRUE, probs = c(25,50,75,95),
  phi.contour = FALSE, phi.levels = NULL,
  phi.contour.col = 8, phi.contour.lty = 2, phi.contour.lwd = 1, phi.contour.labcex = 0.75,
  pt.pch = 16, pt.col = adjustcolor(1, 0.25), pt.cex = 0.5, pt.bg = 4,
  xlab = expression(italic("L")[infinity]), ylab = expression(italic("K")),
  xlim = NULL, ylim = NULL
){

  res <- res$bootRaw

  # Called internally
  add_phiprime <- function(gridsize = 20, ...){
    usr <- par()$usr
    if(par()$`xlog`){usr[1:2] <- 10^usr[1:2]}
    if(par()$`ylog`){usr[3:4] <- 10^usr[3:4]}
    usr <- replace(usr, which(usr < 0), 0)
    Linf <- seq(usr[1], usr[2], length.out = gridsize)
    K <- seq(usr[3], usr[4], length.out = gridsize)
    Linf <- Linf[which(Linf>=0)]
    K <- K[which(K>=0)]
    grd <- expand.grid(
      Linf = Linf,
      K = K
    )
    grd$phiL <- log10(grd$K) + 2 * log10(grd$Linf)

    M <- list(x = Linf, y = K, z = matrix(grd$phiL, nrow = gridsize, ncol = gridsize))
    contour(x = M, add = TRUE, ...)
  }


  zones <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  op <- par(no.readonly = TRUE)
  nf <- layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5), respect = FALSE)
  # layout.show(nf)
  par(cex = 1)

  # histogram data
  xhist = hist(res[,"Linf"], plot=FALSE, breaks = Linf.breaks)
  yhist = hist(res[,"K"], plot=FALSE, breaks = K.breaks)
  top = max(c(xhist$counts, yhist$counts))

  # density estimation
  par(mar=c(3,3,0,0), mgp = c(2,0.5,0), tcl = -0.25)
  if(is.null(H)){
    H <- ks::Hpi(res[,c("Linf", "K")])
  }
  kk <- ks::kde(
    x = res[,c("Linf", "K")], gridsize = gridsize, H = H
  )

  # plot limits (defaults are eval.points of kde)
  if(is.null(xlim)){xlim <- range(kk$eval.points[[1]])}
  if(is.null(ylim)){ylim <- range(kk$eval.points[[2]])}

  # 2d density plot
  image(
    x = kk$eval.points[[1]], y = kk$eval.points[[2]], z = kk$estimate,
    col = if(shading){shading.cols}else{NA},
    xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim
  )
  if(phi.contour){
    if(is.null(phi.levels)){
      add_phiprime(
        col = phi.contour.col, lty = phi.contour.lty, lwd = phi.contour.lwd,
        labcex = phi.contour.labcex
      )
    }else{
      add_phiprime(
        levels = phi.levels,
        col = phi.contour.col, lty = phi.contour.lty, lwd = phi.contour.lwd,
        labcex = phi.contour.labcex
      )
    }
  }
  points(kk$x[,1], kk$x[,2], pch = pt.pch, cex = pt.cex, col = pt.col, bg = pt.bg)
  if(dens.contour){
    plot(kk, type = "slice", add = TRUE, cont = probs)
  }
  box()

  # x histogram
  par(mar=c(0,3,0.5,0))
  plot(xhist$mids, xhist$counts, axes=FALSE, xlab="", ylab="", t="n", ylim=c(0, top), xlim = xlim, yaxs="i", xaxs="i")
  rect(
    xleft = xhist$breaks[-length(xhist$breaks)], ybottom = 0,
    xright = xhist$breaks[-1], ytop = xhist$counts,
    col = 8, border = 1
  )

  # y histogram
  par(mar=c(3,0,0,0.5))
  plot(yhist$counts, yhist$mids, axes=FALSE, xlab="", ylab="", t="n", xlim=c(0, top), ylim = ylim, yaxs="i", xaxs="i")
  rect(
    xleft = 0, ybottom = yhist$breaks[-length(yhist$breaks)],
    xright = yhist$counts, ytop = yhist$breaks[-1],
    col = 8, border = 1
  )

  par(op)
}

