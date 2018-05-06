#' Link/K scatterplot of bootstrapping results
#'
#' @param res list of class `lfqBoot`
#' @param Linf.breaks vector. Breaks for Linf histogram.
#' @param K.breaks vector. Breaks for K histogram.
#' @param gridsize vector. 2 values for defining the resolution of the grid
#' @param H object from `ks::H`
#' @param shading xx
#' @param shading.cols xx
#' @param dens.contour xx
#' @param probs xx
#' @param phi.contour xx
#' @param phi.levels xx
#' @param phi.contour.col xx
#' @param phi.contour.lty xx
#' @param phi.contour.lwd xx
#' @param phi.contour.labcex xx
#' @param pt.pch xx
#' @param pt.col xx
#' @param pt.cex xx
#' @param pt.bg xx
#' @param xlab xx
#' @param ylab xx
#' @param xlim xx
#' @param ylim xx
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
  gridsize = rep(151, 2), H = ks::Hpi(res[,c("Linf", "K")]),
  shading = TRUE, shading.cols = colorRampPalette(c("white", blues9))(50),
  dens.contour = TRUE, probs = c(25,50,75,95),
  phi.contour = FALSE, phi.levels = NULL,
  phi.contour.col = 8, phi.contour.lty = 2, phi.contour.lwd = 1, phi.contour.labcex = 0.75,
  pt.pch = 16, pt.col = rgb(0,0,0,0.25), pt.cex = 0.5, pt.bg = 4,
  xlab=expression(italic("L")[infinity]), ylab=expression(italic("K")),
  xlim = NULL, ylim = NULL
){

  res <- res$bootRaw

  # Called internally
  add_phiprime <- function(gridsize = 20, ...){
    usr <- par()$usr
    Linf = seq(usr[1], usr[2], length.out = gridsize)
    K = seq(usr[3], usr[4], length.out = gridsize)
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


