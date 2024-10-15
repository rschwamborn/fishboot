#' @title Linf/K scatterplot of bootstrapping results
#'
#' @description
#' Pending...
#'
#' @param res Object of class \code{data.frame}, \code{tbl_df}, \code{lfqBoot}
#' or \code{grotagBoot}.
#' @param Linf.breaks,K.breaks Arguments passed to \link[graphics]{hist}
#' function to compute the breakpoints (argument \code{breaks}) for Linf and K
#' histograms respectively.
#' @param gridsize \code{numeric} 2-length vector specifying the resolution of
#' the grid.
#' @param H Plug-in bandwidth object from \code{\link[ks]{Hpi}} (Default:
#' \code{H = ks::Hpi(res[,c("Linf", "K")])})
#' @param shading \code{logical}. Do you want to colour 2D field of density
#' estimates? (Default \code{TRUE})
#' @param shading.cols Colors or color palette used for background shading of 2D
#' field of density estimates. No considered if \code{shading = FALSE}.
#' @param dens.contour \code{logical}. Do you want to add contour lines?
#' (\code{TRUE} by default).
#' @param probs \code{numeric} Density probability cutoffs (in %) to be plotted
#' by contours. By default \code{probs = c(25, 50, 75, 95)} and not considered
#' if \code{dens.contour = FALSE}.
#' @param phi.contour \code{logical}. Do you want to display phi prime isolines?
#' (\code{FALSE} by default)
#' @param phi.levels \code{numeric} vector that controls Phi prime values.
#' Omitted if \code{phi.contour = FALSE} and if \code{NULL} values will be
#' chosen automatically by the \link[graphics]{contour} function.
#' @param phi.contour.col,phi.contour.lty,phi.contour.lwd,phi.contour.labcex Extra
#' arguments used to control the color, line type, line width and labels size of
#' the phi prime contour isolines.
#' @param pt.pch,pt.col,pt.cex,pt.bg Extra arguments to control type, color,
#' size and background color of resampling points.
#' @param xlab,ylab Labels for X and Y axis respectively.
#' @param ... Extra arguments passed to main plot function.
#'
#' @details
#' If NULL, it will be defined as \code{colorRampPalette(c("white", blues9))(1e3)}.
#'
#' @export
#'
#' @examples
#' data(alba_boot) # lfqBoot object
#' LinfK_scatterhist(res = alba_boot)
#'
#' data(bonito_boot) # grotagBoot object
#' LinfK_scatterhist(res = bonito_boot)
LinfK_scatterhist <- function(res,
                              Linf.breaks = "Sturges", K.breaks = "Sturges",
                              gridsize = rep(151, 2), H = NULL,
                              shading = TRUE, shading.cols = NULL,
                              dens.contour = TRUE, probs = c(25, 50, 75, 95),
                              phi.contour = TRUE, phi.levels = NULL,
                              phi.contour.col = 8, phi.contour.lty = 2,
                              phi.contour.lwd = 1, phi.contour.labcex = 0.75,
                              pt.pch = 16, pt.col = adjustcolor(1, 0.25),
                              pt.cex = 0.5, pt.bg = 4,
                              xlab = expression(italic("L")[infinity]),
                              ylab = expression(italic("K")),
                              ...){

  # Extract values of Linf and K from res
  res <- get_LinfK(x = res)

  # Catch original par
  op <- par(no.readonly = TRUE)

  # Define and set order of subplots
  matrix(data = c(2, 0, 1, 3), ncol = 2, byrow = TRUE) |>

    layout(widths = c(4/5, 1/5), heights = c(1/5, 4/5), respect = FALSE)


  # density estimation
  par(mar = c(3, 3, 0, 0), mgp = c(2, 0.5, 0), tcl = -0.25, cex = 1)

  # Define default value for H
  if(is.null(H)) H <- Hpi(x = res)

  kk <- kde(x = res, gridsize = gridsize, H = H)

  # Define default plot limits
  xlim <- list(...)$xlim
  ylim <- list(...)$ylim

  if(is.null(xlim)) xlim <- c(0, max(kk$eval.points[[1]]))
  if(is.null(ylim)) ylim <- c(0, max(kk$eval.points[[2]]))

  if(is.null(shading.cols)) shading.cols <- colorRampPalette(c("white", blues9))(1e3)

  # (I) Main plot: 2D density plot
  image(x = kk$eval.points[[1]],
        y = kk$eval.points[[2]],
        z = kk$estimate,
        col = if(isTRUE(shading)) shading.cols else NA,
        xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)

  # Phi' contour
  if(isTRUE(phi.contour)){
    if(is.null(phi.levels)){
      add_phiprime(col = phi.contour.col, lty = phi.contour.lty,
                   lwd = phi.contour.lwd, labcex = phi.contour.labcex)
    }else{
      add_phiprime(levels = phi.levels,
                   col = phi.contour.col, lty = phi.contour.lty,
                   lwd = phi.contour.lwd, labcex = phi.contour.labcex)
    }
  }

  # Add density contour lines
  if(dens.contour) plot(x = kk, type = "slice", add = TRUE, cont = probs)

  # Add resampling points
  points(x = kk$x[,1], y = kk$x[,2],
         pch = pt.pch, cex = pt.cex, col = pt.col, bg = pt.bg)

  box()

  # histogram data
  xhist <- hist(x = res$Linf, plot = FALSE, breaks = Linf.breaks)
  yhist <- hist(x = res$K, plot = FALSE, breaks = K.breaks)
  top <- max(c(xhist$counts, yhist$counts))


  # (II) Linf histogram
  par(mar = c(0, 3, 0.5, 0))

  # Empty canvas
  plot(x = xhist$mids, y = xhist$counts, type = "n", axes = FALSE,
       xlab = NA, ylab = NA, ylim = c(0, top), xlim = xlim,
       yaxs = "i", xaxs = "i")

  # Drawing histogram bars
  rect(xleft = xhist$breaks[-length(xhist$breaks)], ybottom = 0,
       xright = xhist$breaks[-1], ytop = xhist$counts,
       col = 8, border = 1)


  # (III) K histogram
  par(mar = c(3, 0, 0, 0.5))

  # Empty canvas
  plot(x = yhist$counts, y = yhist$mids, axes = FALSE, type = "n",
       xlab = NA, ylab = NA, xlim = c(0, top), ylim = ylim,
       yaxs = "i", xaxs = "i")

  # Drawing histogram bars
  rect(xleft = 0, ybottom = yhist$breaks[-length(yhist$breaks)],
       xright = yhist$counts, ytop = yhist$breaks[-1],
       col = 8, border = 1)

  # Reverting changes on graphic parameters
  par(op)

  invisible()
}

add_phiprime <- function(gridsize = 20, ...){
  usr <- par()$usr

  if(par()$`xlog`){usr[1:2] <- 10^usr[1:2]}
  if(par()$`ylog`){usr[3:4] <- 10^usr[3:4]}

  usr <- replace(usr, which(usr < 0), 0)

  Linf <- seq(from = usr[1], to = usr[2], length.out = gridsize)
  Linf <- Linf[which(Linf >= 0)]

  K <- seq(from = usr[3], to = usr[4], length.out = gridsize)
  K <- K[which(K >= 0)]

  grd <- expand.grid(Linf = Linf, K = K)
  grd$phiL <- with(grd, log10(K) + 2 * log10(Linf))

  M <- list(x = Linf, y = K,
            z = matrix(data = grd$phiL, nrow = gridsize, ncol = gridsize))

  contour(x = M, add = TRUE, ...)
}

get_LinfK <- function(x){

  # Define functions for extracting Linf and K values from several objects (classes)
  getData <- list("tbl_df" = \(obj) obj,
                  "tbl" = \(obj) obj,
                  "data.frame" = \(obj) obj,
                  "grotagBoot" = \(obj) obj,
                  "lfqBoot" = \(obj) obj$bootRaw)

  # Searching the class of 'x' within the getData definitions
  index <- match(x = class(x), table = names(getData))[1]

  # If there is not a defined way (function) to extract Linf and K variables,
  # return an error msg
  if(is.na(index)){
    sprintf(fmt = "Internal funtion 'get_LinfK' do not know how to extract Linf or K from a '%s' object.",
            class(x)) |> stop()
  }

  # Extracting data an coerce to data.frame
  out <- getData[[index]](x) |> as.data.frame()

  # Set the column names as lowercase
  colnames(out) <- tolower(colnames(out))

  # Indexing Linf and K only and set standard column names
  out[,c("linf", "k")] |> setNames(c("Linf", "K"))
}
