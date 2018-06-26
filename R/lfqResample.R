#' Resampling of length-frequency data
#'
#' @param lfq a length frequency object of the class `lfq`.
#'
#' @description The function resamples the `lfq` data by sample date.
#'   Sampling is done in a non-parametric manner that follows the relative
#'   frequencies of the original data, allowing for individual counts to be
#'   selected more than once (i.e. `replace = TRUE` in \link[base]{sample}),
#'   and resulting in total counts (by sample) equal to the original data.
#'
#' @return a resampled version of the `lfq` class dataset.
#'
#' @examples
#' # load data
#' library(TropFishR)
#' data(alba)
#'
#' # resample lfq data
#' alba_p <- lfqResample(alba)
#'
#' # side-by-side plot
#' op <- par( mfcol = c(2,1), mar=c(4,4,2,1) )
#' plot(lfqRestructure(alba), Fname="rcounts")
#' mtext("original", side=3, line=0.25)
#' plot(lfqRestructure(alba_p), Fname="rcounts")
#' mtext("resampled", side=3, line=0.25)
#' par(op)
#'
#' # relative difference
#' alba_diff <- alba
#' alba_diff$delta <- (alba$catch - alba_p$catch) / alba$catch
#' alba_diff$delta[is.na(alba_diff$delta)] <- 0
#' plot(alba, Fname = "catch", image.col = NA, hist.col=NA, draw = FALSE)
#' with(alba_diff, image(
#'   x=dates, y=midLengths, z=t(delta),
#'   zlim = max(abs(delta))*c(-1,1), col=rev(cm.colors(21)),
#'   add=TRUE
#' ))
#' with(alba_diff, contour(x=dates, y=midLengths, z=t(delta), add=TRUE))
#' box()
#'
#' @importFrom grDevices adjustcolor blues9 colorRampPalette rgb
#'
#' @export
#'
lfqResample <- function(lfq){
  # define bin breaks
  bin.width <- diff(lfq$midLengths) # bin width (should allow for uneven bin sizes)
  bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2) # upper bin limit
  bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2) # lower bin limit
  breaks <- unique(c(bin.lower, bin.upper))
  # copy lfq
  lfqb <- lfq
  # resample with replacement (n = sum(lfq$catch[,i]))
  for(i in seq(length(lfq$dates))){
    # resample with replacement using bin frequencies to inform probability weights
    inds <- sample(x = lfq$midLengths, size = sum(lfq$catch[,i]), prob = lfq$catch[,i], replace = TRUE)
    h <- hist(inds, breaks = breaks, plot = FALSE)
    lfqb$catch[,i] <- h$counts
  }
  return(lfqb)
}

