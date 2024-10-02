#' @title Resampling of length-frequency data
#'
#' @description This function resamples the \code{lfq} data by sampling dates.
#' Sampling is done in a non-parametric way following the relative frequencies
#' of the original data, allowing for individual counts to be selected more than
#' once (i.e. \code{replace = TRUE} in \link[base]{sample}), and resulting in
#' total counts (by sample) equal to the original data.
#'
#' @param lfq A length frequency object of the class \code{lfq}.
#'
#' @return A resampled version of the \code{lfq} class dataset.
#'
#' @export
#'
#' @examples
#' # Load data
#' data("alba", package = "TropFishR")
#'
#' # Resample lfq data
#' alba_p <- lfqResample(lfq = alba)
#'
#' # Side-by-side plot
#' op <- par(mfcol = c(2,1), mar = c(4, 4, 2, 1))
#'
#' # Original
#' plot(x = TropFishR::lfqRestructure(alba), Fname = "rcounts")
#' mtext("original", side=3, line=0.25)
#'
#' # Resampled
#' plot(TropFishR::lfqRestructure(alba_p), Fname = "rcounts")
#' mtext("resampled", side=3, line=0.25)
#'
#' par(op)
lfqResample <- function(lfq){
  # Bin width (should allow for uneven bin sizes)
  bin.width <- diff(x = lfq$midLengths)

  # Upper bin limit
  bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2)

  # Lower bin limit
  bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2)

  breaks <- unique(c(bin.lower, bin.upper))

  # copy lfq
  lfqb <- lfq

  # Resample with replacement (n = sum(lfq$catch[,i]))
  for(i in seq(length(lfq$dates))){
    # Resample with replacement using bin frequencies to inform probability weights
    inds <- sample(x = lfq$midLengths,
                   size = sum(lfq$catch[,i]),
                   prob = lfq$catch[,i],
                   replace = TRUE)

    lfqb$catch[,i] <- hist(x = inds, breaks = breaks, plot = FALSE)$counts
  }

  return(lfqb)
}
