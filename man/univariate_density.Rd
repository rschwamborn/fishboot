% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univariate_density.R
\name{univariate_density}
\alias{univariate_density}
\title{Univariate kernel density estimate plot of VBGF parameter from
bootstrapping results}
\usage{
univariate_density(
  res,
  CI = 95,
  use_hist = FALSE,
  nbreaks = 10,
  mar = c(1.5, 2, 2, 0),
  oma = c(1.5, 0, 0, 0.5),
  mgp = c(2, 0.5, 0),
  tcl = -0.25,
  cex = 1,
  ...
)
}
\arguments{
\item{res}{Object of class \code{lfqBoot}.}

\item{CI}{\code{numeric}. Confidence interval in \% (default: 95).}

\item{use_hist}{\code{logical} Plot histogram in addition to smoothed kernel
density.}

\item{nbreaks}{\code{numeric} vector specifying the number of breaks in the
histogram.}

\item{mar, oma, mgp, tcl, cex, ...}{Additional arguments passed to
\link[graphics]{par}.}
}
\value{
This function returns just the described plot.
}
\description{
This function plots a set of vertical plots with kernel density
distributions for univariate posterior distributions of the VBGF growth
parameters \eqn{L_{inf}}, \eqn{K}, and \eqn{Phi’}. The 95% confidence
interval and the most likely optimum fit estimate are shown for each parameter.
}
\details{
This function used the function \link[ks]{kde} to obtain kernel density
estimates for the VBGF growth parameters \eqn{L_{inf}}, \eqn{K}, and
\eqn{Phi’}. The 95 % confidence interval (i.e., the 95% quantiles for each
posterior distribution) and the most likely optimum fit estimate (i.e., the
mode of each posterior distribution) are then plotted inside each vertical
plot. The input used for plotting is usually the result of a bootstrapped
growth analysis (i.e. a \code{lfqBoot} object generated by \strong{fishboot}
functions such as \link{ELEFAN_SA_boot}, \link{ELEFAN_GA_boot},
\link{grotag_boot}, or \link{grolenage_boot}).
}
\examples{
data(alba_boot)
univariate_density(alba_boot)
}
