% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ELEFAN_GA_boot.R
\name{ELEFAN_GA_boot}
\alias{ELEFAN_GA_boot}
\title{Bootstrapped ELEFAN_GA}
\usage{
ELEFAN_GA_boot(lfq, seasonalised = FALSE, low_par = NULL,
  up_par = NULL, parallel = TRUE, nresamp = 200,
  no_cores = detectCores() - 1, clusterType = "PSOCK",
  outfile = "output.txt", popSize = 60, maxiter = 50, run = 200,
  pmutation = 0.2, pcrossover = 0.8, elitism = base::max(1,
  round(popSize * 0.05)), MA = 5, addl.sqrt = FALSE, agemax = NULL,
  flagging.out = TRUE, resample = TRUE, seed = NULL)
}
\arguments{
\item{lfq}{a length frequency object of the class `lfq`}

\item{seasonalised}{logical; indicating if the seasonalised von Bertalanffy
growth function should be applied (default: FALSE).}

\item{low_par}{a list providing the minimum of the search space in case
of real-valued or permutation encoded optimizations. When set to NULL the
following default values are used:
 \itemize{
  \item \strong{Linf} length infinity in cm (default is calculated from maximum
  length class in the data),
  \item \strong{K} curving coefficient (default: 0.01),
  \item \strong{t_anchor} time point anchoring growth curves in year-length
  coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
  \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
  \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
}}

\item{up_par}{a list providing the maximum of the search space in case of
real-valued or permutation encoded optimizations. When set to NULL the
following default values are used:
 \itemize{
  \item \strong{Linf} length infinity in cm (default is calculated from maximum
  length class in the data),
  \item \strong{K} curving coefficient (default: 1),
  \item \strong{t_anchor} time point anchoring growth curves in year-length
  coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 1),
  \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 1),
  \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 1);
}}

\item{parallel}{logical; should parallelized computing be used. This differs from the
`parallel` argument in \code{\link[TropFishR]{ELEFAN_GA}} in that it is not used
within the \code{\link[GA]{ga}} function for calculation at the population level,
but rather for each resampling. Depending on platform operating system,
the argument `clusterType` can be adjusted (see argument description for
details). (Default: `parallel = TRUE`)}

\item{nresamp}{numeric; the number of permutations to run (Default: `nresamp = 200`)}

\item{no_cores}{numeric (Default: `no_cores = detectCores() - 1`)}

\item{clusterType}{(Default: `clusterType = "PSOCK"`)}

\item{outfile}{character; text file name (Default: `outfile = "output.txt"`) which will
records the progress of the permutation completions.}

\item{popSize}{the population size. (Default: `popSize = 60`)}

\item{maxiter}{the maximum number of iterations to run before the GA search is halted.
(Default: `maxiter = 50`)}

\item{run}{the number of consecutive generations without any improvement
in the best fitness value before the GA is stopped. (Default: `run = 10`)}

\item{pmutation}{numeric. A small fraction of 1.0. The probability of mutation in a
parent chromosome. Usually mutation occurs with a small probability.
(Default: `pmutation = 0.2`)}

\item{pcrossover}{the probability of crossover between pairs of chromosomes.
Typically this is a large value.  (Default: `pcrossover = 0.8`)}

\item{elitism}{the number of best fitness individuals to survive at each generation.
By default the top 5\% individuals will survive at each iteration.}

\item{MA}{number indicating over how many length classes the moving average
should be performed (default: 5, for more information see \code{\link[TropFishR]{lfqRestructure}}}

\item{addl.sqrt}{logical. Should counts be square root transformed prior to restructuring.}

\item{agemax}{numeric. maximum age}

\item{flagging.out}{logical Should flagging out be done}

\item{resample}{logical. Should `lfq` object be resampled (Default: `resample = TRUE`).
When `resample = FALSE`, a `partial bootstrap` is performed, reflecting solution
variation due only to the search algorithm.}

\item{seed}{seed value for random number reproducibility (Default: NULL)}
}
\value{
a list of class `lfqBoot` containing 2 levels: `$bootRaw` - a data.frame of fitted VBGF parameters
(columns) by resampling (rows), `$seed` - a vector of seed values set prior to each resampling
call to `lfqResample`.
}
\description{
`ELEFAN_GA_boot` performs a bootstrapped fitting of
  von Bertalanffy growth function (VBGF) via the \code{\link[TropFishR]{ELEFAN_GA}} function.
  Most of the arguments are simply passed to the function within many permutations
  (resampling) of the original lfq data.
}
\examples{
\donttest{
# load data
library(TropFishR)
data(alba)

# settings (these settings may not be optimal - for demo only)
MA <- 7
low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
popSize <- 60
maxiter <- 50
run <- 10
pmutation <- 0.2
nresamp <- 12


# parallel version
library(parallel)
t1 <- Sys.time()
res <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
  up_par = up_par, low_par = low_par,
  popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
  nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
  seed = 1
)
t2 <- Sys.time()
t2 - t1
res


# non-parallel version
t1 <- Sys.time()
res <- ELEFAN_GA_boot(lfq=alba, seasonalised = FALSE,
  up_par = up_par, low_par = low_par,
  popSize = popSize, maxiter = maxiter, run = run,
  pmutation = pmutation, nresamp = nresamp, MA = MA, parallel = FALSE,
  seed = 1
)
t2 <- Sys.time()
t2 - t1
res

# plot resulting distributions
univariate_density(res, use_hist = TRUE)

# plot scatterhist of Linf and K
LinfK_scatterhist(res)
}

}
