[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/fishboot)](https://cran.r-project.org/package=fishboot) [![CRAN_time_from_release](https://www.r-pkg.org/badges/ago/fishboot)](https://cran.r-project.org/package=fishboot) [![metacran downloads](https://cranlogs.r-pkg.org/badges/fishboot)](https://cran.r-project.org/package=fishboot)
# fishboot :fish: :boot: <a><img src="man/figures/logo.png" align="right" height="150" /></a>

A suite of bootstrap-based models and tools for analyzing fish stocks and aquatic populations. Designed for ecologists and fisheries scientists, it supports data from length-frequency distributions, tag-and-recapture studies, and hard structure readings (e.g., otoliths). See Schwamborn et al., 2019 <doi:10.1016/j.ecolmodel.2018.12.001> for background. The package includes functions for bootstrapped fitting of growth curves and plotting.

**To install** (using `remotes`):
```
# master branch
remotes::install_github(repo = "rschwamborn/fishboot")

# dev branch (unstable)
remotes::install_github(repo = "rschwamborn/fishboot", ref = "dev")
```

**to cite**:
```
citation("fishboot")
```

> Ralf Schwamborn, Tobias K. Mildenberger, Marc H. Taylor, Margit Wilhelm, Wencheng Lau-Medrano (2025).
  fishboot: Bootstrap-Based Methods for the Study of Fish Stocks and
  Aquatic Populations. R package version 1.0.0.
  https://github.com/rschwamborn/fishboot
