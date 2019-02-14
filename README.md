# The Linked Selection Signature of Rapid Adaptation in Temporal Genomic Data

![An image of decaying temporal autocovariance](https://i.imgur.com/tbEMzs3.png)

This file contains the code, manuscript, and some of the smaller intermediate
files necessary to reproduce the manuscript Buffalo and Coop (2019), *The
Linked Selection Signature of Rapid Adaptation in Temporal Genomic Data*.

## Required Software

### msprime

We use [msprime](https://github.com/tskit-dev/msprime) to generate gametes at
mutation-drift-recombination equilibrium. Note that `mspms` must be in the
`$PATH` of the server the simulations are running on. The following version was
used:

```
$ mspms --version
mspms 0.4.0
```

### Required R Packages

```
msr
tidyverse
gsl
parallel
rmarkdown
knitr
```

`msr` is my [MS data parser for R](https://github.com/vsbuffalo/msr).

## Simulation Data

All simulations are meant to be run in the main repository directory. They can
be executed with:

```
Rscript --vanilla scriptname.r
```

## Reproducing the Results

The following files were created by the simulation routines in `r_sims`,

```
expfit-fluct-covs.Rdata
expfit-varyl-varyn-covs.Rdata
ld-res.Rdata
neutral-covs.Rdata
sl-covs.Rdata
```

each by running `Rscript --vanilla scriptname.r` for the corresponding script.

We analyze each of these files in a separate Rmarkdown file, in `analyses/`. In
many case, these files create a processed data file stored in `data/`. The
`Makefile` renders the Rmarkdown files (and as as side effect, creates these
data files).

The processed data created by these Rmarkdown files is then loaded by
`analyses/figures/ms-plot.r` (and `analyses/figures/fluct-plot.r`) for
publication-quality figures. The final figures are then copied over to
`manuscript/images`. This is all done automatically through the `Makefile` in
`analyses/figures/`.

Overall:

1. Ensure you have downloaded the simulation results into `simdata/`. Make sure
   the necessary packages are installed.

2. Run `make` in `analyses/`.

3. Run `make in `analyses/figures/`.

### Session Info

Below is the session info for the server the simulations were run on:

```{r}
sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-redhat-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/R/lib/libRblas.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] TempSigLinkedSel_0.0.0.9000 msr_0.1.1
 [3] reshape2_1.4.3              forcats_0.3.0
 [5] stringr_1.3.1               dplyr_0.7.7
 [7] purrr_0.2.5                 readr_1.1.1
 [9] tidyr_0.8.1                 tibble_1.4.2
[11] ggplot2_3.0.0               tidyverse_1.2.1
[13] devtools_1.13.6

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.18     cellranger_1.1.0 plyr_1.8.4       pillar_1.3.0
 [5] compiler_3.5.1   bindr_0.1.1      tools_3.5.1      digest_0.6.17
 [9] lubridate_1.7.4  jsonlite_1.5     memoise_1.1.0    nlme_3.1-137
[13] gtable_0.2.0     lattice_0.20-35  pkgconfig_2.0.2  rlang_0.2.2
[17] cli_1.0.0        rstudioapi_0.7   commonmark_1.5   haven_1.1.2
[21] bindrcpp_0.2.2   httr_1.3.1       withr_2.1.2      roxygen2_6.1.0
[25] xml2_1.2.0       hms_0.4.2        grid_3.5.1       tidyselect_0.2.5
[29] glue_1.3.0       R6_2.2.2         readxl_1.1.0     modelr_0.1.2
[33] magrittr_1.5     backports_1.1.2  scales_1.0.0     rvest_0.3.2
[37] assertthat_0.2.0 colorspace_1.3-2 stringi_1.2.4    gsl_1.9-10.3
[41] lazyeval_0.2.1   munsell_0.5.0    broom_0.5.0      crayon_1.3.4
```

And the laptop analyses were conducted on:

```{r}
> sessionInfo()
R version 3.5.1 (2018-07-02)
Platform: x86_64-apple-darwin17.6.0 (64-bit)
Running under: macOS  10.14.3

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] parallel  stats     graphics  grDevices utils     datasets  methods
[8] base

other attached packages:
 [1] TempSigLinkedSel_0.0.0.9000 bindrcpp_0.2.2
 [3] msr_0.1.1                   reshape2_1.4.3
 [5] forcats_0.3.0               stringr_1.3.1
 [7] dplyr_0.7.6                 purrr_0.2.5
 [9] readr_1.1.1                 tidyr_0.8.1
[11] tibble_1.4.2                ggplot2_3.0.0
[13] tidyverse_1.2.1             usethis_1.4.0
[15] devtools_2.0.1

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.18      lubridate_1.7.4   lattice_0.20-35   prettyunits_1.0.2
 [5] ps_1.1.0          assertthat_0.2.0  rprojroot_1.3-2   digest_0.6.16
 [9] R6_2.2.2          wesanderson_0.3.6 cellranger_1.1.0  plyr_1.8.4
[13] backports_1.1.2   httr_1.3.1        pillar_1.3.0      rlang_0.2.2
[17] lazyeval_0.2.1    readxl_1.1.0      rstudioapi_0.7    callr_2.0.4
[21] desc_1.2.0        munsell_0.5.0     broom_0.5.0       compiler_3.5.1
[25] modelr_0.1.2      pkgconfig_2.0.2   gsl_1.9-10.3      pkgbuild_1.0.2
[29] tidyselect_0.2.4  crayon_1.3.4      withr_2.1.2       grid_3.5.1
[33] nlme_3.1-137      jsonlite_1.5      gtable_0.2.0      magrittr_1.5
[37] scales_1.0.0      cli_1.0.0         stringi_1.2.4     fs_1.2.6
[41] remotes_2.0.2     xml2_1.2.0        latex2exp_0.4.0   tools_3.5.1
[45] glue_1.3.0        hms_0.4.2         processx_3.2.0    pkgload_1.0.2
[49] colorspace_1.3-2  sessioninfo_1.1.1 rvest_0.3.2       memoise_1.1.0
[53] bindr_0.1.1       haven_1.1.2
```
