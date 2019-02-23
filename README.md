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

## Reproducing the Results

1. The following files were created by the simulation routines in `r_sims` on
   the server, 

   ```
   expfit-fluct-covs.Rdata
   expfit-varyl-varyn-covs.Rdata
   ld-res.Rdata
   neutral-covs.Rdata
   sl-covs.Rdata
   ```

each by running `Rscript --vanilla scriptname.r` for the corresponding script.

2. On the server,  run `make server` in `analyses/`, which generates some
   processed data from the main simulations, filtering out the N=1000 case, and
   then summarizing the complete (varying N) dataset on the server (too
   computationally intensive to run on my laptop). This generates:

    ```
    tempautocov/data/pd5f_varyn.rda
    tempautocov/data/predf_varyn.rda
    tempautocov/simdata/expfit-varyl-covs.Rdata
    ```

3. (If using separate machine for analysis): Copy the files above, and all
   files in `simdata/` over to analysis machine. `touch` each file to update
   timestamp so `make` does not try to regenerate them.

4. Run `make laptop` in `analyses/`; this renders all Rmarkdown files into
   HTML. It also has the side effect of creating all processed data file used
   for graphics stored in `data/`. 

5. Finally, in `analyses/figures`, run `make`. This runs
   `analyses/figures/ms-plot.r` (which sources `analyses/figures/fluct-plot.r`)
   for publication-quality figures. The final figures are then copied by the
   Makefile over to `manuscript/images`. 

### Session Info

Below is the session info for the server the simulations were run on:

```{r}
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
 [1] TempSigLinkedSel_0.0.0.9000 reshape2_1.4.3
 [3] forcats_0.3.0               stringr_1.3.1
 [5] dplyr_0.7.7                 purrr_0.2.5
 [7] readr_1.1.1                 tidyr_0.8.1
 [9] tibble_1.4.2                ggplot2_3.0.0
[11] tidyverse_1.2.1             msr_0.1.1
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
