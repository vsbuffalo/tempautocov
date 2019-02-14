# The Linked Selection Signature of Rapid Adaptation in Temporal Genomic Data

## Required Packages

### MSPrime


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
