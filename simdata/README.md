## Simulation Data

## Contents

Most simulations are run on a large machine. Here's the mapping of simulation
scripts from `r_sims` to data files. Note that `expon-varyl-varyn-covs.R`
filters the main complete simulation file to the N=1000 case, leading to
`expfit-varyl-covs.Rdata` which the main analysis of the paper is done on.

```
# main analysis
expfit-varyl-sims.r -> expfit-varyl-varyn-covs.Rdata
analyses/expon-varyl-varyn-covs.R + expon-varyl-varyn-covs.R -> expfit-varyl-covs.Rdata

# additional analyses
expfit-fluct-sims.r -> expfit-fluct-covs.Rdata
sl-sims.r -> sl-covs.Rdata

# validation
ld-test.r -> ld-res.Rdata
neutral-sims.r -> neutral-covs.Rdata
```

## Running Simulations

All simulations are meant to be run in the main repository directory. They an
be executed with:

```
Rscript --vanilla scriptname.r
```
