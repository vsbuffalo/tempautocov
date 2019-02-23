all: expon-varyl-covs.html neutral-validation.html sl-cov.html expon-fluct-covs.html

# The data set with all Ns of simulation runs is quite large, and our focal
# analysis focuses on N=1000. 
#
# The script expon-varyl-varyn-covs.R filters out the N=1000 case, writing it
# to '../simdata/expfit-varyl-covs.Rdata', as well as the processed
# ../data/pd5f_varyn.rda ../data/predf_varyn.rda files (and their TSV
# versions). This script requires more computational power than a laptop.

# the processed data
varyn_data=../data/pd5f_varyn.rda ../data/predf_varyn.rda 

# RUN ON SERVER:
#$(varyn_data) ../simdata/expfit-varyl-covs.Rdata: ../simdata/expfit-varyl-varyn-covs.Rdata
#	Rscript --vanilla expon-varyl-varyn-covs.R


expon-varyl-covs.html: expon-varyl-covs.Rmd ../simdata/expfit-varyl-covs.Rdata $(varyn_data)
neutral-validation.html: neutral-validation.Rmd ../simdata/neutral-covs.Rdata
sl-cov.html: sl-cov.Rmd ../simdata/sl-covs.Rdata
expon-fluct-covs.html: expon-fluct-covs.Rmd ../simdata/expfit-fluct-covs.Rdata

%.html: %.Rmd
	Rscript --vanilla -e 'library(rmarkdown); render("$<")'

