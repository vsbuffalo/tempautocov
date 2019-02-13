# single-locus-sims.r -- single locus sims
#
library(tidyverse)
load_all()

set.seed(42)

NREPS <- 500
MCCORES <- 15

sl_res <- crossing(N=1000,
                   alpha=c(0.05, 0.1, 0.2, 0.4),
                   rec_frac=c(0, 1e-4, 2e-04, 5e-04, 1e-03, 2e-02, 4.99999900e-01),
                   rep=seq_len(NREPS)) %>% mutate(id=seq_len(nrow(.)))  %>%
 		   # add seeds
                   mutate(seeds = map(id, ~ rand_seed()))
     


sl_res$res <- mclapply(seq_len(nrow(sl_res)),
                       function(i) {
                         t0 <- Sys.time()
                         p <- sl_res[i, ]
                         burnin <- sl_prime_pop(1, N=p$N, rec_frac=p$rec_frac, 
						seeds=unlist(p$seeds))
                         burnin_pop <- burnin$samp
                         gams <- burnin_pop$gametes[[1]][]
		 	 #  only use the two outer most mutations, which are guaranteed to 
			 # span the finite sites recomb point
			 gams <- gams[, c(1, ncol(gams))]
                         # we swap the gamete alleles, as we want the selected
                         # allele to be equal chance of ancestral/derived
                         gams <- swap_gamete_alleles(gams)
                         out <- sl_simpop(gams, ngens=100, N=p$N, 
					  alpha=p$alpha, rec_frac=p$rec_frac)
                         elapsed <- as.numeric(Sys.time() - t0)
                         message(sprintf('completed replicate %d/%d (elapsed: %f)', 
                                         p$id, nrow(sl_res), elapsed))
                         out$elapsed <- elapsed
                         out$theta <- burnin$theta
                         return(out)
                       }, mc.cores=MCCORES)

sl_res <- sl_res %>% mutate(covs=map(res, process_covs), 
                            stats=map(res, 'stats'))

save(sl_res, file='simdata/sl-covs.Rdata')

