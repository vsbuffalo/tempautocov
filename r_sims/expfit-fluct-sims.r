# expfit-flip-sims.r -- exponential fitness function 
# models with flipped fitness function
library(devtools)
load_all()

NREPS <- 100
MCCORES <- 10
NUM_NEUTSITES <- 200

# parameter space
expfit_fluct_res <- crossing(L=c(500),
                             N=c(1000),
                             # last two are the other end of a chromosome ~1.5M long
                             # and two (nearly) unlinked loci
                             #
                             # see fitness notes
                             Va=c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1), 
                             genlen=c(0, 0.005, 0.01, 0.05, 0.1, 0.5, 1.5, 4.5),
                             rep=seq_len(NREPS)) %>% 
                             mutate(r=haldane(genlen),
                                    rho=4*N*genlen) %>%
                             # calculate the target alpha
                             mutate(theta_L=target_theta(L, 2*N), 
                                    alpha=alpha(Va, theta_L)) %>%
			     # add seeds
			     mutate(seeds = map(L, ~ rand_seed()))

expfit_fluct_res$res <- mclapply(seq_len(nrow(expfit_fluct_res)),
              function(i) {
                  t0 <- Sys.time()
	          p <- expfit_fluct_res[i, ]
                  burnin <- prime_pop(1, N=p$N, total_sites=p$L+NUM_NEUTSITES,
				      #seeds=unlist(p$seeds),
				      r=p$genlen)
                  burnin_pop <- burnin$samp
                  gams <- burnin_pop$gametes[[1]]
                  pos <- burnin_pop$positions[[1]]
                  out <- simpop(gams, pos, ngens=50, N=p$N, alpha=p$alpha,
                                sel_sites=p$L, r=p$genlen,
                                wfunc=w_exp_fluct_factory(5, 15))
                  elapsed <- as.numeric(Sys.time() - t0)
                  message(sprintf('completed replicate %d/%d (elapsed: %f)', 
                                  p$id, nrow(expfit_fluct_res), elapsed))
                  out$elapsed <- elapsed
                  out$theta <- burnin$theta
                  return(out)
              }, mc.cores=MCCORES)
	

# save, since this is a lot of simulations
#save(expfit_fluct_res, file='simdata/expfit-flip-sim.Rdata')
expfit_fluct_res <- expfit_fluct_res %>% 
	             mutate(success=map_lgl(res, ~ class(.) != 'try-error')) %>%
		     filter(success)

save(expfit_fluct_res, file='simdata/expfit-fluct-sims')

# calculate all the covariances
expfit_fluct_res <- expfit_fluct_res %>% mutate(covs=map(res, process_covs), 
                                         stats=map(res, 'stats'))

#
expfit_fluct_res <- expfit_fluct_res %>%
    mutate(covs_mat=map(res, process_covs, as_df=FALSE))
   
# save the covariances
expfit_fluct_res <- expfit_fluct_res %>% select(-res)
save(expfit_fluct_res, file='simdata/expfit-fluct-covs.Rdata')

