# expfit-varyl-sims.r -- exponential fitness with fixed Va, varying L
#
# These are to detect if it matters if we pass the genetic length as 
# r to 4Nr, or the recombination fraction. We run both.
set.seed(42)
library(tidyverse)
load_all()

NREPS <- 100
MCCORES <- 15
NUM_NEUTSITES <- 200


# # examples of the Va under stationary distribution
# theta <- 17
# L <- 100
# p <- rbeta(L, theta/L, theta/L)
# Va <- 0.001
# 2*p*(1-p)*alpha(Va, theta, L)

# parameter space
expfit_varyl_res <- crossing(L=c(10, 50, 100, 500),
                             N=c(100, 500, 1000),
                             # last two are the other end of a chromosome ~1.5M long
                             # and two (nearly) unlinked loci
                             #
                             # see fitness notes
                             Va=c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1), 
                             genlen=c(0, 0.005, 0.01, 0.05, 0.1, 0.5, 1.5, 4.5),
                             rep=seq_len(NREPS)) %>% 
                             mutate(r=haldane(genlen),
                                    rho=4*N*genlen) %>%
                             # calculate the target alpha
                             mutate(theta_L=target_theta(L, 2*N), 
                                    alpha=alpha(Va, theta_L)) %>%
			     # add seeds
			     mutate(seeds = map(L, ~ rand_seed()))
      

# neutral parameter space
expfit_varyl_neut_res <- crossing(L=0,
                             N=1000,
                             # last two are the other end of a chromosome ~1.5M long
                             # and two (nearly) unlinked loci
                             #
                             # see fitness notes
                             Va=c(0), 
                             genlen=c(0, 0.005, 0.01, 0.05, 0.1, 0.5, 1.5, 4.5),
                             rep=seq_len(NREPS)) %>% 
                             mutate(r=haldane(genlen),
                                    rho=4*N*genlen) %>%
                             # calculate the target alpha
                             mutate(theta_L=target_theta(L, 2*N), 
                                    alpha=0) %>%
			     # add seeds
			     mutate(seeds = map(L, ~ rand_seed()))

expfit_varyl_res <- rbind(expfit_varyl_res, expfit_varyl_neut_res) %>%  
	             mutate(id=seq_len(nrow(.)))

expfit_varyl_res$res <- mclapply(seq_len(nrow(expfit_varyl_res)),
              function(i) {
                  t0 <- Sys.time()
	          p <- expfit_varyl_res[i, ]
                  burnin <- prime_pop(1, N=p$N, total_sites=p$L+NUM_NEUTSITES,
				      #seeds=unlist(p$seeds),
				      r=p$genlen)
                  burnin_pop <- burnin$samp
                  gams <- burnin_pop$gametes[[1]]
                  pos <- burnin_pop$positions[[1]]
                  out <- simpop(gams, pos, ngens=50, N=p$N, alpha=p$alpha,
                                sel_sites=p$L, r=p$genlen,
                                wfunc=w_exp_factory(5), include_allelic_cov=TRUE)
                  elapsed <- as.numeric(Sys.time() - t0)
                  message(sprintf('completed replicate %d/%d (elapsed: %f)', 
                                  p$id, nrow(expfit_varyl_res), elapsed))
                  out$elapsed <- elapsed
                  out$theta <- burnin$theta
                  return(out)
              }, mc.cores=MCCORES)

# some large L simulations fail, as the some of the stochastic selection of 
# selected loci doesn't produce enough (TODO). For now this is a small concern; 
# we filter these out.
expfit_varyl_res <- expfit_varyl_res %>% 
	             mutate(success=map_lgl(res, ~ class(.) != 'try-error')) %>%
		     filter(success)

# save, since this is a lot of simulations
#save(expfit_varyl_res, file='simdata/expfit-largel-sim.Rdata')

message('sims complete, processing...')
# calculate all the covariances
expfit_varyl_res <- expfit_varyl_res %>% mutate(covs=map(res, process_covs), 
                                                stats=map(res, 'stats'),
						genic_va=map(res, process_genic_va, as_df=TRUE),
						hets=map(res, process_het),
						# extract out allelic cov, since 
						# we need this later but don't want to
						# store res
						allelic_cov=map(res, 'allelic_cov'))

message('sims complete, processing cov mats...')
# covariances in matrix form
expfit_varyl_res <- expfit_varyl_res %>% 
  mutate(covs_mat=map(res, process_covs, as_df=FALSE))

# save the covariances
message('writing sim results...')
save(expfit_varyl_res, file='simdata/expfit-varyl-varyn-covs.Rdata')
message('writing sim complete.')

expfit_varyl_res_orig <- expfit_varyl_res
expfit_varyl_res <- expfit_varyl_res %>% select(-res)

save(expfit_varyl_res, file='simdata/expfit-varyl-varyn-covs-only.Rdata')

