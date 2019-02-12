# neutral-sims.r -- neutral validation 

NREPS <- 100
MCCORES <- 13
NUM_NEUTSITES <- 200

# parameter space
neut_res <- crossing(N=c(100, 500, 1000),
                   rho=c(0, 100, 2000),
                   rep=seq_len(NREPS)) %>% mutate(id=seq_len(nrow(.)))

neut_res$res <- mclapply(seq_len(nrow(neut_res)),
              function(i) {
                  t0 <- Sys.time()
	          p <- neut_res[i, ]
                  burnin <- prime_pop(1, N=p$N, r=p$rho/(4*p$N), 
                                      total_sites=NUM_NEUTSITES) 
                  burnin_pop <- burnin$samp
                  gams <- burnin_pop$gametes[[1]]
                  pos <- burnin_pop$positions[[1]]
                  out <- simpop(gams, pos, ngens=50, N=p$N, alpha=0,
                                sel_sites=0, r=p$rho/(4*p$N), include_rsq=TRUE, 
                                wfunc=w_constant)
                  elapsed <- as.numeric(Sys.time() - t0)
                  message(sprintf('completed replicate %d/%d (elapsed: %f)', 
                                  p$id, nrow(neut_res), elapsed))
                  out$elapsed <- elapsed
                  out$theta <- burnin$theta
                  return(out)
              }, mc.cores=MCCORES)

neut_res <- neut_res %>% mutate(covs=map(res, process_covs), 
                                stats=map(res, 'stats'))

save(neut_res, file='simdata/neutral-covs.Rdata')

