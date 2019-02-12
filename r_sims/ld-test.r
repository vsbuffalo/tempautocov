## ld-test.r --

# we create a population with a high-level of LD: N/2 gametes of 11...1 and N/2 gametes of 00....0
# and then measure the decay.

L <- 100
NREPS <- 50
MCCORES <- 13

# parameter space
ld_res <- crossing(N=c(100, 500, 1000),
                   rho=c(0, 50, 100, 2000),
                   rep=seq_len(NREPS)) %>% mutate(id=seq_len(nrow(.)))

ld_res$res <- glapply(seq_len(nrow(ld_res)),
              function(i) {
                  t0 <- Sys.time()
                  p <- ld_res[i, ]
                  gams <- rbind(matrix(1, nrow=p$N, ncol=L), matrix(0, nrow=p$N, ncol=L))
                  pos <- sort(runif(L))
                  out <- simpop(gams, pos, ngens=50, N=p$N, alpha=0,
                                sel_sites=0, r=p$rho/(4*p$N),
                                include_rsq=TRUE, include_D=TRUE,
                                wfunc=w_constant)
                  elapsed <- as.numeric(Sys.time() - t0)
                  message(sprintf('completed replicate %d/%d (elapsed: %f)',
                                  p$id, nrow(ld_res), elapsed))
                  out$elapsed <- elapsed
                  return(out)
              }, mc.cores=MCCORES)

ld_res_inv_haldane <- ld_res
ld_res_inv_haldane$res <- glapply(seq_len(nrow(ld_res)),
              function(i) {
                  t0 <- Sys.time()
                  p <- ld_res[i, ]
                  gams <- rbind(matrix(1, nrow=p$N, ncol=L), matrix(0, nrow=p$N, ncol=L))
                  pos <- sort(runif(L))
                  out <- simpop(gams, pos, ngens=50, N=p$N, alpha=0,
                                sel_sites=0, r=inv_haldane(p$rho/(4*p$N)),
                                include_rsq=TRUE, include_D=TRUE,
                                wfunc=w_constant)
                  elapsed <- as.numeric(Sys.time() - t0)
                  message(sprintf('completed replicate %d/%d (elapsed: %f)',
                                  p$id, nrow(ld_res), elapsed))
                  out$elapsed <- elapsed
                  return(out)
              }, mc.cores=MCCORES)

save(ld_res, ld_res_inv_haldane, file='simdata/ld-res.Rdata')
