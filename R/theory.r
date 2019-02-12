## theory.r -- theoretic covariances 

ok71 <- function(rho) (10 + rho)/(22 + 13*rho + rho^2)

sved71 <- function(rho) 1/(1 + rho)

triangle_dist_pdf <- function(rd, R) (2/R^2) * (R - rd)


assoc_int <- Vectorize(function(t, R, N, D2fun=c('ok71', 'sved71')) {
  if (R == 0) return(NA)
  D2fun <- match.arg(D2fun)
  funs <- list(ok71=function(rd) {
                triangle_dist_pdf(rd, R)*ok71(haldane(rd)*4*N)*(1-haldane(rd))^t
               }, 
               sved71=function(rd) {
                triangle_dist_pdf(rd, R)*sved71(haldane(rd)*4*N)*(1-haldane(rd))^t
               })
  fun <- funs[[D2fun]]
  out <- integrate(fun, lower=0, upper=R, stop.on.error=FALSE)
  # the integral is divergent for large genlen and large time
  # but it likely settles down to zero. TODO: maybe work out 
  # more explicitly, but not priority.
  if (out$message != 'OK' && t > 10 && R > 1) return(0)
  # other cases:
  if (out$message != 'OK') return(NA)
  out$value
})


analytic_cov <- function(r, gen, va, N, s=1, D2fun='ok71') {
  if (r==0) return(va)
  s^2 * va/2 * assoc_int(gen, r, N, D2fun=D2fun)
}

sl_analytic_cov <- function(r, gen, va, N) {
  if (r==0) return(va)
  va/2 * ok71(4*r*N) * (1-r)^gen
}


#' Calculate the predicted cummulative variance and covariance 
pred_covvar <- function(vas, ntime, R, N) {
  x <- matrix(NA, ncol=ntime, nrow=ntime)
  stopifnot(length(vas) == nrow(x))
  vas_m <- vas[pmin(row(x), col(x))]
  gen <- pmax(row(x), col(x)) - pmin(row(x), col(x))
  covs <- sapply(seq_along(gen), function(i) analytic_cov(R, gen[i], vas_m[i], N))
  dim(covs) <- dim(x)
  # add on the WF variation
  diag(covs) <- diag(covs) + 1/2e3
  labs <- matrix('cov', nrow=nrow(x), ncol=ncol(x))
  diag(labs) <- 'var'
  gen_t <- row(x) 
  gen_s <- col(x) 
  d <- tibble(type=as.vector(labs), vals=as.vector(covs)) 
  d <- d %>% group_by(type) %>% summarize(vals=sum(vals, na.rm=TRUE))
  d
}
