## Multivariate Gaussian

#va_logistic <- function(va0, t, G, tinfl, b) va0 / (1 + exp(G * (t-tinfl))) + b
#cov_form_neut_LD <- function(ngens, r, hets, # known stuff
#                     N, va0, G, tinfl, b, # inferred stuff
#                     ld_int=NULL # optional
#                     ) {
# S <- matrix(NA, ncol=ngens, nrow=ngens) 
# # the minimum generation of each cell, for VA(0)
# T <- pmax(row(S), col(S))
# # get Va for each point
# VA <- va_logistic(va0, T, G, tinfl, b)  / 2
# # LD / recomb bit
# E <- abs(row(S) - col(S))  # elapsed generations between covars
# if (is.null(ld_int)) {
#   A <- assoc_int(E, r, N)
# } else {
#   A <- ld_int
# }
# dim(A) <- dim(S)
# out <- VA * A
# diag(out) <- diag(out) + 1/(2*N)
# out * hets[T]
#}


#freqs2deltas <- function(x, remove_fixations=FALSE) {
#  M <- fixations2NAs(x)
#  if (remove_fixations) {
#    M <- M[, colSums(is.na(M)) == 0]
#  }
#  D <- t(calc_deltas(t(M)))
#  H <- (M[-nrow(M), , drop=FALSE] * (1-M[-nrow(M), , drop=FALSE]))
#  return(list(deltas=D, hets=H))
#}



#tdat <- expfit_varyl_res %>% 
#  filter(Va == 0.1, genlen == 0.1, rep==1)
#tdat <- tdat$res[[1]]$neut_freqs
#F <- freqs2deltas(tdat)
#log_mvnorm <- function(x, sigma) {
#  dim(x) <- c(length(x), 1)
#  k <- nrow(x)
#  -0.5 * t(x) %*% solve(sigma) %*% x - log(sqrt((2*pi)^k * det(sigma)))
#}
#nll_va_logistic_neut_LD <- function(deltas, hets, r, init, ld_int=NULL) {
#  ngens <- nrow(deltas)
#  objfun <- function(par) { 
#       -sum(sapply(seq_len(ngens), function(i) {
#         delts <- deltas[,i, drop=FALSE]
#         hts <- hets[,i, drop=FALSE]
#         keep <- !is.na(delts) & !is.na(hts)
#         nkeep <- sum(keep)
#         if (!nkeep) return(0)
#         #sigma <- cov_form_neut_LD(nkeep, r, hts[keep], par['N'], 0, 0, 0, 0)
#         sigma <- cov_form_neut_LD(nkeep, r, hts[keep], par['N'], par['va0'], 
#                                   par['G'], par['tinfl'], par['b'], 
#                                   ld_int=ld_int)
#         log_mvnorm(delts[keep], sigma)
#    }))
#  }
#  #return(objfun)
#  optim(init, objfun)
#}
#make_objfun <- function(r, deltas, hets) {
#  ngens <- nrow(deltas)
#  objfun <- function(N, va0, G, tinfl, b) { 
#    S <- cov_form_neut_LD(ngens, r, N, va0, G, tinfl, b)
#    out <- apply(deltas/hets, 2, log_mvnorm, sigma=S)
#    sum(out[is.finite(out)])
#  }
#  objfun
#}
#fun <- make_objfun(0.1, F$deltas, F$hets)
#par_init <- list(N=1e3, va0=0.1, G=0.9, tinfl=10, b=0.1)
#ll_va_logistic_neut_LD(F, 0.1, par_init)
###
#neut_cov <- function(ngens, N) {
#  S <- matrix(0, ncol=ngens, nrow=ngens)
#  diag(S) <- 1/(2*N)
#  S
#}
#ndat <- na.exclude(data.frame(d=F$deltas[1,], h=F$hets[1,]))
#optimize(function(N) -sum(dnorm2(ndat$d, 0, sqrt(ndat$h/(2*N)), log=TRUE)), 
#         interval=c(0, 1e7))
