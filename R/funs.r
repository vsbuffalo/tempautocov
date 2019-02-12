# funs.r -- misc functions

γ = 0.5772156649

calc_trait_mu <- function(nloci, N) nloci/(4*N * (γ + log(N)))

segsites <- function(theta, n) theta * (γ + log(n))

# these should be inverse functions of each other
stopifnot(segsites(4e3*calc_trait_mu(100, 1e3), 1e3) == 100)

center_variant <- function(effects, neut_freqs=NULL) {
  # get the middle variant (if neut_freqs is NULL).
  # if neut_freqs is a list of trajectories, then get the middle variant 
  # that's polymorphic the whole time
  neut_sites <- seq_along(effects)[is.na(effects)]
  if (is.null(neut_freqs)) {
    # get the index of the neutral variant (NA) closest to the middle
    # assumes variants are in rank order
    i <- order(abs(neut_sites - length(neut_sites)/2))
    return(i[1])
  } else {
    is_poly <- !apply(neut_freqs == 0 | neut_freqs == 1, 2, any)
    is_poly_neut <- neut_sites[is_poly]
    i <- order(abs(is_poly_neut - length(is_poly_neut)/2))
    return(i[1])
  }
}

center_freqs <- function(effects, neut_freqs) {
  neut_freqs[, center_variant(effects, neut_freqs)]
}

offspring_var <- function(N, mu, sigma2, Ep=1/N) {
  wbar <- exp(sigma2/2)
  vwz <- exp(2*(mu + sigma2)) - 2*exp(mu+sigma2/2)*wbar + wbar^2
  Vp <- vwz/(N*wbar)^2
  N*Ep*(1-Ep) + N*(N-1)*Vp
}

#' Generarlized lapply, for one core uses lapply(), for many cores uses mclapply()
#' 
#' @param X list or vector to apply function over.
#' @param FUN function to apply function to each element.
#' @param ... additional arguments to pass to FUN.
#' @export
glapply <- function(X, FUN, ..., mc.cores=1) {
  if (mc.cores == 1) {
    lapply(X, FUN, ...)
  } else {
    mclapply(X, FUN, ..., mc.cores=mc.cores)
  }
}

#' @export
calc_freqs <- function(gams) {
  colSums(gams)/(nrow(gams))
}

#' @export
#calc_rsq <- function(gams, i=1, j=2, min_freq=0.1) {
calc_rsq <- function(gams, i=1, j=2) {
  p <- colMeans(gams[, c(i, j)])
  # if (p[1] < min_freq || p[2] < min_freq || 
      # p[1] > 1-min_freq || p[2] > 1-min_freq)
    # return(NA)
  rsq <- (mean(gams[, i]*gams[, j] - prod(p))^2 / (p[1] * (1-p[1]) * p[2] * (1-p[2])))
  tibble(rsq=rsq, pi=p[1], pj=p[2])
}

#' Return a logical vector of which frequencies pass thresholds
freq_filter <- function(x, a=0.1) x > a & x < 1-a

#' Return the R^2 for the entire region (first polymorphic site to last)
#'
#' @gams the gamete matrix, nsims x nloci
#'
#' @export
calc_rsq_region <- function(gams, frac=1) {
  calc_rsq(gams, i=1, j=floor(frac*ncol(gams)))
}



#' @export
calc_D <- function(gams, i=1, j=2) {
  p <- colMeans(gams[, c(i, j)])
  tibble(D=mean(gams[, i]*gams[, j] - prod(p)), pi=p[1], pj=p[2])
}


#' Inverse of Haldane's mapping function
#'
#' @param r recombination fraction
#'
#' Takes a recombination fraction, returns the corresponding map distance based
#' on a Poisson model of recombination. The map distance will be larger, due to 
#' double crossovers.
#' @export
inv_haldane <- function(r) -0.5*log(1-2*r)


#' Haldane's mapping function @param d genetic distance in Morgans
#' 
#' Takes a genetic distance and calculates the probability of an odd number of
#' recombination events, e.g. the recombination fraction.
#'
#' @export 
haldane <- function(d) 0.5*(1-exp(-2*d))


#' Find fixations (freq = 1) in the matrix of allele frequencies over time
#'
#' @param mat matrix of allele frequencies, ntime x nloci.
#' @export
fixations <- function(mat) {
  apply(mat == 1, 2, any)
}


#' Find when fixations (freq = 1) occur in the matrix of allele frequencies
#' over time
#'
#' @param mat matrix of allele frequencies, ntime x nloci.
#' @export
when_fix <- function(mat) {
  times <- apply(mat, 2, function(x) min(which(x == 1)))
  times[is.finite(times)]
}

#' Swap the alleles of a gamete randomly
#' ignore is a vector of indices to ignore
swap_gamete_alleles <- function(gams, effects=NULL) {
  if (is.null(effects)) {
    swap <- rbinom(ncol(gams), 1, 0.5)
    return(abs(sweep(gams, 2, swap)))
  }
  
  # if we have effect sizes, we swap neutral alleles, and 
  # polarize by effect size
  polarized <- ifelse(effects < 0, 0, 1)
  swap <- rbinom(sum(is.na(polarized)), 1, 0.5)
  polarized[is.na(polarized)] <- swap
  abs(sweep(gams, 2, polarized))
}


