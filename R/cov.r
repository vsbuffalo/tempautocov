# cov.r -- temporal covariance code

#' Randomly swap the alleles of the supplied frequency matrix
#'
#' @param freqs a matrix of allele frequencies, \code{nloci} by 
#'    \code{ntimepoints}.
#' @export
swap_alleles <- function(freqs) {
  # this returns the transpose of the matrix
  flip <- matrix(0, ncol=ncol(freqs), nrow=nrow(freqs))
  flip[sample(c(TRUE, FALSE), nrow(freqs), replace=TRUE)] <- 1
  abs(flip-freqs)
}

#' Convert Fixations/Losses to NA
#'
#' @params freqs a matrix of allele frequencies
fixations2NAs <- function(freqs) {
  freqs[freqs == 0 | freqs == 1] <- NA
  freqs 
}


#' Calculate the allele frequency change between adjacent generations
#'
#' @params freqs a matrix of allele frequencies, 
#'  \code{nloci} x \code{ntimepoints}.
#' @export
calc_deltas <- function(freqs) {
  freqs[, -1] - freqs[, -ncol(freqs), drop=FALSE]
}


#' Calculate the allele frequency change between arbitrary generations
#'
#' @params freqs a matrix of allele frequencies, 
#'  \code{nloci} x \code{ntimepoints}.
#' @params t0 the first timepoint.
#' @params t0 the second timepoint.
#' @export
calc_deltas2 <- function(freqs, t0=1, t1=ncol(freqs)) {
  freqs[, t1, drop=FALSE] - freqs[, t0, drop=FALSE]
}


#' Calculate the temporal covariance for a matrix of frequencies across time
#'
#' A note, the differences are calculated using the matrix
#' subtraction below as (by row): Δp_{t-1} = p_t - p_{t-1}. Do note that our 
#' notation in the paper is always Δp_t = p_t - p_{t-1}.
#'
#' Covariances between adjacent timepoints (cov(Δp_t, Δp_{t+1})) also 
#' share sampling noise; this leads to negative covariance. This is corrected
#' with the vector N of number of individuals sampled, using a binomial error 
#' model.
#'
#' @param freqs a matrix of allele frequencies, \code{nloci} by 
#'  \code{ntimepoints}
#' @param as_df return the results as a data.frame/tibble
#' @param swap randomly swap the alleles of the frequency matrix
#' @param N a vector of the number of individuals sampled, to the adjacent time point correction.
#' @export
temp_cov <- function(freqs, as_df=FALSE, swap=TRUE, N=NULL, upper_tri=TRUE, 
                     remove_fixed=TRUE, standardize=TRUE) {
  # NOTE: all downstream code work with transpose of the matrix
  # provided by R simulation routines. swap_alleles() returns transpose
  # otherwise, we do it here.
  if (swap) {
    mat <- swap_alleles(freqs)
  } else {
    mat <- freqs
  }

  # remove fixations, converting to NA. This was source of very subtle bug 
  # where fixations were removed in calc_deltas but not hetmat
  if (remove_fixed) {
    mat <- fixations2NAs(mat)
  }

  dmat <- calc_deltas(mat)

  if (FALSE) {  # experimental! Idea is to permute loci.
    N <- NULL  # no correction with permutation
    dmat <- apply(dmat, 2, function(x) sample(x))
  }

  # remove loci with no allele frequency changes
  if (nrow(dmat) > 1) {
    covmat <- cov(dmat, use='pairwise.complete.obs')
  } else {
    covmat <- t(dmat) %*% dmat
  }

  if (!is.null(N)) {
    # TODO: this is *all* an approximation; needs to be worked through
    stopifnot(length(N) == 1) # TODO: support different sizes
    # corrective variation term
    corr_var <- matrix(0, nrow=nrow(covmat), ncol=nrow(covmat))
    binom_var <- (mat[, -ncol(mat)] * (1-mat[, -ncol(mat)]))/N
    off_diag <- abs(row(covmat) - col(covmat)) == 1
    corr_var[off_diag] <- binom_var[off_diag]
    covmat <- covmat + corr_var
  }

  # heterozygosity denominator. Drops the last timepoint, and finds
  #  p_{t} (1 - p_{t})
  hets <- colMeans((mat[, -ncol(mat), drop=FALSE] * 
                   (1-mat[, -ncol(mat), drop=FALSE])), na.rm=TRUE)

  # the heterozygosity is the minimum for the generation that's min(t, s) of
  # the cov(Δp_t, Δp_s)
  hetmat <- hets[pmin(row(covmat), col(covmat))]
  dim(hetmat) <- dim(covmat)
  #attributes(covmat)$nloci <- sum()

  if (standardize)  {
    norm_covmat <- covmat/hetmat
  } else {
    norm_covmat <- covmat
  }

  if (!as_df) return(norm_covmat)

  if (upper_tri) {
    # only look at the upper triangle, as these contain 
    # cov(Δp_t, Δp_s) where s ≥ t
    out_norm_covmat <- norm_covmat[upper.tri(norm_covmat, diag=TRUE)]
    t0 <- row(norm_covmat)[upper.tri(norm_covmat, diag=TRUE)]
    t1 <- col(norm_covmat)[upper.tri(norm_covmat, diag=TRUE)]
    # TODO: this drops time labels, but this may not matter
  } else {
    out_norm_covmat <- norm_covmat
    t0 <- row(norm_covmat)
    t1 <- col(norm_covmat)
  }

  covdf <- tibble(t0=t0, t1=t1, cov=out_norm_covmat)
  colnames(covdf) <- c('t0', 't1', 'cov')
  return(covdf)
}

#' Apply the temporal covariance calculations to simulation results
#' 
#' @param sims a tibble of simulation results
#' @param as_df whether to compress results from matrix to data.frame
#'
#' Simulation results are transposed (e.g. \code{ntimepoints} x \code{nloci})
#' so this function transposes them before passing to \code{temp_cov}.
#'
#' @export
process_covs <- function(sims, as_df=TRUE, remove_fixed=TRUE,
                         standardize=TRUE) {
  # this averages across loci#
  temp_cov(t(sims$neut_freqs), as_df=as_df, swap=TRUE, 
           remove_fixed=remove_fixed, standardize=standardize)
}


process_genic_va <- function(res, as_df=FALSE) {
  if (all(is.na(res$effects))) return(tibble(gen=1:nrow(res$neut_freqs), va=0))
  a <- res$effects[!is.na(res$effects)]
  p <- res$sel_freqs
  het <- 2*p*(1-p)
  va <- sweep(het, 2, a^2, FUN='*')
  if (!as_df) return(rowSums(va))
  return(tibble(gen=1:nrow(res$sel_freqs), genic_va=rowSums(va)))
}

#add_genic_va <- function(stats, res) {
#  a <- res$effects[!is.na(res$effects)]
#  p <- res$sel_freqs
#  het <- 2*p*(1-p)
#  va <- sweep(het, 2, a^2, FUN='*')
#  stats$genic_va <- rowSums(va)
#  stats
#}


## Functions for working with simulated covariance data:

#' Summarize a covariance matrix into components
#' 
#' @param x covariance matrix
#' @param after generations (row/col index) to start sum 
#'
#' @export
sum_covmat <- function(x, after=1, before=ncol(x), bygen=FALSE, use_abs=FALSE) {
  x <- x[after:before, after:before]
  labs <- matrix('cov', nrow=nrow(x), ncol=ncol(x))
  diag(labs) <- 'var'
  if (!bygen) {
    d <- tibble(type=as.vector(labs), vals=as.vector(x)) 
    if (use_abs) 
      d <- d %>% group_by(type) %>% summarize(vals=sum(abs(vals)))
    else 
      d <- d %>% group_by(type) %>% summarize(vals=sum(vals))
  } else {
    gen_t <- row(x) + (after-1)
    gen_s <- col(x) + (after-1)
    if (!use_abs)  {
      d <- tibble(t=as.vector(gen_t), s=as.vector(gen_s),
                  type=as.vector(labs), vals=as.vector(x))  %>% 
      mutate(gen=pmax(t, s)) %>% group_by(type, gen) %>% 
      summarize(var=sum(vals))
    } else {
      d <- tibble(t=as.vector(gen_t), s=as.vector(gen_s),
                  type=as.vector(labs), vals=as.vector(x))  %>% 
      mutate(gen=pmax(t, s)) %>% group_by(type, gen) %>% 
      summarize(var=sum(abs(vals)))
    }
  }
  return(d)
}


#' 
partition_cov <- function(cov) {
  stopifnot(colnames(cov) == c('t0', 't1', 'cov'))
  vars <- filter(cov, t0 == t1) %>% pull(cov) %>% sum()
  covs <- filter(cov, t0 != t1) %>% pull(cov) %>% sum()
  # we double the covariances, since the data here 
  # are only the upper right triangle
  tibble(component=c('var', 'cov'), value=c(vars, 2*covs))
}

cum_cov <- function(cov) {
  stopifnot(colnames(cov) == c('t0', 't1', 'cov'))
  all_covs <- filter(cov, t1 != t0) %>% arrange(t0, t1) %>% mutate(cumcov=cumsum(cov), gen=ifelse(t1 > t0, t1, t0)) %>% select(gen, cumcov)
  all_vars <- filter(cov, t1 == t0) %>% arrange(t0) %>% mutate(cumvar=cumsum(cov), gen=t1) %>% select(gen, cumvar)
  inner_join(all_covs, all_vars, by='gen')
}
