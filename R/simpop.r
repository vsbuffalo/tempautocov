# simpop.r -- key simulation functions for the neutral burnin model
library(tidyverse)
library(reshape2)
library(msr)
library(parallel)

# These are used to calculate α given a particular Va and L
# Initially I thought finite-sites θ was needed, but this did 
# not work as well as assuming SSH = π = θ.
ssh <- function(theta, L=Inf) theta/(1 + 2*theta/L) # Gillespie, p. 54 (not used)
#alpha <- function(Va, theta, L) sqrt(Va/(ssh(theta, L)))
alpha <- function(Va, theta) sqrt(Va/theta)

rand_seed <- function(n=3) {
  x <- runif(3)
  as.integer(x * .Machine$integer.max)
}

w_constant <- function(w, t) rep(1, length(w))

w_truncation_factory <- function(start, p) {
  fun <- function(z, t) {
    # select the top p proportion
    if (t >= start)
      return(ifelse(z >= quantile(z, 1-p), 1, 0))
    else
      return(rep(1, length(z)))
  }
  return(fun)
}

w_exp_fluct_factory <- function(start, flip) {
  fun <- function(z, t) {
    s <- 0
    if (t >= start && t < flip)
      s <- 1
    if (t >= flip)
      s <- -1
    exp(s*z)
  }
  return(fun)
}


w_exp_factory <- function(gen, s_new=1) {
  fun <- function(z, t) {
    s <- 0
    if (t >= gen)
      s <- s_new
    exp(s*z)
  }
  return(fun)
}

w_gss_factory <- function(gen, o=0, Vs=1) {
  fun <- function(z, t) {
    opt <- 0
    if (t >= gen)
      opt <- o
    exp(-(z-opt)^2 / (2*Vs))
  }
  return(fun)
}


gam2diploid <- function(gam) {
  # use the pre-built linear transformation matrix trans_mat to combine rows of
  # the gamete matrix into diploids, adding rows to gene content
  L <- ncol(gam)
  t(sapply(seq(1, nrow(gam), 2), function(i) .colSums(gam[i:(i+1L), ], 2, L)))
}


recombine <- function(g1, g2, r, pos) {
  nbreaks <- rpois(1, r)
  breaks <- integer(length(g1))
  if (nbreaks == 0) {
    return(list(g1, g2)[[sample(1:2, 1)]])
  }
  # sample break point positions 
  break_pos <- runif(nbreaks)
  breaks <- rowSums(sapply(sort(break_pos), function(p) p < pos))
  segs <- as.integer(breaks %% 2 == 0)
  if (rbinom(1, 1, prob=0.5)) {
    segs <- as.integer(!segs)
  }
  as.integer((g1 & segs) | (g2 & !segs))
}


random_mate <- function(gams, pos, parents, r) {
  N <- nrow(gams)/2
  # mapping of individuals to their gametes
  #ind_gam_map <- split(1:nrow(gams), rep(1:N, each=2))

  haps <- Map(function(par_i) {
    # get the index of the parent's two gametes
    #par_i <- ind_gam_map[[i]]
    #g1 <- gams[par_i[1], ]
    #g2 <- gams[par_i[2], ]
    g1 <- gams[2L * par_i - 1L, ]
    g2 <- gams[2L * par_i, ]
    # stopifnot(!(0 %in% dim(g1)))
    # stopifnot(!(0 %in% dim(g2)))
    recombine(g1, g2, r, pos)
  }, parents)

  new_gams <- do.call(rbind, haps)
  new_gams
}

target_theta <- function(nsites, n) {
  # from Wakely, E(S) = \theta \sum_i^{n-1} 1/i \approx \theta(\log(n) + \gamma)
  γ = 0.5772156649
  nsites/(γ + log(n))
}

#' Run msprime to generate gametes for a population
#' 
#' @param nreps number of replicates
#' @param N population size
#' @param r recombination fraction between ends of segment (see mspms doc)
#' @param total_sites the total target number of sites
#' @param inflation factor by which the target number of sites are inflated by
#'
#' We inflate the target number of segregating sites by the inflation factor, as
#' many neutral simulations will not have the required number of sites (since
#' it's stochastic). We could re-run until have the right number of sites, 
#' but this is conditioning on a lower bound. In our simulations we 
#' allocate the number of selected sites, and the rest are neutral.
#'
#' We use 5*total_sites as the number of sites for MS's crossover regime.
#' The idea here is that crossovers over continuous regions may great events that
#' would not in fact change the crossover in basepairs. By having a finite-sites model 
#' for recombination, these are ignored. Here we just need to ensure sufficent resolution 
#' that an event between two mutations is called a recombination, e.g. that these two 
#' sites aren't in one "basepair" bin.

prime_pop <- function(nreps, N, r, total_sites, seeds=NULL, inflation=1.5) {
  rho <- 4*N*r
  # intiate with mspms
  theta <- ceiling(target_theta(inflation*total_sites, 2*N))
  if (is.null(seeds)) { 
    samp <- call_ms(2*N, nreps, t=theta, r=c(rho, 5*total_sites), 
                    ms='mspms') %>% parse_ms()
  } else {
    samp <- call_ms(2*N, nreps, t=theta, r=c(rho, 5*total_sites), 
		    random_seed=seeds,
                    ms='mspms') %>% parse_ms()

  }
  list(samp=samp, theta=theta)
}

#' Simulate a population, using an initial set of specified burnin gametes
#'
#' @param gams
#' @param pos
#' @param ngens how long to run forward simulation for
#' @param N population size (in diploids)
#' @param alpha the effect size (absolute value)
#' @param sel_sites the number of selected sites
#' @param r genetic distance (expected number of breakpoints)
#' @param wfunc the fitness function
#' @param include_rsq include rsq calculation
#' @param include_D include D calculation
#' @param include_ld for neutral site, include average D to all selected sites
#' @param include_gams include the gametes each population (very memory intensive!)
#' @param verbose be verbose
#'
#' The number of sites coming from the gametes are divided into:
#' sel_sites selected loci, and the rest are neutral
simpop <- function(gams, pos, ngens, N, alpha, sel_sites, r, wfunc, 
                   include_rsq=FALSE, include_D=FALSE, 
		   include_ld=FALSE,
		   include_allelic_cov=FALSE,
                   include_gams=FALSE,
		   verbose=FALSE) {
  # input are gametes from coalescent simulation
  # choose selected sites
  sel_idx <- logical(ncol(gams))
  sel_idx[sample(seq_along(sel_idx), sel_sites)] <- TRUE
  effects  <- sample(c(-alpha, alpha), length(sel_idx), replace=TRUE)
  
  summaries <- vector('list', length(ngens))

  for (t in seq_len(ngens)) {
    if (verbose) message(sprintf('generation %d of %d', t, ngens))
    dips <- gam2diploid(gams)
    zmat <- sweep(dips, 2, sel_idx*effects, FUN='*')
    z <- rowSums(zmat)
    w <- wfunc(z, t)
    wbar <- mean(w)
    p <- w/(N*wbar)

    # sample parents of the offspring
    parents <- sample(seq_len(N), 2*N, replace=TRUE, prob=p)
  
    # create the recombined gametes
    new_gams <- random_mate(gams, pos, parents, r)
    stopifnot(dim(new_gams) == dim(gams))
    gams <- new_gams

    # summaries
    par_dist <- table(factor(parents, levels=seq_len(N)))
    stats <- tibble(gen=t, 
                    wbar=wbar, wvar=var(w),
                    zbar=mean(z), zvar=var(z), 
                    kvar=var(par_dist), 
                    kbar=mean(par_dist))
    rsq <- NULL
    D <- NULL
    LD <- NULL
    neut_ld <- NULL
    if (include_rsq) {
      rsq <- calc_rsq(gams, 1, ncol(gams))
      rsq$gen <- t
    }
    if (include_D) {
      D <- calc_D(gams, 1, ncol(gams))
      D$gen <- t
    }
    if (include_ld) {
      ldg <- cov(gams)^2
      neut_ld_raw <- colMeans(ldg[sel_idx, !sel_idx])
      # commented out version isn't averaged
      #neut_ld <- tibble(gen=t, neutsite=seq_along(neut_ld_raw), ld=neut_ld_raw)
      neut_ld <- tibble(gen=t, ld=mean(neut_ld_raw))
    }
    allelic_cov <- NULL
    if (include_allelic_cov) {
      ac <- cov(sweep(gams[, sel_idx], 2, effects[sel_idx], FUN='*'))
      
      allelic_cov <- tibble(ac_var=sum(diag(ac)), 
			    ac_cov=sum(ac[upper.tri(ac)]) + sum(ac[lower.tri(ac)]),
			    gen=t)
    }

    all_gams <- NULL
    if (include_gams) {
      all_gams <- gams
    }

    sel_freqs <- calc_freqs(gams[, sel_idx, drop=FALSE])
    neut_freqs <- calc_freqs(gams[, !sel_idx, drop=FALSE])
    summaries[[t]] <- list(stats=stats, 
                           rsq=rsq,
                           D=D,
			   neut_ld=neut_ld,
			   allelic_cov=allelic_cov,
                           sel_freqs=sel_freqs, 
                           gams=all_gams,
                           neut_freqs=neut_freqs)
  }
  effects <- sel_idx*effects
  effects[effects == 0] <- NA
  stats <- bind_rows(map(summaries, 'stats'))
  rsq <- bind_rows(map(summaries, 'rsq'))
  D <- bind_rows(map(summaries, 'D'))
  neut_ld <- bind_rows(map(summaries, 'neut_ld'))
  all_gams <- map(summaries, 'gams')
  allelic_cov <- bind_rows(map(summaries, 'allelic_cov'))
  sel_freqs <- do.call(rbind, map(summaries, 'sel_freqs'))
  neut_freqs <- do.call(rbind, map(summaries, 'neut_freqs'))
  geno <- t(apply(dips, 2, function(x) table(factor(x, levels=0:2))))
  list(stats=stats, rsq=rsq, D=D, neut_ld=neut_ld, sel_freqs=sel_freqs, geno=geno,
       neut_freqs=neut_freqs, effects=effects, allelic_cov=allelic_cov, pos=pos,
       gams=all_gams)
}

het <- function(x) {
  # x is a matrix of ntimepoints x nloci
  rowSums(2*x*(1-x))
}


#add_genic_va <- function(stats, res) {
#  a <- res$effects[!is.na(res$effects)]
#  p <- res$sel_freqs
#  va <- sweep(het(p), 2, a^2, FUN='*')
#  stats$genic_va <- rowSums(va)
#  stats
#}

sum_sitehet <- function(x) rowSums(2 * x * (1-x))

process_het <- function(x) {
 neut_freqs <- x$neut_freqs
 sel_freqs <- x$sel_freqs
 if (is.null(sel_freqs)) {
   return(tibble(gen=1:nrow(neut_freqs), 
        neut_ssh=sum_sitehet(neut_freqs),
        nneut=ncol(neut_freqs)))

 }
 tibble(gen=1:nrow(neut_freqs), 
        neut_ssh=sum_sitehet(neut_freqs), sel_ssh=sum_sitehet(sel_freqs), 
        nneut=ncol(neut_freqs), nsel=ncol(sel_freqs))
}


## functions for emulating finite sampling
add_binomial_noise <- function(x, size) {
  out <- rbinom(length(x), size, x)/size
  dim(out) <- dim(x)
  out
}

# process covariances with finite sampling
process_sampled_covs <- function(neut_freqs, as_df=TRUE, remove_fixed=TRUE,
                         standardize=TRUE, sample_size=NULL, use_conditional_variance=FALSE) {
  # this averages across loci#
  temp_cov(t(neut_freqs), as_df=as_df, swap=TRUE,
           remove_fixed=remove_fixed, standardize=standardize, 
           sample_size=sample_size, use_conditional_variance=use_conditional_variance)
}


