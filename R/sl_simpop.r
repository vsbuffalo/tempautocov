# sl-simpop.r 

sl_recombine <- function(g1, g2, rec_frac) {
  # sample break point positions 
  xover <- rbinom(1, 1, rec_frac)
  seg_state <- rbinom(1, 1, 0.5)
  segs <- c(seg_state, ifelse(xover, !seg_state, seg_state))
  as.integer((g1 & segs) | (g2 & !segs))
}


sl_random_mate <- function(gams, parents, rec_frac) {
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
    sl_recombine(g1, g2, rec_frac)
  }, parents)

  new_gams <- do.call(rbind, haps)
  new_gams
}


sl_prime_pop <- function(nreps, N, rec_frac, seeds=NULL) { 
  rho <- 4*N*rec_frac
  # intiate with mspms
  stopifnot(nreps == 1)  # not supported currently
  msprime <- FALSE
  if (msprime) {
    enough_sites <- FALSE
    tries <- 0
    while(!enough_sites) {
      # theta sufficiently high that two sites are placed on tree
      if  (is.null(seeds))
      	samp <- call_ms(2*N, nreps, t=100, r=c(rho, 2), ms='mspms') %>% parse_ms()
      else
      	samp <- call_ms(2*N, nreps, .t=100, .r=c(rho, 2), .seeds=seeds, strict=TRUE,
			      ms='mspms') %>% parse_ms()
      enough_sites <- all(samp$segsites >= 2)
      spanned_break <- samp$positions[[1]][1] < 0.5 & samp$positions[[1]][length(samp$positions[[1]])] > 0.5
      if (enough_sites && spanned_break) break
      tries <- tries + 1
      if (tries > 10) stop('increase theta!')
      warning('theta not sufficiently high to generate two sites!')
    }

  } else {
    if (is.null(seeds)) {
      samp <- call_ms(2*N, nreps, t=100, r=c(rho, 2), 
                      ms='ms') %>% parse_ms()
    } else {
      samp <- call_ms(2*N, nreps, .t=100, .r=c(rho, 2), 
                      .seeds=seeds, strict=TRUE, ms='ms') %>% parse_ms()
      spanned_break <- samp$positions[[1]][1] < 0.5 & samp$positions[[1]][length(samp$positions[[1]])] > 0.5
      stopifnot(spanned_break)
    }

  }
  list(samp=samp)
}


sl_simpop <- function(gams, ngens, N, alpha, rec_frac, sel_init_gen=5,
                      include_rsq=TRUE, include_D=FALSE, verbose=FALSE) {
  summaries <- vector('list', length(ngens))

  # swap alleles, so the gametes from MS which are ancestral/derived
  # are randomized
  #gams <- sweep(gams, 2, rbinom(2, 1, 0.5), function(x, y) abs(x-y)) # 

  for (t in seq_len(ngens)) {
    if (verbose) message(sprintf('generation %d of %d', t, ngens))
    dips <- gam2diploid(gams)

    if (t >= sel_init_gen) {
      z <- dips[, 2]*alpha  # second column is selected site
      w <- 1 + z
    } else {
      z <- 0
      w <- rep(1, N)
    }
    wbar <- mean(w)
    p <- w/(N*wbar)
    # if (!all(is.finite(p))) browser()

    # sample parents of the offspring
    parents <- sample(seq_len(N), 2*N, replace=TRUE, prob=p)
  
    # create the recombined gametes
    new_gams <- sl_random_mate(gams, parents, rec_frac)
    stopifnot(dim(new_gams) == dim(gams))
    gams <- new_gams

    # summaries
    par_dist <- table(factor(parents, levels=seq_len(N)))
    stats <- tibble(gen=t, 
                    wbar=wbar, wvar=var(w),
                    #zvar=var(z),
                    kvar=var(par_dist), 
                    kbar=mean(par_dist))
    rsq <- NULL
    D <- NULL
    if (include_rsq) {
      rsq <- calc_rsq(gams, 1, ncol(gams))
      rsq$gen <- t
    }
    if (include_D) {
      D <- calc_D(gams, 1, ncol(gams))
      D$gen <- t

    }
    sel_freqs <- calc_freqs(gams[, 2, drop=FALSE])
    neut_freqs <- calc_freqs(gams[, 1, drop=FALSE])
    summaries[[t]] <- list(stats=stats, 
                           rsq=rsq,
                           D=D,
                           sel_freqs=sel_freqs, 
                           neut_freqs=neut_freqs,
                           pos=c(0, 1))
  }
  stats <- bind_rows(map(summaries, 'stats'))
  rsq <- bind_rows(map(summaries, 'rsq'))
  D <- bind_rows(map(summaries, 'D'))
  sel_freqs <- do.call(rbind, map(summaries, 'sel_freqs'))
  neut_freqs <- do.call(rbind, map(summaries, 'neut_freqs'))
  geno <- t(apply(dips, 2, function(x) table(factor(x, levels=0:2))))
  list(stats=stats, rsq=rsq, D=D, sel_freqs=sel_freqs, geno=geno,
       neut_freqs=neut_freqs)
}


