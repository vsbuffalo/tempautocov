## ldcalcfuns.r -- these are mostly for calculating various LDs to try to solve
## the mystery of the E(D1n D2n) dynamics

offdiag <- function(x) {
  diag(x) <- 0
  return(x)
}

calc_D_prod <- function(x) {
  gamlist <- x$gams
  effects <- x$effects
  neutsites <- which(is.na(effects))
  selsites <- which(!is.na(effects))
  out <- lapply(seq_along(gamlist), function(g) {
     gam <- gamlist[[g]]
     dprodsums <- lapply(neutsites, function(ns) {
       mat <- gam[, c(ns, selsites)]
       # multiply by effect sizes
       matef <- sweep(mat, 2, c(1, effects[!is.na(effects)]), FUN='*')
       covmat <- cov(matef) 
       # this is the LD * effect matrix for all sites, including sel/sel sites.
       # we just need to keep the first column minus the first row (D with neutral 
       # site itself)
       Ds <- covmat[-1, 1]
       Dprod <- outer(Ds, Ds)
       tibble(gen = g, neutsite=ns, diag=sum(diag(Dprod)), offdiag=sum(offdiag(Dprod)))
     })
    do.call(rbind, dprodsums) 
  })
  do.call(rbind, out) 
}

get_rec_pos <- function(x, r) {
  selsites <- which(!is.na(x$effects))
  r*x$pos[selsites]
}


calc_Ds <- function(x, gen, effects, pos=NULL, usesign=FALSE, unify=TRUE, polarize_swap=TRUE) {
  # for single gen (included as col here for easier later processing)
  if (polarize_swap) {
    x <- swap_gamete_alleles(x, effects)
  }
  mat <- cov(x)
  # labels are always ordered!
  sitelabs <- paste('D', pmin(row(mat), col(mat)), pmax(row(mat), col(mat)), sep='_')
  # make effect columns for pairwise D
  if (usesign) {
    z <- ifelse(sign(effects) > 0, '+', '-')
    sitetype <- ifelse(is.na(z), 'N', z)
  } else {
    sitetype <- ifelse(is.na(effects), 'N', 'S')
    if (unify)  # make NS = SN for easier grouping
      sitetype <- ifelse(sitetype == "NS", "SN", sitetype)
  }
  labs <- outer(sitetype, sitetype, FUN=paste0)
  vals <- mat[upper.tri(mat, diag=FALSE)] 
  labsmat <- labs
  dim(labsmat) <- dim(mat)
  # distances between sites
  dists <- NA
  if (!is.null(pos)) {
    dists <- abs(outer(pos, pos, FUN='-'))
    dists <- dists[upper.tri(dists, diag=FALSE)]
  }
  out <- tibble(gen=gen, val=c(vals), 
                D=c(labs[upper.tri(labs, diag=FALSE)]),
                dists=dists)
  return(out)
}

threeloci_LD <- function(x, effects) {
  # this is too slow -- 
  neut_sites <- which(is.na(effects))
  sel_sites <- which(!is.na(effects))
  out <- do.call(rbind, lapply(sel_sites, function(i) {
      combs <- cbind(i, t(combn(sel_sites, 2)))
      do.call(rbind, lapply(1:nrow(combs), function(cols) {
         gams <- x[, combs[cols, ]]
         D123 = mean(apply(sweep(gams, 2, colMeans(gams), '-'), 1, prod))
         tibble(D3=paste0('D', paste0(cols, collapse='')), D3val=D123)
      }))
  }))
  out
}

calc_3D <- function(pop, polarize_swap=TRUE) {
  x <- pop$gams
  effects <- pop$effects
  out <- lapply(1:length(x), function(g) {
    if (polarize_swap) {
      y <- swap_gamete_alleles(x[[g]], effects)
      threeloci_LD(y, effects)
    }
  })
  # three way -- E((IA-pA)(IB - pB)(IC - pC))
  tibble(gen=gen, D123=D123)
}

before_any_fix <- function(x) {
  max(cumsum(rowSums(!(x$sel_freqs == 0 | x$sel_freqs == 1)) > 1))
}


calc_all_Ds <- function(x, pos, spread=FALSE, usesign=FALSE) {
  # calc Ds across all generations
  effects <- x$effects
  max_gen <- before_any_fix(x)
  out <- map2_df(x$gams[1:max_gen], 1:max_gen, function(gam, g) calc_Ds(gam, g, effects, x$pos, usesign=usesign))
  if (spread) {
    return(spread(out, D, val))
  }
  return(out)
}

zerodiag <- function(x) {
  diag(x)  <- 0
  return(x)
}

# calculate the sum of all Ds
calc_D_sum <- function(x) {
  # x is the sim res
  neut_sites <- which(is.na(x$effects))
  tmp <- lapply(x$gams, function(g) {
                  # covar between first neutral site and all non-neutral sites
                  d <- cov(g)[-neut_sites[1], neut_sites[1]]
                  sum(zerodiag(d %o% d))
  })
  dd <- as.data.frame(do.call(rbind, tmp))
  dd$gen <- 1:nrow(dd)
  dd
}


effect_signs <- function(x) {
  e <- x$effects
  e <- e[!is.na(e)]
  paste0(c('+', '-')[(e/abs(e) > 0) + 1L], collapse=' ')
}

simapop <- function(r, alpha, N=1e3, freqs=NULL, nsites=2) {
  if (is.null(freqs)) {
    samp <- prime_pop(1, N, r, nsites+10)
    pos <- samp$samp$positions[[1]]
    gams <- samp$samp$gametes[[1]]
  } else {
    gams  <- do.call(cbind, lapply(freqs, function(x) rbinom(N*2, 1, x)))
    pos <- seq(0, 1, length.out=length(freqs))
  }
  #browser()
  pop <- simpop(gams, pos, 40, N, alpha, nsites, r, 
                w_exp_factory(1), include_gams=TRUE, include_allelic_cov=TRUE)
  pop
}


Dsum <- function(x, effects) {
  neut_sites <- which(is.na(effects))
  out <- x %>% gather(D, val, -gen) %>%
    separate(D, into=c('D', 'site1', 'site2', 'site3'), convert=TRUE) %>%
    filter(site1 %in% neut_sites | site2 %in% neut_sites) %>% 
    filter(is.na(site3)) # ditch three sites
  out %>% group_by(gen) %>% summarize(D=sum(val))
}


matrix2tibble <- function(x, upper_tri=TRUE, diag=FALSE, rowlab='row', collab='col', vallab='val') {
  if (upper_tri) {
    i <- upper.tri(x, diag=diag)
    vals <- x[i]
    out <- tibble::tibble(vals=vals, row=row(x)[i], col=col(x)[i])
    colnames(out) <- c(vallab, rowlab, collab)
    return(out)
  } else {
    stop('not implemented yet') 
  }
}

calc_barton_otto_delD <- function(sel_freqs, recpos, s, N=1e3, ignore_fixed=TRUE, finite_only=TRUE) {
  if (ncol(sel_freqs) != length(recpos)) {
    browser()
  }
  delD <- lapply(1:nrow(sel_freqs), function(g) {
    p <- sel_freqs[g, ]
    if (ignore_fixed) {
      fixed <- p == 1 | p == 0
      p[fixed] <- NA
    }
    # matrix of heterozygosities of pairs of selected alleles
    qp2 <- outer(p*(1-p), p*(1-p))
    r <- abs(outer(recpos, recpos, FUN='-'))
    delD <- -2 * qp2 * s^2 * (1-r) / (2*N*r^3)
    if (finite_only) {
      delD[!is.finite(delD)] <- NA
    }
    out <- matrix2tibble(delD, vallab='barton_otto')
    out$gen <- g
    out
  })
  bind_rows(delD)
}

calc_D_diff <- function(x, return_raw=FALSE, polarize_swap=TRUE) {
  # return raw are all D matrices, no diff taken. This is for debugging.
  gamlist <- x$gams
  effects <- x$effects
  neutsites <- which(is.na(effects))
  selsites <- which(!is.na(effects))
  out <- lapply(seq_along(gamlist), function(g) {
     gam <- gamlist[[g]]
     if (polarize_swap) {
       gam <- swap_gamete_alleles(gam, effects)
     }
     mat <- gam[, selsites]
     Dmat <- cov(mat) 
     Dmat
   })
  Ds <- do.call(abind::abind, list(out, along=3)) 
  if (return_raw) return(Ds)
  diffmat <- aperm(apply(Ds, 1:2, diff), c(2, 3, 1))
  z <- ifelse(sign(effects[selsites]) > 0, '+', '-')
  effectsmat <- outer(z,z, paste0)
  bind_rows(lapply(seq_len(dim(diffmat)[3]), function(i) {
           m <- diffmat[, , i]
           vals <- m[upper.tri(m)]
           tibble(gen=i, D=paste0('D', row(m)[upper.tri(m)], col(m)[upper.tri(m)]), 
                  val=vals, sign=effectsmat[upper.tri(effectsmat)])
  }))
}


ave_dprodsum <- function(x) x %>% group_by(gen) %>% 
  summarize(diag=sum(diag), offdiag=sum(offdiag))


Cijp <- function(rij, Cij, pi, pj, si, sj)  (1-rij) * Cij * (1 + (1-2*pi)*si + (1-2*pj)*sj)
dpip <- function(si, sj, pi, Cij) si * pi * (1-pi) + sj * Cij
dpnp  <- function(si, sj, Cin, Cjn) si*Cin + sj*Cjn
Cinp <- function(rin, Cin, pi, si, sj, Cijn) (1-rin)*(Cin*(1 + (1 - 2*pi)*si) + sj*Cijn)
Cijnp <- function(rin, rjn, Cijn, pi, pj, si, sj, Cin, Cjn, Cij)  (1-rin)*(1-rjn)*(Cijn * (1 + (1-2*pi)*si + (1-2*pj)*sj) - 2*Cij*(si*Cin + sj*Cjn))



run_ld_dyns <- function(nsteps, pi0, pj0, pn0, rin, rjn, rij=rin + rjn, Cin0=ok71(4e3*rin), Cjn0=ok71(4e3*rjn), 
                        si=0.01, sj=si, Cijn0=0, Cij0=0, dCij0=NULL) {
  state <- matrix(NA, nrow=nsteps, ncol=14)
  pialt <- pi0; pjalt <- pj0; Cijalt <- Cij0; Cinalt <- Cin0; Cjnalt <- Cjn0
  state[1, ] <- c(0, pn0, pi0, pj0, Cin0, Cjn0, Cij0, Cijn0, pialt, pjalt, Cij0, Cin0, Cjn0, Cijn0, NA)
  pn <- pn0; pi <- pi0; pj <- pj0; Cin <- Cin0; Cjn <- Cjn0; Cij <- Cij0; Cijn <- Cijn0
  for (i in seq_len(nsteps-1)) {
#    stopifnot(pi >= 0 & pi <= 1)
#    stopifnot(pj >= 0 & pj <= 1)
#    stopifnot(pn >= 0 & pn <= 1)
    # update freqs, e.g. as they'd change in viability selection
    pi <- dpip(si, sj, pi, Cij) + pi
    pj <- dpip(sj, si, pj, Cij) + pj
    pialt <- dpip(si, sj, pialt, 0) + pialt
    pjalt <- dpip(sj, si, pjalt, 0) + pjalt
    if (pi %in% c(0, 1) || pj %in% c(0, 1)) break
    pn <- dpnp(si, sj, Cin, Cjn) + pn
    # now changes in LD after meiosis 
    Cij <- Cijp(rij, Cij, pi, pj, si, sj)
    Cin <- Cinp(rin, Cin, pi, si, sj, Cijn)
    Cjn <- Cinp(rjn, Cjn, pj, sj, si, Cijn)
    Cijn <- Cijnp(rin, rjn, Cijn, pi, pj, si, sj, Cin, Cjn, Cij)
    # alts
    Cijalt <- Cijp(rij, Cijalt, pialt, pjalt, si, sj)
    Cinalt <- Cinp(rin, Cinalt, pialt, si, sj, 0)
    Cjnalt <- Cinp(rjn, Cjnalt, pjalt, sj, si, 0)
    if (!is.null(dCij0)) {
      Cij <- Cij - 2 * si * sj * pj * (1-pj) * pi * (1-pi) * (1-rij) / (2e3 * rij^3)
    }
    state[i,] <- c(i, pn, pi, pj, Cin, Cjn, Cij, Cijn, pialt, pjalt, Cijalt, Cinalt, Cjnalt, Cijalt)
  }
  colnames(state) <- c('gen', 'pn', 'pi', 'pj', 'Cin', 'Cjn', 'Cij', 'Cijn', 'pialt', 'pjalt', 
                       'Cijalt', 'Cinalt', 'Cjnalt', 'Cijnalt')
  as_tibble(state)
}



