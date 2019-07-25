## ld.r -- empirical neutral LD calculations
# these functions, assuming LD levels determined solely by drift and mutation,
# calculate the averge LD in a region between all selected sites with all neutral
# sites


calc_assoc_sum2 <- function(x, genlen, N, t=0) {
  # x is results list from simulations
  sel_sites <- x$pos[!is.na(x$effects)]
  neut_sites <- x$pos[is.na(x$effects)]
  r <- haldane(abs(genlen*outer(sel_sites, neut_sites, FUN='-')))
  sapply(t, function(tt) mean(ok71(4*N*r)*(1-r)^tt))
}

calc_assoc_sum <- function(pos, genlen, N, t=0, prune=100) {
  r <- haldane(abs(genlen*outer(pos, pos, FUN='-')))
  sapply(t, function(tt) mean(ok71(4*N*r)*(1-r)^tt))
}
