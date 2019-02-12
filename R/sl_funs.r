## sl-funs.r -- functions for single locus analysis

get_sel_freqs <- function(res) {
  tibble(gen=1:nrow(res$sel_freqs), freq=res$sel_freqs[,1])
}


get_neut_freqs <- function(res) {
  tibble(gen=1:nrow(res$neut_freqs), freq=res$neut_freqs[,1])
}

reshape_neut_freqs <- function(x) {
  out <- x %>% unnest(neut_freqs) %>% select(-id) %>%
    spread(rep, freq) 
  t(as.matrix(out[,-1]))
}

fixes <- function(res) {
  any(res$neut_freqs == 0 | res$neut_freqs== 1)
}

sl_va <- function(p, s, h=1/2) {
  q <- 1-p
  a11 <- 1 + s; a12 <- 1+h*s; a22 <- 1    
  alpha1 <- p * a11 + q*a12
  alpha2 <- p * a12 + q*a22    
  2*p*q * (alpha2 - alpha1)^2  
}

 
sel_freq <- Vectorize(function(s, p0, t, h=1/2, last=FALSE) {
  if (t == 0) return(p0)
  w11 <- (1+s)
  w12 <- (1+h*s)
  freqs <- numeric(t+1)
  freqs[1] <- p0
  for (i in seq(2, t+1)) {
    p <- freqs[i-1]
    p_prime <- (w11*p^2 + w12*p*(1-p))/wbar(p, s, h)
    freqs[i] <- p_prime
  }
  if (!last) return(freqs)
  return(freqs[length(freqs)])
}, c('s', 'p0', 't'))



sl_var_traj <- function(p0, t, s, h=1/2) {
  p <- sel_freq(s, p0, t)[,1]
  tibble::tibble(gen=seq_along(p), theory_va=sl_va(p, s, h))
}
 
wbar <- function(p, s, h=1/2) (1+s)*p^2 + 2*(1+h*s)*p*(1-p) + (1-p)^2

Va <- Vectorize(function(p, s, h=1/2){
  wb = wbar(p, s, h)
  (1+s-wb)^2 * p^2 + (1+h*s - wb)^2 * 2*p*(1-p) + (1 - wb)^2 *(1-p)^2
}, c('p', 's'))
   
