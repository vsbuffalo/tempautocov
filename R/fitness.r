## fitness.r -- check math of exponential and other fitness models

trait2offspring <- function(N=1000, sigma2=1) {
  # sigma2 is trait variation, for the trait that goes 
  # into fitness function
  z <- rnorm(N, 0, sqrt(sigma2))
  wz <- exp(z)
  wbar <- mean(wz)
  p <- wz/(N*wbar)
  stopifnot(abs(sum(p) - 1) < 1e-5)
  k <- rmultinom(1, N, p)
  k
}

## analytic
offspring_var <- function(N, mu, sigma2) {
  wbar <- exp(sigma2/2)
  vwz <- exp(2*(mu + sigma2)) - 2*exp(mu+sigma2/2)*wbar + wbar^2
  Vp <- vwz/(N*wbar)^2
  Ep <- 1/N
  N*Ep*(1-Ep) + N*(N-1)*Vp
}

## simulation
sim_trait2offspring <- function(nreps, N=1000, sigma2=1) {
  replicate(nreps, {
              d <- trait2offspring(sigma2=sigma2)
              var(d)
  }, simplify=TRUE) 
}


### fitness approximations
mult_fit <- function(s, nreps=50, N=1000, L=1000, trait_mu=3e-4) {
  theta <- 4*N*trait_mu
  bind_rows(replicate(nreps, {
     theta_per_L <- theta/L
     g <- map_dbl(rbeta(L, theta_per_L, theta_per_L), ~ rbinom(1, 2, .))
     tibble(mult=prod(1 + s*g), approx=1 + 2*sum(s*g))
  }, simplify=FALSE))
}


### Trunctation Selection
trunc_Va <- function(prop, h2, Va0, gens) {
  # from Bulmer, p. 154
  z <- qnorm(prop)
  cc <- dnorm(z)/prop  # truncated standard normal
  Vas <- numeric(gens)
  Vas[1] <- Va0
  for (t in 2:gens) {
    Vas[t] <- 0.5*(1 - h2 * cc * (cc-z))*Vas[t-1] + 0.5*Va0
  }
  Vas
}

