## mom.r -- method of moments 

SIM_SSH_COLS <- c('gen', 'neut_ssh', 'sel_ssh', 'nneut', 'nsel')

make_neutld_df <- function(cov, R, N, ssh=NULL, pos=NULL, trange=NULL) {
 if (!is.null(trange)) {
    cov <- cov[trange[1]:trange[2], trange[1]:trange[2]]
    if (!is.null(ssh)) {
      if (is.data.frame(ssh) && 
          colnames(ssh) == SIM_SSH_COLS) {
        # if a dataframe with my sim columns, grab appropriate data
        ssh <- ssh$neut_ssh
      }
      ssh <- ssh[trange[1]:(trange[2]+1)]
      stopifnot(nrow(cov) == length(ssh)-1L)
    }
  }
  # knock off last timepoint of SSH, since we lose this in covs
  #ssh <- ssh[-length(ssh)]
  # for LD decay:
  E <- abs(row(cov) - col(cov))  # elapsed generations between covar
  # for VA:
  T <- pmin(row(cov), col(cov))  # the minimum of col, row indices
  t <- row(cov)
  t <- t[upper.tri(t, diag=TRUE)]
  s <- col(cov)
  s <- s[upper.tri(s, diag=TRUE)]

  # association integral
  # if the positions aren't given (=NULL), then we use theory
  if (is.null(pos)) {
    A <- assoc_int(E, R, N)
  } else {
    A <- calc_assoc_sum(pos, R, N, t=T)
  }
  dim(A) <- dim(cov)
  # get upper tri
  assoc <- A[upper.tri(A, diag=TRUE)]


  # ssh decay. 
  if (!is.null(ssh)) {
    sshr <- ssh/ssh[1]
    ssh_mat <- sshr[T+1] 
    dim(ssh_mat) <- dim(cov)
    utri_ssh_mat <- ssh_mat[upper.tri(ssh_mat, diag=TRUE)] 
    assoc <- utri_ssh_mat * assoc
  }

  assoc <- assoc / 2 # pick up Va/2 factor

  # empirical covs, upper triangle with diag
  covs <- cov[upper.tri(cov, diag=TRUE)] 

  # diagonal drift component
  I <- diag(nrow(cov))
  B <- I[upper.tri(I, diag=TRUE)] 
  tibble(cov = covs, assoc=assoc, B=B, min_gen=T[upper.tri(T, diag=TRUE)], t=t, s=s)
}


fit_mom <- function(df) {
  fit <- lm(cov ~ assoc + B + 0, data=df)
  fit
}

fit_mom_nls <- function(df, ...) {
  fit <- nls(cov ~ assoc*SSlogis(gen, va0, tmid, k) + B + 0, data=df, ...)
  fit
}

fit2params <- function(x) {
  cfs <- coef(x)  
  F <- cfs[2]
  va0 <- cfs[1]
  tibble(N_est=1/(2*F), va0_est=va0)
}


## Fractions of var due to LS
total_var <- function(freqs, start=1, end=nrow(freqs)) {
  freqs <- freqs[start:end, ]
  freqs <- fixations2NAs(freqs)
  #dmat <- calc_deltas(mat)
  p0 <- freqs[1, ]
  t <- (nrow(freqs)-1)
  out <- (var(freqs[nrow(freqs), ] - freqs[1, ], na.rm=TRUE) / (t*mean(p0*(1-p0), na.rm=TRUE)))
  out
}

G <- function(covs_mat, start=1, end=nrow(covs_mat), use_abs=FALSE) {
  if (use_abs) {
    covs_mat <- abs(covs_mat[start:end, start:end])
  } else { 
    covs_mat <- covs_mat[start:end, start:end]
  }
  #sum(offdiag(covs_mat))/sum(covs_mat)
  (sum(covs_mat) - sum(diag(covs_mat)))/sum(covs_mat)
}

Gprime <- function(std_var, Nhat) {
  #std_var = var(pt - p0) / (t p0(1-p0)) = 1/2N
  #1 - (1/Nhat) * 2*std_var
  1 - 1/(Nhat * 2*std_var)
}


# non linear stuff -- experimental 

va_logistic_sse_grad <- function(par, cov, assoc, t, B) {
  if (is.list(par)) {
    par <- unlist(par)  # for debugging 
  }
  # parameters
  #va0 <- par['va0']; g <- par['g']; tinfl <- par['tinfl']; F <- par['F']
  va0 <- par[1]; g <- par[2]; tinfl <- par[3]; F <- par[4]
  Sum <- sum
  Power <- function(b, x) b^x
  E <- exp(1)
  out <- c(
     Sum((2*va0*Power(assoc,2))/Power(1 + Power(E,-(g*(-tinfl + t))),2) + 
         (2*F*assoc*B)/(1 + Power(E,-(g*(-tinfl + t)))) - 
         (2*assoc*cov)/(1 + Power(E,-(g*(-tinfl + t))))),

     Sum((-2*tinfl*Power(va0,2)*Power(assoc,2))/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),3)) - 
         (2*F*tinfl*va0*assoc*B)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2)) + 
         (2*tinfl*va0*assoc*cov)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2)) + 
         (2*Power(va0,2)*Power(assoc,2)*t)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),3)) + 
         (2*F*va0*assoc*B*t)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2)) - 
         (2*va0*assoc*cov*t)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2))),

     Sum((-2*g*Power(va0,2)*Power(assoc,2))/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),3)) - 
         (2*F*g*va0*assoc*B)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2)) + 
         (2*g*va0*assoc*cov)/
         (Power(E,g*(-tinfl + t))*Power(1 + Power(E,-(g*(-tinfl + t))),2))),

           Sum((2*va0*assoc*B)/(1 + Power(E,-(g*(-tinfl + t)))) + 
               2*F*Power(B,2) - 2*B*cov)
           )
  out
}

va_logistic_sse <- function(par, cov, assoc, t, B) {
  # parameters
  va0 <- par['va0']; g <- par['g']; tinfl <- par['tinfl']; F <- par['F']
  # expression
  yhat <- assoc * va0/(1 + exp(-g * (t - tinfl))) + F*B
  sse <- sum((cov - yhat)^2)
  sse
}

optim2logistic <- function(res, t=1:50, as_df=FALSE) {
  va <- va_logistic(res$par['va0'], t, res$par['g'], res$par['tinfl'], 0)
  if (as_df) return (tibble(t=t, va=va))
  va
}


# fit_mom_nls <- function(data, use_grad=TRUE, ..., max_gens=100) {
#   grid <- tidyr::crossing(va0=0.01,
#                           g=seq(0, 1, length.out=10),
#                           tinfl=seq(0, max_gens, length.out=10),
#                           F=1/seq(1, 1e6, length.out=5))
#   out <- apply(grid, 1, va_logistic_sse, 
#                cov=data$cov, assoc=data$assoc, t=data$t, B=data$B)
#   i <- which.min(out)
#   init_params <- grid[i, ]
#   if (use_grad) {
#     out <- optim(par=init_params, fn=va_logistic_sse, gr=grad, 
#                  cov=data$cov, assoc=data$assoc, t=data$min_gen, B=data$B, ...)
#   } else {
#     out <- optim(par=init_params, fn=va_logistic_sse, 
#                  cov=data$cov, assoc=data$assoc, t=data$min_gen, B=data$B, ...)
#   }
#   out
# }

