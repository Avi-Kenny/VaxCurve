#' Logit function
#' 
#' @param x Numeric input
#' @return Numeric output
logit <- function(x) {
  log(x/(1-x))
}



#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {
  1 / (1+exp(-x))
}



#' Derivative of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_expit <- function(x) {
  exp(x) / ((1+exp(x))^2)
}



#' Derivative of logit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_logit <- function(x) {
  1 / (x-x^2)
}



#' Probability of sampling
#' 
#' @param y Vector `y` of dataset returned by generate_data()
#' @param w1 Vector `w1` of dataset returned by generate_data()
#' @param w2 Vector `w2` of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
Pi <- function(delta_star, w1, w2) {
  
  pi <- function(w1,w2) { expit(w1+w2) }
  return(delta_star + (1-delta_star)*pi(w1,w2))
  
}



#' Return IP weights
#' 
#' @param dat Dataset returned by generate_data()
#' @param scale One of c("none", "sum 1", "mean 1")
#' @return A sum-to-one vector of weights
#' @notes
#'   - Pi is accessed globally
wts <- function(dat, scale) {
  
  if (attr(dat, "sampling")=="iid") {
    weights <- rep(1, nrow(dat))
  } else if (attr(dat, "sampling")=="two-phase") {
    weights <- 1 / Pi(dat$delta_star, dat$w1, dat$w2)
  }
  
  if (scale=="sum 1") {
    weights <- weights / sum(weights)
  }
  
  if (scale=="mean 1") {
    weights <- weights / mean(weights)
  }
  
  return(weights)
  
}



#' Return IP weight stabilization factor
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return IP weight stabilization factor "s"
#' @notes
#'   - wts is accessed globally (and therefore so is Pi)
ss <- function(dat_orig) {
  
  if (attr(dat_orig, "sampling")=="iid") {
    return(1)
  } else if (attr(dat_orig, "sampling")=="two-phase") {
    n_orig <- nrow(dat_orig)
    dat_orig %<>% filter(!is.na(a))
    weights <- wts(dat_orig, scale="none")
    s <- sum(weights)/n_orig
    return(s)
  }
  
}



#' Construct conditional survival estimator S_n
#' 
#' @param dat Dataset returned by generate_data(); accepts either full data or
#'     truncated data
#' @param type Currently only "Cox PH" is implemented
#' @param csf Logical; if TRUE, estimate the conditional survival
#'     function of the censoring distribution instead
#' @return Conditional density estimator function
construct_S_n <- function(dat, type, csf=FALSE) {
  
  if (type=="Cox PH") {
    
    weights <- wts(dat, scale="mean 1")
    
    if (csf) dat$delta_star <- 1-dat$delta_star
    
    # Fit Cox model
    model <- coxph(
      Surv(y_star, delta_star)~w1+w2+a,
      data = dat,
      weights = weights
    )
    coeffs <- model$coefficients
    
    # Get cumulative hazard estimate
    bh <- basehaz(model, centered=FALSE)
    
    return(memoise(Vectorize(function(t, w1, w2, a){
      
      index <- which.min(abs(bh$time-t))
      H_0 <- bh$hazard[index]
      lin <- coeffs[[1]]*w1 + coeffs[[2]]*w2 + coeffs[[3]]*a
      return(exp(-1*H_0*exp(lin)))
      
    })))
    
  }
  
}



#' Construct g-computation estimator function of theta_0
#' 
#' @param dat Dataset returned by generate_data(); must be full data
#' @param mu_n A regression estimator returned by construct_mu_n()
#' @return G-computation estimator of theta_0
construct_gcomp <- function(dat_orig, S_n) {
  
  return(memoise(Vectorize(function(a) {
    1 - mean(S_n(C$t_e, dat_orig$w1, dat_orig$w2, a))
  })))
  
}



#' Construct derivative estimator theta'_n
#' 
#' @param gcomp_n G-comp estimator of theta_0 returned by construct_gcomp()
construct_deriv_theta_n <- function(gcomp_n) {
  
  deriv_theta_n <- function(a) {
    
    # Set derivative appx x-coordinates
    width <- 0.05
    p1 <- a - width/2
    p2 <- a + width/2
    if (p1<0) {
      p2 <- p2 - p1
      p1 <- 0
    }
    if (p2>1) {
      p1 <- p1 - p2 + 1
      p2 <- 1
    }
    
    return( (gcomp_n(p2)-gcomp_n(p1))/width )
    
  }
  
  return(memoise(Vectorize(deriv_theta_n)))

}



#' Construct tau_n Chernoff scale factor function
#' 
#' @param deriv_theta_n A derivative estimator returned by
#'     deriv_theta_n()
#' @param gamma_n Nuisance function estimator returned by construct_gamma_n()
#' @param f_a_n Density estimator returned by construct_f_a_n()
#' @return Chernoff scale factor estimator function
construct_tau_n <- function(deriv_theta_n, gamma_n, f_a_n) {
  
  return(memoise(Vectorize(function(x){
    (4*deriv_theta_n(x)*f_a_n(x)*gamma_n(x))^(1/3)
  })))

}



#' Construct gamma_n nuisance estimator function
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param type Type of regression; currently only c("linear")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, type, omega_n, f_aIw_n) {
  
  # Construct weights
  s <- ss(dat_orig)
  dat_orig$weights <- wts(dat_orig, scale="none")
  
  # Construct pseudo-outcomes
  dat_orig %<>% mutate(
    po = ifelse(
      is.na(a),
      0,
      ( (omega_n(w1,w2,y_star,delta_star,a))*(weights/s)) / f_aIw_n(a,w1,w2) )^2
  )
  
  # Run regression
  if (type=="linear") {
    model <- lm(po~a, data=filter(dat_orig,!is.na(a)), weights=weights)
    coeff <- as.numeric(model$coefficients)
    
    return(memoise(Vectorize(function(x){
      coeff[1] + coeff[2]*x
    })))
  }
  
}



#' Construct estimator of marginal density f_A
#' 
#' @param dat Dataset returned by generate_data(); must be full data
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return Marginal density estimator function
construct_f_a_n <- function(dat, f_aIw_n) {
  
  return(memoise(Vectorize(function(a) {
    mean(f_aIw_n(a,dat$w1,dat$w2))
  })))

}



#' Construct estimator of conditional density f_{A|W}
#' 
#' @param dat Dataset returned by generate_data(); accepts either full data or
#'     truncated data
#' @param type One of c("parametric", "binning")
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of A is [0,1]
construct_f_aIw_n <- function(dat, type) {
  
  dat_trunc <- filter(dat, !is.na(a))
  n_trunc <- nrow(dat_trunc)
  weights <- wts(dat_trunc, scale="none")
  
  if (type=="parametric") {
    
    # Set up weighted likelihood
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        shape1 <- par[1] + par[2]*dat_trunc$w1[i]
        shape2 <- par[3] + par[4]*dat_trunc$w2[i]
        loglik <- dbeta(dat_trunc$a[i], shape1=shape1, shape2=shape2, log=TRUE)
        return(loglik*weights[i])
      }))
      
      return(-1*sum_loglik)
      
    }
    
    # Run optimizer
    opt <- optim(par=c(a1=0.5, a2=0.1, a3=0.5, a4=0.1), fn=wlik)
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: optim() did not converge")
    }
    
    f_aIw_n <- function(a, w1, w2){
      shape1 <- opt$par[1] + opt$par[2]*w1
      shape2 <- opt$par[3] + opt$par[4]*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    }
    
  }
  
  if (type=="binning") {
    
    # k is fixed for now; later choose via cross-validation
    k <- 10
    alphas <- seq(0, 1, length.out=k+1)
    
    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # par[1] through par[k-1] are the hazard components for the bins 1 to k-1
    # par[k] and par[k+1] are the coefficients for W1 and W2
    dens <- Vectorize(function(a, w1, w2, par) {
      bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
      hz <- sapply(c(1:(ifelse(bin==k,k-1,bin))), function(j) {
        expit(par[j] + par[k]*w1 + par[k+1]*w2)
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      dens <- k*p1*p2
      return(dens)
    }, vectorize.args=c("a","w1","w2"))
    
    # Set up weighted likelihood
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        lik <- dens(a=dat_trunc$a[i], w1=dat_trunc$w1[i], w2=dat_trunc$w2[i],
                    par)
        return(weights[i]*log(lik))
      }))
      
      return(-1*sum_loglik)
      
    }
    
    # Run optimizer
    opt <- optim(par=rep(0,k+1), fn=wlik, method="CG")
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: optim() did not converge")
    }
    
    f_aIw_n <- function(a, w1, w2){
      
      bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
      par <- opt$par
      hz <- sapply(c(1:(k-1)), function(j) {
        expit(par[j] + par[k]*w1 + par[k+1]*w2)
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      
      return(k*p1*p2)
      
    }
    
  }
  
  return(memoise(Vectorize(f_aIw_n)))
  
}



#' Construct density ratio estimator g_n
#' 
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
construct_g_n <- function(f_aIw_n, f_a_n) {
  
  return(memoise(Vectorize(function(a,w1,w2) {
    f_aIw_n(a,w1,w2) / f_a_n(a)
  })))
  
}



#' Construct estimator of nuisance influence function omega_n
#' 
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param Sc_n Conditional censoring survival function estimator returned by
#'     construct_S_n
#' @param m Number of partitions used for the integral approximation
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(S_n, Sc_n, m=400) {
  
  return(memoise(Vectorize(function(w1,w2,y_star,delta_star,a) {
    
    i <- c(1:m)
    k <- min(y_star,C$t_e)
    # integral <- sum(
    #   ( S_n(t=(i*k)/m, w1,w2,a) - S_n(t=((i-1)*k)/m, w1,w2,a) ) *
    #   ( (S_n(t=(i*k)/m, w1,w2,a))^2 * Sc_n(t=(i*k)/m, w1,w2,a) )^-1
    # )
    # S_n(t=C$t_e, w1,w2,a) * (
    #   ( delta_star * as.integer(y_star<=C$t_e) ) /
    #   ( S_n(t=y_star, w1,w2,a) * Sc_n(t=y_star, w1,w2,a) ) +
    #   integral
    # )
    integral <- sum(
      ( S_n(t=round((i*k)/m), w1,w2,a) - S_n(t=round(((i-1)*k)/m), w1,w2,a) ) *
        ( (S_n(t=round((i*k)/m), w1,w2,a))^2 * Sc_n(t=round((i*k)/m), w1,w2,a) )^-1
    )
    S_n(t=C$t_e, w1,w2,a) * (
      ( delta_star * as.integer(y_star<=C$t_e) ) /
        ( S_n(t=round(y_star), w1,w2,a) * Sc_n(t=round(y_star), w1,w2,a) ) +
        integral
    )
    
  })))
  
}



#' Construct nuisance estimator eta_n
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return Estimator function of nuisance eta_0
construct_eta_n <- function(dat_orig, mu_n) {
  
  s <- ss(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- wts(dat, scale="none")
  
  return(memoise(Vectorize(function(x,w1,w2) {
    sum(weights * as.integer(dat$a<=x) * mu_n(dat$a,w1,w2)) / (n_orig*s)
  })))
  
}



#' Construct nuisance estimator theta_naive_n
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return The "naive" estimator of theta_0
#' @notes This is a "simple" estimator of theta_n rather than the Grenander type
construct_theta_naive_n <- function(dat_orig, mu_n) {
  
  return(memoise(Vectorize(function(a) {
    mean(mu_n(a,dat_orig$w1,dat_orig$w2))
  })))
  
}



#' lambda estimator
#' 
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return Value of lambda
lambda <- function(k, G, dat_orig) {
  
  if (attr(dat_orig, "sampling")=="iid") {
    
    return( mean((G(dat_orig$a))^k) )
    
  } else if (attr(dat_orig, "sampling")=="two-phase") {
    
    n_orig <- nrow(dat_orig)
    s_0 <- ss(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights_0 <- wts(dat, scale="none")
    lambda <- (1/n_orig) * sum(
      (weights_0/s_0) * (G(dat$a))^k
    )
    return(lambda)
    
  }
  
}



#' Construct Gamma_n primitive one-step estimator
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param S_n A conditional survival function returned by construct_S_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @return Gamma_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Gamma_n <- function(dat_orig, omega_n, S_n, g_n) {
  
  dat <- dat_orig
  
  s <- ss(dat)
  n_orig <- nrow(dat)
  dat %<>% filter(!is.na(a))
  n <- nrow(dat)
  weights <- wts(dat, scale="none")
  
  i_long <- rep(c(1:n), each=n)
  j_long <- rep(c(1:n), times=n)
  a_i_long <- dat$a[i_long]
  w1_i_long <- dat$w1[i_long]
  w1_j_long <- dat$w1[j_long]
  w2_i_long <- dat$w2[i_long]
  w2_j_long <- dat$w2[j_long]
  delta_star_i_long <- dat$delta_star[i_long]
  delta_star_j_long <- dat$delta_star[j_long]
  
  subpiece_1a <- ( 1 + (
    omega_n(dat$w1,dat$w2,dat$y_star,dat$delta_star,dat$a) /
    g_n(dat$a,dat$w1,dat$w2)
  ) ) * ( weights / s )
  subpiece_2a <- S_n(C$t_e, w1_j_long, w2_j_long, a_i_long) / (
    s^2 * Pi(delta_star_i_long, w1_i_long, w2_i_long) *
      Pi(delta_star_j_long, w1_j_long, w2_j_long)
  )
  
  return(
    memoise(Vectorize(function(x) {
      
      subpiece_1b <- as.integer(dat$a<=x)
      piece_1 <- sum(subpiece_1a*subpiece_1b) * (1/n_orig)
      
      subpiece_2b <- as.integer(a_i_long<=x)
      piece_2 <- sum(subpiece_2a*subpiece_2b) * (1/(n_orig^2))
      
      return(piece_1-piece_2)
      
    }))
  )
  
}



#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param type One of c("ecdf", "inverse").
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
#'   - This accesses wts() globally
construct_Phi_n <- function (dat_orig, type="ecdf") {
  
  dat <- dat_orig
  s <- ss(dat)
  dat <- cbind(dat, wts=wts(dat, scale="none"))
  dat %<>% arrange(a)
  n_orig <- nrow(dat)
  dat %<>% filter(!is.na(a))
  vals_x <- unique(dat$a)
  
  vals_y <- c()
  for (j in 1:length(vals_x)) {
    indices <- which(dat$a==vals_x[j])
    wts_j <- dat$wts[indices]
    new_y_val <- sum(wts_j) / (n_orig*s)
    vals_y <- c(vals_y, new_y_val)
  }
  vals_y <- cumsum(vals_y)
  
  if (type=="ecdf") {
    rval <- approxfun(vals_x, vals_y, method="constant", yleft=0,
                      yright=1, f=0, ties="ordered")
  } else if (type=="inverse") {
    rval <- approxfun(vals_y, vals_x, method="constant", yleft=min(vals_x),
                      yright=max(vals_x), f=1, ties="ordered")
  }
  return(rval)
  
}



#' Hypothesis test based on logistic regression
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_wald <- function(dat, alt_type="incr", params) {
  
  model <- glm(y~w1+w2+a, data=dat, family="binomial")
  one_sided_p <- pnorm(summary(model)$coefficients["a",3], lower.tail=F)
  reject <- as.integer(one_sided_p<0.05)
  
  return (reject)
  
}



#' Westling 2020 test of the causal null
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_causalnull <- function(dat, alt_type="incr", params) {
  
  if (attr(dat, "sampling")=="iid") {
    
    Y <- dat$y
    A <- dat$a
    W <- subset(dat, select=c(w1,w2))
    
    test_obj <- causalNullTest(
      Y = Y,
      A = A,
      W = W,
      p = 2
      # control = list()
    )
    
  } else if (attr(dat, "sampling")=="two-phase") {
    # Multiple imputation step for two-phase sampled data
    # !!!!! TO DO !!!!!
    stop("Not yet implemented for two-phase sampling")
  }
  
  # !!!!! Right now this is a two-sided test, being compared to a one-sided test
  two_sided_p <- test_obj$test$p.val
  reject <- as.integer(two_sided_p<0.05)
  
  return (reject)
  
}
