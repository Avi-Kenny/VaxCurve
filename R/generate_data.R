#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Dose-response "relationship strength" parameter; should be
#'     zero or negative.
#' @param distr_S Distribution of S, possibly dependent on covariates; one of
#'     c("Unif(0,1)", "Unif(0.3,0.7)", "N(0.5,0.01)", "N(0.5,0.04)",
#'     "N(0.3+0.4x2,0.04)"). The Normals are truncated to lie in [0,1].
#' @param edge Propensity mechanism for point mass at edge. The distribution of
#'     S will be (1-pi(x))*distr_S. One of the following: "none" (no point mass
#'     at the edge), "expit 0.2" (expit propensity model, prob=0.2), "expit 0.5"
#'     (expit propensity model, prob=0.5), "Complex 0.2" (non-logit propensity
#'     model, prob=0.2).
#' @param surv_true True form of the survival function; one c("Cox PH",
#'     "Non PH", "Complex", "exp")
#' @param sc_params Weibull parameters for the survival and censoring
#'     distributions; a list of the form list(lmbd=1,v=1,lmbd2=1,v2=1)
#' @param sampling Two-phase sampling mechanism; one of c("iid",
#'     "two-phase (6%)", "two-phase (72%)", "two-phase (70% random)", two-phase
#'     (6% random))
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @param wts_type One of c("true", "estimated")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_S, edge, surv_true, sc_params,
                          sampling, dir, wts_type="true") {
  
  # Sample baseline covariates
  x <- data.frame(
    x1 = sample(round(seq(0,1,0.1),1), size=n, replace=T),
    x2 = rbinom(n, size=1, prob=0.5)
  )
  
  # Sample S (before point mass)
  if (distr_S=="Unif(0,1)") {
    s <- runif(n)
  } else if (distr_S=="Unif(0.3,0.7)") {
    s <- runif(n, min=0.3, max=0.7)
  } else if (distr_S=="N(0.5,0.01)") {
    s <- rtruncnorm(n, a=0, b=1, mean=0.5, sd=0.1)
  } else if (distr_S=="N(0.5,0.04)") {
    s <- rtruncnorm(n, a=0, b=1, mean=0.5, sd=0.2)
  } else if (distr_S=="N(0.3+0.4x2,0.04)") {
    s <- rtruncnorm(n, a=0, b=1, mean=0.3+0.4*x$x2, sd=0.2)
  } else if (distr_S=="Tri UP") {
    s <- sqrt(runif(n))
  } else if (distr_S=="Tri DN") {
    s <- 1 - sqrt(runif(n))
  } else {
    stop("distr_S incorrectly specified")
  }
  
  # Adjust S for point mass at the edge
  if (edge=="expit 0.2") {
    edge_probs <- expit(x$x1+x$x2-2.5)
  } else if (edge=="expit 0.4") {
    edge_probs <- expit(x$x1+x$x2-1.4)
  } else if (edge=="Complex 0.2") {
    edge_probs <- 1.7*x$x2*pmax(0,1-4*abs(x$x1-0.5))
  } else if (edge=="none") {
    edge_probs <- 0
  }
  edge_val <- rbinom(n, size=1, prob=edge_probs)
  s <- (1-edge_val)*s
  
  # Generate event times
  {
    # Survival times (Weibull)
    U <- runif(n)
    H_0_inv <- function(t) { ((1/sc_params$lmbd)*t)^(1/sc_params$v) }
    if (surv_true=="Cox PH") {
      if (dir=="decr") {
        lin <- C$alpha_1*x$x1 + C$alpha_2*x$x2 + alpha_3*s - 1.7
      } else {
        lin <- C$alpha_1*x$x1 + C$alpha_2*x$x2 + alpha_3*(1-s) - 1.7
      }
    } else if (surv_true=="Complex") {
      if (dir=="decr") {
        lin <- alpha_3*expit(10*(2*s-1)) - 0.5
      } else {
        lin <- alpha_3*expit(10*(2*(1-s)-1)) - 0.5
      }
    } else if (surv_true=="exp") {
      H_0_inv <- function(t) { ((1/sc_params$lmbd)*t) }
      lin <- 0
    }
    if (surv_true=="Non PH") {
      if (dir=="decr") {
        t <- -1 * (1.1*expit(10-50*s))^-1 * log(U) - 1
      } else {
        # !!!!!
      }
    } else {
      t <- H_0_inv(-1*log(U)*exp(-1*lin))
    }
    
    # Censoring times (Weibull)
    U <- runif(n)
    H_0_inv2 <- function(t) { ((1/sc_params$lmbd2)*t)^(1/sc_params$v2) }
    if (surv_true %in% c("Cox PH", "Complex", "Non PH")) {
      lin <- C$alpha_1*x$x1 + C$alpha_2*x$x2 - 1
    } else if (surv_true=="exp") {
      H_0_inv2 <- function(t) { ((1/sc_params$lmbd2)*t) }
      lin <- 0
    }
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))
    
    # Generate survival variables
    y <- pmin(t,c)
    delta <- as.integer(y==t)
    
  }
  
  # Conduct sampling
  if (sampling=="iid") {
    dat_orig <- list(id=c(1:n), x=x, s=s, z=rep(1,n), y=y, delta=delta)
  } else {
    z <- rbinom(n, size=1, prob=Pi(sampling,delta,y,x))
    dat_orig <- list(id=c(1:n), x=x, s=ifelse(z==1,s,NA), z=z, y=y, delta=delta)
  }
  
  # Set up function to calculate true marginalized risk values over C$points
  # These are Monte Carlo approximations
  {
    m <- 10^5 # 10^6
    x1 <- sample(round(seq(0,1,0.1),1), size=m, replace=T)
    x2 <- rbinom(m, size=1, prob=0.5)
    
    lin <- function(x1,x2,s) {
      if (surv_true=="Cox PH") {
        if (dir=="decr") {
          C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*s - 1.7
        } else {
          C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*(1-s) - 1.7
        }
      } else if (surv_true=="Complex") {
        if (dir=="decr") {
          alpha_3*expit(10*(2*s-1)) - 0.5
        } else {
          alpha_3*expit(10*(2*(1-s)-1)) - 0.5
        }
      }
    }
    
    if (surv_true %in% c("Cox PH", "Complex")) {
      Q_0 <- function(t, x1, x2, s) {
        exp( -1 * sc_params$lmbd * (t^sc_params$v) * exp(lin(x1,x2,s)) )
      }
    } else if (surv_true=="Non PH") {
      Q_0 <- function(t, x1, x2, s) {
        if (dir=="decr") {
          exp(-1*1.1*(t+1)*expit(10-50*s))
        } else {
          # !!!!!
        }
      }
    } else if (surv_true=="exp") {
      Q_0 <- function(t, x1, x2, s) { exp(-1*sc_params$lmbd*t) }
    }
    
    r_M0_f <- Vectorize(function(s) { 1 - mean(Q_0(C$t_0,x1,x2,s)) })
    
    # Note: Uncomment this to return true Theta_true
    #       Only holds if S and X are independent
    if (T) {
      
      if (distr_S=="Unif(0,1)") {
        s <- runif(m)
      } else if (distr_S=="Unif(0.3,0.7)") {
        s <- runif(m, min=0.3, max=0.7)
      } else if (distr_S=="N(0.5,0.01)") {
        s <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.1)
      } else if (distr_S=="N(0.5,0.04)") {
        s <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.2)
      } else if (distr_S=="N(0.3+0.4x2,0.04)") {
        s <- rtruncnorm(m, a=0, b=1, mean=0.3+0.4*x2, sd=0.2)
      } else if (distr_S=="Tri UP") {
        s <- sqrt(runif(m))
      } else if (distr_S=="Tri DN") {
        s <- 1 - sqrt(runif(m))
      }
      
      Gamma_true_f <- Vectorize(function(u) {
        mean( as.integer(s<=u) * (1-Q_0(C$t_0,x1,x2,s)) )
      })
      attr(dat_orig, "Gamma_true") <- Gamma_true_f(C$points)
      
    }
    
  }
  
  # Add attributes to dataframe
  attr(dat_orig, "r_M0") <- r_M0_f(C$points)
  attr(dat_orig, "sampling") <- sampling
  
  # Add (stabilized) inverse weights
  wts_str <- wts(dat_orig, type=wts_type, return_strata=T)
  dat_orig$weights <- wts_str$weights
  dat_orig$strata <- wts_str$strata
  
  if (F) {
    dat_orig$weights_true <- wts(dat_orig, type="true", return_strata=T)$weights
  } # DEBUG: Return true weights in addition to estimated weights
  
  return(dat_orig)
  
}
