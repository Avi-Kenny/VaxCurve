#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Dose-response "relationship strength" parameter; should be
#'     zero or negative.
#' @param distr_A Distribution of A, possibly dependent on covariates; one of
#'     c("Unif(0,1)", "Unif(0.3,0.7)", "N(0.5,0.01)", "N(0.5,0.04)",
#'     "N(0.4+0.2w1+0.1w2,0.01)"). The Normals are truncated to lie in [0,1].
#' @param edge Propensity mechanism for point mass at edge. The distribution of
#'     A will be (1-pi(w))*distr_A. One of the following: "none" (no point mass
#'     at the edge), "expit 0.2" (expit propensity model, prob=0.2), "expit 0.5"
#'     (expit propensity model, prob=0.5), "complex 0.2" (non-logit propensity
#'     model, prob=0.2).
#' @param surv_true True form of the survival function; one c("Cox PH",
#'     "complex")
#' @param sc_params Weibull parameters for the survival and censoring
#'     distributions; a list of the form list(lmbd=1,v=1,lmbd2=1,v2=1)
#' @param sampling Two-phase sampling mechanism; one of c("iid",
#'     "two-phase (6%)", "two-phase (72%)", "two-phase (70% random)", two-phase
#'     (6% random))
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_A, edge, surv_true, sc_params,
                          sampling, dir) {
  
  # Sample baseline covariates
  w <- data.frame(
    w1 = sample(round(seq(0,1,0.1),1), size=n, replace=T),
    w2 = rbinom(n, size=1, prob=0.5)
  )
  
  # Sample A (before point mass)
  if (distr_A=="Unif(0,1)") {
    a <- runif(n)
  } else if (distr_A=="Unif(0.3,0.7)") {
    a <- runif(n, min=0.3, max=0.7)
  } else if (distr_A=="N(0.5,0.01)") {
    a <- rtruncnorm(n, a=0, b=1, mean=0.5, sd=0.1)
  } else if (distr_A=="N(0.5,0.04)") {
    a <- rtruncnorm(n, a=0, b=1, mean=0.5, sd=0.2)
  } else if (distr_A=="N(0.4+0.2w1+0.1w2,0.01)") {
    a <- rtruncnorm(n, a=0, b=1, mean=0.4+0.2*w$w1+0.1*w$w2, sd=0.1)
  } else {
    stop("distr_A incorrectly specified")
  }
  
  # Adjust A for point mass at the edge
  if (edge=="expit 0.2") {
    edge_probs <- expit(w$w1+w$w2-2.5)
  } else if (edge=="expit 0.5") {
    edge_probs <- expit(w$w1+w$w2-1.0)
  } else if (edge=="complex 0.2") {
    edge_probs <- 1.7*w$w2*pmax(0,1-4*abs(w$w1-0.5))
  } else if (edge=="none") {
    edge_probs <- 0
  }
  edge_val <- rbinom(n, size=1, prob=edge_probs)
  a <- (1-edge_val)*a
  
  # Generate event times
  {
    # Generate survival times (Weibull)
    U <- runif(n)
    H_0_inv <- function(t) {
      ((1/sc_params$lmbd)*t)^(1/sc_params$v)
    }
    if (surv_true=="Cox PH") {
      if (dir=="decr") {
        lin <- C$alpha_1*w$w1 + C$alpha_2*w$w2 + alpha_3*a - 1.7
      } else {
        lin <- C$alpha_1*w$w1 + C$alpha_2*w$w2 + alpha_3*(1-a) - 1.7
      }
    } else if (surv_true=="complex") {
      if (dir=="decr") {
        lin <- C$alpha_1*pmax(0,2-8*abs(w$w1-0.5)) + 2.5*alpha_3*w$w2*a +
          0.7*alpha_3*(1-w$w2)*a - 1.3
      } else {
        lin <- C$alpha_1*pmax(0,2-8*abs(w$w1-0.5)) + 2.5*alpha_3*w$w2*(1-a) +
          0.7*alpha_3*(1-w$w2)*(1-a) - 1.3
      }
    }
    t <- H_0_inv(-1*log(U)*exp(-1*lin))
    
    # Generate censoring times (Weibull)
    U <- runif(n)
    H_0_inv2 <- function(t) {
      ((1/sc_params$lmbd2)*t)^(1/sc_params$v2)
    }
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w$w1 + C$alpha_2*w$w2 - 1
    } else if (surv_true=="complex") {
      lin <- C$alpha_1*pmax(0,2-8*abs(w$w1-0.5)) - 0.35
    }
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))
    
    # Generate survival variables
    y_star <- pmin(t,c)
    delta_star <- as.integer(y_star==t)
    
  }
  
  # Conduct sampling
  if (sampling=="iid") {
    dat_orig <- list(id=c(1:n), w=w, a=a, delta=rep(1,n),
                     y_star=y_star, delta_star=delta_star)
  } else {
    delta <- rbinom(n, size=1, prob=Pi(sampling,delta_star,y_star,w))
    dat_orig <- list(id=c(1:n), w=w, a=ifelse(delta==1,a,NA), delta=delta,
                     y_star=y_star, delta_star=delta_star)
  }
  
  # Set up function to calculate true regression values over C$points
  # These are Monte Carlo approximations
  {
    m <- 10^6 # 10^5
    w1 <- sample(round(seq(0,1,0.1),1), size=m, replace=T)
    w2 <- rbinom(m, size=1, prob=0.5)
    
    lin <- function(w1,w2,a) {
      if (surv_true=="Cox PH") {
        if (dir=="decr") {
          C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1.7
        } else {
          C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*(1-a) - 1.7
        }
      } else if (surv_true=="complex") {
        if (dir=="decr") {
          C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 2.5*alpha_3*w2*a +
            0.7*alpha_3*(1-w2)*a - 1.3
        } else {
          C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 2.5*alpha_3*w2*(1-a) +
            0.7*alpha_3*(1-w2)*(1-a) - 1.3
        }
      }
    }
    
    S_0 <- function(t, w1, w2, a) {
      exp( -1 * sc_params$lmbd * (t^sc_params$v) * exp(lin(w1,w2,a)) )
    }
    
    theta_true_f <- Vectorize(function(a) {
      return(1 - mean(S_0(C$t_e,w1,w2,a)))
    })
    
    a <- runif(m)
    Theta_true_f <- Vectorize(function(x) {
      return(mean( as.integer(a<=x) * (1-S_0(C$t_e,w1,w2,a)) ))
    })
    attr(dat_orig, "Theta_true") <- Theta_true_f(C$points)
    
  }
  
  # Add attributes to dataframe
  attr(dat_orig, "theta_true") <- theta_true_f(C$points)
  attr(dat_orig, "sampling") <- sampling
  
  # Add (stabilized) inverse weights
  dat_orig$weights <- wts(dat_orig)
  
  return(dat_orig)
  
}
