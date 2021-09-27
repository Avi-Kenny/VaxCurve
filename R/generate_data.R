#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Dose-response "relationship strength" parameter
#' @param distr_A Distribution of A, possibly dependent on covariates; one of
#'     c("Unif(0,1)", "Beta(1.5+w1,1.5+w2)", "N(0.5,0.01)", "N(0.5,0.04)",
#'     "N(0.4+0.2w1+0.1w2,0.01)"). The Normals are truncated to lie in [0,1].
#' @param edge Propensity mechanism for point mass at edge. The distribution of
#'     A will be (1-pi(w))*distr_A. One of the following: "none" (no point mass
#'     at the edge), "expit" (expit propensity model), "complex" (non-logit
#'     propensity model). Both "expit" and "complex" put roughly 10% of the mass
#'     on the edge.s
#' @param surv_true True form of the survival function; one c("Cox PH",
#'     "complex")
#' @param sc_params Weibull parameters for the survival and censoring
#'     distributions; a list of the form list(lmbd=1,v=1,lmbd2=1,v2=1)
#' @param sampling Two-phase sampling mechanism; one of c("iid",
#'     "two-phase (6%)", "two-phase (72%)")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_A, edge, surv_true, sc_params,
                          sampling) {
  
  # Sample baseline covariates
  w1 <- sample(seq(0,1,0.1), size=n, replace=T) # w1 <- runif(n)
  w2 <- rbinom(n, size=1, prob=0.5)
  
  # Sample A (before point mass)
  if (distr_A=="Unif(0,1)") {
    a <- runif(n)
  } else if (distr_A=="Beta(1.5+w1,1.5+w2)") {
    a <- rbeta(n, shape1=1.5+w1, shape2=1.5+w2)
  } else if (distr_A=="N(0.5,0.01)") {
    a <- rnorm(n, mean=0.5, sd=0.1)
    a <- ifelse(a>1,1,ifelse(a<0,0,a))
  } else if (distr_A=="N(0.5,0.04)") {
    a <- rnorm(n, mean=0.5, sd=0.2)
    a <- ifelse(a>1,1,ifelse(a<0,0,a))
  } else if (distr_A=="N(0.4+0.2w1+0.1w2,0.01)") {
    a <- rnorm(n, mean=0.4+0.2*w1+0.1*w2, sd=0.1)
    a <- ifelse(a>1,1,ifelse(a<0,0,a))
  } else {
    stop("distr_A incorrectly specified")
  }
  
  # Round A values to speed computation
  a <- round(a, -log10(C$appx$a))
  
  # Adjust A for point mass at the edge
  if (edge=="expit") {
    edge_probs <- expit(w1+w2-3.3)
  } else if (edge=="expit2") {
    edge_probs <- expit(w1+w2-1)
  } else if (edge=="complex") {
    edge_probs <- 0.84*w2*pmax(0,1-4*abs(w1-0.5))
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
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1.7
    } else if (surv_true=="complex") {
      lin <- C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 1.2*alpha_3*w2*a - 1
    }
    t <- H_0_inv(-1*log(U)*exp(-1*lin))
    
    # Generate censoring times (Weibull)
    U <- runif(n)
    H_0_inv2 <- function(t) {
      ((1/sc_params$lmbd2)*t)^(1/sc_params$v2)
    }
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 - 1
    } else if (surv_true=="complex") {
      lin <- C$alpha_1*pmax(0,2-8*abs(w1-0.5)) - 0.35
    }
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))
    # c <- pmin(c,C$t_e)
    
    # Generate survival variables
    y_star <- pmin(t,c)
    delta_star <- as.integer(y_star==t)
    
  }
  
  # Conduct sampling
  if (sampling=="iid") {
    dat_orig <- data.frame(w1=w1, w2=w2, a=a, delta=1, y_star=y_star,
                      delta_star=delta_star)
  } else {
    delta <- rbinom(n, size=1, prob=Pi(sampling,delta_star,y_star,w1,w2))
    dat_orig <- data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA), delta=delta,
                           y_star=y_star, delta_star=delta_star)
  }
  
  # Set up function to calculate true regression values over C$points
  # These are Monte Carlo approximations
  {
    m <- 10^5
    w1 <- runif(m)
    w2 <- rbinom(m, size=1, prob=0.5)
    
    theta_true_f <- Vectorize(function(a) {
      
      lin <- function(w1,w2,a) {
        if (surv_true=="Cox PH") {
          C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1.7
        } else if (surv_true=="complex") {
          C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 1.2*alpha_3*w2*a - 1
        }
      }
      
      S_0 <- function(t, w1, w2, a) {
        exp( -1 * sc_params$lmbd * (t^sc_params$v) * exp(lin(w1,w2,a)) )
      }
      
      return(1 - mean(S_0(C$t_e, w1, w2, a)))
      
    })
  }
  
  # Add attributes to dataframe
  attr(dat_orig, "theta_true") <- theta_true_f(C$points)
  attr(dat_orig, "sampling") <- sampling
  
  return(dat_orig)
  
}
