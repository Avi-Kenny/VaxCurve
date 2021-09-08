#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Dose-response "relationship strength" parameter
#' @param distr_A Distribution of A, possibly dependent on covariates; one of
#'     c("Unif(0,1)", "Beta(1.5+w1,1.5+w2)")
#' @param edge Propensity mechanism for point mass at edge. The distribution of
#'     A will be (1-pi(w))*distr_A. One of the following: "none" (no point mass
#'     at the edge), "expit" (expit propensity model), "complex" (non-logit
#'     propensity model). Both "expit" and "complex" put roughly 10% of the mass
#'     on the edge.s
#' @param surv_true True form of the survival function; one c("Cox PH",
#'     "complex")
#' @param sampling Two-phase sampling mechanism; one of c("iid","two-phase")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_A, edge, surv_true, sampling) {
  
  # Sample baseline covariates
  w1 <- runif(n)
  w2 <- rbinom(n, size=1, prob=0.5)
  
  # Sample A (before point mass)
  if (distr_A=="Unif(0,1)") {
    a <- runif(n)
  } else if (distr_A=="Beta(1.5+w1,1.5+w2)") {
    a <- rbeta(n, shape1=1.5+w1, shape2=1.5+w2)
  } else {
    stop("distr_A incorrectly specified")
  }
  
  # Adjust A for point mass at the edge
  if (edge=="expit") {
    edge_probs <- expit(w1+w2-3.3)
  } else if (edge=="complex") {
    edge_probs <- w2*as.integer(abs(w1-0.5)<0.11)
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
      ((1/C$lambda)*t)^(1/C$v)
    }
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
    } else if (surv_true=="complex") {
      lin <- as.numeric(abs(w1-0.5)<0.2) + as.numeric(abs(a-0.5)<0.2)
    }
    t <- H_0_inv(-1*log(U)*exp(-1*lin))
    
    # Generate censoring times (Weibull)
    # t_study_end <- 300
    U <- runif(n)
    H_0_inv2 <- function(t) {
      ((1/C$lambda2)*t)^(1/C$v2)
    }
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
    } else if (surv_true=="complex") {
      lin <- as.numeric(abs(w1-0.5)<0.2) + alpha_3*w2*a
    }
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))
    # c <- pmin(c,t_study_end)
    
    # Generate survival variables
    y_star <- pmin(t,c)
    delta_star <- as.integer(y_star==t)
    
  }
  
  # IID sampling
  if (sampling=="iid") {
    dat <- data.frame(w1=w1, w2=w2, a=a, delta=1, y_star=y_star,
                      delta_star=delta_star)
  }
  
  # Two-phase sampling
  if (sampling=="two-phase") {
    delta <- rbinom(n, size=1, prob=Pi("two-phase",delta_star,w1,w2))
    dat <- data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA), delta=delta,
                      y_star=y_star, delta_star=delta_star)
  }
  
  if (!(sampling %in% c("iid", "two-phase"))) {
    stop("`sampling` incorrectly specified")
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
          C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
        } else if (surv_true=="complex") {
          as.numeric(abs(w1-0.5)<0.2) + alpha_3*w2*a
        }
      }
      
      S_0 <- function(t, w1, w2, a) {
        exp( -1 * C$lambda * (t^C$v) * exp(lin(w1,w2,a)) )
      }
      
      return(1 - mean(S_0(C$t_e, w1, w2, a)))
      
    })
  }
  
  # Add attributes to dataframe
  attr(dat, "theta_true") <- theta_true_f(C$points)
  attr(dat, "sampling") <- sampling
  
  return(dat)
  
}
