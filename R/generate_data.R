#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Dose-response "relationship strength" parameter
#' @param distr_A Distribution of A; possibly dependent on covariates. One of
#'     c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)", "Beta(0.8+0.9*w1,0.8+0.4*w2)")
#' @param surv_true True form of the survival function; one c("CoxPH","Complex")
#' @param sampling A list. One of c("iid","two-phase")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_A, surv_true, sampling) {
  
  # Sample baseline covariates
  w1 <- runif(n)
  w2 <- rbinom(n, size=1, prob=0.5)
  
  # Sample A
  if (distr_A=="Unif(0,1)") {
    a <- runif(n)
  } else if (distr_A=="Beta(0.9,1.1+0.4*w2)") {
    shape1 <- 0.9
    shape2 <- 1.1 + 0.4*w2
    a <- rbeta(n, shape1=shape1, shape2=shape2)
  } else if (distr_A=="Beta(0.8+0.9*w1,0.8+0.4*w2)") {
    shape1 <- 0.8 + 0.9*w1
    shape2 <- 0.8 + 0.4*w2
    a <- rbeta(n, shape1=shape1, shape2=shape2)
  } else {
    stop("distr_A incorrectly specified")
  }
  
  # Generate event times
  {
    # Generate survival times (Weibull)
    U <- runif(n)
    H_0_inv <- function(t) {
      ((1/C$lambda)*t)^(1/C$v)
    }
    if (surv_true=="CoxPH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
    } else if (surv_true=="Complex") {
      lin <- ( C$alpha_1*w1 + (C$alpha_2*w2 * alpha_3*a) ) * w1
    }
    t <- H_0_inv(-1*log(U)*exp(-1*lin))
    
    # Generate censoring times (Weibull)
    t_study_end <- 300
    U <- runif(n)
    H_0_inv2 <- function(t) {
      ((1/C$lambda2)*t)^(1/C$v2)
    }
    if (surv_true=="CoxPH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
    } else if (surv_true=="Complex") {
      lin <- ( C$alpha_1*w1 + (C$alpha_2*w2 * alpha_3*a) ) * w1
    }
    c <- H_0_inv2(-1*log(U)*exp(-1*lin))
    c <- pmin(c,t_study_end)
    
    # Generate survival variables
    y_star <- pmin(t,c)
    delta_star <- as.integer(y_star==t)
    
  }
  
  
  # IID sampling
  if (sampling=="iid") {
    dat <- data.frame(w1=w1, w2=w2, a=a, y_star=y_star, delta_star=delta_star)
  }
  
  # Two-phase sampling
  if (sampling=="two-phase") {
    Pi <- Pi(delta_star,w1,w2)
    delta <- rbinom(n, size=1, prob=Pi)
    dat <- data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA),
                      y_star=y_star, delta_star=delta_star)
  }

  # Set up function to calculate true regression values over C$points
  # These are Monte Carlo approximations
  # Values depend on surv_true and alpha_3
  {
    m <- 10^5
    w1 <- runif(m)
    w2 <- rbinom(m, size=1, prob=0.5)
    
    theta_true_f <- Vectorize(function(a) {
      if (surv_true=="CoxPH") {
        lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
      } else if (surv_true=="Complex") {
        lin <- ( C$alpha_1*w1 + (C$alpha_2*w2 * alpha_3*a) ) * w1
      }
      
      return(mean(
        1 - exp( -1 * C$lambda * (C$t_e^C$v) * exp(lin) )
      ))
    })
  }
  
  # Add attributes to dataframe
  attr(dat, "theta_true") <- theta_true_f(C$points)
  attr(dat, "sampling") <- sampling
  
  return(dat)
  
}
