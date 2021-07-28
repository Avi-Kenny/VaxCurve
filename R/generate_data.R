#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Height of the dose-response curve
#' @param distr_A Distribution of A; possibly dependent on covariates. One of
#'     c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)", "Beta(0.8+0.9*w1,0.8+0.4*w2)")
#' @param reg_true True functional form of the regression; one of the following:
#'     - "Logistic": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*A)
#'     - "GAM": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*sqrt(A))
#'     - "Complex": E[Y|W,A]=expit(a0+a1*sin(2*pi*W1)+a2*W2+a3*sqrt(A)+a4*W1*W2)
#' @param sampling A list. One of c("iid","two-phase")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, distr_A, reg_true, sampling) {
  
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
  
  # Compute true regression function (i.e. P(Y=1|W,A))
  if (reg_true=="Logistic") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a)
  } else if (reg_true=="GAM") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*sqrt(a))
  } else if (reg_true=="Complex") {
    probs <- expit(C$alpha_0 + C$alpha_1*sin(2*pi*w1) + C$alpha_2*w2 +
                     alpha_3*sqrt(a) + C$alpha_4*w1*w2)
  }
  
  # Sample outcome
  y <- rbinom(n, size=1, prob=probs)
  
  # IID sampling
  if (sampling=="iid") {
    dat <- data.frame(w1=w1, w2=w2, a=a, y=y)
  }
  
  # Two-phase sampling
  if (sampling=="two-phase") {
    pi <- Pi(y,w1,w2)
    delta <- rbinom(n, size=1, prob=pi)
    dat <- data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA), y=y)
  }

  # Set up function to calculate true regression values over C$points
  # These are Monte Carlo approximations
  # Values depend on reg_true and alpha_3
  {
    m <- 10^5
    w1 <- runif(m)
    w2 <- rbinom(m, size=1, prob=0.5)
    
    theta_true_f <- Vectorize(function(a) {
      if (reg_true=="Logistic") {
        return(mean(
          expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a)
        ))
      } else if (reg_true=="GAM") {
        return(mean(
          expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*sqrt(a))
        ))
      } else if (reg_true=="Complex") {
        return(mean(
          expit(C$alpha_0 + C$alpha_1*sin(2*pi*w1) + C$alpha_2*w2 +
                alpha_3*sqrt(a) + C$alpha_4*w1*w2)
        ))
      }
    })
  }
  
  # Add attributes to dataframe
  attr(dat, "theta_true") <- theta_true_f(C$points)
  attr(dat, "sampling") <- sampling
  
  return(dat)
  
}
