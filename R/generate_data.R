#' Generate data
#' 
#' @param n Sample size
#' @param alpha_3 Height of the dose-response curve
#' @param mono_form Functional form of dose-response curve; should be a monotone
#'     increasing function with domain [0,1] such that f(0)=0 and f(1)=1;
#'     choices are c("identity", "square", "sqrt", "step_0.2", "step_0.8")
#' @param sampling A list. One of c("iid","two-phase")
#' @return A dataframe representing the study population
generate_data <- function(n, alpha_3, mono_form, sampling="iid") {
  
  # Fix parameters
  alpha_0 <- -1.5
  alpha_1 <- 0.3
  alpha_2 <- 0.7
  
  # Sample baseline covariates
  w1 <- rnorm(n)
  w2 <- rbinom(n, size=1, prob=0.5)
  shape_2 <- ifelse(w2==1, 1.5, 1.1)
  a <- rbeta(n, shape1=0.9, shape2=shape_2)
  
  # Set transformation function
  if (mono_form=="identity") {
    mono_f <- function(x) {x}
  } else if (mono_form=="square") {
    mono_f <- function(x) {x^2}
  } else if (mono_form=="sqrt") {
    mono_f <- function(x) {sqrt(x)}
  } else if (mono_form=="step_0.2") {
    mono_f <- function(x) {as.numeric(x>0.2)}
  } else if (mono_form=="step_0.8") {
    mono_f <- function(x) {as.numeric(x>0.8)}
  } else {
    stop("mono_form incorrectly specified")
  }
  
  # Sample outcome
  probs <- expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*mono_f(a))
  y <- rbinom(n, size=1, prob=probs)
  
  # IID sampling
  if (sampling=="iid") {
    return (data.frame(w1=w1, w2=w2, a=a, y=y))
  }
  
  # Two-phase sampling
  if (sampling=="two-phase") {
    pi <- expit(2*y-w2)
    delta <- rbinom(n, size=1, prob=pi)
    return (data.frame(w1=w1, w2=w2, a=ifelse(delta==1,a,NA), y=y))
  }
  
}
