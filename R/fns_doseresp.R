#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {1/(1+exp(-x))}



#' Integral of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
intexpit <- function(x) {log(1+exp(x))}



#' Hypothesis test based on logistic regression
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_wald <- function(dat, alt_type="incr", params) {
  
  model <- glm(y~w1+w2+a, data=dat, family="binomial")
  one_sided_p <- pnorm(summary(model)$coefficients["a",3], lower.tail=F)
  reject <- as.numeric(one_sided_p<0.05)
  
  return (reject)
  
}
