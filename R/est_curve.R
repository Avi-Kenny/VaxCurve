#' Estimate the causal dose-response curve at a point
#' 
#' @param dat Data returned by generate_data()
#' @param estimator Which estimator to use; one of c("G-comp (logistic)",
#'     "G-comp (GAM)","Grenander")
#' @param params A list, containing the following:
#'   - `x` !!!!!
#'   - `x` !!!!!
#' @param points A vector representing the points to estimate
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)
#' @notes
#'   - Note

est_curve <- function(dat, estimator, params, points) {
  
  # !!!!! Params currently unused
  
  if (estimator=="G-comp (logistic)") {
    
    res <- list()
    
    model <- glm(y~w1+w2+a, data=dat, family="binomial")
    coefs <- as.numeric(summary(model)$coefficients[,1])
    vcov <- vcov(model)
    
    for (p in 1:length(points)) {
      
      a <- points[p]
      
      gcomp_i <- apply(
        X = dat,
        MARGIN = 1,
        FUN = function(r) {
          M_i <- c(1,r[["w1"]],r[["w2"]],a)
          lin_i <- as.numeric(M_i %*% coefs)
          vlin_i <- as.numeric(M_i %*% vcov %*% M_i)
          est_i <- expit(lin_i)
          var_i <- (deriv_expit(lin_i))^2 * vlin_i
          return(c(est_i,var_i))
        }
      )
      
      est <- mean(gcomp_i[1,])
      se <- sqrt(mean(gcomp_i[2,]))
      
      res[[p]] <- list(point=a, est=est, se=se)
      
    }
    
    return (res)
    
  } else if (estimator=="logistic GAM") {
    # !!!!!
  } else if (estimator=="Grenander") {
    # !!!!!
  } else {
    stop("estimator incorrectly specified")
  }
  
}
