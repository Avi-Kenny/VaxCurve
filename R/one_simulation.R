######################.
##### Estimation #####
######################.

if (cfg$which_sim=="estimation") {
  
  #' Run a single simulation (estimation)
  #'
  #' @return A list with the estimates and standard errors of the causal dose-
  #'     response curve evaluated at the midpoint and the endpoint of the domain

  one_simulation <- function() {
    
    # # !!!!!
    # L <- list(n=20,alpha_3=0.7,distr_A="Unif(0,1)",mono_form="sqrt",sampling="iid",estimator="G-comp (logistic)")
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$mono_form, L$sampling)
    
    # Obtain estimates
    ests <- est_curve(dat, estimator=L$estimator, params=999, points=c(0.5,1))
    
    # Return results ("mp"=midpoint, "ep"=endpoint)
    return (list(
      "theta_mp" = attr(dat, "theta_mp"),
      "theta_ep" = attr(dat, "theta_ep"),
      "est_mp" = ests[[1]]$est,
      "se_mp" = ests[[1]]$se,
      "est_ep" = ests[[2]]$est,
      "se_ep" = ests[[2]]$se
    ))
    
  }
  
}



##############################.
##### Hypothesis testing #####
##############################.

if (cfg$which_sim=="testing") {
  
  #' Run a single simulation (testing)
  #'
  #' @return !!!!!
  
  one_simulation <- function() {
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$mono_form, L$sampling)
    
    # Perform hypothesis test
    reject <- do.call(L$test$type, list(dat, "incr", L$test$params))
    
    # Return results
    return (list("reject"=reject))
    
  }
  
}
