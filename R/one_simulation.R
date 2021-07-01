######################.
##### Estimation #####
######################.

if (cfg$which_sim=="estimation") {
  
  #' Run a single simulation (estimation)
  #'
  #' @return A list with the estimates and standard errors of the causal dose-
  #'     response curve evaluated at the midpoint and the endpoint of the domain

  one_simulation <- function() {
    
    # # !!!!! For testing !!!!!
    # L <- list(n=500, alpha_3=0.7, distr_A="Unif(0,1)", mono_form="identity",
    #   sampling="iid", estimator=list(
    #     # est="G-comp (logistic)",params=list(boot_reps=5)))
    #     est="Generalized Grenander",params=list(ci_type="logit")))
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$mono_form, L$sampling)
    
    # Obtain estimates
    ests <- est_curve(
      dat = dat,
      estimator = L$estimator$est,
      params = L$estimator$params,
      points = c(0.5,1)
    )
    
    # Return results ("mp"=midpoint, "ep"=endpoint)
    return (list(
      "theta_mp" = attr(dat, "theta_mp"),
      "theta_ep" = attr(dat, "theta_ep"),
      "est_mp" = ests[[1]]$est,
      "ci_lo_mp" = ests[[1]]$ci_lo,
      "ci_hi_mp" = ests[[1]]$ci_hi,
      "est_ep" = ests[[2]]$est,
      "ci_lo_ep" = ests[[2]]$ci_lo,
      "ci_hi_ep" = ests[[2]]$ci_hi
    ))
    
  }
  
}



##############################.
##### Hypothesis testing #####
##############################.

if (cfg$which_sim=="testing") {
  
  #' Run a single simulation (testing)
  #'
  #' @return A list that gives the result of the hypothesis test

  one_simulation <- function() {
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$mono_form, L$sampling)
    
    # Perform hypothesis test
    reject <- use_method(L$test$type, list(dat, "incr", L$test$params))
    
    # Return results
    return (list("reject"=reject))
    
  }
  
}
