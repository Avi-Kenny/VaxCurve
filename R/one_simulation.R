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
    # L <- list(n=1000, alpha_3=0.7, distr_A="Beta(0.9,1.1+0.4*w2)",
    #           reg_true="Logistic", sampling="two-phase",
    #           estimator=list(
    #             # est="Generalized Grenander",params=list(ci_type="logit")))
    #             est="G-comp",params=list(boot_reps=5)))
    # C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3)
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$reg_true, L$sampling)
    
    # Obtain estimates
    ests <- est_curve(
      dat = dat,
      estimator = L$estimator$est,
      params = L$estimator$params,
      points = C$points
    )
    
    # Return results
    theta_true <- attr(dat, "theta_true")
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=1)
      res_list[paste0("theta_",m)] <- theta_true[i]
      res_list[paste0("est_",m)] <- ests[[i]]$est
      res_list[paste0("ci_lo_",m)] <- ests[[i]]$ci_lo
      res_list[paste0("ci_hi_",m)] <- ests[[i]]$ci_hi
    }
    return(res_list)

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
    
    # # !!!!! For testing !!!!!
    # set.seed(1)
    # L <- list(n=1000, alpha_3=0.7, distr_A="Unif(0,1)", reg_true="Logistic",
    #   sampling="iid", test=list(type="test_2",params=list(boot_reps=2)))
    # C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3)
    # dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$reg_true, L$sampling)
    # reject <- test_2(dat, "incr", L$test$params)
    
    # Generate dataset
    dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$reg_true, L$sampling)
    
    # Perform hypothesis test
    reject <- use_method(L$test$type, list(
      dat = dat,
      alt_type = "incr",
      params = L$test$params
    ))
    
    # Return results
    return (list("reject"=reject))
    
  }
  
}
