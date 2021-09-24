######################.
##### Estimation #####
######################.

if (cfg$which_sim=="estimation") {
  
  #' Run a single simulation (estimation)
  #'
  #' @return A list with the estimates and CI limits of the causal dose-response
  #'     curve evaluated at the midpoint and the endpoint of the domain
  
  one_simulation <- function() {
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge,
                              L$surv_true, L$sampling)
    
    # Obtain estimates
    ests <- est_curve(
      dat_orig = dat_orig,
      estimator = L$estimator$est,
      params = L$estimator$params,
      points = C$points
    )
    
    # Return results
    theta_true <- attr(dat_orig, "theta_true")
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=1)
      res_list[paste0("theta_",m)] <- theta_true[i]
      res_list[paste0("est_",m)] <- ests[[i]]$est
      res_list[paste0("ci_lo_",m)] <- ests[[i]]$ci_lo
      res_list[paste0("ci_hi_",m)] <- ests[[i]]$ci_hi
    }
    
    # Return extra results
    res_list[["ex_S_n"]] <- ests$ex_S_n
    res_list[["ex_gamma_n"]] <- ests$ex_gamma_n
    res_list[["ex_deriv_theta_n"]] <- ests$ex_deriv_theta_n
    res_list[["ex_tau_n"]] <- ests$ex_tau_n
    res_list[[".complex"]] <- ests$timestamps # dat_orig
    
    return(res_list)
    
  }
  
}



########################################.
##### Estimation (edge value only) #####
########################################.

if (cfg$which_sim=="edge") {
  
  #' Run a single simulation (estimation; only the portions related to the edge)
  #'
  #' @return A list with edge estimates

  one_simulation <- function() {
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge,
                              L$surv_true, L$sampling)
    
    # Prep
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig, C$appx)
    
    # Construct component functions
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=L$estimator$params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=L$estimator$params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    pi_n <- construct_pi_n(dat_orig, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat_orig, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat_orig, pi_n, S_n, omega_n,
                                   theta_os_n_est)
    
    # Return results
    return(list(
      theta_true = attr(dat_orig, "theta_true")[1],
      theta_est = theta_os_n_est,
      ci_lo = theta_os_n_est - 1.96*sqrt(sigma2_os_n_est/n_orig),
      ci_hi = theta_os_n_est + 1.96*sqrt(sigma2_os_n_est/n_orig)
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge,
                              L$surv_true, L$sampling)
    
    # Perform hypothesis test
    test_results <- use_method(L$test$type, list(
      dat_orig = dat_orig,
      alt_type = "incr",
      params = L$test$params
    ))
    
    # Return results
    return (list(
      "reject" = test_results$reject,
      "p_val" = test_results$p_val
    ))
    
  }
  
}
