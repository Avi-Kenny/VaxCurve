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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Obtain estimates
    ests <- est_curve(
      dat_orig = dat_orig,
      estimator = L$estimator$est,
      params = L$estimator$params,
      points = C$points,
      dir = L$dir,
      # return_extra = "Theta_os_n",
      which = L$which
    )
    
    # Return results
    theta_true <- attr(dat_orig, "theta_true")
    # Theta_true <- attr(dat_orig, "Theta_true") # !!!!!
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=1)
      # res_list[paste0("Theta_",m)] <- Theta_true[i] # !!!!!
      res_list[paste0("theta_",m)] <- theta_true[i]
      res_list[paste0("est_",m)] <- ests$est[i]
      # res_list[paste0("esT_",m)] <- ests$Theta_os_n(C$points[i]) # !!!!!
      res_list[paste0("ci_lo_",m)] <- ests$ci_lo[i]
      res_list[paste0("ci_hi_",m)] <- ests$ci_hi[i]
    }
    
    # # Return extra results
    # res_list[[".complex"]] <- list(
    #   Theta_os_n = ests$Theta_os_n
    # )
    
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
    
    # !!!!! Update function calls (htab-->superfunc)
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat, C$appx)
    
    # Construct component functions
    S_n <- construct_S_n(dat, vlist$S_n, type=L$estimator$params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=L$estimator$params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                 type=params$omega_n_type)
    pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat, pi_n, S_n, omega_n, theta_os_n_est)
    
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    msg <- "Direction of monotonicity does not align with test type"
    if (L$dir=="incr" && L$test$alt_type=="decr") { stop(msg) }
    if (L$dir=="decr" && L$test$alt_type=="incr") { stop(msg) }
    
    # Perform hypothesis test
    test_results <- use_method(L$test$type, list(
      dat_orig = dat_orig,
      alt_type = L$test$alt_type,
      params = L$test$params,
      test_stat_only = L$test$test_stat_only
    ))
    
    # Return results
    return (list(
      "reject" = test_results$reject,
      "p_val" = test_results$p_val,
      "beta_n" = test_results$beta_n,
      "sd_n" = test_results$sd_n,
      "psi1psi2" = test_results$psi1psi2, # !!!!!
      "Psi_1_est" = test_results$Psi_1_est, # !!!!!
      "Psi_2_est" = test_results$Psi_2_est # !!!!!
      # "Gamma_n_5" = test_results$Gamma_n_5, # !!!!!
      # "Gamma_true_5" = attr(dat_orig, "Theta_true")[26], # !!!!!
      # "theta_true_5" = attr(dat_orig, "theta_true")[26], # !!!!!
      # "Gamma_var_n" = test_results$Gamma_var_n, # !!!!!
      # "lambda_2" = test_results$lambda_2, # !!!!!
      # "lambda_3" = test_results$lambda_3, # !!!!!
      # "g_n" = test_results$g_n, # !!!!!
      # "S_n" = test_results$S_n, # !!!!!
      # "Sc_n" = test_results$Sc_n, # !!!!!
      # "omega_n1" = test_results$omega_n1, # !!!!!
      # "omega_n0" = test_results$omega_n0, # !!!!!
      # "gcomp_n" = test_results$gcomp_n, # !!!!!
      # "eta_n" = test_results$eta_n, # !!!!!
      # "xi_n" = test_results$xi_n, # !!!!!
      # "infl_fn_1" = test_results$infl_fn_1, # !!!!!
      # "infl_fn_Gamma" = test_results$infl_fn_Gamma, # !!!!!
      # "infl_fn_2" = test_results$infl_fn_2, # !!!!!
      # "partial_est" = test_results$partial_est, # !!!!!
      # "partial_var" = test_results$partial_var, # !!!!!
      # "Phi_n_5" = test_results$Phi_n_5 # !!!!!
    ))
    
  }
  
}



##############################.
##### Hypothesis testing #####
##############################.

if (cfg$which_sim=="infl_fn_1 (temp)") {
  
  one_simulation <- function() {
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=L$test$params$ecdf_type)
    
    # Gamma_0 true function
    Gamma_0 <- Vectorize(function(a) {
      Theta_true <- attr(dat_orig,"Theta_true")
      grid <- seq(0,1,0.02)
      index <- which.min(abs(a-seq(0,1,0.02)))
      return(Theta_true[index])
    })
    
    # Psi_4
    if (F) {
      Psi_4_n <- (1/n_orig) * sum( weights * ( Phi_n(dat$a) * Gamma_0(dat$a) ))
      construct_infl_fn_4 <- function(dat, Gamma_0, Phi_n) {
        n_orig <- sum(dat$weights)
        weights_j <- dat$weights
        a_j <- dat$a
        fnc <- function(a_i) {
          piece_1 <- (1/n_orig) * sum(
            weights_j * (as.integer(a_i<=a_j)-2*Phi_n(a_j)) * Gamma_0(a_j)
          )
          piece_2 <- Phi_n(a_i) * Gamma_0(a_i)
          return(piece_1+piece_2)
        }
        return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
      }
      infl_fn_Psi_4 <- construct_infl_fn_4(dat, Gamma_0, Phi_n)
      Psi_4_var <- (1/n_orig^2) * sum((weights*infl_fn_Psi_4(dat$a))^2)
    }
    
    # Psi_5
    if (F) {
      Psi_5_n <- (1/n_orig) * sum(weights*Phi_n(dat$a))
      construct_infl_fn_5 <- function(dat, Phi_n) {
        n_orig <- sum(dat$weights)
        weights_j <- dat$weights
        a_j <- dat$a
        fnc <- function(a_i) {
          piece_1 <- (1/n_orig) * sum(
            weights_j * (as.integer(a_i<=a_j)-2*Phi_n(a_j))
          )
          piece_2 <- Phi_n(a_i)
          return(piece_1+piece_2)
        }
        return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
      }
      infl_fn_Psi_5 <- construct_infl_fn_5(dat, Phi_n)
      Psi_5_var <- (1/n_orig^2) * sum((weights*infl_fn_Psi_5(dat$a))^2)
    }
    
    # infl_fn_1
    if (T) {
      lambda_2 <- 1/3 # lambda(dat,2,Phi_n)
      lambda_3 <- 1/4 # lambda(dat,3,Phi_n)
      Phi_0 <- function(x) { x }
      Psi_1_n <- (1/n_orig) * sum( weights * (
        # (lambda_2*(Phi_n(dat$a))^2 - lambda_3*Phi_n(dat$a)) * Gamma_0(dat$a)
        (lambda_2*(Phi_0(dat$a))^2 - lambda_3*Phi_0(dat$a)) * Gamma_0(dat$a)
      ))
      # infl_fn_1 <- construct_infl_fn_1(dat, Gamma_0, Phi_n, lambda_2, lambda_3)
      infl_fn_1 <- construct_infl_fn_1(dat, Gamma_0, Phi_0, lambda_2, lambda_3)
      Psi_1_var <- (1/n_orig^2) * sum((weights*infl_fn_1(dat$a))^2)
    }
    
    # Return results
    return (list(
      # "partial_est" = partial_est,
      # "partial_sd" = sqrt(partial_var),
      # "rho_5" = rho_n(0.5),
      # "lambda_2" = lambda_2,
      # "lambda_3" = lambda_3,
      "Psi_1_n" = Psi_1_n,
      "Psi_1_sd" = sqrt(Psi_1_var)
      # "Psi_4_n" = Psi_4_n,
      # "Psi_4_sd" = sqrt(Psi_4_var),
      # "Psi_5_n" = Psi_5_n,
      # "Psi_5_sd" = sqrt(Psi_5_var)
    ))
    
  }
  
}
