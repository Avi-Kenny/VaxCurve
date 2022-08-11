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
      return_extra = "deriv_theta_n" # "Theta_os_n"
    )

    # Return results
    theta_true <- attr(dat_orig, "theta_true")
    Gamma_true <- attr(dat_orig, "Gamma_true")
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=1)
      res_list[paste0("theta_",m)] <- theta_true[i]
      res_list[paste0("est_",m)] <- ests$est[i]
      res_list[paste0("ci_lo_",m)] <- ests$ci_lo[i]
      res_list[paste0("ci_hi_",m)] <- ests$ci_hi[i]
      if (F) {
        res_list[paste0("Gamma_",m)] <- Gamma_true[i]
        res_list[paste0("estG_",m)] <- ests$ests_Gamma[i]
      } # DEBUG: return Gamma estimates
    }
    
    if (F) {
      res_list$g_n_star <- ests$g_n_star
      res_list$eta_ss_n3 <- ests$eta_ss_n3
      res_list$eta_ss_n5 <- ests$eta_ss_n5
      res_list$alpha_star_n3 <- ests$alpha_star_n3
      res_list$alpha_star_n5 <- ests$alpha_star_n5
      res_list$gcomp_n <- ests$gcomp_n
    } # DEBUG: return extras
    
    # # Return extra results
    # res_list[[".complex"]] <- list(Theta_os_n = ests$Theta_os_n)

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
    
    res <- list(
      "reject" = test_results$reject,
      "p_val" = test_results$p_val,
      "beta_n" = test_results$beta_n,
      "sd_n" = test_results$sd_n,
      "var_n" = test_results$var_n
    )
    
    # Debugging
    if (F) {
      res$if1_mean = test_results$if1_mean
      res$if2_mean = test_results$if2_mean
      res$r_1n = test_results$r_1n
      res$r_2n = test_results$r_2n
      res$Psi_1_var_est = test_results$Psi_1_var_est
      res$sum12_est = test_results$sum12_est
      res$sum12_var_est = test_results$sum12_var_est
      res$Psi_1_est = test_results$Psi_1_est
      res$Psi_2_est = test_results$Psi_2_est
      res$Psi_G_est = test_results$Psi_G_est
      res$Psi_1_var_est = test_results$Psi_1_var_est
      res$Psi_2_var_est = test_results$Psi_2_var_est
      res$Psi_12_covar = test_results$Psi_12_covar
      res$p_val_Psi_1 = test_results$p_val_Psi_1
      res$p_val_Psi_2 = test_results$p_val_Psi_2
      res$p_val_Psi_G = test_results$p_val_Psi_G
      res$p_val_sum12 = test_results$p_val_sum12
      res$p_val_sum12b = test_results$p_val_sum12b
      res$p_val_alt = test_results$p_val_alt
      res$reject_Psi_1 = test_results$reject_Psi_1
      res$reject_Psi_2 = test_results$reject_Psi_2
      res$reject_Psi_G = test_results$reject_Psi_G
      res$reject_sum12 = test_results$reject_sum12
      res$reject_sum12b = test_results$reject_sum12b
      res$reject_alt = test_results$reject_alt
    }
    
    return(res)
    
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    srvSL <- construct_S_n(dat, vlist$S_n, type=L$estimator$params$S_n_type)
    S_n <- srvSL$srv
    Sc_n <- srvSL$cens
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



###########################################.
##### Estimation (Cox model variance) #####
###########################################.

if (cfg$which_sim=="Cox") {
  
  #' Run a single simulation (Cox model variance)
  #'
  #' @return A list
  
  one_simulation <- function() {
    
    # Debugging
    if (F) {
      dat_orig <- generate_data(n=1000, -2, "Unif(0,1)", "none", surv_true="Cox PH", list(lmbd=1e-3,v=1.5,lmbd2=5e-5,v2=1.5), "two-phase (50%)", "decr", "estimated")
    }
    
    # Generate dataset
    dat_orig <- generate_data(
      L$n, L$alpha_3, L$distr_A, L$edge, surv_true="Cox PH",
      L$sc_params, L$sampling, L$dir, L$wts_type
    )
    
    # Round data values and construct dat
    dat_orig$a <- round(dat_orig$a,2)
    dat_orig$weights <- round(dat_orig$weights,3)
    dat_orig$y_star <- round(dat_orig$y_star,0)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    
    # Calculate variance estimates
    a <- 0.5
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_e, points=a,
                       se_beta=T, se_bshz=T, se_surv=T, se_marg=T)
    
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_e, points=0.5, se_marg=T, verbose=T)
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_e, points=c(0.2,0.5,0.8), se_marg=T, verbose=T)
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_e, points=seq(0.1,0.9,0.1), se_marg=T, verbose=T)
    
    # Calculate certain true values
    if (F) {
      z_0 <- c(0.3,1,0.5) # Needs to be consistent with the value in cox_var()
      H_0_true <- function(t) {
        L$sc_params$lmbd*exp(-1.7) * t^L$sc_params$v
      }
      true_lp <- sum(c(C$alpha_1,C$alpha_2,L$alpha_3)*z_0)
      true_surv <- exp(-1*exp(true_lp)*H_0_true(C$t_e))
      true_marg <- 1-attr(dat_orig, "theta_true")[26] # Corresponds to A=0.5
    } # DEBUG: intermediate objects
    
    # Construct simulation results object
    # This needs to line up with res_cox based on the se_* flags
    sim_res <- list(
      true_w1 = C$alpha_1,
      est_w1 = res_cox$beta_n[1],
      se_w1 = sqrt(res_cox$var_est_beta[1]),
      true_w2 = C$alpha_2,
      est_w2 = res_cox$beta_n[2],
      se_w2 = sqrt(res_cox$var_est_beta[2]),
      true_a = L$alpha_3,
      est_a = res_cox$beta_n[3],
      se_a = sqrt(res_cox$var_est_beta[3]),
      true_bshz = H_0_true(C$t_e),
      est_bshz = res_cox$est_bshz,
      se_est_bshz = sqrt(res_cox$var_est_bshz),
      true_surv = true_surv,
      est_surv = res_cox$est_surv,
      se_est_surv = sqrt(res_cox$var_est_surv),
      true_marg = true_marg,
      est_marg = res_cox$est_marg,
      se_est_marg = sqrt(res_cox$var_est_marg)
    )
    
    # Debugging
    if (F) {
      
      # # Fix a covariate vector and time of interest
      # z_0 <- c(0.3,1,0.5) # c(W1,W2,A)
      # 
      # # Calculate the cumulative hazard via predict()
      # newdata <- data.frame(y_star=C$t_e, delta_star=1, w1=z_0[1],
      #                       w2=z_0[2], a=z_0[3])
      # pred <- predict(res_cox$model, newdata=newdata, type="expected", se.fit=T)
      
      
      # Add additional results
      # sim_res$true_cmhz = exp(true_lp) * H_0_true(C$t_e) # Should roughly equal pred$se.fit
      # sim_res$est_cmhz = exp(sum(res_cox$beta_n*z_0))*est_bshz # Should equal pred$fit
      # sim_res$est_surv = exp(-1*exp(sum(res_cox$beta_n*z_0))*est_bshz)
      # sim_res$se_cmhz_MC = sqrt(res_cox$var_cmhz_est)
      # sim_res$se_surv_MC = sqrt(res_cox$var_surv_est)
      # sim_res$se_bshz_Cox = pred$se.fit
      
      # Debugging the Breslow estimator
      # sim_res$se_bshz_MC = sqrt(res_cox$var_bshz_est)
      # sim_res$se_bshz_MC2 = sqrt(res_cox$var_bshz_est2) # New derivation
      
      # sim_res$se_w1_Cox = as.numeric(sqrt(diag(vcov(res_cox$model)))[1])
      # sim_res$se_w1_alt = as.numeric(sqrt(diag(res_cox$I_tilde_inv)/L$n)[1])
      
    }
    
    return(sim_res)
    
  }
  
}
