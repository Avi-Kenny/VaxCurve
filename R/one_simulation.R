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
      dir = L$dir
      # return_extra = "Theta_os_n"
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
      "var_n" = test_results$var_n
      # "if1_mean" = test_results$if1_mean, # !!!!!
      # "if2_mean" = test_results$if2_mean, # !!!!!
      # "r_1n" = test_results$r_1n, # !!!!!
      # "r_2n" = test_results$r_2n, # !!!!!
      # "Psi_1_var_est" = test_results$Psi_1_var_est, # !!!!!
      # "sum12_est" = test_results$sum12_est, # !!!!!
      # "sum12_var_est" = test_results$sum12_var_est, # !!!!!
      # "Psi_1_est" = test_results$Psi_1_est, # !!!!!
      # "Psi_2_est" = test_results$Psi_2_est, # !!!!!
      # "Psi_G_est" = test_results$Psi_G_est, # !!!!!
      # "Psi_1_var_est" = test_results$Psi_1_var_est, # !!!!!
      # "Psi_2_var_est" = test_results$Psi_2_var_est, # !!!!!
      # "Psi_12_covar" = test_results$Psi_12_covar, # !!!!!
      # "p_val_Psi_1" = test_results$p_val_Psi_1, # !!!!!
      # "p_val_Psi_2" = test_results$p_val_Psi_2, # !!!!!
      # "p_val_Psi_G" = test_results$p_val_Psi_G, # !!!!!
      # "p_val_sum12" = test_results$p_val_sum12, # !!!!!
      # "p_val_sum12b" = test_results$p_val_sum12b, # !!!!!
      # "p_val_alt" = test_results$p_val_alt, # !!!!!
      # "reject_Psi_1" = test_results$reject_Psi_1, # !!!!!
      # "reject_Psi_2" = test_results$reject_Psi_2, # !!!!!
      # "reject_Psi_G" = test_results$reject_Psi_G # !!!!!
      # "reject_sum12" = test_results$reject_sum12, # !!!!!
      # "reject_sum12b" = test_results$reject_sum12b, # !!!!!
      # "reject_alt" = test_results$reject_alt # !!!!!
    ))
    
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
    
    # Fix a covariate vector and time of interest
    z_0 <- c(0.3,1,0.5) # c(W1,W2,A)
    
    # Generate dataset
    # dat_orig <- generate_data(n=100, -2, "Unif(0,1)", "none", surv_true="Cox PH", list(lmbd=1e-3,v=1.5,lmbd2=5e-5,v2=1.5), "iid", "decr")
    dat_orig <- generate_data(
      L$n, L$alpha_3, L$distr_A, L$edge, surv_true="Cox PH",
      L$sc_params, L$sampling, L$dir
    )
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    
    # Fit a Cox model
    df <- data.frame(y_star=dat$y_star, delta_star=dat$delta_star,
                     w1=dat$w$w1, w2=dat$w$w2, a=dat$a)
    wts2 <- dat$weights * (length(dat$weights)/sum(dat$weights)) # !!!!! Is this necessary?
    model <- coxph(Surv(y_star,delta_star)~w1+w2+a, data=df, weights=wts2)
    cfs <- model$coefficients
    
    # Calculate the cumulative hazard via predict()
    newdata <- data.frame(y_star=C$t_e, delta_star=1, w1=z_0[1],
                          w2=z_0[2], a=z_0[3])
    pred <- predict(model, newdata=newdata, type="expected", se.fit=T)
    
    # Calculate certain true values
    H_0_true <- function(t) {
      L$sc_params$lmbd*exp(-1.7) * t^L$sc_params$v # !!!!!
    }
    true_lp <- sum(c(C$alpha_1,C$alpha_2,L$alpha_3)*z_0)
    
    # Calculate Fisher information and variance estimate
    res <- cox_var(dat=dat, theta_n=as.numeric(cfs), z_0=z_0, t=C$t_e)
    
    # Calculate the Breslow estimator
    bh <- basehaz(model, centered=FALSE)
    index <- max(which((bh$time<C$t_e)==T))
    est_bshz <- bh$hazard[index]
    
    # Calculate marginalized survival for A=0.5
    est_marg <- (1/L$n) * sum(sapply(c(1:L$n), function(i) {
      exp(-1*exp(sum(cfs*c(dat$w$w1[i],dat$w$w2[i],0.5)))*est_bshz)
    }))
    
    # S_0 <- function(t, w1, w2, a) {
    #   exp( -1 * sc_params$lmbd * (t^sc_params$v) * exp(lin(w1,w2,a)) )
    # }
    # theta_true_f <- Vectorize(function(a) {
    #   return(1 - mean(S_0(C$t_e,w1,w2,a)))
    # })
    
    # Return results
    return(list(
      true_bshz = H_0_true(C$t_e),
      true_cmhz = exp(true_lp) * H_0_true(C$t_e),
      true_surv = exp(-1*exp(true_lp)*H_0_true(C$t_e)),
      true_marg = 1-attr(dat_orig, "theta_true")[26], # Corresponds to A=0.5
      est_bshz = est_bshz,
      est_cmhz = exp(sum(cfs*z_0))*est_bshz, # Should equal pred$fit
      est_surv = exp(-1*exp(sum(cfs*z_0))*est_bshz),
      est_marg = est_marg,
      se_bshz_MC = sqrt(res$var_bshz_est), # Should roughly equal pred$se.fit
      se_cmhz_MC = sqrt(res$var_cmhz_est),
      se_surv_MC = sqrt(res$var_surv_est),
      se_marg_MC = 999
      # se_bshz_Cox = pred$se.fit,
      # est_w1_Cox = as.numeric(cfs[1]),
      # se_w1_Cox = as.numeric(sqrt(diag(vcov(model)))[1]),
      # se_w1_MC = as.numeric(sqrt(diag(res$I_tilde_inv)/L$n)[1])
    ))
    
  }
  
}
