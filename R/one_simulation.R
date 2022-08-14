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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir, L$wts_type)
    
    # Obtain estimates
    ests <- est_curve(
      dat_orig = dat_orig,
      estimator = L$estimator$est,
      params = L$estimator$params,
      points = C$points,
      dir = L$dir
      # return_extra = "deriv_r_Mn" # "Theta_os_n"
    )
    
    # Return results
    r_M0 <- attr(dat_orig, "r_M0")
    Gamma_true <- attr(dat_orig, "Gamma_true")
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=2)
      res_list[paste0("r_M0_",m)] <- r_M0[i]
      res_list[paste0("r_Mn_",m)] <- ests$est[i]
      res_list[paste0("ci_lo_",m)] <- ests$ci_lo[i]
      res_list[paste0("ci_hi_",m)] <- ests$ci_hi[i]
      if (F) {
        res_list[paste0("Gamma_",m)] <- Gamma_true[i]
        res_list[paste0("estG_",m)] <- ests$ests_Gamma[i]
        res_list[paste0("Phi_",m)] <- C$points[i] # Only works for Unif(0,1)
        res_list[paste0("estP_",m)] <- ests$ests_Phi[i]
      } # DEBUG: return Gamma/Phi estimates
    }
    
    if (F) {
      res_list$g_n_star <- ests$g_n_star
      res_list$eta_n3 <- ests$eta_n3
      res_list$eta_n5 <- ests$eta_n5
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir, L$wts_type)
    
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir, L$wts_type)
    
    # Prep
    n_orig <- length(dat_orig$z)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=L$estimator$params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n,
                                 type=params$omega_n_type)
    pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, pi_n, Q_n, omega_n)
    sigma2_edge_est <- sigma2_edge(dat, pi_n, Q_n, omega_n, r_Mn_edge_est)
    
    # Return results
    return(list(
      r_M0 = attr(dat_orig, "r_M0")[1],
      r_Mn = r_Mn_edge_est,
      ci_lo = r_Mn_edge_est - 1.96*sqrt(sigma2_edge_est/n_orig),
      ci_hi = r_Mn_edge_est + 1.96*sqrt(sigma2_edge_est/n_orig)
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
      L$n, L$alpha_3, L$distr_S, L$edge, surv_true="Cox PH",
      L$sc_params, L$sampling, L$dir, L$wts_type
    )
    
    # Round data values and construct dat
    dat_orig$s <- round(dat_orig$s,2)
    dat_orig$weights <- round(dat_orig$weights,3)
    dat_orig$y <- round(dat_orig$y,0)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    
    # Calculate variance estimates
    s <- 0.5
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_0, points=s,
                       se_beta=T, se_bshz=T, se_surv=T, se_marg=T)
    
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_0, points=0.5, se_marg=T, verbose=T)
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_0, points=c(0.2,0.5,0.8), se_marg=T, verbose=T)
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_0, points=seq(0.1,0.9,0.1), se_marg=T, verbose=T)
    
    # Calculate certain true values
    if (F) {
      z_0 <- c(0.3,1,0.5) # Needs to be consistent with the value in cox_var()
      H_0_true <- function(t) {
        L$sc_params$lmbd*exp(-1.7) * t^L$sc_params$v
      }
      true_lp <- sum(c(C$alpha_1,C$alpha_2,L$alpha_3)*z_0)
      true_surv <- exp(-1*exp(true_lp)*H_0_true(C$t_0))
      true_marg <- 1-attr(dat_orig, "r_M0")[26] # Corresponds to A=0.5
    } # DEBUG: intermediate objects
    
    # Construct simulation results object
    # This needs to line up with res_cox based on the se_* flags
    sim_res <- list(
      true_x1 = C$alpha_1,
      est_x1 = res_cox$beta_n[1],
      se_x1 = sqrt(res_cox$var_est_beta[1]),
      true_x2 = C$alpha_2,
      est_x2 = res_cox$beta_n[2],
      se_x2 = sqrt(res_cox$var_est_beta[2]),
      true_s = L$alpha_3,
      est_s = res_cox$beta_n[3],
      se_s = sqrt(res_cox$var_est_beta[3]),
      true_bshz = H_0_true(C$t_0),
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
      # newdata <- data.frame(y=C$t_0, delta=1, x1=z_0[1],
      #                       x2=z_0[2], s=z_0[3])
      # pred <- predict(res_cox$model, newdata=newdata, type="expected", se.fit=T)
      
      
      # Add additional results
      # sim_res$true_cmhz = exp(true_lp) * H_0_true(C$t_0) # Should roughly equal pred$se.fit
      # sim_res$est_cmhz = exp(sum(res_cox$beta_n*z_0))*est_bshz # Should equal pred$fit
      # sim_res$est_surv = exp(-1*exp(sum(res_cox$beta_n*z_0))*est_bshz)
      # sim_res$se_cmhz_MC = sqrt(res_cox$var_cmhz_est)
      # sim_res$se_surv_MC = sqrt(res_cox$var_surv_est)
      # sim_res$se_bshz_Cox = pred$se.fit
      
      # Debugging the Breslow estimator
      # sim_res$se_bshz_MC = sqrt(res_cox$var_bshz_est)
      # sim_res$se_bshz_MC2 = sqrt(res_cox$var_bshz_est2) # New derivation
      
      # sim_res$se_x1_Cox = as.numeric(sqrt(diag(vcov(res_cox$model)))[1])
      # sim_res$se_x1_alt = as.numeric(sqrt(diag(res_cox$I_tilde_inv)/L$n)[1])
      
    }
    
    return(sim_res)
    
  }
  
}



#####################.
##### Debugging #####
#####################.

if (cfg$which_sim=="debugging") {
  
  #' Run a single simulation (estimation)
  #'
  #' @return A list with the estimates and CI limits of the causal dose-response
  #'     curve evaluated at the midpoint and the endpoint of the domain
  
  one_simulation <- function() {
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir, L$wts_type)
    
    points <- C$points
    params <- L$estimator$params
    
    # Obtain estimates (Phi_n only)
    {
      # Set default params
      .default_params <- list(
        Q_n_type="Super Learner", g_n_type="binning", deriv_type="linear",
        ecdf_type="linear (mid)", gamma_type="Super Learner",
        omega_n_type="estimated", boot_reps=1000, ci_type="trunc", cf_folds=1, m=5,
        edge_corr="none", marg="Gamma_star2", lod_shift="none", n_bins=5,
        convex_type="GCM", f_sIx_n_bins=15
      )
      for (i in c(1:length(.default_params))) {
        if (is.null(params[[names(.default_params)[i]]])) {
          params[[names(.default_params)[i]]] <- .default_params[[i]]
        }
      }
      p <- params
      
      # Rescale S to lie in [0,1] and round values
      s_min <- min(dat_orig$s,na.rm=T)
      s_max <- max(dat_orig$s,na.rm=T)
      s_shift <- -1 * s_min
      s_scale <- 1/(s_max-s_min)
      dat_orig$s <- (dat_orig$s+s_shift)*s_scale
      dat_orig <- round_dat(dat_orig)
      
      # Obtain minimum value (excluding edge point mass)
      if (p$edge_corr=="min") { a_min2 <- min(dat_orig$s[dat_orig$s!=0],na.rm=T) }
      
      # Rescale points and remove points outside the range of S
      points_orig <- points
      na_head <- sum(round(points,-log10(C$appx$s))<round(s_min,-log10(C$appx$s)))
      points <- round((points+s_shift)*s_scale, -log10(C$appx$s))
      na_tail <- sum(points>1)
      if (na_head>0) {
        points <- points[-c(1:na_head)]
      }
      if (na_tail>0) {
        points <- points[-c((length(points)-na_tail+1):length(points))]
      }
      
      dat <- ss(dat_orig, which(dat_orig$z==1))
      dat2 <- ss(dat, which(dat$s!=0))
      Phi_n <- construct_Phi_n(dat2, type=p$ecdf_type)
      ests_Phi <- Phi_n(points)
      
    }
    
    # Construct new estimator #2 (manual summation calc)
    n_orig <- sum(dat2$weights)
    Phi_n2 <- Vectorize(function(x) {
      (1/n_orig) * sum(dat2$weights*as.integer(dat2$s<=x))
    })
    ests_Phi2 <- Phi_n2(points)
    
    # Construct new estimator #3 (manual summation calc with original n-value)
    n_orig <- sum(dat$weights)
    Phi_n3 <- Vectorize(function(x) {
      (1/n_orig) * sum(dat2$weights*as.integer(dat2$s<=x))
    })
    ests_Phi3 <- Phi_n3(points)
    
    # Construct new estimator #4 (not excluding edge points)
    Phi_n4 <- construct_Phi_n(dat, type=p$ecdf_type)
    ests_Phi4 <- Phi_n4(points)
    
    # Return results
    res_list <- list()
    for (i in 1:length(C$points)) {
      m <- format(C$points[i], nsmall=1)
      res_list[paste0("Phi_",m)] <- C$points[i] # Only works for Unif(0,1)
      res_list[paste0("estP_",m)] <- ests_Phi[i]
      res_list[paste0("estP2_",m)] <- ests_Phi2[i]
      res_list[paste0("estP3_",m)] <- ests_Phi3[i]
      res_list[paste0("estP4_",m)] <- ests_Phi4[i]
    }
    
    return(res_list)
    
  }
  
}
