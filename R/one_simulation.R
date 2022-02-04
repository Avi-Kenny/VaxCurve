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



#####################################.
##### Hypothesis testing (temp) #####
#####################################.

if (cfg$which_sim=="infl_fn_G (temp)") {
  
  one_simulation <- function() {
    
    params <- L$test$params # list(g_n_type="true", S_n_type="true", omega_n_type="true")
    
    # dat_orig <- generate_data(400, 0, "Unif(0,1)", "none", "Cox PH",
    #                           list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
    #                           "two-phase (50% random)", "decr")
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Prep
    {
      .default_params <- list(
        var="asymptotic", ecdf_type="step", g_n_type="binning", boot_reps=200,
        S_n_type="Super Learner", omega_n_type="estimated", cf_folds=1
      )
      for (i in c(1:length(.default_params))) {
        if (is.null(params[[names(.default_params)[i]]])) {
          params[[names(.default_params)[i]]] <- .default_params[[i]]
        }
      }
      
      a_lims <- c(min(dat_orig$a,na.rm=T),max(dat_orig$a,na.rm=T))
      a_shift <- -1 * a_lims[1]
      a_scale <- 1/(a_lims[2]-a_lims[1])
      dat_orig$a <- (dat_orig$a+a_shift)*a_scale
      
      dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
      for (i in c(1:length(dat_orig$w))) {
        rnd <- 8
        tol <- C$appx$w_tol
        n_unique <- tol + 1
        while(n_unique>tol) {
          rnd <- rnd - 1
          n_unique <- length(unique(round(dat_orig$w[,i],rnd)))
        }
        dat_orig$w[,i] <- round(dat_orig$w[,i],rnd)
      }
      
    }
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat, C$appx)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    lambda_2 <- lambda(dat,2,Phi_n)
    lambda_3 <- lambda(dat,3,Phi_n)
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                 type=params$omega_n_type)
    n_orig <- length(dat_orig$delta)
    z_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$a!=0))
    g_n_star <- construct_g_n_star(f_aIw_n, f_a_n, z_n)
    eta_ss_n <- construct_eta_ss_n(dat, S_n, z_n, vals=NA)
    gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
    alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, z_n, vals=NA)
    q_n <- construct_q_n(which="q_n", type="GAM", dat, dat_orig,
                         omega_n=omega_n, g_n_star=g_n_star, z_n=z_n,
                         gcomp_n=gcomp_n, alpha_star_n=alpha_star_n)
    Gamma_os_n <- construct_Gamma_os_n_star2(dat, dat_orig, omega_n,
                                             g_n_star, eta_ss_n, z_n, q_n,
                                             gcomp_n, alpha_star_n, vals=NA)
    infl_fn_Gamma <- construct_infl_fn_Gamma2(omega_n, g_n_star, gcomp_n, z_n,
                                              alpha_star_n, q_n, eta_ss_n,
                                              Gamma_os_n_star=Gamma_os_n)
    
    # Compute variance estimate and test statistic
    Theta_true <- attr(dat_orig,"Theta_true")
    Gamma_0 <- Vectorize(function(x) { # This only works for Unif(0,1)
      Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
    })
    xx <- 0.7
    Psi_G_est <- Gamma_os_n(xx)
    Psi_G_var_est <- (1/n_orig^2) * sum((
      infl_fn_Gamma(rep(xx, n_orig), dat_orig$w, dat_orig$y_star,
                    dat_orig$delta_star, dat_orig$a, dat_orig$weights)
    )^2)
    test_stat_Psi_G <- (Psi_G_est-Gamma_0(xx))^2/Psi_G_var_est
    p_val_Psi_G <- pchisq(test_stat_Psi_G, df=1, lower.tail=FALSE)
    
    # Estimates from old influence function
    eta_n <- construct_eta_n(dat, vlist$AW_grid, S_n)
    infl_fn_Gamma_old <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
                                             eta_n, Gamma_os_n)
    Psi_G_var_est2 <- (1/n_orig^2) * sum((
      dat$weights * infl_fn_Gamma_old(rep(xx,length(dat$a)), dat$w, dat$y_star,
                                      dat$delta_star, dat$a)
    )^2)
    test_stat_Psi_G2 <- (Psi_G_est-Gamma_0(xx))^2/Psi_G_var_est2
    p_val_Psi_G2 <- pchisq(test_stat_Psi_G2, df=1, lower.tail=FALSE)
    
    # Return results
    return (list(
      Psi_G_est = Psi_G_est-Gamma_0(xx),
      Psi_G_var_est = Psi_G_var_est,
      p_val_Psi_G = p_val_Psi_G,
      reject_Psi_G = as.integer(p_val_Psi_G<0.05),
      Psi_G_var_est2 = Psi_G_var_est2,
      p_val_Psi_G2 = p_val_Psi_G2,
      reject_Psi_G2 = as.integer(p_val_Psi_G2<0.05)
    ))
    
  }
  
}
