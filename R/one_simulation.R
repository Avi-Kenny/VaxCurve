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
      res_list$g_n <- ests$g_n
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
    
    # Parse results object
    res <- list()
    # for (i in c(1:length(test_results))) {
    for (i in c(1:(length(test_results)-1))) { # !!!!!
      r <- test_results[[i]]
      res[[paste0("type_",i)]] <- r$type
      res[[paste0("reject_",i)]] <- as.integer(r$p_val<0.05)
      res[[paste0("p_val_",i)]] <- r$p_val
      res[[paste0("beta_n_",i)]] <- r$beta_n
      res[[paste0("var_n_",i)]] <- r$var_n
    }
    
    # !!!!! Debugging
    if (T) {
      res$Theta_0.2 <- test_results$extras$Theta_0.2
      res$Theta_0.5 <- test_results$extras$Theta_0.5
      res$etastar_0.2 <- test_results$extras$etastar_0.2
      res$etastar_0.5 <- test_results$extras$etastar_0.5
    }
    
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
    } # DEBUG
    
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
    g_sn <- construct_g_sn(dat, vlist$W_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    sigma2_edge_est <- sigma2_edge(dat, g_sn, Q_n, omega_n, r_Mn_edge_est)
    
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
  
  one_simulation <- function() {
    
    # # !!!!!
    # L <- list(n=500, alpha_3=0, distr_S="Unif(0,1)", edge="none",
    #           surv_true="Cox PH",
    #           sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
    #           sampling="two-phase (50%)", dir="decr", wts_type="true")
    # C <- list(alpha_1=0.5, alpha_2=0.7, t_0=200,
    #           appx=list(t_0=1, x_tol=25, s=0.01))
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir, L$wts_type)
    
    # Set default params
    params <- L$test$params
    .default_params <- list(
      type="simple", ecdf_type="step", g_n_type="binning", boot_reps=200,
      Q_n_type="Super Learner", omega_n_type="estimated", q_n_type="new",
      cf_folds=1, f_sIx_n_bins=15
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
    
    # Setup
    n_orig <- length(dat_orig$z)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type, print_coeffs=T)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
    f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                 k=p$f_sIx_n_bins, edge_corr="none",
                                 s_scale=s_scale, s_shift=s_shift)
    f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
    q_tilde_n <- construct_q_tilde_n(type=p$q_n_type, f_n_srv, f_sIx_n,
                                     omega_n)
    
    # Old etastar
    etastar_n <- construct_etastar_n(Q_n, vals=NA, tmp="old etastar")
    Theta_os_n <- construct_Theta_os_n(dat, dat_orig, omega_n, f_sIx_n,
                                       q_tilde_n, etastar_n)
    
    # New etastar
    etastar_n2 <- construct_etastar_n(Q_n, vals=NA, tmp="new etastar")
    Theta_os_n2 <- construct_Theta_os_n(dat, dat_orig, omega_n, f_sIx_n,
                                       q_tilde_n, etastar_n2)
    
    if (F) {
      
      construct_etastar_n2 <- function(Q_n, vals=NA) {
        fnc <- function(u,x) {
          u <- round(u,-log10(C$appx$s))
          if (u==0) {
            return(0)
          } else {
            smp <- round(runif(10^4),2)
            res <- mean(sapply(smp, function(s) {
              In(s<=u) * (1 - Q_n(C$t_0, x, s))
            }))
            return(res)
            # s_seq <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
            # integral <- C$appx$s * sum(sapply(s_seq, function(s) {
            #   Q_n(C$t_0, x, s)
            # }))
            # return(u-integral)
          }
        }
        return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
      }
      etastar_n2 <- construct_etastar_n2(Q_n, vals=NA)
      
      construct_etastar_n3 <- function(Q_n, vals=NA) {
        fnc <- function(u,x) {
          u <- round(u,-log10(C$appx$s))
          if (u==0) {
            return(0)
          } else {
            smp <- round(seq(0,1,C$appx$s),2)
            res <- mean(sapply(smp, function(s) {
              In(s<=u) * (1 - Q_n(C$t_0, x, s))
            }))
            return(res)
            # s_seq <- round(seq(C$appx$s,u,C$appx$s),-log10(C$appx$s))
            # integral <- C$appx$s * sum(sapply(s_seq, function(s) {
            #   Q_n(C$t_0, x, s)
            # }))
            # return(u-integral)
          }
        }
        return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
      }
      etastar_n3 <- construct_etastar_n3(Q_n, vals=NA)
      
    }
    
    # Decompose Theta_n
    piece_1 <- (dat$weights*omega_n(dat$x,dat$s,dat$y,dat$delta)) /
      f_sIx_n(dat$s,dat$x)
    piece_2 <- (1-dat_orig$weights)
    
    u <- 0.5
    Th_cmp_1 <- (1/(n_orig)) * sum(piece_1*In(dat$s<=u))
    Th_cmp_2 <- (1/n_orig) * sum(
      piece_2 * q_tilde_n(dat_orig$x, dat_orig$y, dat_orig$delta, u)
    )
    Th_cmp_3 <- (1/n_orig) * sum( etastar_n(rep(u,length(dat_orig$z)),dat_orig$x) )
    Theta_manual <- Th_cmp_1 + Th_cmp_2 + Th_cmp_3
    
    # Compare to eta_n
    p_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$s!=0))
    eta_n <- construct_eta_n(dat, Q_n, p_n, vals=NA)
    Th_cmp_3c <- (1/n_orig) * sum( etastar_n2(rep(u,length(dat_orig$z)),dat_orig$x) )
    
    if (F) {
      # # From est_curve.R
      # f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
      # g_n <- construct_g_n(f_sIx_n, f_s_n)
      # p_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$s!=0))
      # eta_n <- construct_eta_n(dat, Q_n, p_n, vals=NA)
      # gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$S_grid, Q_n)
      # alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, p_n, vals=NA)
      # f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
      # q_n <- construct_q_n(type=p$q_n_type, dat, dat_orig, omega_n=omega_n, g_n=g_n,
      #                      p_n=p_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
      #                      Q_n=Q_n, Qc_n=Qc_n, f_n_srv=f_n_srv)
      # Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n,
      #                                    q_n, gcomp_n, alpha_star_n,
      #                                    vals=vlist$S_grid)
    } # DEBUG
    
    return(list(
      Theta_0.5 = Theta_os_n(0.5),
      Theta2_0.5 = Theta_os_n2(0.5),
      Theta_manual = Theta_manual,
      Th_cmp_1 = Th_cmp_1,
      Th_cmp_2 = Th_cmp_2,
      Th_cmp_3 = Th_cmp_3,
      Th_cmp_3c = Th_cmp_3c,
      etastar1 = etastar_n(u=0.5, x=c(0,0)),
      etastar2 = etastar_n2(u=0.5, x=c(0,0))
    ))
    
    if (F) {
      # return(list(
      #   Theta_0.0 = Theta_os_n(0),
      #   Gamma_0.0 = Gamma_os_n(0),
      #   Theta_0.2 = Theta_os_n(0.2),
      #   Gamma_0.2 = Gamma_os_n(0.2),
      #   Theta_0.5 = Theta_os_n(0.5),
      #   Gamma_0.5 = Gamma_os_n(0.5),
      #   Theta_0.7 = Theta_os_n(0.7),
      #   Gamma_0.7 = Gamma_os_n(0.7),
      #   Theta_1.0 = Theta_os_n(1),
      #   Gamma_1.0 = Gamma_os_n(1)
      # ))
    } # DEBUG
    
  }
  
}
