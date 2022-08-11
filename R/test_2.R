#' Testing approach 2: regression slope
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param alt_type Type of alternative hypothesis; either "incr", "decr", or
#'     "two-tailed"; currently unused and defaults to "two-tailed".
#' @param params A list, containing the following:
#'   - `S_n_type` S_n_type Type of survival function estimator; currently only
#'     c("Cox PH")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `var` Variance estimation; one of c("boot","mixed boot")
#'   - `boot_reps` Number of bootstrap replicates to run
#' @param test_stat_only Boolean; if TRUE, only calculate the test statistic and
#'     do not estimate its variance
#' @return Binary; is null rejected (1) or not (0)
test_2 <- function(dat_orig, alt_type="two-tailed", params,
                   test_stat_only=F, return_extras=F) {
  
  # Set default params
  .default_params <- list(
    var="monte carlo", ecdf_type="step", g_n_type="binning", boot_reps=200,
    S_n_type="Super Learner", omega_n_type="estimated", cf_folds=1,
    f_aIw_n_bins=15
  )
  for (i in c(1:length(.default_params))) {
    if (is.null(params[[names(.default_params)[i]]])) {
      params[[names(.default_params)[i]]] <- .default_params[[i]]
    }
  }
  p <- params
  
  # Rescale A to lie in [0,1] and round values
  a_min <- min(dat_orig$a,na.rm=T)
  a_max <- max(dat_orig$a,na.rm=T)
  a_shift <- -1 * a_min
  a_scale <- 1/(a_max-a_min)
  if (F) {
    a_shift <- 0
    a_scale <- 1
  } # DEBUG
  dat_orig$a <- (dat_orig$a+a_shift)*a_scale
  dat_orig <- round_dat(dat_orig)
  
  if (p$var=="asymptotic") {
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=p$ecdf_type)
    lambda_2 <- lambda(dat,2,Phi_n)
    lambda_3 <- lambda(dat,3,Phi_n)
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
                                 k=p$f_aIw_n_bins)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type)
    S_n <- srvSL$srv
    Sc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type=p$omega_n_type)
    
    # Construct regular Gamma_0 estimator
    if (p$cf_folds==1) {
      # Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
      
      # !!!!! New functions
      n_orig <- length(dat_orig$delta)
      z_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$a!=0))
      g_n_star <- construct_g_n_star(f_aIw_n, f_a_n, z_n)
      eta_ss_n <- construct_eta_ss_n(dat, S_n, z_n, vals=NA)
      gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
      alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, z_n, vals=NA)
      # Gamma_os_n <- construct_Gamma_os_n_star(dat, omega_n, g_n_star,
      #                                              eta_ss_n, z_n, gcomp_n,
      #                                              alpha_star_n, vals=NA)
      q_n <- construct_q_n(type="new", dat, dat_orig, omega_n=omega_n,
                           g_n_star=g_n_star, z_n=z_n, gcomp_n=gcomp_n,
                           alpha_star_n=alpha_star_n)
      Gamma_os_n <- construct_Gamma_os_n_star2(dat, dat_orig, omega_n, g_n_star,
                                               eta_ss_n, z_n, q_n, gcomp_n,
                                               alpha_star_n, vals=NA)
    }
    
    # Construct cross-fitted Gamma_0 estimator
    if (p$cf_folds>1) {
      # Gamma_os_n <- construct_Gamma_cf(dat_orig, p, vlist)
    }
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * ( lambda_2*(Phi_n(dat$a))^2 - lambda_3*Phi_n(dat$a) ) *
        Gamma_os_n(round(dat$a,-log10(C$appx$a)))
    )
    
    if (test_stat_only) {
      
      p_val <- NA
      sd_n <- NA
      
    } else {
      
      # Construct influence functions
      # eta_n <- construct_eta_n(dat, vlist$AW_grid, S_n)
      infl_fn_1 <- construct_infl_fn_1(dat, Gamma_os_n, Phi_n, lambda_2,
                                       lambda_3)
      # infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
      #                                          eta_n, Gamma_os_n)
      infl_fn_Gamma <- construct_infl_fn_Gamma2(omega_n, g_n_star, gcomp_n, z_n,
                                                alpha_star_n, q_n, eta_ss_n,
                                                Gamma_os_n_star=Gamma_os_n)
      infl_fn_2 <- construct_infl_fn_2(dat, Phi_n, infl_fn_Gamma, lambda_2,
                                       lambda_3)
      
      # Estimate variance
      var_n <- beta_n_var_hat(dat, dat_orig, infl_fn_1, infl_fn_2) / n_orig
      sd_n <- sqrt(var_n)
      
      # Calculate P-value
      if (alt_type=="incr") {
        p_val <- pnorm(beta_n, mean=0, sd=sd_n, lower.tail=FALSE)
      } else if (alt_type=="decr") {
        p_val <- pnorm(beta_n, mean=0, sd=sd_n)
      } else if (alt_type=="two-tailed") {
        test_stat_chisq <- beta_n^2/var_n
        p_val <- pchisq(test_stat_chisq, df=1, lower.tail=FALSE)
        
        if (F) {
          test_stat_chisq_alt <- (Psi_1_est+Psi_2_est)^2/var_n
          p_val_alt <- pchisq(test_stat_chisq_alt, df=1, lower.tail=FALSE)
        } # DEBUG: alternate test stat
      }
      
    }
    
  }
  
  if (p$var=="monte carlo") {
    
    # Setup
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    # f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
    #                              k=p$f_aIw_n_bins)
    # srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type)
    # S_n <- srvSL$srv
    # Sc_n <- srvSL$cens
    # omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type=p$omega_n_type)
    # etastar_n <- construct_etastar_n(S_n)
    # q_star_n <- construct_q_n(which="q_star_n", type="Super Learner", dat, # ?????
    #                           dat_orig, omega_n=omega_n, f_aIw_n=f_aIw_n) # ?????
    Theta_os_n <- construct_Theta_os_n2(dat, dat_orig, omega_n, f_aIw_n,
                                        q_star_n, etastar_n, vals=NA)
    # infl_fn_Theta <- construct_infl_fn_Theta(omega_n, f_aIw_n, q_star_n,
    #                                          etastar_n, Theta_os_n)
    
    # Construct pieces needed for hypothesis test
    u_mc <- round(seq(C$appx$a,1,C$appx$a),-log10(C$appx$a))
    m <- length(u_mc)
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    
    # Compute test statistic and variance estimate
    beta_n <- mean((lambda_2*u_mc^2-lambda_3*u_mc)*Theta_os_n(u_mc))
    var_n <- 0
    for (i in c(1:n_orig)) {
      a_m <- rep(dat_orig$a[i],m)
      y_star_m <- rep(dat_orig$y_star[i],m)
      delta_star_m <- rep(dat_orig$delta_star[i],m)
      weight_m <- rep(dat_orig$weight[i],m)
      w_m <- as.data.frame(
        matrix(rep(dat_orig$w[i,],m), ncol=length(dat_orig$w[i,]), byrow=T)
      )
      var_n <- var_n + (sum(
        infl_fn_Theta(x=u_mc, w_m, y_star_m, delta_star_m, a_m, weight_m) *
          (lambda_2*u_mc^2-lambda_3*u_mc)
      ))^2
    }
    var_n <- var_n/(n_orig^2*m^2)
    sd_n <- sqrt(var_n)
    
    # Compute P-value
    test_stat_chisq <- beta_n^2/var_n
    p_val <- pchisq(test_stat_chisq, df=1, lower.tail=FALSE)
    
  }
  
  if (p$var=="boot") {
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat_orig_b <- ss(dat_orig, indices)
      
      # Calculate components on bootstrapped dataset
      {
        # Prep
        n_orig <- length(dat_orig_b$delta)
        dat <- ss(dat_orig_b, which(dat_orig_b$delta==1))
        weights <- dat$weights
        
        # Construct dataframes of values to pre-compute functions on
        vlist <- create_val_list(dat_orig)
        
        # Construct component functions
        Phi_n <- construct_Phi_n(dat, type=p$ecdf_type)
        lambda_2 <- lambda(dat,2,Phi_n)
        lambda_3 <- lambda(dat,3,Phi_n)
        f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
                                     k=p$f_aIw_n_bins)
        f_a_n <- construct_f_a_n(dat_orig_b, vlist$A_grid, f_aIw_n)
        g_n <- construct_g_n(f_aIw_n, f_a_n)
        srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type)
        S_n <- srvSL$srv
        Sc_n <- srvSL$cens
        omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                     type=p$omega_n_type)
        Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
        
      }
      
      # Compute the test statistic
      beta_n <- (1/n_orig) * sum(
        weights * ( lambda_2*(Phi_n(dat$a))^2 - lambda_3*Phi_n(dat$a) ) *
          Gamma_os_n(round(dat$a,-log10(C$appx$a)))
      )
      
      if (F) {
        
        # Psi_1
        Theta_true <- attr(dat_orig,"Theta_true")
        Gamma_0 <- Vectorize(function(x) {
          Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
        })
        Psi_1_est <- (1/n_orig) * sum(weights*(
          lambda_2*(Phi_n(dat$a))^2*Gamma_0(dat$a) -
            lambda_3*Phi_n(dat$a)*Gamma_0(dat$a)
        ))
        
        # Psi_2
        Phi_0 <- function(x) {x}
        a_mc <- runif(10^6)
        Psi_2_est <- mean(
          (1/3)*(Phi_0(a_mc))^2*Gamma_os_n(round(a_mc,-log10(C$appx$a))) -
            (1/4)*Phi_0(a_mc)*Gamma_os_n(round(a_mc,-log10(C$appx$a)))
        )
        
      } # DEBUG
      
      return(list(beta_n=beta_n, Psi_1_est=Psi_1_est, Psi_2_est=Psi_2_est))
      
    }
    
    # Run bootstrap
    boot_results <- sapply(c(1:p$boot_reps), function(i) {
      indices <- sample(c(1:length(dat_orig$a)), replace=T)
      bootstat(dat_orig,indices)
    })
    
    # Calculate beta_n, SD, variance
    stats <- bootstat(dat_orig, c(1:length(dat_orig$delta)))
    beta_n <- stats$beta_n
    sd_n <- sd(as.numeric(boot_results["beta_n",]))
    var_n <- var(as.numeric(boot_results["beta_n",]))
    
    Psi_1_est <- stats$Psi_1_est
    Psi_2_est <- stats$Psi_2_est
    Psi_1_var_est <- var(as.numeric(boot_results["Psi_1_est",]))
    Psi_2_var_est <- var(as.numeric(boot_results["Psi_2_est",]))
    Psi_12_covar <- cov(
      as.numeric(boot_results["Psi_1_est",]),
      as.numeric(boot_results["Psi_2_est",])
    )
    
    # Calculate P-values
    test_stat_chisq <- beta_n^2/var_n
    p_val <- pchisq(test_stat_chisq, df=1, lower.tail=FALSE)
    test_stat_Psi_1 <- Psi_1_est^2/Psi_1_var_est
    p_val_Psi_1 <- pchisq(test_stat_Psi_1, df=1, lower.tail=FALSE)
    test_stat_Psi_2 <- Psi_2_est^2/Psi_2_var_est
    p_val_Psi_2 <- pchisq(test_stat_Psi_2, df=1, lower.tail=FALSE)
    
  }
  
  if (p$var=="mixed boot") {
    
    # # !!!!! Update all of this if needed
    # 
    # # Pre-calculate non-bootstrapped pieces
    # {
    #   dat_0_orig <- dat
    #   n_orig <- length(dat_0_orig$delta)
    #   dat_0 <- ss(dat_0_orig, which(dat_0_orig$delta==1))
    #   n_0 <- length(dat_0$a)
    #   weights_0 <- wts(dat_0)
    #   
    #   G_0 <- construct_Phi_n(dat_0, type=p$ecdf_type)
    #   srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type)
    #   S_0 <- srvSL$srv
    #   Sc_0 <- srvSL$cens
    #   omega_0 <- construct_omega_n(S_0, Sc_0)
    #   f_aIw_n <- construct_f_aIw_n(dat_0, type=p$g_n_type, k=p$f_aIw_n_bins)
    #   f_a_n <- construct_f_a_n(dat_0_orig, f_aIw_n=f_aIw_n)
    #   g_0 <- construct_g_n(f_aIw_n, f_a_n)
    #   Gamma_0 <- construct_Gamma_os_n(dat_0, vlist$A_grid, omega_0, S_0, g_0)
    #   lambda_2 <- lambda(dat_0_orig, k=2, G_0)
    #   lambda_3 <- lambda(dat_0_orig, k=3, G_0)
    #   eta_0 <- construct_eta_n(dat_0, vlist$AW_grid, S_0)
    #   gcomp_0 <- construct_gcomp_n(dat_0_orig, vlist$A_grid, S_0)
    # 
    #   beta_0 <- (1/n_orig) * sum(
    #     weights_0 *
    #       (
    #         lambda(dat_0_orig,2,G_0)*(G_0(dat_0$a))^2 -
    #           lambda(dat_0_orig,3,G_0)*G_0(dat_0$a)
    #       ) *
    #       (Gamma_0(dat_0$a))
    #   )
    #   
    #   piece_3 <- -2*beta_0
    #   
    # }
    # 
    # # Define the statistic to bootstrap
    # bootstat <- function(dat_orig,indices) {
    #   
    #   dat_b_orig <- dat_orig[indices,]
    #   dat_b <- ss(dat_b_orig, which(dat_b_orig$delta==1))
    #   n_b <- length(dat_b$a)
    #   weights_b <- wts(dat_b)
    #   Phi_n <- construct_Phi_n(dat_b, type=p$ecdf_type)
    #   
    #   piece_1 <- (1/n_orig) * sum(
    #     weights_b *
    #       (
    #         lambda(dat_b_orig,k=2,G_0)*(G_0(dat_b$a))^2 -
    #           lambda(dat_b_orig,k=3,G_0)*G_0(dat_b$a)
    #       ) *
    #       (Gamma_0(dat_b$a))
    #   )
    #   
    #   piece_2 <- (1/n_orig) * sum(
    #     weights_0 *
    #       (
    #         lambda(dat_0_orig,2,Phi_n)*(Phi_n(dat_0$a))^2 -
    #           lambda(dat_0_orig,3,Phi_n)*Phi_n(dat_0$a)
    #       ) *
    #       Gamma_0(dat_0$a)
    #   )
    #   
    #   index_0 <- rep(c(1:n_0), each=n_b)
    #   index_b <- rep(c(1:n_b), times=n_0)
    #   a_0_long <- dat_0$a[index_0]
    #   a_b_long <- dat_b$a[index_b]
    #   y_b_long <- dat_b$y[index_b]
    #   w1_b_long <- dat_b$w1[index_b]
    #   w2_b_long <- dat_b$w2[index_b]
    #   weights_0_long <- weights_0[index_0]
    #   weights_b_long <- weights_b[index_b]
    #   
    #   piece_4 <- (1/n_orig)^2 * sum(
    #     weights_0_long * weights_b_long *
    #       (
    #         (
    #           (as.integer(a_b_long<=a_0_long)) *
    #             (
    #               ( y_b_long - mu_0(a_b_long, w1_b_long, w2_b_long) ) /
    #                 g_0(a_b_long, w1_b_long, w2_b_long) +
    #                 gcomp_0(a_b_long)
    #             )
    #         ) +
    #           eta_0(a_0_long, w1_b_long, w2_b_long) -
    #           (2*Gamma_0(a_0_long))
    #       ) *
    #       (
    #         lambda(dat_0_orig,2,G_0)*(G_0(a_0_long))^2 -
    #           lambda(dat_0_orig,3,G_0)*G_0(a_0_long)
    #       )
    #   )
    #   
    #   return (beta_0+piece_1+piece_2+piece_3+piece_4)
    #   
    # }
    # 
    # # Run bootstrap
    # boot_obj <- boot(data=dat_0_orig, statistic=bootstat, R=p$boot_reps)
    # 
    # # Calculate critical value (for a one-sided test)
    # crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  res <- list(
    reject = as.integer(p_val<0.05),
    p_val = p_val,
    beta_n = beta_n,
    sd_n = sd_n,
    var_n = var_n
  )
  
  # Return debugging components (var="Monte Carlo")
  if (return_extras) {
    res$Theta_os_n <- Theta_os_n
  }
  
  if (F) {
    # res$if1_mean <- if1_mean
    # res$if2_mean <- if2_mean
    # res$r_1n <- Psi_1_est - if1_mean
    # res$r_2n <- Psi_2_est - if2_mean
    # res$Psi_1_var_est <- Psi_1_var_est
    # res$sum12_est <- Psi_1_est+Psi_2_est
    # res$sum12_var_est <- sum12_var_est
    # res$Psi_1_est <- Psi_1_est
    # res$Psi_2_est <- Psi_2_est
    # res$Psi_G_est <- Psi_G_est-Gamma_0(xx)
    # res$Psi_1_var_est <- Psi_1_var_est
    # res$Psi_2_var_est <- Psi_2_var_est
    # res$Psi_12_covar <- Psi_12_covar
    # res$p_val_Psi_1 <- p_val_Psi_1
    # res$p_val_Psi_2 <- p_val_Psi_2
    # res$p_val_Psi_G <- p_val_Psi_G
    # res$p_val_sum12 <- p_val_sum12
    # res$p_val_sum12b <- p_val_sum12b
    # res$p_val_alt <- p_val_alt
    # res$reject_Psi_1 <- as.integer(p_val_Psi_1<0.05)
    # res$reject_Psi_2 <- as.integer(p_val_Psi_2<0.05)
    # res$reject_Psi_G <- as.integer(p_val_Psi_G<0.05)
    # res$reject_sum12 <- as.integer(p_val_sum12<0.05)
    # res$reject_sum12b <- as.integer(p_val_sum12b<0.05)
    # res$reject_alt <- as.integer(p_val_alt<0.05)
  } # DEBUG: return debugging components (var="asymptotic")
  
  return(res)
  
}
