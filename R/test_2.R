#' Testing approach 2: regression slope
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr";
#'     currently unused
#' @param params A list, containing the following:
#'   - `S_n_type` S_n_type Type of survival function estimator; currently only
#'     c("Cox PH")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `var` Variance estimation; one of c("boot","mixed boot")
#'   - `boot_reps` Number of bootstrap replicates to run
#' @param return_sd Boolean; if TRUE, return the standard deviation estimate
#'     instead of the test result
#' @return Binary; is null rejected (1) or not (0)
test_2 <- function(dat_orig, alt_type="incr", params, return_sd=FALSE) {
  
  if (params$var=="asymptotic") {
    
    # Prep
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig, C$appx)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat_orig)
    lambda_2 <- lambda(dat_orig,2,Phi_n)
    lambda_3 <- lambda(dat_orig,3,Phi_n)
    
    # Construct component functions
    f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid, type=params$g_n_type,
                                 k=15)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    
    # Construct regular Gamma_0 estimator
    if (params$cf_folds==1) {
      Gamma_n <- construct_Gamma_n(dat_orig, vlist$A_grid, omega_n, S_n, g_n)
    }
    
    # Construct cross-fitted Gamma_0 estimator
    if (params$cf_folds>1) {
      Gamma_n <- construct_Gamma_cf(dat_orig, params, vlist)
    }
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * (
        lambda_2*(Phi_n(dat$a))^2 -
          lambda_3*Phi_n(dat$a)
      ) *
        (Gamma_n(dat$a))
    )
    
    # Use known values of nuisances to estimate variance more accurately
    if (params$est_known_nuis==FALSE) {
      lambda_2 <- 1/3
      lambda_3 <- 1/4
      rho_n <- function(x) { 0 }
    }
    
    # Construct additional component functions
    gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
    eta_n <- construct_eta_n(dat_orig, vlist$AW_grid, S_n)
    xi_n <- construct_xi_n(Phi_n, lambda_2, lambda_3)
    # rho_n <- construct_rho_n(dat_orig, Phi_n)
    infl_fn_1 <- construct_infl_fn_1(dat_orig, Gamma_n, Phi_n, xi_n, rho_n,
                                     lambda_2, lambda_3)
    infl_fn_Gamma <- construct_infl_fn_Gamma(dat_orig, omega_n, g_n, gcomp_n,
                                             eta_n, Gamma_n)
    infl_fn_2 <- construct_infl_fn_2(dat_orig, Phi_n, infl_fn_Gamma,
                                     lambda_2, lambda_3)
    
    # Estimate variance
    var_n <- beta_n_var_hat(dat_orig, infl_fn_1, infl_fn_2) / nrow(dat_orig)
    sd_n <- sqrt(var_n)
    
    if (return_sd) return(sd_n)
    
    # Calculate P-value (for a one-sided test)
    p_val <- pnorm(0, mean=beta_n, sd=sd_n)

  }
  
  if (params$var=="boot") {
    
    # !!!!! Update all of this
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat <- dat_orig[indices,]
      n_orig <- nrow(dat)
      dat_0 <- dat %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0)
      G_0 <- construct_Phi_n(dat_orig)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      S_0 <- construct_S_n(dat, type=params$S_n_type)
      Sc_0 <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
      omega_0 <- construct_omega_n(S_0, Sc_0)
      Gamma_0 <- construct_Gamma_n(dat_orig, omega_0, S_0, g_0)
      
      beta_0 <- (1/n_orig) * sum(
        weights_0 *
          (
            lambda(dat_orig,2,G_0)*(G_0(dat_0$a))^2 -
              lambda(dat_orig,3,G_0)*G_0(dat_0$a)
          ) *
          (Gamma_0(dat_0$a))
      )
      
      return (beta_0)
      
    }
    
    # Run bootstrap
    boot_obj <- boot(data=dat, statistic=bootstat, R=params$boot_reps)
    
    # Calculate beta_n
    beta_0 <- bootstat(dat, c(1:nrow(dat)))
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  if (params$var=="mixed boot") {
    
    # !!!!! Update all of this
    
    # Pre-calculate non-bootstrapped pieces
    {
      dat_0_orig <- dat
      n_orig <- nrow(dat_0_orig)
      dat_0 <- dat_0_orig %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0)
      
      G_0 <- construct_Phi_n(dat_0_orig)
      S_0 <- construct_S_n(dat_0, type=params$S_n_type)
      Sc_0 <- construct_S_n(dat_0, type=params$S_n_type, csf=TRUE)
      omega_0 <- construct_omega_n(S_0, Sc_0)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_0_orig, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      Gamma_0 <- construct_Gamma_n(dat_0_orig, omega_0, S_0, g_0)
      lambda_2 <- lambda(dat_0_orig, k=2, G_0)
      lambda_3 <- lambda(dat_0_orig, k=3, G_0)
      eta_0 <- construct_eta_n(dat_0_orig, S_0)
      gcomp_0 <- construct_gcomp_n(dat_0_orig, S_0)

      beta_0 <- (1/n_orig) * sum(
        weights_0 *
          (
            lambda(dat_0_orig,2,G_0)*(G_0(dat_0$a))^2 -
              lambda(dat_0_orig,3,G_0)*G_0(dat_0$a)
          ) *
          (Gamma_0(dat_0$a))
      )
      
      piece_3 <- -2*beta_0
      
    }
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat_b_orig <- dat_orig[indices,]
      dat_b <- dat_b_orig %>% filter(!is.na(a))
      n_b <- nrow(dat_b)
      weights_b <- wts(dat_b)
      Phi_n <- construct_Phi_n(dat_b_orig)
      
      piece_1 <- (1/n_orig) * sum(
        weights_b *
          (
            lambda(dat_b_orig,k=2,G_0)*(G_0(dat_b$a))^2 -
              lambda(dat_b_orig,k=3,G_0)*G_0(dat_b$a)
          ) *
          (Gamma_0(dat_b$a))
      )
      
      piece_2 <- (1/n_orig) * sum(
        weights_0 *
          (
            lambda(dat_0_orig,2,Phi_n)*(Phi_n(dat_0$a))^2 -
              lambda(dat_0_orig,3,Phi_n)*Phi_n(dat_0$a)
          ) *
          Gamma_0(dat_0$a)
      )
      
      index_0 <- rep(c(1:n_0), each=n_b)
      index_b <- rep(c(1:n_b), times=n_0)
      a_0_long <- dat_0$a[index_0]
      a_b_long <- dat_b$a[index_b]
      y_b_long <- dat_b$y[index_b]
      w1_b_long <- dat_b$w1[index_b]
      w2_b_long <- dat_b$w2[index_b]
      weights_0_long <- weights_0[index_0]
      weights_b_long <- weights_b[index_b]
      
      piece_4 <- (1/n_orig)^2 * sum(
        weights_0_long * weights_b_long *
          (
            (
              (as.integer(a_b_long<=a_0_long)) *
                (
                  ( y_b_long - mu_0(a_b_long, w1_b_long, w2_b_long) ) /
                    g_0(a_b_long, w1_b_long, w2_b_long) +
                    gcomp_0(a_b_long)
                )
            ) +
              eta_0(a_0_long, w1_b_long, w2_b_long) -
              (2*Gamma_0(a_0_long))
          ) *
          (
            lambda(dat_0_orig,2,G_0)*(G_0(a_0_long))^2 -
              lambda(dat_0_orig,3,G_0)*G_0(a_0_long)
          )
      )
      
      return (beta_0+piece_1+piece_2+piece_3+piece_4)
      
    }

    # Run bootstrap
    boot_obj <- boot(data=dat_0_orig, statistic=bootstat, R=params$boot_reps)
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  return(list(
    reject = as.integer(p_val<0.05),
    p_val = p_val
  ))
  
}
