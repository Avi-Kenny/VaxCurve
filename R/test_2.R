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
    s <- stab(dat_orig)
    n_orig <- nrow(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- wts(dat, scale="none")
    
    # Construct vals
    vals_A <- data.frame(a=dat$a)
    vals_AW <- data.frame(a=dat$a, w1=dat$w1, w2=dat$w2)
    vals_eta_n <- expand.grid(x=seq(0,1,0.01), w1=seq(0,1,0.01), w2=c(0,1))
    vals_A_grid <- data.frame(a=seq(0,1,0.001))
    vals_S_n <- expand.grid(t=c(0:C$t_e), w1=seq(0,1,0.1), w2=c(0,1),
                            a=seq(0,1,0.1))
    
    # Construct component functions
    G_n <- construct_Phi_n(dat_orig)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(vals_AW, f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vals=dat, S_n, Sc_n)
    Gamma_n <- construct_Gamma_n(dat_orig, vals_A, omega_n, S_n, g_n)
    gcomp_n <- construct_gcomp_n(dat_orig, vals_A_grid, S_n)
    eta_n <- construct_eta_n(dat_orig, vals_eta_n, S_n)
    xi_n <- construct_xi_n(Phi_n=G_n, lambda_2, lambda_3)
    rho_n <- construct_rho_n(dat_orig, Phi_n=G_n)
    lambda_2 <- lambda(2,G_n,dat_orig)
    lambda_3 <- lambda(3,G_n,dat_orig)
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      (weights/s) * (
        lambda_2*(G_n(dat$a))^2 -
        lambda_3*G_n(dat$a)
      ) *
      (v("Gamma_n",dat$a))
    )
    
    # Estimate the test statistic variance
    if (params$est_known_nuis==FALSE) {
      lambda_2 <- 1/3
      lambda_3 <- 1/4
      rho_n <- function(x) { 0 }
    }
    infl_fn_1 <- construct_infl_fn_1(dat_orig, Gamma_n, Phi_n=G_n, xi_n, rho_n,
                                     lambda_2, lambda_3)
    infl_fn_Gamma <- construct_infl_fn_Gamma(dat_orig, omega_n, g_n, gcomp_n,
                                             eta_n, Gamma_n)
    infl_fn_2 <- construct_infl_fn_2(dat_orig, Phi_n=G_n, infl_fn_Gamma,
                                     lambda_2, lambda_3)
    
    var_hat <- beta_n_var_hat(dat_orig, infl_fn_1, infl_fn_2) / nrow(dat_orig)
    sd_hat <- sqrt(var_hat)
    
    if (return_sd) return(sd_hat)
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=beta_n, sd=sd_hat)
    
  }
  
  if (params$var=="boot") {
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat <- dat_orig[indices,]
      
      s_0 <- stab(dat_orig)
      n_orig <- nrow(dat)
      dat_0 <- dat %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0, scale="none")
      G_0 <- construct_Phi_n(dat_orig)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_orig, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      S_0 <- construct_S_n(dat, type=params$S_n_type)
      Sc_0 <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
      omega_0 <- construct_omega_n(S_0, Sc_0)
      Gamma_0 <- construct_Gamma_n(dat_orig, omega_0, S_0, g_0)
      
      beta_0 <- (1/n_orig) * sum(
        (weights_0/s_0) *
          (
            lambda(2, G_0, dat_orig)*(G_0(dat_0$a))^2 -
            lambda(3, G_0, dat_orig)*G_0(dat_0$a)
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
    
    # Pre-calculate non-bootstrapped pieces
    {
      dat_0_orig <- dat
      s_0 <- stab(dat_0_orig)
      n_orig <- nrow(dat_0_orig)
      dat_0 <- dat_0_orig %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0, scale="none")
      
      G_0 <- construct_Phi_n(dat_0_orig)
      S_0 <- construct_S_n(dat_0, type=params$S_n_type)
      Sc_0 <- construct_S_n(dat_0, type=params$S_n_type, csf=TRUE)
      omega_0 <- construct_omega_n(S_0, Sc_0)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_0_orig, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      Gamma_0 <- construct_Gamma_n(dat_0_orig, omega_0, S_0, g_0)
      lambda_2 <- lambda(k=2, G_0, dat_0_orig)
      lambda_3 <- lambda(k=3, G_0, dat_0_orig)
      eta_0 <- construct_eta_n(dat_0_orig, S_0)
      gcomp_0 <- construct_gcomp_n(dat_0_orig, S_0)

      beta_0 <- (1/n_orig) * sum(
        (weights_0/s_0) *
          (
            lambda(2, G_0, dat_0_orig)*(G_0(dat_0$a))^2 -
              lambda(3, G_0, dat_0_orig)*G_0(dat_0$a)
          ) *
          (Gamma_0(dat_0$a))
      )
      
      piece_3 <- -2*beta_0
      
    }
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat_b_orig <- dat_orig[indices,]
      s_b <- stab(dat_b_orig)
      dat_b <- dat_b_orig %>% filter(!is.na(a))
      n_b <- nrow(dat_b)
      weights_b <- wts(dat_b, scale="none")
      G_n <- construct_Phi_n(dat_b_orig)
      
      piece_1 <- (1/n_orig) * sum(
        (weights_b/s_b) *
          (
            lambda(k=2, G_0, dat_b_orig)*(G_0(dat_b$a))^2 -
              lambda(k=3, G_0, dat_b_orig)*G_0(dat_b$a)
          ) *
          (Gamma_0(dat_b$a))
      )
      
      piece_2 <- (1/n_orig) * sum(
        (weights_0/s_0) *
          (
            lambda(2, G_n, dat_0_orig)*(G_n(dat_0$a))^2 -
              lambda(3, G_n, dat_0_orig)*G_n(dat_0$a)
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
        (weights_0_long/s_0) * (weights_b_long/s_b) *
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
            lambda(2, G_0, dat_0_orig)*(G_0(a_0_long))^2 -
              lambda(3, G_0, dat_0_orig)*G_0(a_0_long)
          )
      )
      
      return (beta_0+piece_1+piece_2+piece_3+piece_4)
      
    }

    # Run bootstrap
    boot_obj <- boot(data=dat_0_orig, statistic=bootstat, R=params$boot_reps)
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  return(as.integer(crit_val>0))
  
}
