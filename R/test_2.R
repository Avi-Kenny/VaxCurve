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
#' @return Binary; is null rejected (1) or not (0)
test_2 <- function(dat_orig, alt_type="incr", params) {
  
  if (params$var=="asymptotic") {
    
    # Construct component functions
    ss <- ss(dat_orig)
    n_orig <- nrow(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights_n <- wts(dat, scale="none")
    G_n <- construct_Phi_n(dat_orig)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    Gamma_n <- construct_Gamma_n(dat_orig, omega_n, S_n, g_n)
    gcomp_n <- construct_gcomp(dat_orig, S_n=S_0)
    eta_n <- construct_eta_n(dat_orig, S_n=S_0)
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      (weights/ss) * (
        lambda(2, G_n, dat)*(G_n(dat$a))^2 -
        lambda(3, G_n, dat)*G_n(dat$a)
      ) *
      (Gamma_n(dat$a))
    )
    
    # Estimate the test statistic variance
    infl_fn_1 <- construct_infl_fn_1(dat_orig, Gamma_n, Phi_n=G_n)
    infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n, eta_n,
                                             Gamma_n)
    infl_fn_2 <- construct_infl_fn_2(dat_orig, Phi_n=G_n, infl_fn_Gamma)
    var_hat <- beta_n_var_hat(dat_orig, infl_fn_1, infl_fn_2) / nrow(dat_orig)
    sd_hat <- sqrt(var_hat)
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=beta_n, sd=sd_hat)
    
  }
  
  if (params$var=="boot") {
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat <- dat_orig[indices,]
      
      ss_0 <- ss(dat)
      n_orig <- nrow(dat)
      dat_0 <- dat %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0, scale="none")
      G_0 <- construct_Phi_n(dat)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_orig, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      S_0 <- construct_S_n(dat, type=params$S_n_type)
      Sc_0 <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
      omega_0 <- construct_omega_n(S_0, Sc_0)
      Gamma_0 <- construct_Gamma_n(dat_orig, omega_0, S_0, g_0)
      
      beta_0 <- (1/n_orig) * sum(
        (weights_0/ss_0) *
          (
            lambda(2, G_0, dat)*(G_0(dat_0$a))^2 -
              lambda(3, G_0, dat)*G_0(dat_0$a)
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
      ss_0 <- ss(dat_0_orig)
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
      lambda_2 <- lambda(k=2, G_0, dat_0)
      lambda_3 <- lambda(k=3, G_0, dat_0)
      eta_0 <- construct_eta_n(dat_0_orig, S_0)
      gcomp_0 <- construct_gcomp(dat_0_orig, S_0)

      beta_0 <- (1/n_orig) * sum(
        (weights_0/ss_0) *
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
      s_b <- ss(dat_b_orig)
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
        (weights_0/ss_0) *
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
        (weights_0_long/ss_0) * (weights_b_long/s_b) *
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
