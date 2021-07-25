#' Testing approach 2: regression slope
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr";
#'     currently unused
#' @param mu_n_type Type of regression estimator; one of One of c("Logistic",
#'     "GAM", "Random forest")
#' @param g_n_type Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#' @param params A list, containing the following:
#'   - `mu_n_type` Type of regression estimator; one of One of c("Logistic",
#'     "GAM", "Random forest")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `var` Variance estimation; one of c("boot","mixed boot")
#'   - `boot_reps` Number of bootstrap replicates to run
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - Note

test_2 <- function(dat, alt_type="incr", params) {
  
  if (params$var=="boot") {
    
    # Define the statistic to bootstrap
    bootstat <- function(dat_orig,indices) {
      
      dat <- dat_orig[indices,]
      
      s_0 <- ss(dat)
      n_orig <- nrow(dat)
      dat_0 <- dat %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0, scale="none")
      G_0 <- construct_Phi_n(dat)
      mu_0 <- construct_mu_n(dat_0, type=params$mu_n_type)
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_0, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      Gamma_0 <- construct_Gamma_n(dat, mu_0, g_0)
      
      beta_0 <- (1/n_orig) * sum(
        (weights_0/s_0) *
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
    
  } else if (params$var=="mixed boot") {
    
    # Pre-calculate non-bootstrapped pieces
    {
      dat_0_orig <- dat
      s_0 <- ss(dat_0_orig)
      n_orig <- nrow(dat_0_orig)
      dat_0 <- dat_0_orig %>% filter(!is.na(a))
      n_0 <- nrow(dat_0)
      weights_0 <- wts(dat_0, scale="none")
      
      G_0 <- construct_Phi_n(dat_0_orig)
      mu_0 <- construct_mu_n(dat_0, type=params$mu_n_type)
      
      f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_0, f_aIw_n=f_aIw_n)
      g_0 <- construct_g_n(f_aIw_n, f_a_n)
      Gamma_0 <- construct_Gamma_n(dat_0_orig, mu_0, g_0)
      lambda_2 <- lambda(k=2, G_0, dat_0)
      lambda_3 <- lambda(k=3, G_0, dat_0)
      eta_0 <- construct_eta_n(dat_0_orig, mu_0)
      theta_naive_0 <- construct_theta_naive_n(dat_0_orig, mu_0)
      
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
                    theta_naive_0(a_b_long)
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
    
  } else {
    stop("Invalid specification for params$var")
  }
  
  return(as.integer(crit_val>0))
  
}
