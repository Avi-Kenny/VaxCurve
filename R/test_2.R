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
                   test_stat_only=FALSE) {
  
  if (params$var=="asymptotic") {
    
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
    f_aIw_n <- 999 # construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
    f_a_n <- 999 # construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- function(a,w) {1} # construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                 type=params$omega_n_type)
    
    # Construct regular Gamma_0 estimator
    if (params$cf_folds==1) {
      Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    }
    
    # Construct cross-fitted Gamma_0 estimator
    if (params$cf_folds>1) {
      # Gamma_os_n <- construct_Gamma_cf(dat_orig, params, vlist)
    }
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * (
        lambda_2*(Phi_n(dat$a))^2 -
          lambda_3*Phi_n(dat$a)
      ) *
        Gamma_os_n(round(dat$a,-log10(C$appx$a)))
    )
    
    if (test_stat_only) {
      
      p_val <- NA
      sd_n <- NA
      
    } else {
      
      # Construct additional component functions
      gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
      eta_n <- construct_eta_n(dat, vlist$AW_grid, S_n)
      infl_fn_1 <- construct_infl_fn_1(dat, Gamma_os_n, Phi_n,
                                       lambda_2, lambda_3)
      infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
                                               eta_n, Gamma_os_n)
      infl_fn_2 <- construct_infl_fn_2(dat, Phi_n, infl_fn_Gamma,
                                       lambda_2, lambda_3)
      
      # Estimate variance
      var_n <- beta_n_var_hat(dat, infl_fn_1, infl_fn_2) / n_orig
      sd_n <- sqrt(var_n)
      
      # !!!!!
      if (!is.null(L$temp_which)) {
        
        if (L$temp_which=="Gamma") {
          
          # x <- 0.3
          # Theta_true <- attr(dat_orig,"Theta_true")
          # Gamma_0 <- Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
          # infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
          #                                          eta_n, Gamma_os_n)
          # Gamma_var_est <- (1/n_orig^2) * sum((
          #   weights*infl_fn_Gamma(rep(x,length(dat$a)),dat$w,dat$y_star,
          #                         dat$delta_star,dat$a)
          # )^2)
          # Gamma_est <- Gamma_os_n(x)
          # test_stat <- (Gamma_est-Gamma_0)^2/Gamma_var_est
          # p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
          
        } else if (L$temp_which=="Psi_1") {
          
          Theta_true <- attr(dat_orig,"Theta_true")
          Gamma_0 <- Vectorize(function(x) {
            Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
          })
          infl_fn_1 <- construct_infl_fn_1(dat, Gamma_0, Phi_n,
                                           lambda_2, lambda_3)
          Psi_1_var_est <- (1/n_orig^2) * sum((weights*infl_fn_1(dat$a))^2)
          Psi_1_est <- (1/n_orig) * sum(weights*(
            lambda_2*(Phi_n(dat$a))^2*Gamma_0(dat$a) -
              lambda_3*Phi_n(dat$a)*Gamma_0(dat$a)
          ))
          test_stat <- Psi_1_est^2/Psi_1_var_est
          p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
          
        } else if (L$temp_which=="Psi_2") {
          
          # Phi_0 <- function(x) {x}
          # infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
          #                                          eta_n, Gamma_os_n)
          # infl_fn_2 <- construct_infl_fn_2(dat, Phi_0, infl_fn_Gamma, 1/3, 1/4)
          # Psi_2_var_est <- (1/n_orig^2) * sum((
          #   weights*infl_fn_2(dat$w,dat$y_star,dat$delta_star,dat$a)
          # )^2)
          # a_mc <- runif(10^6)
          # Psi_2_est <- mean(
          #   (1/3)*(Phi_0(a_mc))^2*Gamma_os_n(round(a_mc,-log10(C$appx$a))) -
          #     (1/4)*Phi_0(a_mc)*Gamma_os_n(round(a_mc,-log10(C$appx$a)))
          # )
          # test_stat <- Psi_2_est^2/Psi_2_var_est
          # p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
          
        } else if (L$temp_which=="Psi_1+Psi_2") {
          
          # # Psi_1
          # Theta_true <- attr(dat_orig,"Theta_true")
          # Gamma_0 <- Vectorize(function(x) {
          #   Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
          # })
          # infl_fn_1 <- construct_infl_fn_1(dat, Gamma_0, Phi_n,
          #                                  lambda_2, lambda_3)
          # Psi_1_est <- (1/n_orig) * sum(weights*(
          #   lambda_2*(Phi_n(dat$a))^2*Gamma_0(dat$a) -
          #     lambda_3*Phi_n(dat$a)*Gamma_0(dat$a)
          # ))
          # 
          # # Psi_2
          # Phi_0 <- function(x) {x}
          # infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
          #                                          eta_n, Gamma_os_n)
          # infl_fn_2 <- construct_infl_fn_2(dat, Phi_0, infl_fn_Gamma, 1/3, 1/4)
          # a_mc <- runif(10^6)
          # Psi_2_est <- mean(
          #   (1/3)*(Phi_0(a_mc))^2*Gamma_os_n(round(a_mc,-log10(C$appx$a))) -
          #     (1/4)*Phi_0(a_mc)*Gamma_os_n(round(a_mc,-log10(C$appx$a)))
          # )
          # 
          # # Combined
          # infl_fn_beta <- function(a,w,y_star,delta_star) {
          #   infl_fn_1(a) + infl_fn_2(w,y_star,delta_star,a)
          # }
          # beta_var_est <- (1/n_orig^2) * sum((
          #   weights*infl_fn_beta(dat$a,dat$w,dat$y_star,dat$delta_star)
          # )^2)
          # test_stat <- (Psi_1_est+Psi_2_est)^2/beta_var_est
          # p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
          
        } else if (L$temp_which=="beta_n") {
          
          # infl_fn_beta <- function(a,w,y_star,delta_star) {
          #   infl_fn_1(a) + infl_fn_2(w,y_star,delta_star,a)
          # }
          # beta_var_est <- (1/n_orig^2) * sum((
          #   weights*infl_fn_beta(dat$a,dat$w,dat$y_star,dat$delta_star)
          # )^2)
          # test_stat <- beta_n^2/beta_var_est
          # p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
          
      }
        
        # Calculate P-value
        if (alt_type=="incr") {
          p_val <- pnorm(beta_n, mean=0, sd=sd_n, lower.tail=FALSE)
        } else if (alt_type=="decr") {
          p_val <- pnorm(beta_n, mean=0, sd=sd_n)
        } else if (alt_type=="two-tailed") {
          test_stat_chisq <- beta_n^2/var_n
          p_val <- pchisq(test_stat_chisq, df=1, lower.tail=FALSE)
        }
        
      }
      
    }
    
  }
  
  if (params$var=="boot") {
    
    # # !!!!! Update all of this if needed
    # 
    # # Define the statistic to bootstrap
    # bootstat <- function(dat_orig,indices) {
    #   
    #   # !!!!! Change dat to dat_orig_b
    #   
    #   dat <- dat_orig[indices,]
    #   n_orig <- length(dat$delta)
    #   dat_0 <- ss(dat_orig, which(dat_orig$delta==1))
    #   n_0 <- length(dat_0$a)
    #   weights_0 <- wts(dat_0) # !!!!!
    #   G_0 <- construct_Phi_n(dat, type=params$ecdf_type)
    #   f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type, k=15)
    #   f_a_n <- construct_f_a_n(dat_orig, f_aIw_n=f_aIw_n)
    #   g_0 <- construct_g_n(f_aIw_n, f_a_n)
    #   S_0 <- construct_S_n(dat, type=params$S_n_type)
    #   Sc_0 <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    #   omega_0 <- construct_omega_n(S_0, Sc_0)
    #   Gamma_0 <- construct_Gamma_os_n(dat, vlist$A_grid, omega_0, S_0, g_0)
    #   
    #   beta_0 <- (1/n_orig) * sum(
    #     weights_0 *
    #       (
    #         lambda(dat_orig,2,G_0)*(G_0(dat_0$a))^2 -
    #           lambda(dat_orig,3,G_0)*G_0(dat_0$a)
    #       ) *
    #       (Gamma_0(dat_0$a))
    #   )
    #   
    #   return (beta_0)
    #   
    # }
    # 
    # # Run bootstrap
    # boot_obj <- boot(data=dat, statistic=bootstat, R=params$boot_reps)
    # 
    # # Calculate beta_n
    # beta_0 <- bootstat(dat, c(1:length(dat$delta)))
    # 
    # # Calculate critical value (for a one-sided test)
    # crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  if (params$var=="mixed boot") {
    
    # # !!!!! Update all of this if needed
    # 
    # # Pre-calculate non-bootstrapped pieces
    # {
    #   dat_0_orig <- dat
    #   n_orig <- length(dat_0_orig$delta)
    #   dat_0 <- ss(dat_0_orig, which(dat_0_orig$delta==1))
    #   n_0 <- length(dat_0$a)
    #   weights_0 <- wts(dat_0) # !!!!!
    #   
    #   G_0 <- construct_Phi_n(dat_0, type=params$ecdf_type)
    #   S_0 <- construct_S_n(dat_0, type=params$S_n_type)
    #   Sc_0 <- construct_S_n(dat_0, type=params$S_n_type, csf=TRUE)
    #   omega_0 <- construct_omega_n(S_0, Sc_0)
    #   f_aIw_n <- construct_f_aIw_n(dat_0, type=params$g_n_type, k=15)
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
    #   weights_b <- wts(dat_b) # !!!!!
    #   Phi_n <- construct_Phi_n(dat_b, type=params$ecdf_type)
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
    # boot_obj <- boot(data=dat_0_orig, statistic=bootstat, R=params$boot_reps)
    # 
    # # Calculate critical value (for a one-sided test)
    # crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    
  }
  
  return(list(
    reject = as.integer(p_val<0.05),
    p_val = p_val,
    beta_n = beta_n,
    sd_n = sd_n
    # Gamma_n_5 = Gamma_os_n(0.5), # !!!!!
    # Gamma_var_n = Gamma_var_n, # !!!!!
    # lambda_2 = lambda_2_, # !!!!!
    # lambda_3 = lambda_3_, # !!!!!
    # g_n = g_n(0.5,c(1,1)), # !!!!!
    # S_n = S_n(100,c(1,1),0.5), # !!!!!
    # Sc_n = Sc_n(100,c(1,1),0.5), # !!!!!
    # omega_n1 = omega_n(c(1,1),0.5,80,1), # !!!!!
    # omega_n0 = omega_n(c(1,1),0.5,80,0), # !!!!!
    # gcomp_n = gcomp_n(0.5), # !!!!!
    # eta_n = eta_n(0.5,c(1,1)), # !!!!!
    # xi_n = xi_n(0.5,0.5), # !!!!!
    # infl_fn_1 = infl_fn_1(0.5), # !!!!!
    # infl_fn_Gamma = infl_fn_Gamma(0.5,c(1,1),80,1,0.5), # !!!!!
    # infl_fn_2 = infl_fn_2(c(1,1),80,1,0.5), # !!!!!
    # partial_est = partial_est, # !!!!!
    # partial_var = partial_var, # !!!!!
    # Phi_n_5 = Phi_n(0.5) # !!!!!
  ))
  
}
