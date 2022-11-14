#' Testing approach 2: regression slope
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param alt_type Type of alternative hypothesis; either "incr", "decr", or
#'     "two-tailed"; currently unused and defaults to "two-tailed".
#' @param params A list, containing the following:
#'   - `Q_n_type` Q_n_type Type of survival function estimator; currently only
#'     c("Cox PH")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `var` Variance estimation; one of c("boot","mixed boot")
#'   - `boot_reps` Number of bootstrap replicates to run
#' @param test_stat_only Boolean; if TRUE, only calculate the test statistic and
#'     do not estimate its variance
#' @return Binary; is null rejected (1) or not (0)
test_2 <- function(dat_orig, alt_type="two-tailed", params, test_stat_only=F) {
  
  # ..p <- list(coeffs_srv=c(NA,NA,NA), coeffs_cens=c(NA,NA,NA)) # !!!!!
  
  # Set default params
  .default_params <- list(
    type="simple (with constant)", ecdf_type="step", g_n_type="binning",
    Q_n_type="Super Learner", omega_n_type="estimated", q_n_type="new",
    cf_folds=1, f_sIx_n_bins=15, boot_reps=200
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
  etastar_n <- construct_etastar_n(Q_n, vals=NA)
  f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                               k=p$f_sIx_n_bins, s_scale=s_scale,
                               s_shift=s_shift)
  f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
  q_tilde_n <- construct_q_tilde_n(type=p$q_n_type, f_n_srv, f_sIx_n,
                                   omega_n)
  Theta_os_n <- construct_Theta_os_n(dat, dat_orig, omega_n, f_sIx_n,
                                     q_tilde_n, etastar_n)
  infl_fn_Theta <- construct_infl_fn_Theta(omega_n, f_sIx_n, q_tilde_n,
                                           etastar_n, Theta_os_n)
  
  # Function to compute P values
  compute_p_val <- function(alt_type, beta_n, var_n) {
    if (alt_type=="incr") {
      p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n), lower.tail=FALSE)
    } else if (alt_type=="decr") {
      p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n))
    } else if (alt_type=="two-tailed") {
      p_val <- pchisq(beta_n^2/var_n, df=1, lower.tail=FALSE)
    }
    return(p_val)
  }
  
  # Set up container to hold results
  res <- list()
  
  if ("simple" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    m <- length(u_mc)
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    
    # Compute test statistic and variance estimate
    beta_n <- mean((lambda_2*u_mc^2-lambda_3*u_mc)*Theta_os_n(u_mc))
    var_n <- 0
    for (i in c(1:n_orig)) {
      s_m <- rep(dat_orig$s[i],m)
      y_m <- rep(dat_orig$y[i],m)
      delta_m <- rep(dat_orig$delta[i],m)
      weight_m <- rep(dat_orig$weight[i],m)
      x_m <- as.data.frame(
        matrix(rep(dat_orig$x[i,],m), ncol=length(dat_orig$x[i,]), byrow=T)
      )
      var_n <- var_n + (sum(
        infl_fn_Theta(u=u_mc, x_m, y_m, delta_m, s_m, weight_m) *
          (lambda_2*u_mc^2-lambda_3*u_mc)
      ))^2
    }
    var_n <- var_n/(n_orig^2*m^2)
    
    res[[length(res)+1]] <- list(
      type = "simple",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    
  }
  
  if ("simple (with constant)" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    infl_fn_beta_n <- construct_infl_fn_beta_n(infl_fn_Theta)
    
    # Compute test statistic and variance estimate
    u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    m <- length(u_mc)
    lambda_1 <- mean(u_mc) # ~1/2
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    
    beta_n <- mean((
      (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
        (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
    ) * Theta_os_n(u_mc))
    var_n <- 0
    for (i in c(1:n_orig)) {
      var_n <- var_n + (
        infl_fn_beta_n(dat_orig$s[i], dat_orig$y[i], dat_orig$delta[i],
                       dat_orig$weight[i], as.numeric(dat_orig$x[i,]))
      )^2
    }
    var_n <- var_n/n_orig^2
    
    res[[length(res)+1]] <- list(
      type = "simple (with constant)",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )

  }
  
  if ("S-weighted (with constant)" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    s <- dat$s
    lmd_1 <- (1/n_orig) * sum(dat$weights*s)
    lmd_2 <- (1/n_orig) * sum(dat$weights*s^2)
    lmd_3 <- (1/n_orig) * sum(dat$weights*s^3)
    xi_0 <- (1/n_orig) * sum(dat$weights*Theta_os_n(s))
    xi_1 <- (1/n_orig) * sum(dat$weights*s*Theta_os_n(s))
    xi_2 <- (1/n_orig) * sum(dat$weights*s^2*Theta_os_n(s))
    h1 <- 2*lmd_1*xi_2 - lmd_3*xi_0 - lmd_2*xi_1
    h2 <- 2*lmd_2*xi_0 - lmd_1*xi_1 - xi_2
    h3 <- xi_1 - lmd_1*xi_0
    h4 <- lmd_2^2 - lmd_1*lmd_3
    h5 <- lmd_3 - lmd_1*lmd_2
    h6 <- lmd_1^2 - lmd_2
    
    # Compute test statistic
    piece_1 <- lmd_3-lmd_1*lmd_2
    piece_2 <- lmd_1^2-lmd_2
    beta_n <- (1/n_orig) * sum(
      dat$weights * (piece_1*(s-lmd_1)+piece_2*(s^2-lmd_2)) * Theta_os_n(s)
    )

    # Influence function psi_1 (!!!!! move this)
    const <- 3*(lmd_1*lmd_3-lmd_2^2)*xi_0 +
      (3*lmd_1*lmd_2-2*lmd_3)*xi_1 +
      (2*lmd_2-3*lmd_1^2)*xi_2
    infl_fn_psi_1 <- function(s_i, weight_i) {
      if (weight_i==0) {
        return(const)
      } else {
        return(
          weight_i *
            (h1*s_i+h2*s_i^2+h3*s_i^3+(h4+h5*s_i+h6*s_i^2)*Theta_os_n(s_i)) +
            const
        )
      }
    }
    infl_fn_psi_1 <- construct_superfunc(infl_fn_psi_1, vec=c(1,1))
    
    # Influence function psi_2 (!!!!! move this)
    infl_fn_psi_2 <- function(x_i, y_i, delta_i, s_i, weight_i) {
      (1/n_orig) * sum(unlist(lapply(c(1:length(dat$s)), function(j) {
        s_j <- dat$s[j]
        dat$weights[j] * (piece_1*(s_j-lmd_1)+piece_2*(s_j^2-lmd_2)) * 
          infl_fn_Theta(u=s_j, x_i, y_i, delta_i, s_i, weight_i)
      })))
    }
    infl_fn_psi_2 <- construct_superfunc(infl_fn_psi_2, vec=c(2,1,1,1,1))
    
    # Compute variance estimate
    var_n <- 0
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      x_i <- as.numeric(dat_orig$x[i,])
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      var_n <- var_n + (
        infl_fn_psi_1(s_i, weight_i) +
          infl_fn_psi_2(x_i, y_i, delta_i, s_i, weight_i)
      )^2
    }
    var_n <- var_n/(n_orig^2)

    res[[length(res)+1]] <- list(
      type = "S-weighted (with constant)",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )

  }
  
  if ("edge" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n,
                                                     r_Mn_edge_est, val=0)
    infl_fn_beta_en <- construct_infl_fn_beta_en(infl_fn_Theta,
                                                 infl_fn_r_Mn_edge)
    
    # Compute test statistic and variance estimate
    beta_en <- Theta_os_n(1) - r_Mn_edge_est
    var_n <- 0
    var_edge <- 0 # !!!!! TEMP
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      x_i <- as.numeric(dat_orig$x[i,])
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      var_n <- var_n + (infl_fn_beta_en(x_i, y_i, delta_i, s_i, weight_i))^2
      var_edge <- var_edge + (infl_fn_r_Mn_edge(weight_i, s_i, x_i, y_i, delta_i))^2 # !!!!! TEMP
    }
    var_n <- var_n/(n_orig^2)
    var_edge <- var_edge/(n_orig^2) # !!!!! TEMP
    
    res[[length(res)+1]] <- list(
      type = "edge",
      p_val = compute_p_val(alt_type, beta_en, var_n),
      beta_n = beta_en,
      var_n = var_n
    )

  }
  
  if ("combined" %in% p$type) {
    
    # Construct pieces needed for beta_n
    u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    m <- length(u_mc)
    lambda_1 <- mean(u_mc) # ~1/2
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    beta_n <- mean((
      (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
        (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
    ) * Theta_os_n(u_mc))
    infl_fn_beta_n <- construct_infl_fn_beta_n(infl_fn_Theta)

    # Construct pieces needed for beta_en
    g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n,
                                                     r_Mn_edge_est, val=0)
    infl_fn_beta_en <- construct_infl_fn_beta_en(infl_fn_Theta,
                                                 infl_fn_r_Mn_edge)
    beta_en <- Theta_os_n(1) - r_Mn_edge_est

    # Calculate variance components
    sigma2_bn <- 0
    sigma2_ben <- 0
    cov_n <- 0
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      x_i <- as.numeric(dat_orig$x[i,])
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      
      if_bn <- infl_fn_beta_n(s_i, y_i, delta_i, weight_i, x_i)
      if_ben <- infl_fn_beta_en(x_i, y_i, delta_i, s_i, weight_i)
      
      sigma2_bn <- sigma2_bn + if_bn^2
      sigma2_ben <- sigma2_ben + if_ben^2
      cov_n <- cov_n + if_bn*if_ben
    }
    sigma_bn <- sqrt(sigma2_bn/n_orig)
    sigma_ben <- sqrt(sigma2_ben/n_orig)
    cov_n <- cov_n/n_orig
    rho_n <- cov_n/(sigma_bn*sigma_ben)
    
    # Calculate combined test statistic
    beta_star_n <- n_orig * (
      beta_n^2/sigma_bn^2 +
        (sigma_bn*beta_en-rho_n*sigma_ben*beta_n)^2 /
        (sigma_bn^2*sigma_ben^2*(1-rho_n^2))
    )
    
    res[[length(res)+1]] <- list(
      type = "combined",
      p_val = pchisq(beta_star_n, df=2, lower.tail=FALSE),
      beta_n = beta_star_n,
      var_n = 999
    )
    
  }
  
  if ("combined 2" %in% p$type) {
    
    # Construct pieces needed for beta_n
    u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    m <- length(u_mc)
    lambda_1 <- mean(u_mc) # ~1/2
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    beta_n <- mean((
      (lambda_1*lambda_2-lambda_3)*(u_mc-lambda_1) +
        (lambda_2-lambda_1^2)*(u_mc^2-lambda_2)
    ) * Theta_os_n(u_mc))
    infl_fn_beta_n <- construct_infl_fn_beta_n(infl_fn_Theta)
    
    # Construct pieces needed for beta_en
    g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    infl_fn_r_Mn_edge <- construct_infl_fn_r_Mn_edge(Q_n, g_sn, omega_n,
                                                     r_Mn_edge_est, val=0)
    infl_fn_beta_en <- construct_infl_fn_beta_en(infl_fn_Theta,
                                                 infl_fn_r_Mn_edge)
    beta_en <- Theta_os_n(1) - r_Mn_edge_est
    
    # Calculate variance components
    sigma2_bn <- 0
    sigma2_ben <- 0
    cov_n <- 0
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      x_i <- as.numeric(dat_orig$x[i,])
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      
      if_bn <- infl_fn_beta_n(s_i, y_i, delta_i, weight_i, x_i)
      if_ben <- infl_fn_beta_en(x_i, y_i, delta_i, s_i, weight_i)
      
      sigma2_bn <- sigma2_bn + if_bn^2
      sigma2_ben <- sigma2_ben + if_ben^2
      cov_n <- cov_n + if_bn*if_ben
    }
    sigma_bn <- sqrt(sigma2_bn/n_orig)
    sigma_ben <- sqrt(sigma2_ben/n_orig)
    cov_n <- cov_n/n_orig
    rho_n <- cov_n/(sigma_bn*sigma_ben)
    
    # Calculate combined test statistic
    beta_star_n <- sqrt(n_orig/(2+2*rho_n)) *
      (beta_n/sigma_bn + beta_en/sigma_ben)
    
    res[[length(res)+1]] <- list(
      type = "combined 2",
      p_val = compute_p_val(alt_type, beta_star_n, 1),
      beta_n = beta_star_n,
      var_n = 1
    )
    
  }
  
  if ("complex" %in% p$type) {} # Archived for now
  if ("asymptotic, Gamma_n" %in% p$type) {} # Archived for now
  if ("boot" %in% p$type) {} # !!!!! Archived for now
  if ("mixed boot" %in% p$type) {} # !!!!! Archived for now
  if ("debug" %in% p$type) {} # !!!!! Archived for now
  
  if (T) {
    res$extras <- list(
      Theta_1.0 = Theta_os_n(1),
      r_Mn_0.0 = r_Mn_edge_est,
      var_edge = var_edge,
      sd_edge = sqrt(var_edge),
      beta_n = beta_n,
      beta_en = beta_en,
      sigma_bn = sigma_bn,
      sigma_ben = sigma_ben,
      rho_n = rho_n
    )
    
    # res$extras <- list(
    #   Theta_0.1 = Theta_os_n(0.1),
    #   Theta_0.4 = Theta_os_n(0.4),
    #   Theta_0.8 = Theta_os_n(0.8),
    #   etastar_0.1 = mean(sapply(c(1:n_orig), function(i) {
    #     etastar_n(u=0.1, x=as.numeric(dat_orig$x[i,]))
    #   })),
    #   etastar_0.4 = mean(sapply(c(1:n_orig), function(i) {
    #     etastar_n(u=0.4, x=as.numeric(dat_orig$x[i,]))
    #   })),
    #   etastar_0.8 = mean(sapply(c(1:n_orig), function(i) {
    #     etastar_n(u=0.8, x=as.numeric(dat_orig$x[i,]))
    #   }))
    # )
  } # DEBUG: return debugging components
  
  return(res)
  
}
