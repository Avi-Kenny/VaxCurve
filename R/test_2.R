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
test_2 <- function(dat_orig, alt_type="two-tailed", params,
                   test_stat_only=F, return_extras=F) {
  
  # Set default params
  chk(0, "Testing: START")
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
  chk(1)
  srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type, print_coeffs=T)
  Q_n <- srvSL$srv
  Qc_n <- srvSL$cens
  chk(2)
  omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
  etastar_n <- construct_etastar_n(Q_n, vals=NA)
  chk(3)
  f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                               k=p$f_sIx_n_bins, edge_corr="none",
                               s_scale=s_scale, s_shift=s_shift)
  f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
  chk(4)
  q_tilde_n <- construct_q_tilde_n(type=p$q_n_type, f_n_srv, f_sIx_n,
                                   omega_n)
  chk(5)
  Theta_os_n <- construct_Theta_os_n(dat, dat_orig, omega_n, f_sIx_n,
                                     q_tilde_n, etastar_n)
  chk(6)
  infl_fn_Theta <- construct_infl_fn_Theta(omega_n, f_sIx_n, q_tilde_n,
                                           etastar_n, Theta_os_n)
  chk(7)

  # Function to compute P values
  compute_p_val <- function(alt_type, beta_n, var_n) {
    if (alt_type=="incr") {
      p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n))
    } else if (alt_type=="decr") {
      p_val <- pnorm(beta_n, mean=0, sd=sqrt(var_n), lower.tail=FALSE)
    } else if (alt_type=="two-tailed") {
      p_val <- pchisq(beta_n^2/var_n, df=1, lower.tail=FALSE)
    }
    return(p_val)
  }
  
  # Set up container to hold results
  res <- list()
  
  if ("debug" %in% p$type) {
    
    # Compute test statistic and variance estimate
    beta_n <- Theta_os_n(0.5) - 0.495
    var_n <- 0
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      x_i <- as.numeric(dat_orig$x[i,])
      var_n <- var_n + (
        infl_fn_Theta(u=0.5, x_i, y_i, delta_i, s_i, weight_i)
      )^2
    }
    var_n <- var_n/(n_orig^2)
    
    res[[length(res)+1]] <- list(
      type = "debug",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    
  }
  
  if ("simple" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    chk(8, "simple: START")
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

    # !!!!! Temporary
    beta_n <- -1 * beta_n
    
    res[[length(res)+1]] <- list(
      type = "simple",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    chk(9, "simple: END")
    
  }
  
  if ("simple (with constant)" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    chk(10, "simple (with constant): START")
    u_mc <- round(seq(C$appx$s,1,C$appx$s),-log10(C$appx$s))
    m <- length(u_mc)
    lambda_1 <- mean(u_mc) # ~1/2
    lambda_2 <- mean((u_mc)^2) # ~1/3
    lambda_3 <- mean((u_mc)^3) # ~1/4
    
    # Compute test statistic and variance estimate
    beta_n <- mean((
      (lambda_3-lambda_1*lambda_2)*(u_mc-lambda_1) +
        (lambda_1^2-lambda_2)*(u_mc^2-lambda_2)
    ) * Theta_os_n(u_mc))
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
        infl_fn_Theta(u=u_mc, x_m, y_m, delta_m, s_m, weight_m) * (
          (lambda_3-lambda_1*lambda_2)*(u_mc-lambda_1) +
          (lambda_1^2-lambda_2)*(u_mc^2-lambda_2)
        )
      ))^2
    }
    var_n <- var_n/(n_orig^2*m^2)
    
    res[[length(res)+1]] <- list(
      type = "simple (with constant)",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    chk(11, "simple (with constant): END")
    
  }
  
  if ("S-weighted (with constant)" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    chk(12, "S-weighted (with constant): START")
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
      # (1/n_orig) * sum(unlist(lapply(dat$s, function(s_j) {
      #   (piece_1*(s_j-lmd_1)+piece_2*(s_j^2-lmd_2)) * 
      #     infl_fn_Theta(u=s_j, x_i, y_i, delta_i, s_i, weight_i)
      # })))
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
    chk(13, "S-weighted (with constant): END")
    
  }
  
  if ("edge" %in% p$type) {
    
    # Construct pieces needed for hypothesis test
    g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
    r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    infl_fn_r_Mn_edge <- function(weight, s, x, y, delta) { # !!!!! Move this later
      if (weight==0) {
        return(0)
      } else {
        val <- weight * (
          Q_n(C$t_0,x,s=0) - (
            (In(s==0)/g_sn(x, 0)) * omega_n(x,s=0,y,delta)
          ) -
            (1-r_Mn_edge_est)
        )
        return(val)
      }
    }
    
    # Compute test statistic and variance estimate
    beta_n <- Theta_os_n(1) - r_Mn_edge_est
    var_n <- 0
    for (i in c(1:n_orig)) {
      s_i <- dat_orig$s[i]
      x_i <- as.numeric(dat_orig$x[i,])
      y_i <- dat_orig$y[i]
      delta_i <- dat_orig$delta[i]
      weight_i <- dat_orig$weight[i]
      var_n <- var_n + (
        infl_fn_Theta(u=1, x_i, y_i, delta_i, s_i, weight_i) -
          infl_fn_r_Mn_edge(weight_i, s_i, x_i, y_i, delta_i)
      )^2
    }
    var_n <- var_n/(n_orig^2)
    
    res[[length(res)+1]] <- list(
      type = "edge",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    
  }
  
  if ("combined" %in% p$type) {
    
    # # Construct pieces needed for hypothesis test
    # g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
    # r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
    # infl_fn_r_Mn_edge <- function(weight, s, x, y, delta) { # !!!!! Move this later
    #   if (weight==0) {
    #     return(0)
    #   } else {
    #     val <- weight * (
    #       Q_n(C$t_0,x,s=0) - (
    #         (In(s==0)/g_sn(x, 0)) * omega_n(x,s=0,y,delta)
    #       ) -
    #         (1-r_Mn_edge_est)
    #     )
    #     return(val)
    #   }
    # }
    
    # # Compute test statistic and variance estimate
    # beta_n <- Theta_os_n(1) - r_Mn_edge_est
    # var_n <- 0
    # for (i in c(1:n_orig)) {
    #   s_i <- dat_orig$s[i]
    #   x_i <- as.numeric(dat_orig$x[i,])
    #   y_i <- dat_orig$y[i]
    #   delta_i <- dat_orig$delta[i]
    #   weight_i <- dat_orig$weight[i]
    #   var_n <- var_n + (
    #     infl_fn_Theta(u=1, x_i, y_i, delta_i, s_i, weight_i) -
    #       infl_fn_r_Mn_edge(weight_i, s_i, x_i, y_i, delta_i)
    #   )^2
    # }
    # var_n <- var_n/(n_orig^2)
    
    res[[length(res)+1]] <- list(
      type = "combined",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    
  }
  
  if ("complex" %in% p$type) {
    
    # Compute test statistic and variance estimate
    lambda_2n <- (1/n_orig) * sum(dat$weights*dat$s^2)
    lambda_3n <- (1/n_orig) * sum(dat$weights*dat$s^3)
    beta_n <- (1/n_orig) * sum(dat$weights*(
      (lambda_2n*dat$s^2-lambda_3n*dat$s) * Theta_os_n(dat$s)
    ))
    
    # Compute variance estimate
    xi_1n <- (1/n_orig)*sum(dat$weights*dat$s*Theta_os_n(dat$s))
    xi_2n <- (1/n_orig)*sum(dat$weights*dat$s^2*Theta_os_n(dat$s))
    piece_2 <- 2*(lambda_3n*xi_1n-lambda_2n*xi_2n)
    var_n <- 0
    n_dat <- length(dat$z)
    dat_s <- dat$s
    dat_s2 <- dat_s^2
    for (i in c(1:n_orig)) {
      x_i <- as.data.frame(
        matrix(rep(dat_orig$x[i,],n_dat), ncol=length(dat_orig$x[i,]), byrow=T)
      )
      y_i <- rep(dat_orig$y[i],n_dat)
      delta_i <- rep(dat_orig$delta[i],n_dat)
      s_i <- rep(dat_orig$s[i],n_dat)
      wt_i <- dat_orig$weights[i]
      if (wt_i==0) {
        piece_1 <- 0
      } else {
        s_i <- dat_orig$s[i]
        piece_1 <- wt_i * ( xi_2n*s_i^2 - xi_1n*s_i^3 +
                             (lambda_2n*s_i^2-lambda_3n*s_i)*Theta_os_n(s_i) )
      }
      piece_3 <- (1/n_orig) * sum(
        dat$weights * (lambda_2n*dat_s2-lambda_3n*dat_s) *
          infl_fn_Theta(dat_s,x_i,y_i,delta_i,s_i,wt_i)
      )
      var_n <- var_n + (piece_1+piece_2+piece_3)^2
    }
    var_n <- var_n/(n_orig^2)
    
    # !!!!! Temporary
    beta_n <- -1 * beta_n
    
    res[[length(res)+1]] <- list(
      type = "complex",
      p_val = compute_p_val(alt_type, beta_n, var_n),
      beta_n = beta_n,
      var_n = var_n
    )
    
  }
  
  if ("asymptotic, Gamma_n" %in% p$type) {} # Archived for now
  if ("boot" %in% p$type) {} # !!!!! Archived for now
  if ("mixed boot" %in% p$type) {} # !!!!! Archived for now
  
  # Return debugging components (var="Monte Carlo")
  if (return_extras) { res$Theta_os_n <- Theta_os_n }
  
  if (T) {
    res$extras <- list(
      Theta_0.1 = Theta_os_n(0.1),
      Theta_0.4 = Theta_os_n(0.4),
      Theta_0.8 = Theta_os_n(0.8),
      etastar_0.1 = mean(sapply(c(1:n_orig), function(i) {
        etastar_n(u=0.1, x=as.numeric(dat_orig$x[i,]))
      })),
      etastar_0.4 = mean(sapply(c(1:n_orig), function(i) {
        etastar_n(u=0.4, x=as.numeric(dat_orig$x[i,]))
      })),
      etastar_0.8 = mean(sapply(c(1:n_orig), function(i) {
        etastar_n(u=0.8, x=as.numeric(dat_orig$x[i,]))
      }))
    )
  } # DEBUG: return debugging components
  
  return(res)
  
}
