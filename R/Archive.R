
# fns_doseresp: construct_rho_n
if (F) {
  
  construct_rho_n <- function(dat, Phi_n, vals=NA) {
    
    n_orig <- sum(dat$weights)
    
    fnc <- function(s) {
      (1/n_orig) * sum(
        dat$weights * (Phi_n(dat$s)^3) * (In(s<=dat$s) - Phi_n(dat$s))
      )
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
    
  }
  
}

# fns_doseresp: construct_xi_n
if (F) {
  
  construct_xi_n <- function(Phi_n, lambda_2, lambda_3, vals=NA) {
    
    fnc <- function(s_i,s_j) {
      (2*In(s_i<=s_j) - 3*Phi_n(s_j))*Phi_n(s_j)*lambda_2 +
        (2*Phi_n(s_j) - In(s_i<=s_j))*lambda_3
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=c(1,1), vals=vals))
    
  }

}

# fns_doseresp: beta_n_var_hat
if (F) {
  
  beta_n_var_hat <- function(dat, dat_orig, infl_fn_1, infl_fn_2) {
    
    n_orig <- sum(dat$weights)
    
    b_sum <- 0
    for (i in c(1:length(dat_orig$s))) {
      if (dat_orig$z[i]==0) {
        b_sum <- b_sum + (infl_fn_2(dat_orig$x[i,], dat_orig$y[i],
                                    dat_orig$delta[i], dat_orig$s[i],
                                    dat_orig$weights[i]))^2
      } else {
        b_sum <- b_sum + (
          dat_orig$weights[i]*infl_fn_1(dat_orig$s[i]) +
            infl_fn_2(dat_orig$x[i,], dat_orig$y[i], dat_orig$delta[i],
                      dat_orig$s[i], dat_orig$weights[i])
        )^2
      }
    }
    
    return( (1/n_orig)*b_sum )
    
  }
  
}

# fns_doseresp: construct_infl_fn_2
if (F) {
  
  construct_infl_fn_2 <- function(dat, Phi_n, infl_fn_Gamma, lambda_2, lambda_3) {
    
    n_orig <- sum(dat$weights)
    weights_j <- dat$weights
    s_j <- dat$s
    len <- length(s_j)
    
    fnc <- function(x,y,delta,s,wt) {
      x_long <- as.data.frame(
        matrix(rep(x,len), ncol=length(x), byrow=T)
      )
      (1/n_orig) * sum(weights_j * (
        ( lambda_2*(Phi_n(s_j)^2) - lambda_3*Phi_n(s_j) ) *
          infl_fn_Gamma(s_j, x_long, rep(y,len), rep(delta,len),
                        rep(s,len), rep(wt,len))
      ))
    }
    
    return(construct_superfunc(fnc, vec=c(2,1,1,1,1)))
    
  }
  
}

# fns_doseresp: construct_infl_fn_Gamma
if (F) {
  
  construct_infl_fn_Gamma <- function(omega_n, g_n, r_tilde_Mn, p_n, Gamma_tilde_n,
                                      q_n, eta_n, Gamma_os_n) {
    
    fnc <- function(u,x,y,delta,s,wt) {
      if (wt==0) {
        piece_1 <- 0
        piece_2 <- 0
        piece_3 <- 0
      } else {
        piece_1 <- In(s!=0 & s<=u)
        piece_2 <- omega_n(x,s,y,delta)/(p_n*g_n(s,x)) + r_tilde_Mn(s)/p_n
        piece_3 <- In(s!=0)*Gamma_tilde_n(u)
      }
      wt*(piece_1*piece_2-piece_3/p_n) +
        (1-wt)*q_n(x,y,delta,u) +
        eta_n(u,x) -
        Gamma_os_n(round(u,-log10(C$appx$s)))
    }
    
    return(construct_superfunc(fnc, vec=c(1,2,1,1,1,1)))
    
  }
  
}

# fns_doseresp: construct_infl_fn_1
if (F) {
  
  construct_infl_fn_1 <- function(dat, Gamma_os_n, Phi_n, lambda_2,
                                  lambda_3, vals=NA) {
    
    n_orig <- sum(dat$weights)
    weights_j <- dat$weights
    s_j <- dat$s
    
    rho_n <- function(s,u) {
      (1/n_orig) * sum(
        weights_j * (Phi_n(s))^u * Gamma_os_n(round(s,-log10(C$appx$s)))
      )
    }
    
    rho_1 <- rho_n(s_j,1)
    rho_2 <- rho_n(s_j,2)
    piece_01 <- Gamma_os_n(round(s_j,-log10(C$appx$s)))
    piece_20 <- (Phi_n(s_j))^2
    piece_10 <- Phi_n(s_j)
    piece_11 <- Phi_n(s_j) * Gamma_os_n(round(s_j,-log10(C$appx$s)))
    
    fnc <- function(s_i) {
      
      piece_1 <- (2*(1/n_orig)*sum(weights_j*In(s_i<=s_j)*piece_10)+
                    (Phi_n(s_i))^2-6*lambda_2)*rho_2
      piece_2 <- lambda_2*(
        2*(1/n_orig)*sum(weights_j*In(s_i<=s_j)*piece_11)+
          (Phi_n(s_i))^2*Gamma_os_n(round(s_i,-log10(C$appx$s)))
      )
      piece_3 <- (3*(1/n_orig)*sum(weights_j*In(s_i<=s_j)*piece_20)+
                    (Phi_n(s_i))^3-6*lambda_3)*rho_1
      piece_4 <- lambda_3*((1/n_orig)*sum(weights_j*In(s_i<=s_j)*
                                            piece_01)+Phi_n(s_i)*
                             Gamma_os_n(round(s_i,-log10(C$appx$s))))
      
      return(piece_1+piece_2-piece_3-piece_4)
      
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=vals))
    
  }

}

# test2.R: "debugging" test statistic (test based on Theta_n alone)
if (F) {
  
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

# test2.R: "complex" test statistic
if (F) {
  
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
  
}

# test2.R: Test statistic based on Gamma_n
if (F) {
  
  # !!!!! Update all of this if needed
  
  # Prep
  n_orig <- length(dat_orig$z)
  dat <- ss(dat_orig, which(dat_orig$z==1))
  weights <- dat$weights
  
  # Construct dataframes of values to pre-compute functions on
  vlist <- create_val_list(dat_orig)
  
  # Construct component functions
  Phi_n <- construct_Phi_n(dat, type=p$ecdf_type)
  lambda_2 <- lambda(dat,2,Phi_n)
  lambda_3 <- lambda(dat,3,Phi_n)
  f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                               k=p$f_sIx_n_bins)
  f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
  g_n <- construct_g_n(f_sIx_n, f_s_n)
  srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type)
  Q_n <- srvSL$srv
  Qc_n <- srvSL$cens
  omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
  
  # Construct regular Gamma_0 estimator
  if (p$cf_folds==1) {
    # Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n, q_n, gcomp_n, alpha_star_n)
    
    # !!!!! New functions
    n_orig <- length(dat_orig$z)
    p_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$s!=0))
    eta_n <- construct_eta_n(dat, Q_n, p_n, vals=NA)
    gcomp_n <- construct_gcomp_n(dat_orig, vlist$S_grid, Q_n)
    alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, p_n, vals=NA)
    f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
    q_n <- construct_q_n(type=p$q_n_type, dat, dat_orig, omega_n=omega_n, g_n=g_n,
                         p_n=p_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
                         Q_n=Q_n, Qc_n=Qc_n, f_n_srv=f_n_srv)
    Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n,
                                       p_n, q_n, gcomp_n, alpha_star_n)
  }
  
  # Construct cross-fitted Gamma_0 estimator
  if (p$cf_folds>1) {
    # Gamma_os_n <- construct_Gamma_cf(dat_orig, p, vlist)
  }
  
  # Compute the test statistic
  beta_n <- (1/n_orig) * sum(
    weights * ( lambda_2*(Phi_n(dat$s))^2 - lambda_3*Phi_n(dat$s) ) *
      Gamma_os_n(round(dat$s,-log10(C$appx$s)))
  )
  
  if (test_stat_only) {
    
    p_val <- NA
    sd_n <- NA
    
  } else {
    
    # Construct influence functions
    infl_fn_1 <- construct_infl_fn_1(dat, Gamma_os_n, Phi_n, lambda_2,
                                     lambda_3)
    infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n, p_n,
                                             alpha_star_n, q_n, eta_n,
                                             Gamma_os_n)
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

# test2.R: (Regular) bootstrap test statistic
if (F) {
  
  # !!!!! Update all of this if needed
  
  # Define the statistic to bootstrap
  bootstat <- function(dat_orig,indices) {

    dat_orig_b <- ss(dat_orig, indices)

    # Calculate components on bootstrapped dataset
    {
      # Prep
      n_orig <- length(dat_orig_b$z)
      dat <- ss(dat_orig_b, which(dat_orig_b$z==1))
      weights <- dat$weights

      # Construct dataframes of values to pre-compute functions on
      vlist <- create_val_list(dat_orig)

      # Construct component functions
      Phi_n <- construct_Phi_n(dat, type=p$ecdf_type)
      lambda_2 <- lambda(dat,2,Phi_n)
      lambda_3 <- lambda(dat,3,Phi_n)
      f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                   k=p$f_sIx_n_bins)
      f_s_n <- construct_f_s_n(dat_orig_b, vlist$S_grid, f_sIx_n)
      g_n <- construct_g_n(f_sIx_n, f_s_n)
      srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type)
      Q_n <- srvSL$srv
      Qc_n <- srvSL$cens
      omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n,
                                   type=p$omega_n_type)
      Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n, q_n, gcomp_n, alpha_star_n)

    }

    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * ( lambda_2*(Phi_n(dat$s))^2 - lambda_3*Phi_n(dat$s) ) *
        Gamma_os_n(round(dat$s,-log10(C$appx$s)))
    )

    if (F) {

      # Psi_1
      Theta_true <- attr(dat_orig,"Theta_true")
      Gamma_0 <- Vectorize(function(u) {
        Theta_true[which.min(abs(u-seq(0,1,0.02)))[1]]
      })
      Psi_1_est <- (1/n_orig) * sum(weights*(
        lambda_2*(Phi_n(dat$s))^2*Gamma_0(dat$s) -
          lambda_3*Phi_n(dat$s)*Gamma_0(dat$s)
      ))

      # Psi_2
      Phi_0 <- function(u) { u }
      s_mc <- runif(10^6)
      Psi_2_est <- mean(
        (1/3)*(Phi_0(s_mc))^2*Gamma_os_n(round(s_mc,-log10(C$appx$s))) -
          (1/4)*Phi_0(s_mc)*Gamma_os_n(round(s_mc,-log10(C$appx$s)))
      )

    } # DEBUG

    return(list(beta_n=beta_n, Psi_1_est=Psi_1_est, Psi_2_est=Psi_2_est))

  }

  # Run bootstrap
  boot_results <- sapply(c(1:p$boot_reps), function(i) {
    indices <- sample(c(1:length(dat_orig$s)), replace=T)
    bootstat(dat_orig,indices)
  })

  # Calculate beta_n, SD, variance
  stats <- bootstat(dat_orig, c(1:length(dat_orig$z)))
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

# test2.R: Mixed bootstrap test statistic
if (F) {
  
  # !!!!! Update all of this if needed

  # Pre-calculate non-bootstrapped pieces
  {
    dat_0_orig <- dat
    n_orig <- length(dat_0_orig$z)
    dat_0 <- ss(dat_0_orig, which(dat_0_orig$z==1))
    n_0 <- length(dat_0$s)
    weights_0 <- wts(dat_0)

    G_0 <- construct_Phi_n(dat_0, type=p$ecdf_type)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type)
    Qc_0 <- srvSL$srv
    Qc_0 <- srvSL$cens
    omega_0 <- construct_omega_n(Qc_0, Qc_0)
    f_sIx_n <- construct_f_sIx_n(dat_0, type=p$g_n_type, k=p$f_sIx_n_bins)
    f_s_n <- construct_f_s_n(dat_0_orig, f_sIx_n=f_sIx_n)
    g_0 <- construct_g_n(f_sIx_n, f_s_n)
    Gamma_0 <- construct_Gamma_os_n(dat_0, vlist$S_grid, omega_0, Qc_0, g_0)
    lambda_2 <- lambda(dat_0_orig, k=2, G_0)
    lambda_3 <- lambda(dat_0_orig, k=3, G_0)
    # eta_0 <- construct_eta_n(dat_0, vlist$SX_grid, Qc_0) # ARCHIVED
    gcomp_0 <- construct_gcomp_n(dat_0_orig, vlist$S_grid, Qc_0)

    beta_0 <- (1/n_orig) * sum(
      weights_0 *
        (
          lambda(dat_0_orig,2,G_0)*(G_0(dat_0$s))^2 -
            lambda(dat_0_orig,3,G_0)*G_0(dat_0$s)
        ) *
        (Gamma_0(dat_0$s))
    )

    piece_3 <- -2*beta_0

  }

  # Define the statistic to bootstrap
  bootstat <- function(dat_orig,indices) {

    dat_b_orig <- dat_orig[indices,]
    dat_b <- ss(dat_b_orig, which(dat_b_orig$z==1))
    n_b <- length(dat_b$s)
    weights_b <- wts(dat_b)
    Phi_n <- construct_Phi_n(dat_b, type=p$ecdf_type)

    piece_1 <- (1/n_orig) * sum(
      weights_b *
        (
          lambda(dat_b_orig,k=2,G_0)*(G_0(dat_b$s))^2 -
            lambda(dat_b_orig,k=3,G_0)*G_0(dat_b$s)
        ) *
        (Gamma_0(dat_b$s))
    )

    piece_2 <- (1/n_orig) * sum(
      weights_0 *
        (
          lambda(dat_0_orig,2,Phi_n)*(Phi_n(dat_0$s))^2 -
            lambda(dat_0_orig,3,Phi_n)*Phi_n(dat_0$s)
        ) *
        Gamma_0(dat_0$s)
    )

    index_0 <- rep(c(1:n_0), each=n_b)
    index_b <- rep(c(1:n_b), times=n_0)
    s_0_long <- dat_0$s[index_0]
    s_b_long <- dat_b$s[index_b]
    y_b_long <- dat_b$y[index_b]
    x1_b_long <- dat_b$x1[index_b]
    x2_b_long <- dat_b$x2[index_b]
    weights_0_long <- weights_0[index_0]
    weights_b_long <- weights_b[index_b]

    piece_4 <- (1/n_orig)^2 * sum(
      weights_0_long * weights_b_long *
        (
          (
            (as.integer(s_b_long<=s_0_long)) *
              (
                ( y_b_long - mu_0(s_b_long, x1_b_long, x2_b_long) ) /
                  g_0(s_b_long, x1_b_long, x2_b_long) +
                  gcomp_0(s_b_long)
              )
          ) +
            eta_0(s_0_long, x1_b_long, x_b_long) -
            (2*Gamma_0(s_0_long))
        ) *
        (
          lambda(dat_0_orig,2,G_0)*(G_0(s_0_long))^2 -
            lambda(dat_0_orig,3,G_0)*G_0(s_0_long)
        )
    )

    return (beta_0+piece_1+piece_2+piece_3+piece_4)

  }

  # Run bootstrap
  boot_obj <- boot(data=dat_0_orig, statistic=bootstat, R=p$boot_reps)

  # Calculate critical value (for a one-sided test)
  crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
  
}

# construct_Gamma_os_n_star
if (F) {
  
  #' Construct Gamma_os_n_star primitive one-step estimator
  #'
  #' @param x !!!!!
  construct_Gamma_os_n_star <- function(dat, omega_n, g_n_star, eta_n, p_n,
                                        gcomp_n, alpha_star_n, vals=NA) {
    
    weights_i <- dat$weights
    n_orig <- sum(weights_i)
    s_i <- dat$s
    x_i <- dat$x
    piece_1 <- omega_n(dat$x,dat$s,dat$y,dat$delta) /
      g_n_star(dat$s,dat$x)
    piece_2 <- In(s_i!=0)
    piece_3 <- gcomp_n(s_i)
    
    # Remove large intermediate objects
    rm(dat,omega_n,g_n_star,gcomp_n)
    
    fnc <- function(u) {
      (1/n_orig) * sum(weights_i * (
        piece_2*In(s_i<=u)*piece_1 +
          eta_n(rep(u,nrow(x_i)),x_i) +
          (piece_2/p_n)*(In(s_i<=u)*piece_3-alpha_star_n(u))
      ))
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
    
  }
  
}

# construct_Gamma_os_n
if (F) {
  
  #' Construct Gamma_os_n primitive one-step estimator
  #' 
  #' @param dat Subsample of dataset returned by ss() for which z==1
  #' @param vals List of values to pre-compute function on; passed to
  #'     construct_superfunc()
  #' @param omega_n A nuisance influence function returned by construct_omega_n()
  #' @param Q_n A conditional survival function returned by construct_Q_n()
  #' @param g_n A density ratio estimator function returned by construct_g_n()
  #' @param type One of c("one-step", "plug-in")
  #' @return Gamma_os_n estimator
  #' @notes This is a generalization of the one-step estimator from Westling &
  #'     Carone 2020
  construct_Gamma_os_n <- function(dat, vals=NA, omega_n, Q_n, g_n,
                                   type="one-step") {
    
    n_orig <- sum(dat$weights)
    n <- length(dat$s)
    i_long <- rep(c(1:n), each=n)
    j_long <- rep(c(1:n), times=n)
    s_i_short <- dat$s
    s_i_long <- dat$s[i_long]
    x_j_long <- dat$x[j_long,]
    delta_i_long <- dat$delta[i_long]
    delta_j_long <- dat$delta[j_long]
    weights_i_long <- dat$weights[i_long]
    weights_j_long <- dat$weights[j_long]
    
    if (type=="one-step") {
      
      subpiece_1a <- dat$weights * ( 1 + (
        omega_n(dat$x,dat$s,dat$y,dat$delta) /
          g_n(dat$s,dat$x)
      ) )
      subpiece_2a <- (weights_i_long*weights_j_long) *
        Q_n(rep(C$t_0, length(s_i_long)),x_j_long,s_i_long)
      
      # Remove large intermediate objects
      rm(dat,delta_i_long,delta_j_long,i_long,j_long,omega_n,Q_n,x_j_long,
         weights_i_long,weights_j_long)
      
      fnc <- function(s) {
        
        subpiece_1b <- In(round(s_i_short,-log10(C$appx$s))<=
                            round(s,-log10(C$appx$s)))
        piece_1 <- (1/n_orig) * sum(subpiece_1a*subpiece_1b)
        
        subpiece_2b <- In(round(s_i_long,-log10(C$appx$s))<=
                            round(s,-log10(C$appx$s)))
        piece_2 <- (1/(n_orig^2)) * sum(subpiece_2a*subpiece_2b)
        
        return(piece_1-piece_2)
        
      }
      
    }
    
    if (type=="plug-in") {
      
      piece <- weights_i_long*weights_j_long *
        (1 - Q_n(rep(C$t_0, length(s_i_long)),x_j_long,s_i_long))
      
      # Remove large intermediate objects
      rm(dat,delta_i_long,delta_j_long,i_long,j_long,omega_n,Q_n,
         x_i_long,x_j_long,weights_i_long,weights_j_long)
      
      fnc <- function(s) {
        (1/n_orig^2) * sum(piece * In(s_i_long<=s))
      }
      
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
    
  }
  
}

# construct_eta_n
if (F) {
  
  #' Construct nuisance estimator eta_n
  #' 
  #' @param dat Subsample of dataset returned by ss() for which z==1
  #' @param vals List of values to pre-compute function on; passed to
  #'     construct_superfunc()
  #' @param Q_n Conditional survival function estimator returned by construct_Q_n
  #' @return Estimator function of nuisance eta_0
  construct_eta_n <- function(dat, vals=NA, Q_n) {
    
    n_orig <- sum(dat$weights)
    
    fnc <- function(u,x) {
      
      x_long <- as.data.frame(
        matrix(rep(x,length(dat$s)), ncol=length(x), byrow=T)
      )
      
      return(
        (1/n_orig) * sum(
          dat$weights * In(dat$s<=u) *
            (1-Q_n(rep(C$t_0,length(dat$s)),x_long,dat$s))
        )
      )
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
    
  }
  
}

# construct_q_n(type=="SuperLearner")
if (F) {
  
  if (type=="Super Learner") {
    
    # Create grid of x-values and container for regression predictions
    x_grid <- round(seq(0.1,1,0.1),1) # Try 0.01, 0.02, or 0.05 !!!!!
    preds <- list()
    
    # Set up objects
    a <- dat$a
    w <- dat$w
    y_star <- dat$y_star
    delta_star <- dat$delta_star
    X <- cbind(w, y_star=y_star, delta_star=delta_star)
    newX <- distinct(cbind(
      dat_orig$w,
      y_star = dat_orig$y_star,
      delta_star = dat_orig$delta_star
    ))
    
    # Set library
    if (type=="Super Learner") {
      SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
                      "SL.svm")
    } else if (type=="GAM") {
      SL.library <- c("SL.gam")
    }
    
    for (i in c(1:length(x_grid))) {
      
      # Create pseudo-outcomes
      x <- x_grid[i]
      if (which=="q_n") {
        po <- (
          (In(a!=0 & a<=x)*omega_n(w,a,y_star,delta_star))/g_n_star(a,w)
        ) + (
          (In(a!=0)/z_n) * (In(a<=x)*gcomp_n(a) - alpha_star_n(x))
        )
      } else if (which=="q_star_n") {
        po <- (In(a<=x)*omega_n(w,a,y_star,delta_star))/f_aIw_n(a,w)
      }
      
      # Fit SuperLearner regression
      model_sl <- SuperLearner(Y=po, X=X, newX=newX, family="gaussian",
                               SL.library=SL.library, verbose=F)
      preds[[i]] <- as.numeric(model_sl$SL.predict)
      rm(model_sl)
      
    }
    
    # Construct function
    newX$index <- c(1:nrow(newX))
    fnc <- function(w, y_star, delta_star, x) {
      
      if (x==0) { return(0) } else {
        
        # Choose which regression to use based on `x`
        pred <- preds[[which.min(abs(x-x_grid))]]
        
        # Dynamically filter to select index
        # !!!!! Test if this works for factors
        # !!!!! Modify construct_S_n to follow this paradigm instead
        cond <- paste0("y_star==",y_star," & delta_star==",delta_star,"")
        for (i in c(1:length(w))) {
          cond <- paste0(cond," & w",i,"==",w[i])
        }
        index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
        if (length(index)!=1) {
          stop(paste0("Error in construct_q_n; ", "w=(", paste(w,collapse=","),
                      "), y_star=",y_star,", delta_star=",delta_star,", x=",x))
        }
        
        # Return prediction
        return(pred[index])
      }
      
    }
    
    # Remove large intermediate objects
    rm(dat,dat_orig,omega_n,g_n_star,gcomp_n,alpha_star_n,f_aIw_n)
    
    # !!!!! Plot regression predictions
    if (F) {
      
      # Omit `newX = newX` and set `newX <- X` to test
      
      # Generate predictions
      sfnc <- construct_superfunc(fnc, aux=NA, vec=c(2,1,1,0), vals=NA)
      pred_y <- sfnc(dat$w, dat$y_star, dat$delta_star, x)
      
      # Plot pseudo-outcomes against predictions
      plot_data <- data.frame(x=po, y=pred_y, w1=dat$w$w1, w2=dat$w$w2,
                              y_star=y_star, delta_star=delta_star)
      ggplot(plot_data, aes(x=x, y=y, color=factor(delta_star))) +
        geom_point() +
        lims(x=c(-1.6,1.6), y=c(-1.6,1.6)) +
        # labs(title="SL.xgboost") +
        geom_abline(slope=1, intercept=0, color="grey")
      
      # Calculate MSE
      mean((pred_y-po)^2)
      
      # See SL weights (need to remove `rm(model_sl)`)
      coef(model_sl)
      
    }
    
  }
  
}

# Cox gcomp archive
if (F) {
  
  # Estimates the integral \int_a^b h(x)dG(x) based on a Riemann sum approximation
  # m is the number of rectangles in the Riemann sum
  int_appx <- function(h, G, a, b, m=500) {
    i <- c(1:m)
    x1 <- (i*(b-a))/m
    x2 <- ((i-1)*(b-a))/m
    sum((G(a+x1)-G(a+x2))*h(a+x1))
  }
  
  # Estimates the integral \int_a^b h(x)dG(x)
  # G should be a dataframe with columns `x` and `y`, where x is a vector of step
  #     function jump points and y=G(x)
  int_step <- function(h, G, a, b) {
    x1 <- G$x[-1]
    y1 <- G$y[-1]
    y2 <- G$y[-length(G$y)]
    sum(as.integer(x1>=a & x1<=b)*(y1-y2)*h(x1))
  }
  
  # TEMP: Influence function components of Breslow estiamtor
  infl_fn_2a_Bres <- function(z_i,delta_i,t_i) {
    pc_2 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
      Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_n))*Lambda_n(T_[j])) *
        Q_n(Z_[,j],D_[j],T_[j])
    }))
    pc_3 <- l_tilde(z_i,delta_i,t_i)
    return(sum(pc_2*pc_3))
  }
  infl_fn_2b_Bres <- function(z_i,delta_i,t_i) {
    pc_1 <- Q_n(z_i,delta_i,t_i)
    return(pc_1)
  }
  
  # Calculate variance of Breslow estimator (Lin formula)
  piece_1 <- (1/N) * Reduce("+", lapply(t_ev, function(t_i) {
    (as.integer(t_i<=t)*m_n(t_i)) / S_0n(t_i)
  }))
  piece_2 <- (1/N) * sum(sapply(t_ev, function(t_i) {
    as.integer(t_i<=t) / (S_0n(t_i))^2
  }))
  var_bshz_est2 <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1)+piece_2)/N
  
  # TEMP: variance component estimate
  var_cmhz_est2a <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_2a_Bres(Z_[,i],D_[i],T_[i]))^2
  }))
  var_cmhz_est2b <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2
    # (N*infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2 # !!!!!
  }))
  # var_bshz_est_a <- as.numeric(N*(t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  var_bshz_est_a <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  var_bshz_est_b <- as.numeric(piece_2)/N
  
  # Estimate information matrix using score; should be consistent
  I_tilde2 <- matrix(0, nrow=d, ncol=d)
  for (i in c(1:N)) {
    score <- l_star(Z_[,i],D_[i],T_[i])
    I_tilde2 <- I_tilde + score %*% t(score)
  }
  I_tilde_inv2 <- solve(I_tilde2)
  return(list(
    I_tilde_inv = I_tilde_inv,
    I_tilde_inv2 = I_tilde_inv2,
    var_cmhz_est = var_cmhz_est
  ))
  
  # Calculate the baseline cumulative hazard manually
  cuml_haz_n <- function(t, w, a) {
    return(base_hz_n(t)*exp(cfs[["w1"]]*w[1]+cfs[["w2"]]*w[2]+cfs[["a"]]*a))
  }
  
  # Calculate the cumulative hazard via predict()
  newdata <- data.frame(y_star=14, delta_star=1, a=0, w1=0, w2=0)
  predict(model,newdata=newdata, type="expected", se.fit=F)
  x=survfit(model, newdata=newdata, stype=2)
  
  # !!!!! Testing new int function
  h <- function(x) { sqrt(x) }
  int_appx(h, G=base_hz_n, 0, 100)
  int_step(h, G=data.frame(x=bh$time, y=bh$hazard), 0, 100)
  
  # !!!!! Use true baseline hazard function
  H_0_true <- function(t) {
    # L$sc_params$lmbd * t^L$sc_params$v
    L$sc_params$lmbd*exp(-1.7) * t^L$sc_params$v # !!!!!
  }
  H_0_true(c(10,50,100,200))
  base_hz_n(c(10,50,100,200))
  
  t_ev <- sort(dat$y_star[dat$delta_star==1])
  Lambda_n2 <- Vectorize(function(x) {
    sum(sapply(t_ev, function(t_i) {
      as.integer(t_i<=x) / (sum(exp(lin)*as.integer(dat$y_star>=t_i)))
    }))
  })
  H_0_true(c(10,50,100,200,400,800))
  Lambda_n(c(10,50,100,200,400,800))
  Lambda_n2(c(10,50,100,200,400,800))
  base_hz_n(c(10,50,100,200,400,800))
  
  
}

# Old pieces related to cox_var()
if (F) {
  
  # Beta estimates that don't account for weights
  var_est_betas <- (1/N^2) * Reduce("+", lapply(c(1:n), function(i) {
    (WT[i] * l_tilde(Z_[,i],Ds_[i],T_[i]))^2
  })) %>% as.numeric()
  res$se_w1_MC <- sqrt(var_est_betas[1])
  res$se_w2_MC <- sqrt(var_est_betas[2])
  res$se_a_MC <- sqrt(var_est_betas[3])
  
  # Calculate information using score
  I_tilde2 <- (1/N) * Reduce("+", lapply(c(1:N), function(i) {
    WT[i]*l_star(Z_[,i],Ds_[i],T_[i]) %*% t(WT[i]*l_star(Z_[,i],Ds_[i],T_[i]))
  }))
  
  # Variance of cumulative hazard estimator
  # !!!!! Adapt to two-phase sampling
  var_cmhz_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (omega_n(Z_[,i],Ds_[i],T_[i],z_0))^2
  }))
  
  # Influence function of cumulative hazard estiamtor
  # !!!!! Adapt to two-phase sampling
  infl_fn_3 <- (function() {
    pc_1 <- exp(sum(beta_n*z_0))
    pc_3 <- z_0 * Lambda_n(t)
    return(function(z_i,ds_i,t_i) {
      pc_2 <- Q_n(z_i,ds_i,t_i)
      pc_4 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
        Z_[,j] * (Ds_[j]-exp(sum(Z_[,j]*beta_n))*Lambda_n(T_[j])) *
          Q_n(Z_[,j],Ds_[j],T_[j])
      }))
      pc_5 <- l_tilde(z_i,ds_i,t_i)
      return(pc_1*(pc_2+sum((pc_3-pc_4)*pc_5)))
    })
  })()
  
  # Variance estimate (cumulative hazard) using influence function
  # !!!!! Adapt to two-phase sampling
  var_cmhz_est2 <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_3(Z_[,i],Ds_[i],T_[i]))^2
  }))
  
  # Variance estimate (survival) using influence function
  # !!!!! Adapt to two-phase sampling
  surv_est <- exp(-exp(sum(z_0*beta_n))*Lambda_n(t))
  var_surv_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (surv_est*infl_fn_3(Z_[,i],Ds_[i],T_[i]))^2
  }))
  
  res$var_cmhz_est <- var_cmhz_est # !!!!! var_cmhz_est or var_cmhz_est2 ?????
  res$var_surv_est <- var_surv_est
  
  # Calculate component estimator Q_n
  Q_n <- memoise(function(z_i,ds_i,t_i) {
    piece_1 <- (ds_i*as.integer(t_i<=t)) / S_0n(t_i)
    piece_2 <- exp(sum(z_i*beta_n))
    piece_3 <- (1/N) * sum(sapply(i_ev, function(j) {
      WT[j] * as.integer(T_[j]<=min(t,t_i)) / (S_0n(T_[j]))^2
    }))
    return(piece_1-piece_2*piece_3)
  })
  
  # Influence function of Breslow estiamtor
  infl_fn_2 <- function(z_i,ds_i,t_i) {
    pc_1 <- Q_n(z_i,ds_i,t_i)
    pc_2 <- (1/N) * Reduce("+", lapply(c(1:n), function(j) {
      WT[j] * Z_[,j] * (Ds_[j]-exp(sum(Z_[,j]*beta_n))*Lambda_n(T_[j])) *
        Q_n(Z_[,j],Ds_[j],T_[j])
    }))
    pc_3 <- l_tilde(z_i,ds_i,t_i)
    return(pc_1-sum(pc_2*pc_3))
  }
  
  # Influence function of Breslow estimator (new derivation)
  infl_fn_2b <- function(z_i,ds_i,t_i) {
    pc_1 <- Q_n(z_i,ds_i,t_i)
    pc_2 <- (1/N^2) * Reduce("+", lapply(c(1:n), function(i) {
      WT[i] * Z_[,i] * exp(sum(Z_[,i]*beta_n)) *
        sum(sapply(t_ev, function(t_j) {
          as.integer(t_j<=min(t,T_[i])) / (S_0n(t_j))^2
        }))
    }))
    pc_3 <- l_tilde(z_i,ds_i,t_i)
    return(pc_1-sum(pc_2*pc_3))
  }
  
  # Variance estimate of Breslow estimator using influence function
  var_est_bshz <- (1/N^2) * sum(sapply(c(1:n), function(i) {
    (WT[i] * infl_fn_2(Z_[,i],Ds_[i],T_[i]))^2
  }))
  
  # Variance estimate of Breslow estimator using influence function (new derivation)
  var_est_bshz2 <- (1/N^2) * sum(sapply(c(1:n), function(i) {
    (WT[i] * infl_fn_2b(Z_[,i],Ds_[i],T_[i]))^2
  }))
  
}

# Comparing old gamma_n estimators
if (F) {
  
  # OLD
  ggplot(df, aes(x=a, y=po)) +
    geom_point(alpha=0.3) +
    # ylim(c(0,1)) +
    geom_line(data=data.frame(
      a = x_grid,
      po = reg(x_grid)
    ),
    color="forestgreen") +
    labs(x=paste0("a (bw=",round(best$bw,4),")"))
  
  # NEW
  reg2 <- construct_superfunc(reg, aux=NA, vec=c(2,1), vals=NA)
  ggplot(
    data.frame(x=po, y=reg2(dat$w,dat$a), z=dat$a),
    aes(x=x, y=y, color=z)
  ) +
    geom_point() +
    geom_abline(slope=1, intercept=0, color="grey")
  
}

# Debugging hypothesis test (test_2.R, after "# Estimate variance")
if (F) {
  
  # Psi_1
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
  test_stat_Psi_1 <- Psi_1_est^2/Psi_1_var_est
  p_val_Psi_1 <- pchisq(test_stat_Psi_1, df=1, lower.tail=FALSE)
  
  # Psi_2
  Phi_0 <- function(x) {x}
  # infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
  #                                          eta_n, Gamma_os_n)
  infl_fn_Gamma <- construct_infl_fn_Gamma2(omega_n, g_n_star, gcomp_n,
                                            z_n, alpha_star_n, q_n,
                                            eta_ss_n,
                                            Gamma_os_n_star=Gamma_os_n)
  infl_fn_2 <- construct_infl_fn_2(dat, Phi_0, infl_fn_Gamma, 1/3, 1/4)
  # Psi_2_var_est <- (1/n_orig^2) * sum((
  #   weights*infl_fn_2(dat$w,dat$y_star,dat$delta_star,dat$a)
  # )^2)
  Psi_2_var_est <- (1/n_orig^2) * sum((
    infl_fn_2(dat_orig$w, dat_orig$y_star, dat_orig$delta_star,
              dat_orig$a, dat_orig$weights)
  )^2)
  a_mc <- runif(10^6)
  Psi_2_est <- mean(
    (1/3)*(Phi_0(a_mc))^2*Gamma_os_n(round(a_mc,-log10(C$appx$a))) -
      (1/4)*Phi_0(a_mc)*Gamma_os_n(round(a_mc,-log10(C$appx$a)))
  )
  test_stat_Psi_2 <- Psi_2_est^2/Psi_2_var_est
  p_val_Psi_2 <- pchisq(test_stat_Psi_2, df=1, lower.tail=FALSE)
  
  # Psi_G
  xx <- 0.3
  Psi_G_est <- Gamma_os_n(xx)
  Psi_G_var_est <- (1/n_orig^2) * sum((
    infl_fn_Gamma(rep(xx, n_orig), dat_orig$w, dat_orig$y_star,
                  dat_orig$delta_star, dat_orig$a, dat_orig$weights)
  )^2)
  test_stat_Psi_G <- (Psi_G_est-Gamma_0(xx))^2/Psi_G_var_est
  p_val_Psi_G <- pchisq(test_stat_Psi_G, df=1, lower.tail=FALSE)
  
  # Covariance
  Psi_12_covar <- 999 # !!!!! Below is wrong bc infl_fn_2 needs dat_orig
  # Psi_12_covar <- (1/n_orig^2) * sum(
  #   weights *
  #     infl_fn_1(dat$a) *
  #     infl_fn_2(dat$w,dat$y_star,dat$delta_star,dat$a,dat$weights)
  # )
  
  # Combined (Psi_1+Psi_2)
  # infl_fn_sum12 <- function(a,w,y_star,delta_star) {
  #   infl_fn_1(a) + infl_fn_2(w,y_star,delta_star,a)
  # }
  # sum12_var_est <- (1/n_orig^2) * sum((
  #   weights*infl_fn_sum12(dat$a,dat$w,dat$y_star,dat$delta_star)
  # )^2)
  # test_stat_sum12 <- (Psi_1_est+Psi_2_est)^2/sum12_var_est
  # p_val_sum12 <- pchisq(test_stat_sum12, df=1, lower.tail=FALSE)
  # test_stat_sum12b <- beta_n^2/sum12_var_est
  # p_val_sum12b <- pchisq(test_stat_sum12b, df=1, lower.tail=FALSE)
  
  # Infl fn means
  if1_mean <- (1/n_orig) * sum(weights * infl_fn_1(dat$a))
  if2_mean <- (1/n_orig) * sum(infl_fn_2(
    dat_orig$w, dat_orig$y_star, dat_orig$delta_star,
    dat_orig$a, dat_orig$weights
  ))
  
}

# Comparing derivative estimators
if (F) {
  
  deriv_theta_0 <- Vectorize(function(x) {
    theta_0 <- attr(dat_orig, "theta_true")
    index <- which.min(abs(x-seq(0,1,0.02)))
    if (index==1) {
      return((theta_0[2]-theta_0[1])/0.02)
    } else if (index==51) {
      return((theta_0[51]-theta_0[50])/0.02)
    } else {
      return((theta_0[round(index+1)]-theta_0[round(index-1)])/0.04)
    }
  })
  deriv_linear <- construct_deriv_theta_n(theta_n, type="linear", dir=dir)
  deriv_line <- construct_deriv_theta_n(theta_n, type="line", dir=dir)
  deriv_spline <- construct_deriv_theta_n(theta_n, type="spline", dir=dir)
  deriv_mspl <- construct_deriv_theta_n(theta_n, type="m-spline", dir=dir)
  
  grid <- round(seq(0,1,0.02), 2)
  # df_plot1 <- data.frame(
  #   x = rep(grid,2),
  #   y = c(attr(dat_orig, "theta_true"),
  #         theta_n(grid)),
  #   which = rep(c("theta_0","theta_n"), each=length(grid))
  # )
  # ggplot(df_plot1, aes(x=x, y=y, color=which)) +
  #   geom_line()
  df_plot2 <- data.frame(
    x = rep(grid,5),
    y = c(deriv_theta_0(grid),
          deriv_linear(grid),
          deriv_line(grid),
          deriv_spline(grid),
          deriv_mspl(grid)),
    which = rep(c("true","linear", "line", "spline", "mspline"),
                each=length(grid))
  )
  ggplot(df_plot2, aes(x=x, y=y, color=which)) +
    geom_line() +
    ylim(-2.5,0)
  
}

# Old code from est_curve
if (F) {
  
  # LOD shift
  a_min <- min(dat_orig$a,na.rm=T)
  lod12 <- a_min
  if (params$lod_shift=="3/4") {
    lod34 <- log10(1.5*(10^lod12))
    dat_orig$a[dat_orig$a==lod12] <- lod34
    a_min <- lod34
  } else if (params$lod_shift=="1") {
    lod <- log10(2*(10^lod12))
    dat_orig$a[dat_orig$a==lod12] <- lod
    a_min <- lod
  }
  
  # Phi_n_inv_notrans
  Phi_n_inv_notrans <- construct_Phi_n(dat, which="inverse",
                                       type=params$ecdf_type)
  
}

# Construct Psi_n
if (F) {
  
  # Phased this out for alternate method based on CIR paper
  if (params$marg=="Theta") {
    if (dir=="incr") {
      Psi_n <- Theta_os_n
    } else {
      Psi_n <- Vectorize(function(x) { -1 * Theta_os_n(x) })
    }
  } else if (params$marg=="Gamma") {
    if (dir=="incr") {
      Psi_n <- Vectorize(function(x) {
        Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    } else {
      Psi_n <- Vectorize(function(x) {
        -1 * Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    }
  } else if (params$marg %in% c("Gamma_star", "Gamma_star2")) {
    if (dir=="incr") {
      Psi_n <- Vectorize(function(x) {
        Gamma_os_n_star(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    } else {
      Psi_n <- Vectorize(function(x) {
        -1 * Gamma_os_n_star(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    }
  }
  
}

# Testing that the two epectations are equal
if (F) {
  
  df_res <- data.frame(
    "n" = integer(),
    "rep" = integer(),
    "E1a" = double(),
    "E2a" = double(),
    "E1b" = double(),
    "E2b" = double(),
    "E1c" = double(),
    "E2c" = double()
  )
  for (nnn in c(1280)) { # c(40,80,160,320)
    for (i in c(1:20)) {
      
      dat_orig <- generate_data(n=nnn, -2, "Unif(0,1)", "none", surv_true="Cox PH", list(lmbd=1e-3,v=1.5,lmbd2=5e-5,v2=1.5), "iid", "decr")
      dat <- ss(dat_orig, which(dat_orig$delta==1))
      
      model <- coxph(
        formula = formula(paste0("Surv(y_star,delta_star)~",
                                 paste(names(dat$w),collapse="+"),"+a")),
        data = cbind(y_star=dat$y_star, delta_star=dat$delta_star,
                     dat$w, a=dat$a),
        weights = dat$weights * (length(dat$weights)/sum(dat$weights))
      )
      theta_n <- as.numeric(model$coefficients)
      WT <- dat$weights
      ST <- dat$strata
      t <- C$t_e
      N <- round(sum(WT))
      n <- length(WT)
      Z_ <- t(as.matrix(cbind(dat$w,a=dat$a)))
      T_ <- dat$y_star
      Ds_ <- dat$delta_star
      lin <- as.numeric(t(theta_n)%*%Z_)
      S_0n <- function(x) {
        (1/N) * sum(WT*as.integer(T_>=x)*exp(lin))
      }
      d <- dim(Z_)[1]
      i_ev <- which(Ds_==1)
      t_ev <- T_[which(Ds_==1)]
      Lambda_n <- Vectorize(memoise(function(x) {
        (1/N) * sum(sapply(i_ev, function(i) {
          (WT[i] * as.integer(T_[i]<=x)) / S_0n(T_[i])
        }))
      }))
      Q_n <- memoise(function(z_i,ds_i,t_i) {
        piece_1 <- (ds_i*as.integer(t_i<=t)) / S_0n(t_i)
        piece_2 <- exp(sum(z_i*theta_n))
        piece_3 <- (1/N) * sum(sapply(i_ev, function(j) {
          WT[j] * as.integer(T_[j]<=min(t,t_i)) / (S_0n(T_[j]))^2
        }))
        return(piece_1-piece_2*piece_3)
      })
      E1 <- (1/N^2) * Reduce("+", lapply(c(1:n), function(i) {
        WT[i] * Z_[,i] * exp(sum(Z_[,i]*theta_n)) *
          sum(sapply(t_ev, function(t_j) {
            as.integer(t_j<=min(t,T_[i])) / (S_0n(t_j))^2
          }))
      }))
      E2 <- (1/N) * Reduce("+", lapply(c(1:n), function(j) {
        WT[j] * Z_[,j] * (Ds_[j]-exp(sum(Z_[,j]*theta_n))*Lambda_n(T_[j])) *
          Q_n(Z_[,j],Ds_[j],T_[j])
      }))
      
      df_res[nrow(df_res)+1,] <- c(nnn,i,E1[1],E2[1],E1[2],E2[2],E1[3],E2[3])
      
      print(paste0("Done with n=", nnn, ", rep=", i))
      print(Sys.time())
      
    }
  }
  df_plot <- pivot_longer(df_res, cols=c(E1a,E2a,E1b,E2b,E1c,E2c))
  ggplot(
    filter(df_plot, name %in% c("E1a","E2a")),
    aes(x=n, y=value, color=name)
  ) + geom_jitter(width=1, height=0, alpha=0.5)
  
  # ggplot(df_res, aes(x=n, y=diff1, color=n)) + geom_point()
  
}

# kernel2 gamma_n estimator
if (F) {
  
  df_0 <- filter(df, I_star==0)
  df_1 <- filter(df, I_star==1)
  
  ks_0 <- ksmooth(
    x = df_0$a,
    y = df_0$po,
    kernel = "normal",
    bandwidth = 0.2, # !!!!! Select via CV
    x.points = vals$a
  )
  
  ks_1 <- ksmooth(
    x = df_1$a,
    y = df_1$po,
    kernel = "normal",
    bandwidth = 0.2, # !!!!! Select via CV
    x.points = vals$a
  )
  
  reg_0 <- function(a) {
    index <- which.min(abs(a-ks_0$x))
    return(ks_0$y[index])
  }
  
  reg_1 <- function(a) {
    index <- which.min(abs(a-ks_1$x))
    return(ks_1$y[index])
  }
  
  # Run logistic regression
  model <- glm(I_star~a, data=df, family="binomial")
  coeff <- as.numeric(model$coefficients)
  
  prob_1 <- function(a) { expit(coeff[1]+coeff[2]*a) }
  prob_0 <- function(a) { 1 - prob_1(a) }
  
  reg <- function(a) {
    reg_0(a)*prob_0(a) + reg_1(a)*prob_1(a)
  }
  
}

# Debugging gamma_n CV selection of bandwidth
if (F) {
  
  # Create pseudo-outcome dataframe
  delta_prob <- mean(dat_orig$delta)
  I_star <- dat$delta_star*as.integer(dat$y_star<=C$t_e)
  po <- ((dat$weights*omega_n(dat$w,dat$a,dat$y_star,dat$delta_star)) /
           f_aIw_n(dat$a,dat$w))^2
  df <- data.frame(a=dat$a, po=po, I_star=I_star)
  df %<>% filter(is.finite(po))
  
  # Select bandwidth via cross-validation
  {
    # CV prep
    n_folds <- 5
    folds <- sample(cut(c(1:nrow(df)), breaks=n_folds, labels=FALSE))
    # bws <- seq(2*C$appx$a, 0.3, length.out=30)
    bws <- seq(0.02,0.4,0.02)
    # bws <- seq(2*0.001, 0.3, length.out=30)
    
    # Conduct CV
    best <- list(bw=999, sum_sse=999)
    for (bw in bws) {
      sum_sse <- 0
      for (i in c(1:n_folds)) {
        df_train <- df[-which(folds==i),]
        df_test <- df[which(folds==i),]
        ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal", bandwidth=bw)
        reg <- Vectorize(function(a) {
          index <- which.min(abs(a-ks$x))
          return(ks$y[index])
        })
        sum_sse <- sum_sse + sum((reg(df_test$a)-df_test$po)^2, na.rm=T)
      }
      if (sum_sse<best$sum_sse || best$sum_sse==999) {
        best$bw <- bw
        best$sum_sse <- sum_sse
      }
      
    }
    
  }
  
  # Construct optimal function from true data
  ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal",
                bandwidth=best$bw)
  reg <- Vectorize(function(a) {
    index <- which.min(abs(a-ks$x))
    return(ks$y[index])
  })
  
  # Generate plots 1-3
  plot_gamma1 <- ggplot(df, aes(x=a, y=po)) +
    geom_point(alpha=0.3) +
    ylim(c(0,10)) +
    geom_line(data=data.frame(
      a = seq(0,1,0.01),
      po = reg(seq(0,1,0.01))
    ),
    color="forestgreen") +
    labs(x=paste0("a (bw=",round(best$bw,4),")"))
  ggsave(
    filename = paste0("Moderna plots/debug_gamma_",
                      Sys.getenv("SLURM_ARRAY_TASK_ID"),"_ylim(10).pdf"),
    plot=plot_gamma1, device="pdf", width=6, height=4
  )
  plot_gamma2 <- plot_gamma1 + ylim(c(0,1))
  ggsave(
    filename = paste0("Moderna plots/debug_gamma_",
                      Sys.getenv("SLURM_ARRAY_TASK_ID"),"_ylim(1).pdf"),
    plot=plot_gamma2, device="pdf", width=6, height=4
  )
  plot_gamma3 <- plot_gamma1 + ylim(c(0,0.1))
  ggsave(
    filename = paste0("Moderna plots/debug_gamma_",
                      Sys.getenv("SLURM_ARRAY_TASK_ID"),"_ylim(0.1).pdf"),
    plot=plot_gamma3, device="pdf", width=6, height=4
  )
  plot_gamma4 <- plot_gamma1 + ylim(c(0,0.01))
  ggsave(
    filename = paste0("Moderna plots/debug_gamma_",
                      Sys.getenv("SLURM_ARRAY_TASK_ID"),"_ylim(0.01).pdf"),
    plot=plot_gamma4, device="pdf", width=6, height=4
  )
  
  # Generate plot 5
  mult_factor <- Vectorize(function(a) {
    delta_prob*(f_a_delta1_n(a)/f_a_n(a))
  })
  plot_gamma5 <- ggplot(data.frame(
    x = seq(0,1,0.01),
    y = mult_factor(seq(0,1,0.01))
  ), aes(x=x, y=y)) +
    geom_line()
  ggsave(
    filename = paste0("Moderna plots/debug_gamma2_",
                      Sys.getenv("SLURM_ARRAY_TASK_ID"),".pdf"),
    plot=plot_gamma5, device="pdf", width=6, height=4
  )
  
}



# Old infl_fn_1
if (F) {
  
  #' !!!!! document
  #'
  #' @param x x
  #' @return x
  construct_infl_fn_1 <- function(dat, Gamma_os_n, Phi_n, xi_n, rho_n,
                                  lambda_2, lambda_3, vals=NA) {

    n_orig <- sum(dat$weights)
    weights_j <- dat$weights
    a_j <- dat$a

    fnc <- function(a_i) {

      piece_1 <- (lambda_2*(Phi_n(a_i)^2) - lambda_3*Phi_n(a_i)) *
        Gamma_os_n(round(a_i, -log10(C$appx$a)))

      piece_2 <- (1/n_orig) * sum(
        weights_j * (lambda_2*(Phi_n(a_j)^2) - lambda_3*Phi_n(a_j)) *
          Gamma_os_n(round(a_j, -log10(C$appx$a)))
      )

      # piece_1 <- (1/n_orig) * sum(
      #   weights_j * (xi_n(a_i,a_j) - rho_n(a_i)) *
      #     Gamma_os_n(round(a_j, -log10(C$appx$a)))
      # )

      # piece_2 <- (lambda_2*(Phi_n(a)^2) - lambda_3*Phi_n(a)) *
      #   Gamma_os_n(round(a, -log10(C$appx$a)))

      return(piece_1-piece_2)

    }

    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=vals))

  }
  
}

# Alternative (incorrect) rho_n and xi_n
if (F) {
  
  rho_n <- Vectorize(function(k,a,a_i) {
    (1/n_orig) * sum( weights * (
      (Phi_n(a))^k * ( as.integer(a_i<=a) - Phi_n(a) ) * (Phi_n(a))^(3-k)
    ))
  })
  xi_n <- Vectorize(function(a,a_i) {
    (Phi_n(a_i)^2-4*lambda_2)*(Phi_n(a))^2 +
      (3*lambda_3 - (Phi_n(a_i))^3 + 2*lambda_2*as.integer(a_i<=a))*Phi_n(a) -
      lambda_3*as.integer(a_i<=a)
  })
  
  construct_infl_fn_1b <- function(dat, Gamma_0, Phi_n, xi_n, rho_n,
                                   lambda_2, lambda_3, vals=NA) {
    
    n_orig <- sum(dat$weights)
    weights_j <- dat$weights
    a_j <- dat$a
    
    fnc <- function(a_i) {
      
      rho_n(1,a_j,a_i)
      rho_n(2,a_j,a_i)
      
      piece_1 <- (1/n_orig) * sum(
        weights_j * 2*rho_n(rep(1,length(a_j)),a_j,rep(a_i,length(a_j))) -
          3*rho_n(rep(1,length(a_j)),a_j,rep(a_i,length(a_j))) +
          xi_n(a_j,rep(a_i,length(a_j))) *
          Gamma_0(a_j) # Gamma_os_n(round(a_j, -log10(C$appx$a)))
      )
      
      piece_2 <- (lambda_2*(Phi_n(a_i)^2) - lambda_3*Phi_n(a_i)) *
        Gamma_0(a_i) # Gamma_os_n(round(a_j, -log10(C$appx$a)))
      
      return(piece_1+piece_2)
      
    }
    
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=vals))
    
  }
  
}

# Old GAM code inside construct_S_n()
if (F) {
  
  fml <- "Surv(y_star,delta_star)~a"
  for (i in 1:length(dat$w)) {
    fml <- paste0(fml, "+w",i)
  }
  fml <- formula(fml)
  df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star, "a"=dat$a,
              dat$w, "weights"=dat$weights)

  model <- gam(
    formula = fml,
    family = cox.ph(),
    data = df,
    weights = dat$weights
  )

  model <- gam(
    time~s(age,by=sex)+sex+s(nodes)+perfor+rx+obstruct+adhere,
           family=cox.ph(),data=col1,weights=status)

  model <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=100,
                 splitrule="logrank", nsplit=0, case.wt=df$weights,
                 samptype="swor")
  if (return_model) { return(model) }



  newX <- cbind(vals$w, a=vals$a)[which(vals$t==0),]
  pred <- predict(model, newdata=newX)

  fnc <- function(t, w, a) {
    r <- list()
    for (i in 1:length(w)) {
      r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
    }
    r[[length(w)+1]] <- which(abs(a-newX[["a"]])<1e-8)
    row <- Reduce(intersect, r)
    col <- which.min(abs(t-pred$time.interest))
    if (length(row)!=1) { stop("Error in construct_S_n (B)") }
    if (length(col)!=1) { stop("Error in construct_S_n (C)") }
    return(pred$survival[row,col])
  }
  
}

# Old Cox PH code inside construct_S_n()
if (F) {
  
  if (type=="Cox PH") {

    weights_m1 <- dat$weights * (length(dat$weights)/sum(dat$weights))

    # Fit Cox model
    model <- coxph(fml, data=df, weights=weights_m1)
    if (return_model) { return(model) }
    coeffs <- model$coefficients

    # Get cumulative hazard estimate
    bh <- basehaz(model, centered=FALSE)

    # Pre-calculate H_0 vector
    H_0 <- c()
    for (t in 0:C$t_e) {
      index <- which.min(abs(bh$time-t))
      H_0[t+1] <- bh$hazard[index]
    }

    fnc <- function(t, w, a) {

      if (length(w)!=(length(coeffs)-1)) { stop("Error in construct_S_n (A)") }
      lin <- coeffs[["a"]]*a
      for (i in 1:length(w)) {
        lin <- lin + coeffs[[paste0("w",i)]]*w[i]
      }
      return(exp(-1*H_0[t+1]*exp(lin)))
    }

  }

}

# saving functions (analysis_705)
if (F) {
  
  # Construct/save component functions (SECTION 1)
  fns <- c("Phi_n","Phi_n_inv","S_n","Sc_n","f_aIw_n","f_a_n","g_n","omega_n")
  if (cfg2$run_recomp_1) {
    print(paste("Check 1:",Sys.time()))
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
    print(paste("Check 2:",Sys.time()))
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    print(paste("Check 3:",Sys.time()))
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15) # !!!!! Also do k=0 for cross-validated selection of k
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    print(paste("Check 4:",Sys.time()))
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    print(paste("Check 5:",Sys.time()))
    if (cfg2$save_fns) { save_or_load(fns, cfg2$folder, "save") }
  } else {
    save_or_load(fns, cfg2$folder, "load")
  }
  
  # Construct/save component functions (SECTION 2)
  fns <- c("Gamma_os_n")
  if (cfg2$run_recomp_2) {
    print(paste("Check 6:",Sys.time()))
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    print(paste("Check 7:",Sys.time()))
    if (cfg2$save_fns) { save_or_load(fns, cfg2$folder, "save") }
  } else {
    save_or_load(fns, cfg2$folder, "load")
  }
  
  # Construct/save component functions (SECTION 3)
  fns <- c("Psi_n","gcm","dGCM","theta_n_Gr","theta_n")
  if (cfg2$run_recomp_3) {
    print(paste("Check 8:",Sys.time()))
    Psi_n <- Vectorize(function(x) {
      Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
    })
    gcm <- gcmlcm(x=seq(0,1,C$appx$a), y=Psi_n(seq(0,1,C$appx$a)), type="lcm")
    dGCM <- Vectorize(function(x) {
      index <- which(round(x,5)<=gcm$x.knots)[1]-1
      if (index==0) { index <- 1 }
      return(gcm$slope.knots[index])
    })
    theta_n_Gr <- Vectorize(function(x) { min(dGCM(Phi_n(x)),1) })
    theta_n <- theta_n_Gr
    print(paste("Check 9:",Sys.time()))
    if (cfg2$save_fns) { save_or_load(fns, cfg2$folder, "save") }
  } else {
    save_or_load(fns, cfg2$folder, "load")
  }
  
  # Construct/save component functions (SECTION 4)
  fns <- c("f_aIw_delta1_n","f_a_delta1_n","gamma_n","deriv_theta_n","tau_n")
  if (cfg2$run_recomp_4) {
    print(paste("Check 10:",Sys.time()))
    f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type,
                                        k=15, delta1=TRUE) # !!!!! Also do k=0 for cross-validated selection of k
    f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                    f_aIw_delta1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, vlist$A_grid,
                                 type=params$gamma_type, omega_n, f_aIw_n,
                                 f_a_n, f_a_delta1_n)
    deriv_theta_n <- construct_deriv_theta_n(theta_n, type=params$deriv_type,
                                             dir="decr")
    tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
    print(paste("Check 11:",Sys.time()))
    if (cfg2$save_fns) { save_or_load(fns, cfg2$folder, "save") }
  } else {
    save_or_load(fns, cfg2$folder, "load")
  }
  
  # Construct/save component functions (SECTION 5)
  fns <- c("gcomp")
  if (cfg2$run_recomp_5) {
    print(paste("Check 12:",Sys.time()))
    S_n2 <- construct_S_n(dat, vlist$S_n, type="Cox PH")
    gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n2)
    print(paste("Check 13:",Sys.time()))
    if (cfg2$save_fns) { save_or_load(fns, cfg2$folder, "save") }
  } else {
    save_or_load(fns, cfg2$folder, "load")
  }
  
}

# simlog() function
if (F) {
  
  simlog <- (function() {
    time_st <- Sys.time()
    timestamps <- data.frame(msg="Start", time="0")
    return(function (msg=NA) {
      if (is.na(msg)) {
        return(timestamps)
      } else {
        time_elapsed <- format(Sys.time()-time_st)
        timestamps[nrow(timestamps)+1,] <<- list(msg, time_elapsed)
      }
    })
  })()
  simlog("1")
  simlog()
  
}

# Beta density estimator (weighted)
if (F) {
  
  library(kdensity)
  library(ks)
  
  n <- 500
  shape <- c(0.6,0.6) # c(0.3,0.3)
  sample <- rbeta(n, shape1=shape[1], shape2=shape[2]) %>% sort()
  grid <- seq(0,1,0.01)
  
  wts <- c(
    # rep(1,length(sample))
    rep(0.5,length(sample)/2),
    rep(1.5,length(sample)/2)
  )
  wts <- wts / sum(wts)
  kd_gauss2_ests <- density(
    x = sample,
    kernel = "gaussian",
    from = 0,
    to = 1,
    weights = wts
  )
  kd_gauss2 <- function(x) {
    index <- which.min(abs(kd_gauss2_ests$x - x))
    area <- mean(kd_gauss2_ests$y)
    return(kd_gauss2_ests$y[index]/area)
  }
  kd_gauss <- kdensity(
    x = sample,
    start = "gumbel",
    kernel = "gaussian",
    support = c(0,1)
  )
  # kd_beta <- kdensity( # Weights unavailable
  #   x = sample,
  #   start = "gumbel",
  #   kernel = "beta",
  #   support = c(0,1)
  # )
  kd_beta2_ests <- kde.boundary(
    x = sample,
    boundary.kernel = "beta",
    w = wts,
    xmin = 0,
    xmax = 1
  )
  kd_beta2 <- function(x) {
    len <- length(kd_beta2_ests$eval.points)
    k_x <- kd_beta2_ests$eval.points[2:(len-1)]
    k_dens <- kd_beta2_ests$estimate[2:(len-1)]
    index <- which.min(abs(k_x - x))
    return(k_dens[index])
  }
  
  # Graph of true Beta distribution against estimates
  dens_true <- sapply(grid, function(x) {
    dbeta(x, shape1=shape[1], shape2=shape[2])
  })
  dens_est_gauss <- sapply(grid, kd_gauss)
  dens_est_gauss2 <- sapply(grid, kd_gauss2)
  dens_est_beta <- sapply(grid, kd_beta)
  dens_est_beta2 <- sapply(grid, kd_beta2)
  dat_plot <- data.frame(
    x = rep(grid, 5),
    density = c(dens_true, dens_est_gauss, dens_est_gauss2,
                dens_est_beta, dens_est_beta2),
    which = rep(c("True", "Est: Gaussian", "Est: Gaussian 2",
                  "Est: Beta", "Est: Beta 2"), each=length(grid))
  )
  ggplot(dat_plot, aes(x=x, y=density, color=factor(which))) +
    geom_line() +
    labs(color="Which")
  
  # # Histogram of sample
  # ggplot(data.frame(x=sample), aes(x=x)) + geom_histogram()
  
  # # Check to see if densities are properly normalized (should appx equal 1)
  # mean(dens_est_beta[1:100])
  # # mean(dens_est_beta2[1:100])
  # mean(dens_est_3[1:100])
  # mean(dens_est_gauss[1:100])
  
  
}

# Old htab function
if (F) {
  
  #' Construct a vectorized/memoized hash table
  #'
  #' @param fnc Function to evaluate
  #' @param vals List of values to pre-compute function on; passed to
  #' @param round_args If provided, a vector of length equal to the number of
  #'     columns in vals; all computations will be rounded accordingly.
  #' @param check_dupes NOT YET IMPLEMENTED; Logical; if TRUE, the values in
  #'     vals will not be double-computed if there are duplicates (including
  #'     duplicates after rounding)
  #' @return A vectorized/memoized hash table
  create_htab <- function(fnc, vals, round_args=NA, check_dupes=FALSE) {
    
    rnd <- !(is.na(round_args[1]))
    
    if (rnd && length(vals)!=length(round_args)) {
      stop("length(vals) must equal length(round_args)")
    }
    
    # Create and populate hash table (for function values)
    htab <- new.env()
    for (i in 1:nrow(vals)) {
      # !!!!! Check keys in here for dupes; if (check_dupes) {...}
      row <- as.numeric(vals[i,])
      if (rnd) row <- round(row, round_args)
      key <- rlang::hash(row)
      htab[[key]] <- do.call(fnc, as.list(row))
    }
    rm(row)
    rm(key)
    
    # Create hash table (for vectorized evaluations; populated on the fly)
    htab_v <- new.env()
    
    # Create function
    return(function(...) {
      
      # if (max(sapply(list(...), length))==1) {
      #   if (!rnd) {
      #     args <- as.numeric(list(...))
      #   } else {
      #     args <- round(as.numeric(list(...)), round_args)
      #   }
      #   key <- rlang::hash(args)
      #   val <- htab[[key]]
      #   if (is.null(val)) {
      #     stop(paste0("Value corresponding to arguments (",
      #                 paste0(args, collapse=","),
      #                 ") has not been set"))
      #   }
      #   return(val)
      # } else {
      
      # Memoize vectorized evaluations
      if (!rnd) {
        lst <- list(...)
      } else {
        lst <- lapply(
          X = c(1:length(list(...))),
          FUN = function(i) { round(list(...)[[i]], round_args[i]) }
        )
      }
      hsh <- rlang::hash(lst)
      vec <- htab_v[[hsh]]
      if (!is.null(vec)) {
        return(vec)
      } else {
        vec <- do.call("mapply", c(
          FUN = function(...) {
            key <- rlang::hash(as.numeric(list(...)))
            val <- htab[[key]]
            if (is.null(val)) {
              stop(paste0("Value corresponding to arguments (",
                          paste0(as.numeric(list(...)), collapse=","),
                          ") has not been set"))
            }
            return(val)
          },
          lst
        ))
        htab_v[[hsh]] <- vec
        return(vec)
      }
      # }
    })
  }
  
}

# Precision-weighted edge correction
if (F) {
  
  if (params$edge_corr=="weighted") {
    
    sigma2_Gr <- ( 0.52 * (n_orig^(-1/3)) * tau_n(0) )^2
    sigma2_OS <- sigma2_os_n_est / n_orig

    wt <- sigma2_Gr / (sigma2_OS + sigma2_Gr)
    theta_n_weighted <- wt*theta_os_n_est + (1-wt)*theta_n_Gr(0)

    theta_n <- function(x) {
      max(theta_n_weighted, theta_n_Gr(x))
    }

    gren_ests <- sapply(points, theta_n_Gr)
    gren_points <- sapply(points, function(x) {
      as.numeric(gren_ests>theta_n_weighted)
    })
    gren_points[1] <- 0
    
  }  
}

# Sample split CI estimates (est_curve.R)
if (F) {
  
  # !!!!! Need to check; haven't run this in a while
  
  # Construct data splits (discarding "extra rows" at end)
  m <- params$m
  splits <- matrix(NA, nrow=m, ncol=2)
  split_size <- as.integer(n/m)
  splits[,2] <- (1:m)*split_size
  splits[,1] <- ((1:m)*split_size+1)-split_size
  
  # Construct estimate separately for each data split
  ests <- c()
  ses <- c()
  for (point in points) {
    split_ests <- sapply(c(1:m), function(x) {
      dat_split <- dat_orig[c(splits[x,1]:splits[x,2]),]
      theta_n <- construct_fns(dat_orig=dat_split, return_tau_n=F)$theta_n
      return(theta_n(point))
    })
    
    ests <- c(ests, mean(split_ests))
    ses <- c(ses, sd(split_ests)/sqrt(m))
  }
  
  # Construct CIs
  t_quant <- qt(1-(0.05/2), df=(m-1))
  ci_lo <- ests - t_quant*ses
  ci_hi <- ests + t_quant*ses
  
}

# G-comp estimator (est_curve.R)
if (F) {
  
  if (estimator=="G-comp") {
    
    # !!!!! May need to rewrite some of this
    
    # Prep
    dat <- dat_orig %>% filter(!is.na(a))
    
    # Construct component functions
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
    
    # Compute estimates
    ests <- gcomp_n(points)
    
    # Run bootstrap for SEs
    {
      my_stat <- function(dat_orig,indices) {
        d <- dat_orig[indices,]
        S_n <- construct_S_n(d, vlist$S_n, type=params$S_n_type)
        gcomp_n <- construct_gcomp_n(d, vlist$A_grid, S_n)
        return (gcomp_n(points))
      }
      boot_obj <- boot(data=dat_orig, statistic=my_stat, R=params$boot_reps)
    }
    
    # Parse results object
    # t_quant <- qt(1-(0.05/2), df=(params$boot_reps-1))
    res <- list()
    for (p in 1:length(points)) {
      boot_sd <- sd(boot_obj$t[,p])
      res[[p]] <- list(
        point = points[p],
        est = ests[p],
        ci_lo = ests[p] - 1.96*boot_sd,
        ci_hi = ests[p] + 1.96*boot_sd
      )
    }
    
    return (res)
    
  }
  
}

# Alternative conditional density estimators
if (F) {
  
  # !!!!! Testing: haldensify
  {
    library(haldensify)

    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights

    haldensify_fit <- haldensify(
      A = dat$a,
      W = subset(dat, select=c(w1,w2)),
      wts = weights,
      n_bins = 10, # c(10,25)
      grid_type = "equal_range", # c("equal_range", "equal_mass")
      lambda_seq = exp(seq(-1, -10, length = 50))
      # arguments passed to hal9001::fit_hal()
      # max_degree = 5,
      # smoothness_orders = 0,
      # num_knots = NULL,
      # reduce_basis = 0.05
    )

    grid <- seq(0.01,0.99,0.01)
    new_A <- rep(grid, 4)
    new_W <- data.frame(
      w1 = rep(c(0.2,0.8,0.2,0.8), each=length(grid)),
      w2 = rep(c(0,0,1,1), each=length(grid))
    )
    pred_haldensify <- predict(
      haldensify_fit,
      new_A = new_A,
      new_W = rep(c(1,2,3,4), each=length(grid))
    )

  }
  
  # !!!!! Testing: kernel density estimator
  {
    
    library(np)
    
    bw <- npcdensbw(a~w1+w2, data=dat)
    grid <- seq(0.01,0.99,0.01)
    fhat <- npcdens(
      bws = bw,
      eydat = data.frame(
        a = rep(grid, 4)
      ),
      exdat = data.frame(
        w1 = rep(c(0.2,0.8,0.2,0.8), each=length(grid)),
        w2 = rep(c(0,0,1,1), each=length(grid))
      )
      # newdata = data.frame(
      #   a = rep(grid, 4),
      #   w1 = rep(c(0.2,0.8,0.2,0.8), each=length(grid)),
      #   w2 = rep(c(0,0,1,1), each=length(grid))
      # )
    )
    fhat$condens
    
    # > predict(fhat)[1:10]
    # [1] 0.7324757 1.2702254 1.2157468 1.7770840 1.0554559
    # [6] 0.8422052 1.6739409 1.1410945 1.0494324 1.0189126
    
  }
  
  
}

# Wald comparator test
if (F) {
  
  #' Hypothesis test based on logistic regression
  #' 
  #' @param dat Data returned by generate_data_dr()
  #' @param params Unused
  #' @return Binary; is null rejected (1) or not (0)
  test_wald <- function(dat, alt_type="incr", params) {
    
    model <- glm(y~w1+w2+a, data=dat, family="binomial")
    one_sided_p <- pnorm(summary(model)$coefficients["a",3], lower.tail=F)
    reject <- as.integer(one_sided_p<0.05)
    
    return (reject)
    
  }
  
}

# Old conditional variance function
if (F) {
  
  #' Construct conditional variance estimator function sigma^2_n(Y|A,W)
  #' 
  #' @param mu_n Dataset returned by generate_data()
  #' @param mu2_n Currently only "logistic"
  #' @return Conditional variane estimator function
  construct_sigma2_n <- function(mu_n, mu2_n) {
    
    return(memoise(Vectorize(function(a, w1, w2){
      mu2_n(a,w1,w2) - (mu_n(a,w1,w2))^2
    })))
    
  }
  
}

# Old visualizations (Bias, MSE, Coverage in one graph)
if (F) {
  
  # Export: 6" x 4"
  distr_A_ <- "Beta(1.5+w1,1.5+w2)"
  estimand_ <- "Midpoint" # Midpoint Endpoint
  ggplot(
    filter(p_data, distr_A==distr_A_ & estimand==estimand_),
    aes(x=n, y=value, color=Estimator)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_point() +
    geom_line() +
    # geom_bar(stat="identity", position=position_dodge(),
    #          width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
                  scales=list(y=list(
                    Bias = scale_y_continuous(labels = percent_format()),
                    Coverage = scale_y_continuous(labels = percent_format()),
                    MSE = scale_y_continuous()
                  ))) +
    theme(legend.position="bottom") +
    # scale_color_manual(values=m_colors) +
    labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
         y=NULL, x=NULL, color="Estimator")
  
  
  # !!!!! Temp (export 7.5 x 4.5)
  distr_A_ <- "Beta(1.5+w1,1.5+w2)"
  estimand_ <- "Endpoint" # Midpoint Endpoint
  ggplot(
    filter(p_data, distr_A==distr_A_ & estimand==estimand_),
    aes(x=sampling, y=value, fill=Estimator)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_bar(stat="identity", position=position_dodge(),
             width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
                  scales=list(y=list(
                    Bias = scale_y_continuous(labels = percent_format()),
                    Coverage = scale_y_continuous(labels = percent_format()),
                    MSE = scale_y_continuous()
                  ))) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors) +
    labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
         y=NULL, x=NULL, color="Estimator")
  
}

# Testing of regression estimators
if (F) {
  
  ##########################################.
  ##### TESTING: Regression estimators #####
  ##########################################.
  
  if (FALSE) {
    
    # Set levels here
    n <- 5000
    reg_true <- "Logistic" # Logistic GAM Complex
    sampling <- "iid" # iid two-phase
    
    # Generate data
    L <- list(alpha_3=0.7)
    C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3)
    dat <- generate_data(
      n = n,
      alpha_3 = 0.7,
      distr_A = "Unif(0,1)", # Unif(0,1) Beta(1.5+w1,1.5+w2)
      edge = "none",
      reg_true = reg_true,
      sampling = sampling
    )
    
    # True regression function
    mu_0 <- function(a,w1,w2) {
      if (reg_true=="Logistic") {
        expit(-1.5 + 0.3*w1 + 0.7*w2 + 0.7*a)
      } else if (reg_true=="GAM") {
        expit(-1.5 + 0.3*w1 + 0.7*w2 + 0.7*sqrt(a))
      } else if (reg_true=="Complex") {
        expit(-1.5 + 0.3*sin(2*pi*w1) + 0.7*w2 + 0.7*sqrt(a) + -0.3*w1*w2)
      }
    }
    
    # Construct regression functions
    mu_n_logistic <- construct_mu_n(dat=dat, type="Logistic")
    mu_n_gam <- construct_mu_n(dat=dat, type="GAM")
    mu_n_rf <- construct_mu_n(dat=dat, type="Random forest")
    
    # Generate plot data
    grid <- seq(0,1,0.01)
    mu_models <- c("Truth", "Logistic", "Random forest", "GAM")
    n_models <- length(mu_models)
    len <- length(grid)
    plot_data <- data.frame(
      a = rep(grid, 4*n_models),
      mu = c(
        sapply(grid, function(a) { mu_0(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_0(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_0(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_0(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.8, w2=1) })
      ),
      which = rep(mu_models, each=len*4),
      covariates = rep(c(
        rep("W1=0.2, W2=0",len),
        rep("W1=0.8, W2=0",len),
        rep("W1=0.2, W2=1",len),
        rep("W1=0.8, W2=1",len)
      ), n_models)
    )
    ggplot(plot_data, aes(x=a, y=mu, color=factor(which))) +
      geom_line() +
      facet_wrap(~covariates, ncol=2) +
      labs(color="Estimator", title="Estimation of regression: E[Y|W,A]")
    
  }

}

# Regression estimator
if (F) {
  
  #' Construct regression function mu_n
  #' 
  #' @param dat Dataset returned by generate_data(); accepts either full data or
  #'     truncated data
  #' @param type One of c("Logistic", "GAM", "Random forest")
  #' @param moment If moment=k, the regression E[Y^k|A,W] is estimated
  #' @return Regression function
  construct_mu_n <- function(dat, type, moment=1) {
    
    dat %<>% filter(!is.na(a))
    
    if (moment!=1) {
      dat %<>% mutate(y=y^moment)
    }
    
    if (type=="Logistic") {
      
      model <- glm(
        y~w1+w2+a,
        data = dat,
        family = "binomial",
        weights = wts(dat, scale="mean 1")
      )
      coeffs <- as.numeric(summary(model)$coefficients[,1])
      
      return(Vectorize(function(a, w1, w2){
        expit( coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2 + coeffs[4]*a )
      }))
      
    }
    
    if (type=="GAM") {
      
      k <- 4
      knots <- seq(0,1,1/(k+1))[2:(k+1)]
      ns_basis <- ns(dat$a, knots=knots, Boundary.knots=c(0,1))
      formula <- "y ~ w1 + w2"
      for (i in 1:(k+1)) {
        dat[paste0("b_",i)] <- ns_basis[,i]
        formula <- paste0(formula," + b_",i)
      }
      
      dat$wts <- wts(dat, scale="mean 1")
      model <- glm(
        formula,
        data = dat,
        family = "binomial",
        weights = wts
      )
      coeffs <- as.numeric(summary(model)$coefficients[,1])
      
      # Construct natural spline basis in advance over grid
      grid <- seq(0,1,0.01)
      nsb_grid <- ns(grid, knots=knots, Boundary.knots=c(0,1))
      
      return(Vectorize(function(a, w1, w2){
        index <- which.min(abs(a-grid))
        nsb_a <- as.numeric(nsb_grid[index,])
        lin_pred <- coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2
        for (i in 1:(k+1)) {
          lin_pred <- lin_pred + coeffs[i+3]*nsb_a[i]
        }
        return(expit(lin_pred))
      }))
      
    }
    
    if (type=="Random forest") {
      
      model <- ranger(
        y~w1+w2+a,
        data = dat,
        num.trees = 500,
        case.weights = wts(dat, scale="mean 1")
      )
      
      return(memoise(Vectorize(function(a, w1, w2){
        predict(model, data.frame(a=a, w1=w1, w2=w2))$predictions
      })))
    }
    
  }
  
}



# Chernoff realizations
if (F) {
  
  library(twostageTE)
  data(chernoff_realizations)
  
}

# Old deriv_theta_n estimator
if (F) {
  
  fnc <- function(a) {
    
    # Set derivative appx x-coordinates
    width <- 0.1
    p1 <- a - width/2
    p2 <- a + width/2
    if (p1<0) {
      p2 <- p2 - p1
      p1 <- 0
    }
    if (p2>1) {
      p1 <- p1 - p2 + 1
      p2 <- 1
    }
    c(p1,p2)
    
    return( (gcomp_n(p2)-gcomp_n(p1))/width )
    
  }
  
}

# generate_data old regression code
if (F) {
  
  #' @param reg_true True functional form of the regression; one of the following:
  #'     - "Logistic": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*A)
  #'     - "GAM": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*sqrt(A))
  #'     - "Complex": E[Y|W,A]=expit(a0+a1*sin(2*pi*W1)+a2*W2+a3*sqrt(A)+a4*W1*W2)
  
  # Compute true regression function (i.e. P(Y=1|W,A))
  if (reg_true=="Logistic") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a)
  } else if (reg_true=="GAM") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*sqrt(a))
  } else if (reg_true=="Complex") {
    probs <- expit(C$alpha_0 + C$alpha_1*sin(2*pi*w1) + C$alpha_2*w2 +
                     alpha_3*sqrt(a) + C$alpha_4*w1*w2)
  }
  
  # Sample outcome
  y <- rbinom(n, size=1, prob=probs)
  
}

# Old GAM code (using mgcv; slow predictions)
if (F) {
  
  model <- gam(
    y~w1+w2+s(a, bs="cr"),
    # y~w1+w2+s(a, bs="tp"),
    # y~w1+w2+s(a, k=5, bs="cr", fx=TRUE),
    data = dat,
    family = "binomial",
    weights = wts(dat, scale="mean 1")
  )

  return(memoise(Vectorize(function(a, w1, w2){
    predict.gam(model, list(a=a, w1=w1, w2=w2), type="response")
  })))
  
}

# Old construct_f_aIw_n function contents
if (F) {
  
  #' @param type One of c("KDE (Beta)", "Beta"). The former fits a KDE with Beta
  #'     kernels (Chen 1999). The latter fits Beta density via MLE.
  #' @notes
  #'   - We trim the first and last evaluation points (i.e. zero and one) because
  #'     the kde.boundary() function estimates the density as zero there; the
  #'     neighboring points are used as the estimates instead
  
  if (type=="simple") {
    
    dat_0 <- dplyr::filter(dat,w2==0 & !is.na(a))
    dat_1 <- dplyr::filter(dat,w2==1 & !is.na(a))
    
    # print(paste("sim_uid:", L$sim_uid))
    # print(paste("dat contains", nrow(dat), "rows.")) # !!!!!
    # print(paste("dat contains", nrow(filter(dat, !is.na(a))), "complete rows.")) # !!!!!
    # print(paste("dat_0 contains", nrow(dat_0), "rows.")) # !!!!!
    # print(paste("dat_1 contains", nrow(dat_1), "rows.")) # !!!!!
    
    # kd_0 <- kdensity(
    #   x = dat_0$a,
    #   start = "gumbel",
    #   kernel = "gaussian"
    # )
    kde_0 <- density(
      x = dat_0$a,
      kernel = "gaussian",
      weights = wts(dat_0, scale="sum 1"),
      from = 0,
      to = 1
    )
    kd_0 <- function(x) {
      index <- which.min(abs(kde_0$x - x))
      area <- mean(kde_0$y)
      return(kde_0$y[index]/area)
    }
    
    # kd_1 <- kdensity(
    #   x = dat_1$a,
    #   start = "gumbel",
    #   kernel = "gaussian"
    # )
    kde_1 <- density(
      x = dat_1$a,
      kernel = "gaussian",
      weights = wts(dat_1, scale="sum 1"),
      from = 0,
      to = 1
    )
    kd_1 <- function(x) {
      index <- which.min(abs(kde_1$x - x))
      area <- mean(kde_1$y)
      return(kde_1$y[index]/area)
    }
    
    f_aIw_n <- function(a,w1,w2) {
      if (w2==0) { return(kd_0(a)) }
      if (w2==1) { return(kd_1(a)) }
    }
    
    return(memoise(Vectorize(f_aIw_n)))
    
  }
  
}

# Old construct_f_a_n function contents
if (F) {
  
  # Run weighted KDE
  dat %<>% filter(!is.na(a))
  
  # KDE (Beta kernels)
  if (type=="KDE (Beta)") {
    
    kde <- kde.boundary(
      x = dat$a,
      boundary.kernel = "beta",
      w = wts(dat, scale="sum 1"),
      xmin = 0,
      xmax = 1
    )
    f_a_n <- function(x) {
      len <- length(kde$eval.points)
      k_x <- kde$eval.points[2:(len-1)] # Trimming off first and last point
      k_dens <- kde$estimate[2:(len-1)] # Trimming off first and last point
      index <- which.min(abs(k_x - x))
      return(k_dens[index])
    }
    
  }
  
  if (type=="Beta") {
    
    # # !!!!!
    # dat <- list()
    # dat$a <- rbeta(100, shape1=0.3, shape2=0.5)
    # # !!!!!
    
    # !!!!! Testing
    {
      # Beta(1.5+w1,1.5+w2)
      # Generate true marginal distribution of A
      n <- 10000
      beta_samp_1 <- rbeta(n, shape1=0.9, shape2=1.1)
      beta_samp_2 <- rbeta(n, shape1=0.9, shape2=1.5)
      beta_samp <- c(beta_samp_1, beta_samp_2)
      ggplot(data.frame(x=beta_samp_1), aes(x=x)) + geom_histogram(bins=50)
      ggplot(data.frame(x=beta_samp_2), aes(x=x)) + geom_histogram(bins=50)
      ggplot(data.frame(x=beta_samp), aes(x=x)) + geom_histogram(bins=50)
      
    }
    
    dat %<>% filter(!is.na(a))
    n <- length(dat$a)
    
    # !!!!! Comparison
    Rfast::beta.mle(dat$a)
    
    # weights <- wts(dat, scale="sum 1") # !!!!!
    
    # Set the objective function (weighted likelihood)
    #   par[1] is alpha and par[2] is beta
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        loglik <- dbeta(dat$a[i], shape1=par[1], shape2=par[2], log=TRUE)
        wt <- 1 # !!!!! Testing
        # wt <- weights[i]
        return(-1*loglik*wt)
      }))
      
      return(sum_loglik)
      
    }
    optim(par=c(alpha=1, beta=1), fn=wlik)

  }
  
}

# Old GAM specification function
if (F) {
  
  if (mono_form=="identity") {
    mono_f <- function(x) {x}
  } else if (mono_form=="square") {
    mono_f <- function(x) {x^2}
  } else if (mono_form=="sqrt") {
    mono_f <- function(x) {sqrt(x)}
  } else if (mono_form=="step_0.2") {
    mono_f <- function(x) {as.integer(x>0.2)}
  } else if (mono_form=="step_0.8") {
    mono_f <- function(x) {as.integer(x>0.8)}
  } else {
    stop("mono_form incorrectly specified")
  }
  
}

# Old kernel density function
if (F) {
  
  # kde <- density(
  #   x = dat$a,
  #   kernel = "gaussian",
  #   weights = wts(dat),
  #   from = 0,
  #   to = 1
  # )
  # f_a_n <- function(x) {
  #   index <- which.min(abs(kde$x - x))
  #   area <- mean(kde$y)
  #   return(kde$y[index]/area)
  # }
  
  # f_a_n <- kdensity(
  #   x = dat$a,
  #   start = "gumbel",
  #   kernel = "beta", # gaussian
  #   support = c(0,1)
  # )
  
}

# Debugging code in construct_fns()
if (F) {
  mu_n_iid <- construct_mu_n(dat_iid, type="logistic")
  f_aIw_n_iid <- construct_f_aIw_n(dat_iid, type="simple")
  f_a_n_iid <- construct_f_a_n(dat_iid, type=NULL)
  g_n_iid <- construct_g_n(f_aIw_n_iid, f_a_n_iid)
  mu_n_tp <- construct_mu_n(dat_tp, type="logistic")
  f_aIw_n_tp <- construct_f_aIw_n(dat_tp, type="simple")
  f_a_n_tp <- construct_f_a_n(dat_tp, type=NULL)
  g_n_tp <- construct_g_n(f_aIw_n_tp, f_a_n_tp)
  Gamma_n_iid <- construct_Gamma_n(dat_iid, mu_n_iid, g_n_iid)
  Gamma_n_tp <- construct_Gamma_n(dat_tp, mu_n_tp, g_n_tp)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=function(a) {
    Gamma_n_iid(a)
  }) + ylim(0,0.5)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=function(a) {
    Gamma_n_tp(a)
  }) + ylim(0,0.5)
}



# Constructor functions for estimator componetns (hypothesis testing)
{
  
  # Set up component functions
  {
    # Set up G constructor
    construct_G <- function(params, dat) {
      if (params$G == "identity") {
        return(
          function(x) { x }
        )
      }
      if (params$G == "marginal") {
        ecdf_G <- ecdf(dat$a)
        return(
          function(x) { ecdf_G(x) }
        )
      }
    }
    
    # lambda function (returns a constant)
    lambda <- function(k, G, dat) {
      return( mean((G(dat$a, dat))^k) )
    }
    
    # Construct influence function
    # !!!!! Vectorize? Memoise?
    # !!!!! Also option to use one-step or plug-in ?????
    construct_infl_fn <- function(sub_x, x) {
      
      # !!!!!
      
      return(999)
      
    }
    
  }
    
}



# Test empirical CDF and empirical quantile functions
{
  dat <- list(a=c(0.6,0.8,0.9))
  Phi_n <- ecdf(dat$a)
  ggplot(data.frame(x=grid), aes(x=x)) +
    stat_function(fun = Phi_n) + ylim(0,1)
  Phi_n_inv <- function(x) {
    qemp(p=x, obs=dat$a, discrete=T)
  }
  ggplot(data.frame(x=grid), aes(x=x)) +
    stat_function(fun = Phi_n_inv) + ylim(0,1)
}



# Old functions
{
  #' Integral of expit function
  #' 
  #' @param x Numeric input
  #' @return Numeric output
  int_expit <- function(x) {log(1+exp(x))}
}



#' Hypothesis testing approach 2: slope (!!!!! OLD CODE !!!!!)
#'
#' @param dat Data returned by generate_data_dr()
#' @param params A list; `est` is the estimator used for Theta_hat; one of
#'     c("glm", "sm spline")
#' @return Binary; is null rejected (1) or not (0)
test_app2_dr <- function(dat, params) {

  Theta_hat_constr <- function(dat, subtype) {

    if (subtype=="glm") {

      # Run model and extract coefficients
      model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
      coeff <- summary(model)$coefficients
      alpha_0_hat <- coeff["(Intercept)",1]
      alpha_1_hat <- coeff["bmi",1]
      alpha_2_hat <- coeff["sex",1]
      alpha_3_hat <- coeff["antib",1]

      # Create Theta_hat function
      Theta_hat <- Vectorize(function(x) {

        t_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            log( (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] +
                  alpha_2_hat*r[["sex"]] + alpha_3_hat*x)) /
                (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]]))
            )
          }
        )

        return ((alpha_3_hat^-1)*mean(t_i))

      })

      return (Theta_hat)

    }

    if (subtype=="ss") {

      # Run model and extract coefficients
      model <- gam(
        infected ~ bmi + sex + s(antib, fx=FALSE, bs="cr", m=2, pc=0),
        data = dat,
        family = "binomial"
      )

      coeff <- model$coefficients
      alpha_0_hat <- coeff[["(Intercept)"]]
      alpha_1_hat <- coeff[["bmi"]]
      alpha_2_hat <- coeff[["sex"]]
      spline_vals <- as.numeric(predict(
        model,
        newdata = list(bmi=rep(25,101), sex=rep(0,101), antib=seq(0,1,0.01)),
        type = "terms"
      )[,3])

      # Construct theta_hat from GAM formula
      theta_hat <- Vectorize(function(x) {
        E_hat_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            expit(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]] +
                  spline_vals[1:(round(100*round(x,2)+1,0))])
          }
        )
        return (mean(E_hat_i))
      })

      # Construct Theta_hat by integrating theta_hat
      # Approximating integral using a Riemann sum with ten intervals
      Theta_hat <- Vectorize(
        function (a) { a * mean(theta_hat(seq(a/10,a,a/10))) }
      )

      return (Theta_hat)

    }

  }

  # Construct Theta_hat function
  Theta_hat <- Theta_hat_constr(dat, params$subtype)

  # Calculate value of test statistic
  x <- dat$antib
  mu_2n <- mean(x^2)
  mu_3n <- mean(x^3)
  beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))

  # Define the statistic to bootstrap
  bootstat <- function(dat,indices) {
    d <- dat[indices,]
    Theta_hat <- Theta_hat_constr(d, params$subtype)
    x <- dat$antib
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    return (mean((mu_2n*x^2 - mu_3n*x)*(Theta_hat(x)-x)))
  }

  # Run bootstrap
  boot_obj <- boot(data=dat, statistic=bootstat, R=100) # 1000

  # Calculate critical value
  # crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
  # !!!!! This is a one-sided test; either make this two-sided or do the Wald test above as a one-sided test
  crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))

  return(as.numeric(crit_val>0))

}

