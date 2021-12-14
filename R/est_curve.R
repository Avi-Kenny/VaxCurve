#' Estimate the causal dose-response curve at a point
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param estimator Which estimator to use; currently only "Grenander"
#' @param params A list, containing the following:
#'   - `S_n_type` Type of survival function estimator; corresponds to
#'     construct_S_n(type=...)
#'   - `g_n_type` Type of conditional density ratio estimator; corresponds to
#'     construct_f_aIw_n(type...)
#'   - `deriv_type` Type of derivative estimator; corresponds to
#'     construct_deriv_theta_n(type=...)
#'   - `boot_reps` Used for G-comp; number of bootstrap replicates
#'   - `ci_type` One of c("regular", "logit", "sample split", "none"). "regular"
#'     is the standard approach. "logit" transforms the CIs so that the bounds
#'     are in [0,1]. "sample split" is the Banerjee method (and also returns a
#'     different estimator).
#'   - `cf_folds` Number of cross-fitting folds; 1 for no cross-fitting
#'   - `m` If params$ci_type=="sample split", the number of splits
#'   - `edge_corr` One of c("none", "spread", "point", "weighted", "min").
#'     "point" uses the root-n estimator to replace the Grenander-based
#'     estimator only at the leftmost point. "weighted" and "min" both use the
#'     precision-weighted estimator to adjust both the leftmost point and the
#'     rest of the curve. "spread" adds a small amount of noise to the edge
#'     points
#' @param points A vector representing the points to estimate
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @param return_extra A character vector of additional components to return
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)
est_curve <- function(dat_orig, estimator, params, points, dir="decr",
                      return_extra=NULL) {
  
  if (estimator=="Grenander") {
    
    if (params$edge_corr=="spread") {
      noise <- runif(length(dat_orig$a))*0.05
      dat_orig$a <- ifelse(dat_orig$a==0, noise, dat_orig$a)
    }
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat_orig, C$appx) # !!!!! vlist <- create_val_list(dat, C$appx)
    
    # Construct regular Gamma_0 estimator
    if (params$cf_folds==1) {
      
      # Construct component functions
      Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
      Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
      S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
      Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
      f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                   type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
      g_n <- construct_g_n(f_aIw_n, f_a_n)
      omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                   type=params$omega_n_type)
      Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
      
      # Construct one-step edge estimator
      if (params$edge_corr!="none") {
        pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
        theta_os_n_est <- theta_os_n(dat, pi_n, S_n, omega_n)
        sigma2_os_n_est <- sigma2_os_n(dat, pi_n, S_n, omega_n, theta_os_n_est)
      }
      
    }
    
    # Construct cross-fitted Gamma_0 estimator
    if (params$cf_folds>1) {
      
      Gamma_os_n <- construct_Gamma_cf(dat_orig, params, vlist)
      
      # Recompute functions on full dataset
      Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
      Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
      
      # Construct cross-fitted one-step edge estimator
      if (params$edge_corr!="none") {
        # !!!!! TO DO
      }
      
    }
    
    # Construct Psi_n
    if (dir=="incr") {
      Psi_n <- Vectorize(function(x) {
        Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    } else {
      Psi_n <- Vectorize(function(x) {
        -1 * Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
      })
    }
    
    # Compute GCM and extract its derivative
    gcm <- gcmlcm(
      x = seq(0,1,C$appx$a),
      y = Psi_n(seq(0,1,C$appx$a)),
      type = "gcm"
    )
    dGCM <- Vectorize(function(x) {
      # The round deals with a floating point issue
      index <- which(round(x,5)<=gcm$x.knots)[1]-1
      if (index==0) { index <- 1 }
      return(gcm$slope.knots[index])
    })
    
    # Construct Grenander-based theta_n
    if (dir=="incr") {
      theta_n_Gr <- Vectorize(function(x) { min(max(dGCM(Phi_n(x)),0),1) })
    } else {
      theta_n_Gr <- Vectorize(function(x) { min(max(-1 * dGCM(Phi_n(x)),0),1) })
    }
    
    # Recompute functions on full dataset
    if (params$cf_folds>1) {
      S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
      Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
      f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                   type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
      omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                   type=params$omega_n_type)
    }
    
    # Compute new functions
    f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                        type=params$g_n_type, k=15, delta1=TRUE)
    f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                    f_aIw_delta1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, vlist$A_grid,
                                 type=params$gamma_type, omega_n, f_aIw_n,
                                 f_a_n, f_a_delta1_n)
    
    # Edge correction
    if (params$edge_corr %in% c("none", "spread")) {
      
      theta_n <- theta_n_Gr
      
    } else if (params$edge_corr=="point") {
      
      theta_n <- function(x) {
        if(x==0) {
          theta_os_n_est
        } else {
          theta_n_Gr(x)
        }
      }
      
    } else if (params$edge_corr=="min") {
      
      theta_n <- Vectorize(function(x) {
        if (dir=="decr") {
          min(theta_os_n_est, theta_n_Gr(x))
        } else {
          max(theta_os_n_est, theta_n_Gr(x))
        }
      })
      
      gren_ests <- theta_n_Gr(points)
      gren_points <- sapply(c(1:length(points)), function(i) {
        as.numeric(gren_ests[i]>theta_os_n_est)
      })
      gren_points[1] <- 0
      
    }
    
    # Generate estimates for each point
    ests <- theta_n(points)
    
    # Construct variance scale factor
    deriv_theta_n <- construct_deriv_theta_n(theta_n, type=params$deriv_type,
                                             dir=dir)
    tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
    
    # Generate confidence limits
    if (params$ci_type=="none") {
      
      ci_lo <- rep(0, length(ests))
      ci_hi <- rep(0, length(ests))
      
    } else {
      
      # Generate variance scale factor for each point
      tau_ns <- tau_n(points)
      
      # Construct CIs
      # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
      # The Normal approximation would use qnorm(0.975, sd=0.52) instead
      qnt <- 1.00
      # qnt <- qnorm(0.975, sd=0.52)
      n_orig <- length(dat_orig$delta)
      if (params$ci_type=="regular") {
        ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
        ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
      } else if (params$ci_type=="logit") {
        ci_lo <- expit(
          logit(ests) - (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
        )
        ci_hi <- expit(
          logit(ests) + (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
        )
      } else if (params$ci_type=="trunc") {
        ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
        ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
        ci_lo %<>% pmax(0) %>% pmin(1)
        ci_hi %<>% pmax(0) %>% pmin(1)
      }
      
      # Edge correction
      if (params$edge_corr=="point") {
        ci_lo[1] <- ests[1] - 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_hi[1] <- ests[1] + 1.96*sqrt(sigma2_os_n_est/n_orig)
      } else if (params$edge_corr %in% c("weighted","min")) {
        ci_lo2 <- ests - 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_hi2 <- ests + 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_lo <- (1-gren_points)*pmin(ci_lo,ci_lo2) + gren_points*ci_lo
        ci_hi <- (1-gren_points)*pmax(ci_hi,ci_hi2) + gren_points*ci_hi
      }
      
    }
    
    # # Add extra return data
    # res[["Phi_n"]] <- Phi_n
    # res[["Gamma_os_n"]] <- Gamma_os_n
    
  }
  
  if (estimator=="Qbins") {
    
    # !!!!! Handle case in which marginal distribution has mass
    # !!!!! Implement cross-fitted version
    
    # Construct bin cutoffs
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
    cutoffs <- Phi_n_inv(seq(0,1,length.out=params$n_bins+1))
    
    # Function to transform A values to categorical bins values
    transform_a <- function(a) {
      cut(a, breaks=cutoffs, right=F, include.lowest=T)
    }
    
    # Create vlist
    vlist <- create_val_list(dat, C$appx, factor_A=unique(transform_a(dat$a)))
    
    # Construct f_aIw_n BEFORE transforming A values
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                 type=params$g_n_type, k=15)
    
    # Transform A values
    dat$a <- transform_a(dat$a)
    
    # !!!!! Check to see if these need to be modified to handle factor A (including different S_n types)
    # Note: S_n and Sc_n will not work with type="true"
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                 type=params$omega_n_type)
    pi_n <- construct_pi_n(dat, vlist$W_grid, type="generalized",
                           f_aIw_n=f_aIw_n, cutoffs=cutoffs)
    
    # Generate estimates and standard deviations for each point
    ests <- sapply(c(1:length(points)), function(i) {
      a_binned <- transform_a(points[i])
      return(theta_os_n(dat, pi_n, S_n, omega_n, val=a_binned))
    })
    sigma2s <- sapply(c(1:length(points)), function(i) {
      a_binned <- transform_a(points[i])
      return(sigma2_os_n(dat, pi_n, S_n, omega_n, ests[i], val=a_binned))
    })
    
    # Construct CIs
    n_orig <- length(dat_orig$delta)
    ci_lo <- ests - 1.96*sqrt(sigma2s/n_orig)
    ci_hi <- ests + 1.96*sqrt(sigma2s/n_orig)
    # !!!!! Deal with CI truncation
    
  }
  
  if (estimator=="gcomp") {
    
    # For debugging purposes related to bias
    # !!!!! This does not give confidence intervals
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat, C$appx)
    
    # Construct component functions
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$A_grid, S_n)
    
    # Compute estimates and dummy CIs
    ests <- gcomp_n(points)
    ci_lo <- rep(0, length(ests))
    ci_hi <- rep(0, length(ests))
    
    # Add extra return data
    # res[["ex_S_n"]] <- S_n(C$t_e, w=c(0.5,1), a=0.5)
    # res[["ex_gamma_n"]] <- 0
    # res[["ex_deriv_theta_n"]] <- 0
    # res[["ex_tau_n"]] <- 0
    
  }
  
  # Parse and return results
  res <- list(point=points, est=ests, ci_lo=ci_lo, ci_hi=ci_hi)
  if ("gcomp" %in% return_extra) {
    S_n2 <- construct_S_n(dat, vlist$S_n, type="Cox PH")
    res$gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n2)
  }
  if ("f_a_n" %in% return_extra) { res$f_a_n <- f_a_n }
  if ("Phi_n_inv" %in% return_extra) { res$Phi_n_inv <- Phi_n_inv }
  
  return(res)
  
}
