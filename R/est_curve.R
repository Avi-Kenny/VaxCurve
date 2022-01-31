#' Estimate the causal dose-response curve at a point
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param estimator Which estimator to use; currently only "Grenander"
#' @param params A list, containing the following:
#'   - `S_n_type` Type of survival function estimator; corresponds to
#'     construct_S_n()
#'   - `g_n_type` Type of conditional density ratio estimator; corresponds to
#'     construct_f_aIw_n()
#'   - `ecdf_type` Type of CDF estimator; corresponds to construct_Phi_n()
#'   - `deriv_type` Type of derivative estimator; corresponds to
#'     construct_deriv_theta_n()
#'   - `gamma_type` Type of nuisance estimator; corresponds to
#'     construct_gamma_n()
#'   - `omega_n_type` Type of nuisance estimator; corresponds to
#'     construct_omega_n()
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
#'   - `marg` One of c("Gamma", "Theta"); whether or not to transform by the
#'     marginal distribution of A
#' @param points A vector representing the points at which estimates and CIs
#'     should be calculated. This should be a unique increasing sequence.
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @param return_extra A character vector of additional components to return
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)
est_curve <- function(dat_orig, estimator, params, points, dir="decr",
                      return_extra=NULL) {
  
  # Set default params
  .default_params <- list(
    S_n_type="Super Learner", g_n_type="binning", deriv_type="m-spline",
    ecdf_type="linear (mid)", gamma_type="kernel", omega_n_type="estimated",
    boot_reps=1000, ci_type="trunc", cf_folds=1, m=5, edge_corr="none",
    marg="Theta", lod_shift="none"
  )
  for (i in c(1:length(.default_params))) {
    if (is.null(params[[names(.default_params)[i]]])) {
      params[[names(.default_params)[i]]] <- .default_params[[i]]
    }
  }
  
  if (estimator=="Grenander") {
    
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
    
    # Rescale A to lie in [0,1]
    # !!!!! Functionize and refactor w/ test_2.R
    a_lims <- c(min(dat_orig$a,na.rm=T),max(dat_orig$a,na.rm=T))
    a_shift <- -1 * a_lims[1]
    a_scale <- 1/(a_lims[2]-a_lims[1])
    dat_orig$a <- (dat_orig$a+a_shift)*a_scale
    
    # Round values
    # !!!!! Functionize and refactor w/ test_2.R
    dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
    for (i in c(1:length(dat_orig$w))) {
      rnd <- 8
      tol <- C$appx$w_tol
      n_unique <- tol + 1
      while(n_unique>tol) {
        rnd <- rnd - 1
        n_unique <- length(unique(round(dat_orig$w[,i],rnd)))
      }
      dat_orig$w[,i] <- round(dat_orig$w[,i],rnd)
    }
    
    # Obtain minimum value (excluding edge point mass)
    if (params$edge_corr=="min") {
      a_min2 <- min(dat_orig$a[dat_orig$a!=0],na.rm=T)
    }
    
    # Rescale points and remove points outside the range of A
    points_orig <- points
    na_head <- sum(points<a_min)
    points <- round((points+a_shift)*a_scale, -log10(C$appx$a))
    na_tail <- sum(points>1)
    if (na_head>0) {
      points <- points[-c(1:na_head)]
    }
    if (na_tail>0) {
      points <- points[-c((length(points)-na_tail+1):length(points))]
    }
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat_orig, C$appx)
    
    # Construct regular Gamma_0 estimator
    if (params$cf_folds==1) {
      
      # !!!!! New G_n
      dat2 <- ss(dat, which(dat$a!=0))
      G_n <- construct_Phi_n(dat2, type=params$ecdf_type)
      G_n_inv <- construct_Phi_n(dat2, which="inverse", type=params$ecdf_type)
      n_orig <- length(dat_orig$delta)
      z_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$a!=0))
      # !!!!!
      
      # Construct component functions
      Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
      Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
      S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
      Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
      
      # !!!!! New G_n
      gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$A_grid, S_n)
      alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, z_n, vals=NA)
      eta_ss_n <- construct_eta_ss_n(dat, S_n, z_n, vals=NA)
      # !!!!!
      
      f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                   type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
      omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                                   type=params$omega_n_type)
      if (params$marg=="Theta") {
        etastar_n <- construct_etastar_n(S_n)
        Theta_os_n <- construct_Theta_os_n(dat, vlist$A_grid, omega_n,
                                           f_aIw_n, etastar_n)
      } else if (params$marg=="Gamma") {
        g_n <- construct_g_n(f_aIw_n, f_a_n)
        Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
      } else if (params$marg=="Gamma_star") {
        g_n_star <- construct_g_n_star(f_aIw_n, f_a_n, z_n)
        Gamma_os_n_star <- construct_Gamma_os_n_star(dat, omega_n, g_n_star,
                                                     eta_ss_n, z_n, gcomp_n,
                                                     alpha_star_n, vals=NA)
      } else if (params$marg=="Gamma_star2") {
        g_n_star <- construct_g_n_star(f_aIw_n, f_a_n, z_n)
        q_n <- construct_q_n(dat, dat_orig, type="Super Learner", omega_n,
                             g_n_star, z_n, gcomp_n, alpha_star_n, vals=NA)
        Gamma_os_n_star <- construct_Gamma_os_n_star2(dat, dat_orig, omega_n,
                                                       g_n_star, eta_ss_n, z_n,
                                                       q_n, gcomp_n,
                                                       alpha_star_n, vals=NA)
      } else {
        stop("`params$marg` must be one of c('Theta', 'Gamma', 'Gamma_star'")
      }
      
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
    if (params$marg=="Gamma") {
      if (dir=="incr") {
        Psi_n <- Vectorize(function(x) {
          Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
        })
      } else {
        Psi_n <- Vectorize(function(x) {
          -1 * Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
        })
      }
    } else if (params$marg=="Theta") {
      if (dir=="incr") {
        Psi_n <- Theta_os_n
      } else {
        Psi_n <- Vectorize(function(x) { -1 * Theta_os_n(x) })
      }
    # } else if (params$marg=="Gamma_star") {
    } else if (params$marg %in% c("Gamma_star", "Gamma_star2")) {
      if (dir=="incr") {
        Psi_n <- Vectorize(function(x) {
          Gamma_os_n_star(round(G_n_inv(x), -log10(C$appx$a)))
        })
      } else {
        Psi_n <- Vectorize(function(x) {
          -1 * Gamma_os_n_star(round(G_n_inv(x), -log10(C$appx$a)))
        })
      }
    }
    
    # Compute GCM and extract its derivative
    grid <- round(seq(0,1,C$appx$a),-log10(C$appx$a))
    gcm <- gcmlcm(x=grid, y=Psi_n(grid), type="gcm")
    dGCM <- approxfun(
      x = gcm$x.knots[-length(gcm$x.knots)],
      y = gcm$slope.knots,
      method = "constant",
      rule = 2,
      f = 0
    )
    
    # Construct Grenander-based theta_n
    if (params$marg=="Gamma") {
      if (dir=="incr") {
        theta_n_Gr <- Vectorize(function(x) { min(max(dGCM(Phi_n(x)),0),1) })
      } else {
        theta_n_Gr <- Vectorize(function(x) { min(max(-1 * dGCM(Phi_n(x)),0),1) })
      }
    } else if (params$marg=="Theta") {
      if (dir=="incr") {
        theta_n_Gr <- Vectorize(function(x) { min(max(dGCM(x),0),1) })
      } else {
        theta_n_Gr <- Vectorize(function(x) { min(max(-1*dGCM(x),0),1) })
      }
    # } else if (params$marg=="Gamma_star") {
    } else if (params$marg %in% c("Gamma_star", "Gamma_star2")) {
      if (dir=="incr") {
        # !!!!! COnsolidate this by renaming G_n to Phi_n and G_n_inv to Phi_n_inv ?????
        theta_n_Gr <- Vectorize(function(x) { min(max(dGCM(G_n(x)),0),1) })
      } else {
        theta_n_Gr <- Vectorize(function(x) { min(max(-1 * dGCM(G_n(x)),0),1) })
      }
    }
    
    # Recompute functions on full dataset
    if (params$cf_folds>1) {
      # !!!!! Update to incorporate params$marg %in% c("Theta", "Gamma_star")
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
        if(x==0 || x<a_min2) {
          theta_os_n_est
        } else {
          if (dir=="incr") {
            max(theta_os_n_est, theta_n_Gr(x))
          } else {
            min(theta_os_n_est, theta_n_Gr(x))
          }
        }
      })
      
      gren_points <- sapply(c(1:length(points)), function(i) {
        if (dir=="incr") {
          as.numeric(theta_n(points[i])>theta_os_n_est)
        } else {
          as.numeric(theta_n(points[i])<theta_os_n_est)
        }
      })
      
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
        # ci_lo <- (1-gren_points)*pmin(ci_lo,ci_lo2) + gren_points*ci_lo
        # ci_hi <- (1-gren_points)*pmax(ci_hi,ci_hi2) + gren_points*ci_hi
        ci_lo <- (1-gren_points)*ci_lo2 + replace_na(gren_points*ci_lo, 0)
        ci_hi <- (1-gren_points)*ci_hi2 + replace_na(gren_points*ci_hi, 0)
      }
      
    }
    
    # # Add extra return data
    # res[["Phi_n"]] <- Phi_n
    # res[["Gamma_os_n"]] <- Gamma_os_n
    
  }
  
  if (estimator=="Qbins") {
    
    # !!!!! Handle case in which marginal distribution has mass
    # !!!!! Implement cross-fitted version
    # !!!!! Implement isotonized version
    # !!!!! Adapt for `points` as above; make sure NA calues can be handled
    
    # Construct bin cutoffs
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
    cutoffs <- Phi_n_inv(seq(0,1,length.out=params$n_bins+1))
    
    # Function to transform A values to categorical bins values
    transform_a <- function(a) {
      cut(a, breaks=cutoffs, right=F, include.lowest=T)
    }
    
    # Create vlist
    vlist <- create_val_list(dat_orig, C$appx,
                             factor_A=unique(transform_a(dat$a)))
    
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
    vlist <- create_val_list(dat_orig, C$appx)
    
    # Construct component functions
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$A_grid, S_n)
    
    # Compute estimates and dummy CIs
    ests <- gcomp_n(points)
    ci_lo <- rep(0, length(ests))
    ci_hi <- rep(0, length(ests))
    
    # Add extra return data
    # res[["ex_S_n"]] <- S_n(C$t_e, w=c(0.5,1), a=0.5)
    
  }
  
  # Parse and return results
  res <- list(
    point = points_orig,
    est = c(rep(NA,na_head), ests, rep(NA,na_tail)),
    ci_lo = c(rep(NA,na_head), ci_lo, rep(NA,na_tail)),
    ci_hi = c(rep(NA,na_head), ci_hi, rep(NA,na_tail))
  )
  if ("gcomp" %in% return_extra) {
    S_n2 <- construct_S_n(dat, vlist$S_n, type="Cox PH")
    res$gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n2)
  }
  fns_extra <- c("f_a_n", "gamma_n", "deriv_theta_n", "Phi_n_inv", "Theta_os_n",
                 "Psi_n", "omega_n", "f_aIw_n", "etastar_n", "S_n", "gcm",
                 "dGCM")
  for (fn in fns_extra) {
    if (fn %in% return_extra) { res[[fn]] <- eval(as.name(fn)) }
  }
  
  return(res)
  
}
