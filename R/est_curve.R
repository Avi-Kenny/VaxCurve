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
#'     construct_deriv_r_Mn()
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
#'   - `edge_corr` One of c("none", "point", "min"). "point" uses the
#'     root-n estimator to replace the Grenander-based estimator only at the
#'     leftmost point. "min" adjusts both the leftmost point and the rest of the
#'     curve.
#'   - `marg` One of c("Gamma", "Theta"); whether or not to transform by the
#'     marginal distribution of A
#'   - `convex_type` One of c("GCM", "LS"); whether to fit the GCM to the
#'     primitive or the least squares line
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
    S_n_type="Super Learner", g_n_type="binning", deriv_type="linear",
    ecdf_type="linear (mid)", gamma_type="Super Learner", q_n_type="new",
    omega_n_type="estimated", boot_reps=1000, ci_type="trunc", cf_folds=1, m=5,
    edge_corr="none", marg="Gamma_star2", lod_shift="none", n_bins=5,
    convex_type="GCM", f_aIw_n_bins=15
  )
  for (i in c(1:length(.default_params))) {
    if (is.null(params[[names(.default_params)[i]]])) {
      params[[names(.default_params)[i]]] <- .default_params[[i]]
    }
  }
  p <- params
  
  # Rescale A to lie in [0,1] and round values
  a_min <- min(dat_orig$a,na.rm=T)
  a_max <- max(dat_orig$a,na.rm=T)
  a_shift <- -1 * a_min
  a_scale <- 1/(a_max-a_min)
  if (F) {
    a_shift <- 0
    a_scale <- 1
  } # DEBUG
  dat_orig$a <- (dat_orig$a+a_shift)*a_scale
  dat_orig <- round_dat(dat_orig)
  
  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr=="min") { a_min2 <- min(dat_orig$a[dat_orig$a!=0],na.rm=T) }
  
  # Rescale points and remove points outside the range of A
  points_orig <- points
  na_head <- sum(round(points,-log10(C$appx$a))<round(a_min,-log10(C$appx$a)))
  points <- round((points+a_shift)*a_scale, -log10(C$appx$a))
  na_tail <- sum(points>1)
  if (na_head>0) {
    points <- points[-c(1:na_head)]
  }
  if (na_tail>0) {
    points <- points[-c((length(points)-na_tail+1):length(points))]
  }
  
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  
  if (estimator=="Grenander") {
    
    print(paste("Check 0:", Sys.time()))
    vlist <- create_val_list(dat_orig)
    srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type, print_coeffs=T)
    S_n <- srvSL$srv
    Sc_n <- srvSL$cens
    print(paste("Check 1:", Sys.time()))
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type=p$omega_n_type)
    
    # !!!!!
    if (F) {
      
      # omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type=p$omega_n_type)
      system.time({
        for (a in seq(0,1,0.01)) {
          ooo <- omega_n(dat$w,rep(a,length(dat$a)),dat$y_star,dat$delta_star)
        }
      })

    }
    
    print(paste("Check 2:", Sys.time()))
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
                                 k=p$f_aIw_n_bins, edge_corr=p$edge_corr)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    print(paste("Check 3:", Sys.time()))
    
    if (p$marg %in% c("Gamma_star", "Gamma_star2")) {
      dat2 <- ss(dat, which(dat$a!=0))
      Phi_n <- construct_Phi_n(dat2, type=p$ecdf_type)
      print(paste("Check 4:", Sys.time()))
      n_orig <- length(dat_orig$delta)
      z_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$a!=0))
      g_n_star <- construct_g_n_star(f_aIw_n, f_a_n, z_n)
      print(paste("Check 5:", Sys.time()))
      eta_ss_n <- construct_eta_ss_n(dat, S_n, z_n, vals=NA)
      print(paste("Check 6:", Sys.time()))
      gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$A_grid, S_n)
      print(paste("Check 7:", Sys.time()))
      alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, z_n, vals=NA)
      print(paste("Check 8:", Sys.time()))
    }
    
    if (p$marg=="Theta") {
      etastar_n <- construct_etastar_n(S_n)
      Theta_os_n <- construct_Theta_os_n(dat, vlist$A_grid, omega_n, f_aIw_n,
                                         etastar_n)
    } else if (p$marg=="Gamma") {
      Phi_n <- construct_Phi_n(dat, type=p$ecdf_type)
      Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    } else if (p$marg=="Gamma_star") {
      Gamma_os_n_star <- construct_Gamma_os_n_star(dat, omega_n, g_n_star,
                                                   eta_ss_n, z_n, gcomp_n,
                                                   alpha_star_n,
                                                   vals=vlist$A_grid)
    } else if (p$marg=="Gamma_star2") {
      print(paste("Check 9:", Sys.time()))
      q_n <- construct_q_n(type=p$q_n_type, dat, dat_orig, omega_n=omega_n, g_n=g_n,
                           z_n=z_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
                           S_n=S_n, Sc_n=Sc_n)
      print(paste("Check 10:", Sys.time()))
      Gamma_os_n_star <- construct_Gamma_os_n_star2(dat, dat_orig, omega_n,
                                                    g_n, eta_ss_n, z_n, # g_n_star
                                                    q_n, gcomp_n, alpha_star_n,
                                                    vals=vlist$A_grid)
      # !!!!!
      if (F) {
        Gamma_os_n_star(0.5) # Don't profile this line
        Gamma_os_n_star(0.6)
        Gamma_os_n_star(0.7)
      }
      
      print(paste("Check 11:", Sys.time()))
      
      if (F) {
        m <- 10^5
        if (L$distr_A=="Unif(0,1)") {
          a <- runif(m)
        } else if (L$distr_A=="Unif(0.3,0.7)") {
          a <- runif(m, min=0.3, max=0.7)
        } else if (L$distr_A=="N(0.5,0.01)") {
          a <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.1)
        } else if (L$distr_A=="N(0.5,0.04)") {
          a <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.2)
        } else if (L$distr_A=="Tri UP") {
          a <- sqrt(runif(m))
        } else if (L$distr_A=="Tri DN") {
          a <- 1 - sqrt(runif(m))
        }
        a <- (a+a_shift)*a_scale
        Gamma_os_n_star <- Vectorize(function(x) {
          x <- (x+a_shift)*a_scale # x <- round((x+a_shift)*a_scale, -log10(C$appx$a))
          mean( as.integer(a<=x) * (1-exp(-1*L$sc_params$lmbd*C$t_e)) )
        })
      } # DEBUG: True Gamma_os_n_star
      
    } else {
      stop(paste0("`params$marg` must be one of c('Theta', 'Gamma', 'Gamma_s",
                  "tar', 'Gamma_star2'"))
    }
    
    # Construct one-step edge estimator
    if (p$edge_corr!="none") {
      pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
      theta_os_n_est <- theta_os_n(dat, pi_n, S_n, omega_n)
      sigma2_os_n_est <- sigma2_os_n(dat, pi_n, S_n, omega_n, theta_os_n_est)
    }
    
    # !!!!!
    if (F) {
      q_n <- construct_q_n(type="new", dat, dat_orig, omega_n=omega_n, g_n=g_n,
                           z_n=z_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
                           S_n=S_n, Sc_n=Sc_n)
      # q_n(dat_orig$w[1,],dat_orig$y_star[1],dat_orig$delta_star[1],x=0.5)
      q_n(dat_orig$w[1:3,],dat_orig$y_star[1:3],dat_orig$delta_star[1:3],x=0.5)
      q_n(dat_orig$w[4:6,],dat_orig$y_star[4:6],dat_orig$delta_star[4:6],x=0.5)
      
      
      # q_n(dat_orig$w[1:99,],dat_orig$y_star[1:99],dat_orig$delta_star[1:99],x=0.5)
      # q_n(dat_orig$w,dat_orig$y_star,dat_orig$delta_star,x=0.5)
      # Gamma_os_n_star(0.1)
    }
    
    # Compute GCM (or least squares line) and extract its derivative
    print(paste("Check 12:", Sys.time()))
    grid <- sort(unique(dat$a))
    x_vals <- Phi_n(grid)
    indices_to_keep <- !base::duplicated(x_vals)
    x_vals <- x_vals[indices_to_keep]
    if (dir=="incr") {
      if (p$marg=="Gamma") {
        y_vals <- Gamma_os_n(grid[indices_to_keep])
      } else {
        y_vals <- Gamma_os_n_star(grid[indices_to_keep])
      }
    } else {
      if (p$marg=="Gamma") {
        y_vals <- -1 * Gamma_os_n(grid[indices_to_keep])
      } else {
        
        # !!!!! Debugging !!!!!
        if (T) {
          for (aa in grid[indices_to_keep]) {
            print("SPECIAL CHECK")
            print(Sys.time())
            print(paste("aa:",aa))
            bb <- Gamma_os_n_star(aa)
            print(paste("bb:",bb))
          }
        }
        
        y_vals <- -1 * Gamma_os_n_star(grid[indices_to_keep])
      }
    }
    print(paste("Check 13:", Sys.time()))
    if (!any(x_vals==0)) {
      x_vals <- c(0, x_vals)
      y_vals <- c(0, y_vals)
    }
    if (p$convex_type=="GCM") {
      gcm <- gcmlcm(x=x_vals, y=y_vals, type="gcm")
      dGCM <- approxfun(
        x = gcm$x.knots[-length(gcm$x.knots)],
        y = gcm$slope.knots,
        method = "constant",
        rule = 2,
        f = 0
      )
    } else if (p$convex_type=="LS") {
      gcm <- function(x) { 1 } # Ignored
      fit <- cvx.lse.reg(t=x_vals, z=y_vals)
      pred_x <- round(seq(0,1,0.001),3)
      pred_y <- predict(fit, newdata=pred_x)
      dGCM <- Vectorize(function(x) {
        width <- 0.05
        x1 <- x - width/2; x2 <- x + width/2;
        if (x1<0) { x2 <- x2 - x1; x1 <- 0; }
        if (x2>1) { x1 <- x1 - x2 + 1; x2 <- 1; }
        x1 <- round(x1,3); x2 <- round(x2,3);
        ind1 <- which(pred_x==x1); ind2 <- which(pred_x==x2);
        y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
        return((y2-y1)/width)
      })
    }
    print(paste("Check 14:", Sys.time()))
    
    # Construct Grenander-based theta_n
    if (p$marg=="Theta") {
      if (dir=="incr") {
        theta_n_Gr <- Vectorize(function(x) {
          min(max(dGCM(x),0),1) # dGCM(x)
        })
      } else {
        theta_n_Gr <- Vectorize(function(x) {
          min(max(-1*dGCM(x),0),1) # -1*dGCM(x)
        })
      }
    } else if (p$marg %in% c("Gamma", "Gamma_star", "Gamma_star2")) {
      if (dir=="incr") {
        theta_n_Gr <- Vectorize(function(x) {
          min(max(dGCM(Phi_n(x)),0),1) # dGCM(Phi_n(x))
        })
      } else {
        theta_n_Gr <- Vectorize(function(x) {
          min(max(-1 * dGCM(Phi_n(x)),0),1) # -1 * dGCM(Phi_n(x))
        })
      }
    }
    
    # Compute variance component functions
    print(paste("Check 15:", Sys.time()))
    f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
                                        k=p$f_aIw_n_bins, delta1=TRUE)
    f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_delta1_n)
    print(paste("Check 16:", Sys.time()))
    gamma_n <- construct_gamma_n(dat_orig, dat, type=p$gamma_type,
                                 vals=vlist$A_grid, omega_n=omega_n,
                                 f_aIw_n=f_aIw_n, f_a_n=f_a_n,
                                 f_a_delta1_n=f_a_delta1_n)
    if (F) {
      gamma_n <- function(w,a) { S_n(C$t_e,w,a)*(1-S_n(C$t_e,w,a)) }
    } # DEBUG: alternate gamma_n estimator when there is no censoring
    print(paste("Check 17:", Sys.time()))
    pi_star_n <- construct_pi_star_n(dat_orig, vals=NA, type="Super Learner",
                                     f_aIw_n, f_aIw_delta1_n)
    print(paste("Check 18:", Sys.time()))
    
    # Edge correction
    if (p$edge_corr=="none") {
      
      theta_n <- theta_n_Gr
      
    } else if (p$edge_corr=="point") {
      
      theta_n <- function(x) {
        if(x==0) { theta_os_n_est } else { theta_n_Gr(x) }
      }
      
    } else if (p$edge_corr=="min") {
      
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
    print(paste("Check 19:", Sys.time()))
    ests <- theta_n(points)
    if (F) {
      ests_Gamma <- Gamma_os_n_star(points)
      ests_Phi <- Phi_n(points)
    } # DEBUG: return Gamma/Phi estimates
    print(paste("Check 20:", Sys.time()))
    
    # Construct variance scale factor
    deriv_r_Mn <- construct_deriv_r_Mn(theta_n, type=p$deriv_type,
                                             dir=dir)
    
    print(paste("Check 21:", Sys.time()))
    tau_n <- construct_tau_n(deriv_r_Mn, gamma_n, f_a_n, pi_star_n, g_n,
                             dat_orig)
    print(paste("Check 22:", Sys.time()))
    
    # Generate confidence limits
    if (p$ci_type=="none") {
      
      ci_lo <- rep(0, length(ests))
      ci_hi <- rep(0, length(ests))
      
    } else {
      
      # Generate variance scale factor for each point
      tau_ns <- tau_n(points)
      print(paste("Check 23:", Sys.time()))
      
      # Construct CIs
      # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
      # The Normal approximation would use qnorm(0.975, sd=0.52) instead
      qnt <- 1.00
      n_orig <- length(dat_orig$delta)
      if (p$ci_type=="regular") {
        ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
        ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
      } else if (p$ci_type=="logit") {
        ci_lo <- expit(
          logit(ests) - (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
        )
        ci_hi <- expit(
          logit(ests) + (qnt*tau_ns*deriv_logit(ests))/(n_orig^(1/3))
        )
      } else if (p$ci_type=="trunc") {
        ci_lo <- ests - (qnt*tau_ns)/(n_orig^(1/3))
        ci_hi <- ests + (qnt*tau_ns)/(n_orig^(1/3))
        ci_lo %<>% pmax(0) %>% pmin(1)
        ci_hi %<>% pmax(0) %>% pmin(1)
      }
      
      # Edge correction
      if (p$edge_corr=="point") {
        ci_lo[1] <- ests[1] - 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_hi[1] <- ests[1] + 1.96*sqrt(sigma2_os_n_est/n_orig)
      } else if (p$edge_corr=="min") {
        ci_lo2 <- ests - 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_hi2 <- ests + 1.96*sqrt(sigma2_os_n_est/n_orig)
        ci_lo <- (1-gren_points)*ci_lo2 + replace_na(gren_points*ci_lo, 0)
        ci_hi <- (1-gren_points)*ci_hi2 + replace_na(gren_points*ci_hi, 0)
      }
      
    }
    
    # # Add extra return data
    # res[["Phi_n"]] <- Phi_n
    
  }
  
  if (estimator=="Qbins") {
    
    # Construct bin cutoffs
    {
      Phi_n_inv <- construct_Phi_n(dat, which="inverse", type="step")
      cutoffs <- Phi_n_inv(seq(0,1,length.out=p$n_bins+1))
      if (sum(cutoffs==0)>1) {
        a_sub <- dat$a[dat$a!=0]
        wt_sub <- dat$weights[dat$a!=0]
        Phi_n_inv_sub <- construct_Phi_n(
          dat = list(a=a_sub, weights=wt_sub),
          which = "inverse",
          type = p$ecdf_type
        )
        cutoffs <- c(0, Phi_n_inv_sub(seq(0,1,length.out=p$n_bins)))
      }
    }
    
    # Function to transform A values to categorical bins values
    transform_a <- Vectorize(function(a) {
      cut(a, breaks=cutoffs, right=F, include.lowest=T)
    })
    
    # Create vlist
    vlist <- create_val_list(dat_orig, factor_A=unique(transform_a(dat$a)))
    
    # Construct f_aIw_n BEFORE transforming A values
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=p$g_n_type,
                                 k=p$f_aIw_n_bins)
    
    # Transform A values
    dat$a <- transform_a(dat$a)
    
    # Note: S_n and Sc_n will not work with type="true"
    srvSL <- construct_S_n(dat, vlist$S_n, type=p$S_n_type, print_coeffs=T)
    S_n <- srvSL$srv
    Sc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type=p$omega_n_type)
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
    
  }
  
  if (estimator=="Cox gcomp") {
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat_orig)
    
    # Fit Cox model and compute variance
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_e,
                   points=points, se_marg=T) # !!!!! verbose=T
    
    # Compute CIs (logit transformed)
    ests <- 1-res_cox$est_marg
    ses <- sqrt(res_cox$var_est_marg)
    ci_lo <- expit(logit(ests) - 1.96*deriv_logit(ests)*ses)
    ci_hi <- expit(logit(ests) + 1.96*deriv_logit(ests)*ses)
    # ci_lo <- (ests - 1.96*ses) %>% pmax(0) %>% pmin(1)
    # ci_hi <- (ests + 1.96*ses) %>% pmax(0) %>% pmin(1)
    
  }
  
  # !!!!! Update
  if (estimator=="Cox GAM") {
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    
    # Generate spline basis (4 degrees of freedom)
    qnt <- as.numeric(quantile(dat$a, seq(0,1,0.25)))
    qnt <- unique(qnt)
    spl <- list(
      K = qnt[-c(1,length(qnt))],
      B = c(qnt[1],qnt[length(qnt)]),
      L = length(qnt)-1
    )
    ns_basis <- ns(dat$a, knots=spl$K, Boundary.knots=spl$B)
    
    # Fit Cox model
    fml <- "Surv(y_star,delta_star)~b1"
    for (i in 2:spl$L) { fml <- paste0(fml, "+b",i) }
    for (i in 1:length(dat$w)) { fml <- paste0(fml, "+w",i) }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star,
                dat$w, "weights"=dat$weights)
    for (i in 1:spl$L) { df[paste0("b",i)] <- ns_basis[,i] }
    model4 <- coxph(fml, data=df, weights=dat$weights)
    beta_n4 <- coefficients(model4)
    
    # Get Breslow estimator
    bh4 <- basehaz(model4, centered=FALSE)
    index4 <- max(which((bh4$time<C$t_e)==T))
    est_bshz4 <- bh4$hazard[index4]
    
    # Construct conditional survival function
    S_n4 <- function(w,a) {
      lin <- c(as.numeric(ns(a, knots=spl$K, Boundary.knots=spl$B)),w)
      exp(-1*est_bshz4*exp(sum(as.numeric(beta_n4)*lin)))
    }
    
    # Construct marginalized survival function
    r_M4 <- Vectorize(function(a) {
      1 - mean(sapply(c(1:length(dat_orig$a)), function(i) {
        S_n4(as.numeric(dat_orig$w[i,]),a)
      }))
    })
    
    # Compute CIs
    ests <- r_M4(points)
    ci_lo <- rep(NA, length(points))
    ci_hi <- rep(NA, length(points))
    
  }
  
  # Parse and return results
  res <- list(
    point = points_orig,
    est = c(rep(NA,na_head), ests, rep(NA,na_tail)),
    ci_lo = c(rep(NA,na_head), ci_lo, rep(NA,na_tail)),
    ci_hi = c(rep(NA,na_head), ci_hi, rep(NA,na_tail))
  )
  
  if (F) {
    res$ests_Gamma = c(rep(NA,na_head), ests_Gamma, rep(NA,na_tail))
    res$ests_Phi = c(rep(NA,na_head), ests_Phi, rep(NA,na_tail))
  } # DEBUG: return Gamma/Phi estimates
  
  if (estimator=="Grenander") {
    res$tau_ns <- c(rep(NA,na_head), tau_ns, rep(NA,na_tail))
    res$n <- n_orig
  }
  
  if (F) {
    p3 <- round((0.3+a_shift)*a_scale, -log10(C$appx$a))
    p5 <- round((0.5+a_shift)*a_scale, -log10(C$appx$a))
    res$g_n_star <- g_n_star(0.5,c(0,0))
    res$eta_ss_n3 <- eta_ss_n(p3,c(0,0))
    res$eta_ss_n5 <- eta_ss_n(p5,c(0,0))
    res$alpha_star_n3 <- alpha_star_n(p3)
    res$alpha_star_n5 <- alpha_star_n(p5)
    res$gcomp_n <- gcomp_n(0.9)
  } # DEBUG: return additional estimates
  
  if ("gcomp" %in% return_extra) {
    S_n2 <- (construct_S_n(dat, vlist$S_n, type="Cox PH"))$srv
    res$gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n2)
  }
  fns_extra <- c("f_a_n", "gamma_n", "deriv_r_Mn", "Phi_n_inv", "Theta_os_n",
                 "Phi_n", "omega_n", "f_aIw_n", "etastar_n", "S_n", "gcm",
                 "dGCM", "grid", "Gamma_os_n_star",
                 "theta_n_Gr")
  for (fn in fns_extra) {
    if (fn %in% return_extra) { res[[fn]] <- eval(as.name(fn)) }
  }
  
  return(res)
  
}
