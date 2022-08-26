#' Estimate the causal dose-response curve at a point
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param estimator Which estimator to use; currently only "Grenander"
#' @param params A list, containing the following:
#'   - `Q_n_type` Type of survival function estimator; corresponds to
#'     construct_Q_n()
#'   - `g_n_type` Type of conditional density ratio estimator; corresponds to
#'     construct_f_sIx_n()
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
#'   - `convex_type` One of c("GCM", "LS"); whether to fit the greatest convex
#'     minorant (GCM) to the primitive or the least squares (LS) line
#' @param points A vector representing the points at which estimates and CIs
#'     should be calculated. This should be a unique increasing sequence.
#' @param dir Direction of monotonicity; one of c("incr", "decr")
#' @param return_extra A character vector of additional components to return
#' @return A list of point estimates, SEs and CIs
est_curve <- function(dat_orig, estimator, params, points, dir="decr",
                      return_extra=NULL) {
  
  # Set default params
  .default_params <- list(
    Q_n_type="Super Learner", g_n_type="binning", deriv_type="linear",
    ecdf_type="linear (mid)", gamma_type="Super Learner", q_n_type="new",
    omega_n_type="estimated", boot_reps=1000, ci_type="logit", cf_folds=1, m=5,
    edge_corr="none", lod_shift="none", n_bins=5,
    convex_type="GCM", f_sIx_n_bins=15
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
  if (F) {
    s_shift <- 0
    s_scale <- 1
  } # DEBUG
  dat_orig$s <- (dat_orig$s+s_shift)*s_scale
  dat_orig <- round_dat(dat_orig)
  
  # Obtain minimum value (excluding edge point mass)
  if (p$edge_corr=="min") { s_min2 <- min(dat_orig$s[dat_orig$s!=0],na.rm=T) }
  
  # Rescale points and remove points outside the range of A
  points_orig <- points
  # na_head <- sum(round(points,-log10(C$appx$s))<round(s_min,-log10(C$appx$s))) # !!!!! For some sim reps, s=-0.01 was kept
  points <- round((points+s_shift)*s_scale, -log10(C$appx$s))
  na_head <- sum(points<0) # !!!!! New
  na_tail <- sum(points>1)
  if (na_head>0) { points <- points[-c(1:na_head)] }
  len_p <- length(points)
  if (na_tail>0) { points <- points[-c((len_p-na_tail+1):len_p)] }
  
  dat <- ss(dat_orig, which(dat_orig$z==1))
  
  if (estimator=="Grenander") {
    
    vlist <- create_val_list(dat_orig)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type, print_coeffs=T)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
    
    if (F) {
      # omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
      system.time({
        for (s in seq(0,1,0.01)) {
          ooo <- omega_n(dat$x,rep(s,length(dat$s)),dat$y,dat$delta)
        }
      })
    } # DEBUG
    
    f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                 k=p$f_sIx_n_bins, edge_corr=p$edge_corr)
    f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    dat2 <- ss(dat, which(dat$s!=0))
    Phi_n <- construct_Phi_n(dat2, type=p$ecdf_type)
    n_orig <- length(dat_orig$z)
    p_n <- (1/n_orig) * sum(dat$weights * as.integer(dat$s!=0))
    eta_n <- construct_eta_n(dat, Q_n, p_n, vals=NA)
    gcomp_n <- construct_gcomp_n(dat_orig, vals=vlist$S_grid, Q_n)
    alpha_star_n <- construct_alpha_star_n(dat, gcomp_n, p_n, vals=NA)
    f_n_srv <- construct_f_n_srv(Q_n=Q_n, Qc_n=Qc_n)
    q_n <- construct_q_n(type=p$q_n_type, dat, dat_orig, omega_n=omega_n, g_n=g_n,
                         p_n=p_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
                         Q_n=Q_n, Qc_n=Qc_n, f_n_srv=f_n_srv)
    Gamma_os_n <- construct_Gamma_os_n(dat, dat_orig, omega_n, g_n, eta_n, p_n,
                                       q_n, gcomp_n, alpha_star_n,
                                       vals=vlist$S_grid)
    
    if (F) {
      Gamma_os_n(0.5) # Don't profile this line
      Gamma_os_n(0.6)
      Gamma_os_n(0.7)
    } # DEBUG
    
    if (F) {
      m <- 10^5
      if (L$distr_S=="Unif(0,1)") {
        s <- runif(m)
      } else if (L$distr_S=="Unif(0.3,0.7)") {
        s <- runif(m, min=0.3, max=0.7)
      } else if (L$distr_S=="N(0.5,0.01)") {
        s <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.1)
      } else if (L$distr_S=="N(0.5,0.04)") {
        s <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.2)
      } else if (L$distr_S=="Tri UP") {
        s <- sqrt(runif(m))
      } else if (L$distr_S=="Tri DN") {
        s <- 1 - sqrt(runif(m))
      }
      s <- (s+s_shift)*s_scale
      Gamma_os_n <- Vectorize(function(u) {
        u <- (u+s_shift)*s_scale # u <- round((u+s_shift)*s_scale, -log10(C$appx$s))
        mean( as.integer(s<=u) * (1-exp(-1*L$sc_params$lmbd*C$t_0)) )
      })
    } # DEBUG: True Gamma_os_n
    
    # Construct one-step edge estimator
    if (p$edge_corr!="none") {
      g_sn <- construct_g_sn(dat, vlist$X_grid, type="logistic")
      r_Mn_edge_est <- r_Mn_edge(dat, g_sn, Q_n, omega_n)
      sigma2_edge_est <- sigma2_edge(dat, g_sn, Q_n, omega_n, r_Mn_edge_est)
    }
    
    if (F) {
      q_n <- construct_q_n(type="new", dat, dat_orig, omega_n=omega_n, g_n=g_n,
                           p_n=p_n, gcomp_n=gcomp_n, alpha_star_n=alpha_star_n,
                           Q_n=Q_n, Qc_n=Qc_n, f_n_srv=f_n_srv)
      q_n(dat_orig$x[1:3,],dat_orig$y[1:3],dat_orig$delta[1:3],u=0.5)
      q_n(dat_orig$x[4:6,],dat_orig$y[4:6],dat_orig$delta[4:6],u=0.5)
    } # DEBUG
    
    # Compute GCM (or least squares line) and extract its derivative
    grid <- sort(unique(dat$s))
    x_vals <- Phi_n(grid)
    indices_to_keep <- !base::duplicated(x_vals)
    x_vals <- x_vals[indices_to_keep]
    if (dir=="incr") {
      y_vals <- Gamma_os_n(grid[indices_to_keep])
    } else {
      y_vals <- -1 * Gamma_os_n(grid[indices_to_keep])
    }
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
      dGCM <- Vectorize(function(u) {
        width <- 0.05
        u1 <- u - width/2; u2 <- u + width/2;
        if (u1<0) { u2 <- u2 - u1; u1 <- 0; }
        if (u2>1) { u1 <- u1 - u2 + 1; u2 <- 1; }
        u1 <- round(u1,3); u2 <- round(u2,3);
        ind1 <- which(pred_x==u1); ind2 <- which(pred_x==u2);
        y1 <- pred_y[ind1]; y2 <- pred_y[ind2];
        return((y2-y1)/width)
      })
    }
    
    # Construct Grenander-based r_Mn
    if (dir=="incr") {
      r_Mn_Gr <- Vectorize(function(u) { min(max(dGCM(Phi_n(u)),0),1) })
    } else {
      r_Mn_Gr <- Vectorize(function(u) { min(max(-1*dGCM(Phi_n(u)),0),1) })
    }
    
    # Compute variance component functions
    f_sIx_z1_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                        k=p$f_sIx_n_bins, z1=TRUE)
    f_s_z1_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_z1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, type=p$gamma_type,
                                 vals=vlist$S_grid, omega_n=omega_n,
                                 f_sIx_n=f_sIx_n, f_s_n=f_s_n,
                                 f_s_z1_n=f_s_z1_n)
    if (F) {
      gamma_n <- function(x,s) { Q_n(C$t_0,x,s)*(1-Q_n(C$t_0,x,s)) }
    } # DEBUG: alternate gamma_n estimator when there is no censoring
    
    g_zn <- construct_g_zn(dat_orig, vals=NA, type="Super Learner", f_sIx_n,
                           f_sIx_z1_n)
    
    # Edge correction
    if (p$edge_corr=="none") {
      
      r_Mn <- r_Mn_Gr
      
    } else if (p$edge_corr=="point") {
      
      r_Mn <- function(u) {
        if(u==0) { r_Mn_edge_est } else { r_Mn_Gr(u) }
      }
      
    } else if (p$edge_corr=="min") {
      
      r_Mn <- Vectorize(function(u) {
        if(u==0 || u<s_min2) {
          r_Mn_edge_est
        } else {
          if (dir=="incr") {
            max(r_Mn_edge_est, r_Mn_Gr(u))
          } else {
            min(r_Mn_edge_est, r_Mn_Gr(u))
          }
        }
      })
      
      gren_points <- sapply(c(1:length(points)), function(i) {
        if (dir=="incr") {
          as.numeric(r_Mn(points[i])>r_Mn_edge_est)
        } else {
          as.numeric(r_Mn(points[i])<r_Mn_edge_est)
        }
      })
      
    }
    
    # Generate estimates for each point
    ests <- r_Mn(points)
    if (F) {
      ests_Gamma <- Gamma_os_n(points)
      ests_Phi <- Phi_n(points)
    } # DEBUG: return Gamma/Phi estimates
    
    # Construct variance scale factor
    deriv_r_Mn <- construct_deriv_r_Mn(r_Mn, type=p$deriv_type,
                                             dir=dir)
    
    tau_n <- construct_tau_n(deriv_r_Mn, gamma_n, f_s_n, g_zn, g_n,
                             dat_orig)
    
    # Generate confidence limits
    if (p$ci_type=="none") {
      
      ci_lo <- rep(0, length(ests))
      ci_hi <- rep(0, length(ests))
      
    } else {
      
      # Generate variance scale factor for each point
      tau_ns <- tau_n(points)
      
      # Construct CIs
      # The 0.975 quantile of the Chernoff distribution occurs at roughly 1.00
      # The Normal approximation would use qnorm(0.975, sd=0.52) instead
      qnt <- 1.00
      n_orig <- length(dat_orig$z)
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
        ci_lo[1] <- ests[1] - 1.96*sqrt(sigma2_edge_est/n_orig)
        ci_hi[1] <- ests[1] + 1.96*sqrt(sigma2_edge_est/n_orig)
      } else if (p$edge_corr=="min") {
        ci_lo2 <- ests - 1.96*sqrt(sigma2_edge_est/n_orig)
        ci_hi2 <- ests + 1.96*sqrt(sigma2_edge_est/n_orig)
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
        s_sub <- dat$s[dat$s!=0]
        wt_sub <- dat$weights[dat$s!=0]
        Phi_n_inv_sub <- construct_Phi_n(
          dat = list(s=s_sub, weights=wt_sub),
          which = "inverse",
          type = p$ecdf_type
        )
        cutoffs <- c(0, Phi_n_inv_sub(seq(0,1,length.out=p$n_bins)))
      }
    }
    
    # Function to transform S values to categorical bins values
    transform_s <- Vectorize(function(s) {
      cut(s, breaks=cutoffs, right=F, include.lowest=T)
    })
    
    # Create vlist
    vlist <- create_val_list(dat_orig, factor_S=unique(transform_s(dat$s)))
    
    # Construct f_sIx_n BEFORE transforming S values
    f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=p$g_n_type,
                                 k=p$f_sIx_n_bins)
    
    # Transform S values
    dat$s <- transform_s(dat$s)
    
    # Note: Q_n and Qc_n will not work with type="true"
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=p$Q_n_type, print_coeffs=T)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n, type=p$omega_n_type)
    g_sn <- construct_g_sn(dat, vlist$X_grid, type="generalized",
                           f_sIx_n=f_sIx_n, cutoffs=cutoffs)
    
    # Generate estimates and standard deviations for each point
    ests <- sapply(c(1:length(points)), function(i) {
      s_binned <- transform_s(points[i])
      return(r_Mn_edge(dat, g_sn, Q_n, omega_n, val=s_binned))
    })
    sigma2s <- sapply(c(1:length(points)), function(i) {
      a_binned <- transform_s(points[i])
      return(sigma2_edge(dat, g_sn, Q_n, omega_n, ests[i], val=s_binned))
    })
    
    # Construct CIs
    n_orig <- length(dat_orig$z)
    ci_lo <- ests - 1.96*sqrt(sigma2s/n_orig)
    ci_hi <- ests + 1.96*sqrt(sigma2s/n_orig)
    
  }
  
  if (estimator=="Cox gcomp") {
    
    # Setup
    dat <- ss(dat_orig, which(dat_orig$z==1))
    vlist <- create_val_list(dat_orig)
    
    # Fit Cox model and compute variance
    res_cox <- cox_var(dat_orig=dat_orig, dat=dat, t=C$t_0,
                   points=points, se_marg=T) # verbose=T
    
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
    dat <- ss(dat_orig, which(dat_orig$z==1))
    
    # Generate spline basis (4 degrees of freedom)
    qnt <- as.numeric(quantile(dat$s, seq(0,1,0.25)))
    qnt <- unique(qnt)
    spl <- list(
      K = qnt[-c(1,length(qnt))],
      B = c(qnt[1],qnt[length(qnt)]),
      L = length(qnt)-1
    )
    ns_basis <- ns(dat$s, knots=spl$K, Boundary.knots=spl$B)
    
    # Fit Cox model
    fml <- "Surv(y,delta)~b1"
    for (i in 2:spl$L) { fml <- paste0(fml, "+b",i) }
    for (i in 1:length(dat$x)) { fml <- paste0(fml, "+x",i) }
    fml <- formula(fml)
    df <- cbind("y"=dat$y, "delta"=dat$delta,
                dat$x, "weights"=dat$weights)
    for (i in 1:spl$L) { df[paste0("b",i)] <- ns_basis[,i] }
    model4 <- coxph(fml, data=df, weights=dat$weights)
    beta_n4 <- coefficients(model4)
    
    # Get Breslow estimator
    bh4 <- basehaz(model4, centered=FALSE)
    index4 <- max(which((bh4$time<C$t_0)==T))
    est_bshz4 <- bh4$hazard[index4]
    
    # Construct conditional survival function
    Q_n4 <- function(x,s) {
      lin <- c(as.numeric(ns(s, knots=spl$K, Boundary.knots=spl$B)),x)
      exp(-1*est_bshz4*exp(sum(as.numeric(beta_n4)*lin)))
    }
    
    # Construct marginalized survival function
    r_M4 <- Vectorize(function(s) {
      1 - mean(sapply(c(1:length(dat_orig$s)), function(i) {
        Q_n4(as.numeric(dat_orig$x[i,]),s)
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
    p3 <- round((0.3+s_shift)*s_scale, -log10(C$appx$s))
    p5 <- round((0.5+s_shift)*s_scale, -log10(C$appx$s))
    res$eta_n3 <- eta_n(p3,c(0,0))
    res$eta_n5 <- eta_n(p5,c(0,0))
    res$alpha_star_n3 <- alpha_star_n(p3)
    res$alpha_star_n5 <- alpha_star_n(p5)
    res$gcomp_n <- gcomp_n(0.9)
  } # DEBUG: return additional estimates
  
  if ("gcomp" %in% return_extra) {
    Q_n2 <- (construct_Q_n(dat, vlist$Q_n, type="Cox PH"))$srv
    res$gcomp <- construct_gcomp_n(dat_orig, vlist$S_grid, Q_n=Q_n2)
  }
  fns_extra <- c("f_s_n", "gamma_n", "deriv_r_Mn", "Phi_n_inv", "Phi_n",
                 "omega_n", "f_sIx_n", "Q_n", "gcm", "dGCM", "grid",
                 "Gamma_os_n", "r_Mn_Gr")
  for (fn in fns_extra) {
    if (fn %in% return_extra) { res[[fn]] <- eval(as.name(fn)) }
  }
  
  return(res)
  
}
