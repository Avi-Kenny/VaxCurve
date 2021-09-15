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
#'   - `edge_corr` One of c("none", "spread", "point", "weighted", "max").
#'     "point" uses the root-n estimator to replace the Grenander-based
#'     estimator only at the leftmost point. "weighted" and "max" both use the
#'     precision-weighted estimator to adjust both the leftmost point and the
#'     rest of the curve. "spread" adds a small amount of noise to the edge
#'     points
#' @param points A vector representing the points to estimate
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)
est_curve <- function(dat_orig, estimator, params, points) {
  
  if (estimator=="Grenander") {
    
    if (params$edge_corr=="spread") {
      dat_orig$noise <- runif(nrow(dat_orig))*0.05
      dat_orig %<>% mutate(a = ifelse(a==0, noise, a))
    }
    
    # Prep
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig, C$appx)
    
    # Construct regular Gamma_0 estimator
    if (params$cf_folds==1) {
      
      # Construct component functions
      Phi_n <- construct_Phi_n(dat_orig)
      Phi_n_inv <- construct_Phi_n(dat_orig, type="inverse")
      S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
      Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                            csf=TRUE)
      f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                   type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
      g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
      omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
      
      # Construct Gamma_n
      Gamma_n <- construct_Gamma_n(dat_orig, vlist$A_grid, omega_n, S_n, g_n)
      
      # Construct one-step edge estimator
      if (params$edge_corr!="none") {
        pi_n <- construct_pi_n(dat_orig, vlist$W_grid, type="logistic")
        theta_os_n_est <- theta_os_n(dat_orig, pi_n, S_n, omega_n)
        sigma2_os_n_est <- sigma2_os_n(dat_orig, pi_n, S_n, omega_n,
                                       theta_os_n_est)
      }
      
    }
    
    # Construct cross-fitted Gamma_0 estimator
    if (params$cf_folds>1) {
      
      Gamma_n <- construct_Gamma_cf(dat_orig, params, vlist)
      
      # Recompute functions on full dataset
      Phi_n <- construct_Phi_n(dat_orig)
      Phi_n_inv <- construct_Phi_n(dat_orig, type="inverse")
      
      # Construct cross-fitted one-step edge estimator
      if (params$edge_corr!="none") {
        # !!!!! TO DO
      }
      
    }
    
    Psi_n <- Vectorize(function(x) {
      return(Gamma_n(Phi_n_inv(x)))
    })
    gcm <- gcmlcm(x=seq(0,1,0.0001), y=Psi_n(seq(0,1,0.0001)), type="gcm")
    dGCM <- Vectorize(function(x) {
      if (x==0) {
        index <- 1
      } else {
        index <- which(x<=gcm$x.knots)[1]-1
      }
      return(gcm$slope.knots[index])
    })
    
    # Construct Grenander-based theta_n
    theta_n_Gr <- function(x) {
      x_trans <- Phi_n(x)
      return(dGCM(x_trans))
    }
    
    # Recompute functions on full dataset
    if (params$cf_folds>1) {
      S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
      Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                            csf=TRUE)
      f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                   type=params$g_n_type, k=15)
      f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
      omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    }
    
    # Compute new functions
    gcomp_n <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n)
    f_aIw_delta1_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                        type=params$g_n_type, k=15, delta1=TRUE)
    f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                    f_aIw_delta1_n)
    gamma_n <- construct_gamma_n(dat_orig, vlist$A_grid, type="kernel",
                                 omega_n, f_aIw_n, f_a_n, f_a_delta1_n)
    
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
      
    } else if (params$edge_corr=="weighted") {
      
      # sigma2_Gr <- ( 0.52 * (n_orig^(-1/3)) * tau_n(0) )^2
      # sigma2_OS <- sigma2_os_n_est / n_orig
      # 
      # wt <- sigma2_Gr / (sigma2_OS + sigma2_Gr)
      # theta_n_weighted <- wt*theta_os_n_est + (1-wt)*theta_n_Gr(0)
      # 
      # theta_n <- function(x) {
      #   max(theta_n_weighted, theta_n_Gr(x))
      # }
      # 
      # gren_ests <- sapply(points, theta_n_Gr)
      # gren_points <- sapply(points, function(x) {
      #   as.numeric(gren_ests>theta_n_weighted)
      # })
      # gren_points[1] <- 0
      
    } else if (params$edge_corr=="max") {
      
      theta_n <- function(x) {
        max(theta_os_n_est, theta_n_Gr(x))
      }
      
      gren_ests <- sapply(points, theta_n_Gr)
      gren_points <- sapply(points, function(x) {
        as.numeric(gren_ests>theta_os_n_est)
      })
      gren_points[1] <- 0
      
    }
    
    # Construct variance scale factor
    if (params$deriv_type=="linear") {
      deriv_theta_n <- construct_deriv_theta_n(theta_n, type="linear")
    } else if (params$deriv_type=="line") {
      deriv_theta_n <- construct_deriv_theta_n(theta_n, type="line")
    } else if (params$deriv_type=="gcomp") {
      deriv_theta_n <- construct_deriv_theta_n(gcomp_n, type="gcomp")
    } else if (params$deriv_type=="spline") {
      deriv_theta_n <- construct_deriv_theta_n(theta_n, type="spline")
    }
    tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
    
    # Generate estimates for each point
    # The pmax() prevents errors when estimates are negative
    ests <- pmax(theta_n(points),0)
    
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
      # qnt <- 1.00
      qnt <- qnorm(0.975, sd=0.52)
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
      }
      
      # Edge correction
      if (params$edge_corr=="point") {
        ci_lo[1] <- ests[1] - 1.96*sqrt(sigma2_os_n_est/n_orig) # !!!!! logit transform here as well
        ci_hi[1] <- ests[1] + 1.96*sqrt(sigma2_os_n_est/n_orig) # !!!!! logit transform here as well
      } else if (params$edge_corr %in% c("weighted","max")) {
        ci_lo2 <- ests - 1.96*sqrt(sigma2_os_n_est/n_orig) # !!!!! logit transform here as well
        ci_hi2 <- ests + 1.96*sqrt(sigma2_os_n_est/n_orig) # !!!!! logit transform here as well
        ci_lo <- (1-gren_points)*pmin(ci_lo,ci_lo2) + gren_points*ci_lo
        ci_hi <- (1-gren_points)*pmax(ci_hi,ci_hi2) + gren_points*ci_hi
      }
      
    }
    
    # Parse and return results
    res <- list()
    for (p in 1:length(points)) {
      res[[p]] <- list(point=points[p], est=ests[p],
                       ci_lo=ci_lo[p], ci_hi=ci_hi[p])
    }
    
    # Add extra return data
    res[["ex_S_n"]] <- S_n(C$t_e, w1=0.5, w2=1, a=0.5)
    res[["ex_gamma_n"]] <- gamma_n(0.5)
    res[["ex_deriv_theta_n"]] <- deriv_theta_n(0.5)
    res[["ex_tau_n"]] <- tau_n(0.5)
    
    return(res)
    
  }

}
