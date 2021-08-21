#' Estimate the causal dose-response curve at a point
#' 
#' @param dat_orig Data returned by generate_data(); FULL data
#' @param estimator Which estimator to use; one of c("G-comp", "Generalized
#'     Grenander")
#' @param params A list, containing the following:
#'   - `S_n_type` Type of survival function estimator; currently only c("Cox PH")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `boot_reps` Used for G-comp; number of bootstrap replicates
#'   - `ci_type` One of c("regular", "logit", "sample split", "none"). "regular"
#'     is the standard approach. "logit" transforms the CIs so that the bounds
#'     are in [0,1]. "sample split" is the Banerjee method (and also returns a
#'     different estimator).
#'   - `cf_folds` Number of cross-fitting folds
#'   - `m` If params$ci_type=="sample split", the number of splits
#' @param points A vector representing the points to estimate
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)
est_curve <- function(dat_orig, estimator, params, points) {
  
  if (estimator=="G-comp") {
    
    # Compute estimates
    S_n <- construct_S_n(dat_orig, type=params$S_n_type)
    gcomp_n <- construct_gcomp(dat_orig, S_n=S_n)
    ests <- gcomp_n(points)
    
    # Run bootstrap for SEs
    {
      my_stat <- function(dat_orig,indices) {
        d <- dat_orig[indices,]
        S_n <- construct_S_n(d, type=params$S_n_type)
        gcomp_n <- construct_gcomp(d, S_n=S_n)
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
  
  if (estimator=="Generalized Grenander") {
    
    if (is.null(params$cf_folds) || params$cf_folds==1) {
      
      
      
    }
    
    # # Construct data splits (discarding "extra rows" at end)
    # m <- params$m
    # splits <- matrix(NA, nrow=m, ncol=2)
    # split_size <- as.integer(n/m)
    # splits[,2] <- (1:m)*split_size
    # splits[,1] <- ((1:m)*split_size+1)-split_size
    # 
    # # Construct estimate separately for each data split
    # ests <- c()
    # ses <- c()
    # for (point in points) {
    #   split_ests <- sapply(c(1:m), function(x) {
    #     dat_split <- dat_orig[c(splits[x,1]:splits[x,2]),]
    #     theta_n <- construct_fns(dat_orig=dat_split, return_tau_n=F)$theta_n
    #     return(theta_n(point))
    #   })
    #   
    #   ests <- c(ests, mean(split_ests))
    #   ses <- c(ses, sd(split_ests)/sqrt(m))
    # }
    # 
    # # Construct CIs
    # t_quant <- qt(1-(0.05/2), df=(m-1))
    # ci_lo <- ests - t_quant*ses
    # ci_hi <- ests + t_quant*ses
    
    # Construct theta_n and tau_n, given a dataset
    construct_fns <- function(dat_orig, return_tau_n=T) {
      
      # Construct component functions
      grid <- seq(0,1,0.01)
      Phi_n <- construct_Phi_n(dat_orig)
      Phi_n_inv <- construct_Phi_n(dat_orig, type="inverse")
      gcomp_n <- construct_gcomp(dat_orig, S_n=S_n)
      f_aIw_n <- construct_f_aIw_n(dat_orig, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat_orig, f_aIw_n=f_aIw_n)
      g_n <- construct_g_n(f_aIw_n, f_a_n)
      S_n <- construct_S_n(dat_orig, type=params$S_n_type)
      Sc_n <- construct_S_n(dat_orig, type=params$S_n_type, csf=TRUE)
      omega_n <- construct_omega_n(S_n, Sc_n)
      Gamma_n <- construct_Gamma_n(dat_orig, omega_n, S_n, g_n)

      Psi_n <- function(x) { Gamma_n(Phi_n_inv(x)) }
      gcm <- gcmlcm(x=grid, y=Psi_n(grid), type="gcm")
      dGCM <- Vectorize(function(x) {
        if (x==0) {
          index <- 1
        } else {
          index <- which(x<=gcm$x.knots)[1]-1
        }
        return(gcm$slope.knots[index])
      })
      
      # Construct theta_n
      theta_n <- function(x) {
        x_trans <- Phi_n(x)
        return(dGCM(x_trans))
      }
      
      if (return_tau_n==F) {
        
        return(list(theta_n=theta_n))
        
      } else {
        
        # Construct tau_n
        if (return_tau_n==T) {
          deriv_theta_n <- construct_deriv_theta_n(gcomp_n)
          gamma_n <- construct_gamma_n(dat_orig, type="cubic", omega_n, f_aIw_n)
          tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
        }
        
        return(list(theta_n=theta_n, tau_n=tau_n))
        
      }
      
    }
    
    # Construct the sample split estimator
    if (params$ci_type=="sample split") {
      
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
      
    } else {
      
      # Construct theta_n
      fns <- construct_fns(dat=dat_orig, return_tau_n=T)
      theta_n <- fns$theta_n
      
      # Generate estimates for each point
      ests <- sapply(points, theta_n)
      
      if (params$ci_type=="none") {
        
        ci_lo <- rep(0, length(ests))
        ci_hi <- rep(0, length(ests))
        
      } else {
        
        tau_n <- fns$tau_n
        
        # Generate variance scale factor for each point
        tau_ns <- sapply(points, tau_n)
        
        # Construct CIs
        # The 0.975 quantile of the Chernoff distribution occurs at roughly x=1.00
        qnt <- 1.00 # qnorm(0.975, sd=0.52)
        n_orig <- nrow(dat_orig)
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
        
      }
      
    }
    
    # Parse and return results
    res <- list()
    for (p in 1:length(points)) {
      res[[p]] <- list(point=points[p], est=ests[p],
                       ci_lo=ci_lo[p], ci_hi=ci_hi[p])
    }
    return(res)
    
  }
  
}
