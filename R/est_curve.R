#' Estimate the causal dose-response curve at a point
#' 
#' @param dat Data returned by generate_data()
#' @param estimator Which estimator to use; one of c("G-comp", "Generalized
#'     Grenander")
#' @param params A list, containing the following:
#'   - `mu_n_type` Type of regression estimator; one of One of c("Logistic",
#'     "GAM", "Random forest")
#'   - `g_n_type` Type of conditional density ratio estimator; one of
#'     c("parametric", "binning")
#'   - `boot_reps` Used for G-comp; number of bootstrap replicates
#'   - `Phi` Used for Grenander; one of c("marginal","identity") # !!!!! Not yet implemented
#'   - `ci_type` One of c("regular", "logit", "sample split"). "regular" is the
#'     standard approach. "logit" transforms the CIs so that the bounds are in
#'     [0,1]. "sample split" is the Banerjee method (and also returns a
#'     different estimator).
#'   - `m` If params$ci_type=="sample split", the number of splits
#' @param points A vector representing the points to estimate
#' @return A list of lists of the form:
#'     list(list(point=1, est=1, se=1), list(...), ...)

est_curve <- function(dat, estimator, params, points) {
  
  if (estimator=="G-comp") {
    
    # Compute estimates
    mu_n <- construct_mu_n(dat, type=params$mu_n_type)
    gcomp_n <- construct_gcomp(dat, mu_n=mu_n)
    ests <- gcomp_n(points)
    
    # Run bootstrap for SEs
    {
      my_stat <- function(dat,indices) {
        d <- dat[indices,]
        mu_n <- construct_mu_n(d, type=params$mu_n_type)
        gcomp_n <- construct_gcomp(d, mu_n=mu_n)
        return (gcomp_n(points))
      }
      boot_obj <- boot(data=dat, statistic=my_stat, R=params$boot_reps)
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
    
  } else if (estimator=="Generalized Grenander") {
    
    # Construct theta_n and tao_n, given a dataset
    construct_fns <- function(dat, return_tao_n=T) {
      
      # Construct component functions
      # Note: mu_n=mu2_n since the y values are binary
      grid <- seq(0,1,0.01)
      Phi_n <- construct_Phi_n(dat)
      Phi_n_inv <- construct_Phi_n(dat, type="inverse")
      mu_n <- construct_mu_n(dat, type=params$mu_n_type)
      gcomp_n <- construct_gcomp(dat, mu_n=mu_n)
      mu2_n <- mu_n
      sigma2_n <- construct_sigma2_n(mu_n, mu2_n)
      f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type)
      f_a_n <- construct_f_a_n(dat, f_aIw_n=f_aIw_n)
      g_n <- construct_g_n(f_aIw_n, f_a_n)
      Gamma_n <- construct_Gamma_n(dat, mu_n, g_n)
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
      
      if (return_tao_n==F) {
        
        return(list(theta_n=theta_n))
        
      } else {
        
        # Construct tao_n
        if (return_tao_n==T) {
          deriv_theta_n <- construct_deriv_theta_n(gcomp_n)
          tao_n <- function(x) {
            (4*deriv_theta_n(x) * mean(apply(
              X = dat,
              MARGIN = 1,
              FUN = function(r) {
                sigma2_n(x,r[["w1"]],r[["w2"]])/f_aIw_n(x,r[["w1"]],r[["w2"]])
              }
            )))^(1/3)
          }
        }
        
        return(list(theta_n=theta_n, tao_n=tao_n))
        
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
          dat_split <- dat[c(splits[x,1]:splits[x,2]),]
          theta_n <- construct_fns(dat=dat_split, return_tao_n=F)$theta_n
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
      
      # Construct theta_n and tao_n
      fns <- construct_fns(dat=dat, return_tao_n=T)
      theta_n <- fns$theta_n
      tao_n <- fns$tao_n
      
      # Generate estimates and variance scale factors for each point
      ests <- sapply(points, theta_n)
      tao_ns <- sapply(points, tao_n)
      
      # Construct CIs
      # The 0.975 quantile of the Chernoff distribution occurs at roughly x=1.00
      qnt <- 1.00 # qnorm(0.975, sd=0.52)
      n <- nrow(filter(dat, !is.na(a)))
      if (params$ci_type=="regular") {
        ci_lo <- ests - (qnt*tao_ns)/(n^(1/3))
        ci_hi <- ests + (qnt*tao_ns)/(n^(1/3))
      } else if (params$ci_type=="logit") {
        ci_lo <- expit( logit(ests) - (qnt*tao_ns*deriv_logit(ests))/(n^(1/3)) )
        ci_hi <- expit( logit(ests) + (qnt*tao_ns*deriv_logit(ests))/(n^(1/3)) )
      }
      
    }
    
    # Parse and return results
    res <- list()
    for (p in 1:length(points)) {
      res[[p]] <- list(point=points[p], est=ests[p],
                       ci_lo=ci_lo[p], ci_hi=ci_hi[p])
    }
    return(res)
    
  } else {
    stop("estimator incorrectly specified")
  }
  
}
