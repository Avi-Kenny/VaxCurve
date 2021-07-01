#' Estimate the causal dose-response curve at a point
#' 
#' @param dat Data returned by generate_data()
#' @param estimator Which estimator to use; one of c("G-comp (logistic)",
#'     "G-comp (GAM)","Generalized Grenander")
#' @param params A list, containing the following:
#'   - `boot_reps` Used for G-comp (logistic); number of bootstrap replicates
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
  
  n <- nrow(dat)
  
  if (estimator=="G-comp (logistic)") {
    
    gcomp <- function(dat, a) {
      
      dat %<>% filter(!is.na(a))
      weights <- wts(dat)
      
      # Run model
      model <- glm(y~w1+w2+a, data=dat, weights=weights, family="binomial")
      coefs <- as.numeric(summary(model)$coefficients[,1])
      
      # Run gcomp
      gcomp_i <- apply(
        X = dat,
        MARGIN = 1,
        FUN = function(r) {
          M_i <- c(1,r[["w1"]],r[["w2"]],a)
          return(expit(as.numeric(M_i %*% coefs)))
        }
      )
      
      return(sum(weights*gcomp_i))

    }
    
    # Note: this is inefficient, but fine if length(points) is small
    res <- list()
    for (p in 1:length(points)) {
      
      # Calculate estimate
      est <- gcomp(dat, points[p])
      
      # Get bootstrap CI
      my_stat <- function(dat,indices) {
        d <- dat[indices,]
        return (gcomp(d, points[p]))
      }
      boot_obj <- boot(data=dat, statistic=my_stat, R=params$boot_reps)
      ci <- boot.ci(boot_obj, type="norm", conf=0.95)
      ci_lo <- ci$normal[2]
      ci_hi <- ci$normal[3]
      
      res[[p]] <- list(point=points[p], est=est, ci_lo=ci_lo, ci_hi=ci_hi)
      
    }
    
    return (res)
    
  } else if (estimator=="logistic GAM") {
    
    # !!!!!
    
  } else if (estimator=="Generalized Grenander") {
    
    # Construct theta_n and tao_n, given a dataset
    construct_fns <- function(dat, return_tao_n=T) {
      
      # Construct component functions
      # Note: mu_n=mu2_n since the y values are binary
      grid <- seq(0,1,0.01)
      Phi_n <- construct_Phi_n(dat)
      Phi_n_inv <- construct_Phi_n(dat, type="inverse")
      mu_n <- construct_mu_n(dat, type="logistic")
      mu2_n <- mu_n # construct_mu_n(dat, type="logistic", moment=2)
      sigma2_n <- construct_sigma2_n(mu_n, mu2_n)
      f_aIw_n <- construct_f_aIw_n(dat, type="simple")
      f_a_n <- construct_f_a_n(dat, type=NULL)
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
          deriv_theta_n <- construct_deriv_theta_n(dat, theta_n, grid)
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
      if (params$ci_type=="regular") {
        ci_lo <- ests - 1.00*tao_ns/(n^(1/3))
        ci_hi <- ests + (1.00*tao_ns)/(n^(1/3))
      } else if (params$ci_type=="logit") {
        ci_lo <- expit( logit(ests) - (1.00*tao_ns*deriv_logit(ests))/(n^(1/3)) )
        ci_hi <- expit( logit(ests) + (1.00*tao_ns*deriv_logit(ests))/(n^(1/3)) )
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
