#' Logit function
#' 
#' @param x Numeric input
#' @return Numeric output
logit <- function(x) {
  log(x/(1-x))
}



#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {
  1 / (1+exp(-x))
}



#' Derivative of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_expit <- function(x) {
  exp(x) / ((1+exp(x))^2)
}



#' Derivative of logit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_logit <- function(x) {
  1 / (x-x^2)
}



#' Construct a vectorized/memoized hash table
#' 
#' @param fnc Function to evaluate
#' @param vals A data frame of values to run the function on
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



#' Probability of sampling
#' 
#' @param sampling One of c("iid", "two-phase")
#' @param y Vector `y` of dataset returned by generate_data()
#' @param w1 Vector `w1` of dataset returned by generate_data()
#' @param w2 Vector `w2` of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
Pi <- function(sampling, delta_star, w1, w2) {
  
  if (sampling=="iid") {
    
    return(rep(1, length(w1)))
    
  } else if (sampling=="two-phase") {
    
    pi_w <- function(w1,w2) { expit(w1+w2) }
    return(delta_star + (1-delta_star)*pi_w(w1,w2))
    
  }
  
}



#' Return IP weights
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param scale One of c("none", "stabilized", "mean 1")
#' @param s Weight stabilization factor
#' @return A sum-to-one vector of weights
#' @notes
#'   - Pi is accessed globally
#'   - "Mean 1" weights are mean 1 EXCLUDING the zero weights
wts <- function(dat_orig, scale="stabilized") {
  
  sampling <- attr(dat_orig,"sampling")
  
  weights <- dat_orig$delta /
    Pi(sampling, dat_orig$delta_star, dat_orig$w1, dat_orig$w2)
  
  if (scale=="stabilized") {
    if (sampling=="iid") {
      s <- 1
    } else if (sampling=="two-phase") {
      s <- sum(weights) / nrow(dat_orig)
    }
    weights <- weights / s
  } else if (scale=="mean 1") {
    weights <- ( weights / sum(weights) ) * nrow(filter(dat_orig,!is.na(a)))
  }
  
  return(weights)
  
}



#' Construct conditional survival estimator S_n
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param type One of c("true", "Cox PH", "KM", "Random Forest")
#' @param csf Logical; if TRUE, estimate the conditional survival
#'     function of the censoring distribution instead
#' @return Conditional density estimator function
construct_S_n <- function(dat_orig, vals, type, csf=FALSE) {
  
  # Construct weights
  n_orig <- nrow(dat_orig)
  dat_orig$weights_m1 <- wts(dat_orig, scale="mean 1")
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  
  if (csf) { dat$delta_star <- 1 - dat$delta_star }
  
  if (type=="true") {
    
    surv_true <- L$surv_true
    alpha_3 <- L$alpha_3
    if (csf) {
      lambda <- L$lambda2
      v <- C$v2
    } else {
      lambda <- C$lambda
      v <- C$v
    }
    
    fnc <- function(t, w1, w2, a) {
      if (L$surv_true=="Cox PH") {
        lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1
        return(exp(-1*lambda*(t^v)*exp(lin)))
      } else if (L$surv_true=="complex") {
        lin <- C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + alpha_3*w2*a - 1
        return(exp(-1*lambda*(t^v)*exp(lin)))
      }
    }
    
  } else if (type=="Cox PH") {
    
    # Fit Cox model
    model <- coxph(
      Surv(y_star,delta_star)~w1+w2+a,
      data = dat,
      weights = weights_m1
    )
    coeffs <- model$coefficients
    
    # Get cumulative hazard estimate
    bh <- basehaz(model, centered=FALSE)
    
    # Pre-calculate H_0 vector
    H_0 <- c()
    for (t in 0:C$t_e) {
      index <- which.min(abs(bh$time-t))
      H_0[t+1] <- bh$hazard[index]
    }
    
    fnc <- function(t, w1, w2, a) {
      return(exp(-1*H_0[t+1]*exp(coeffs[[1]]*w1+coeffs[[2]]*w2+coeffs[[3]]*a)))
    }
    
  } else if (type=="Random Forest") {
    
    model <- rfsrc(
      Surv(y_star,delta_star)~w1+w2+a,
      data = dat,
      ntree = 500,
      mtry = 2,
      nodesize = 100,
      splitrule = "logrank",
      nsplit = 0,
      case.wt = dat$weights
      # samptype = "swr"
    )
    
    newX <- subset(filter(vals, t==0), select=-c(t))
    pred <- predict(model, newdata=newX)
    
    fnc <- function(t, w1, w2, a) {
      r1 <- which(abs(w1-newX$w1)<1e-10)
      r2 <- which(abs(w2-newX$w2)<1e-10)
      r3 <- which(abs(a-newX$a)<1e-10)
      row <- intersect(r1,intersect(r2,r3))
      col <- which.min(abs(t-pred$time.interest))
      return(pred$survival[row,col])
    }
    
  } else if (type=="Random Forest Ted") {
    
    method <- "survSL.rfsrc"
    
    newX <- subset(filter(vals, t==0), select=-c(t))
    new.times <- unique(vals$t)

    srv <- survSuperLearner(
      time = dat$y_star,
      event = dat$delta_star,
      X = subset(dat, select=c(w1,w2,a)),
      newX = newX,
      new.times = new.times,
      event.SL.library = c(method),
      cens.SL.library = c(method),
      obsWeights = dat$weights,
      control = list(
        initWeightAlg = method,
        max.SL.iter = 10
      )
    )

    fnc <- function(t, w1, w2, a) {
      r1 <- which(abs(w1-newX$w1)<1e-10)
      r2 <- which(abs(w2-newX$w2)<1e-10)
      r3 <- which(abs(a-newX$a)<1e-10)
      row <- intersect(r1,intersect(r2,r3))
      col <- which(t==new.times)
      return(srv$event.SL.predict[row,col])
    }
    
  }
  
  round_args <- c(-log10(C$appx$t_e), -log10(C$appx$w1b), 0, -log10(C$appx$a))
  return (create_htab(fnc, vals, round_args=round_args))
  
}



#' Construct g-computation estimator function of theta_0
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @return G-computation estimator of theta_0
construct_gcomp_n <- function(dat_orig, vals, S_n) {
  
  fnc <- function(a) {
    1 - mean(S_n(C$t_e,dat_orig$w1,dat_orig$w2,a))
  }
  
  return (create_htab(fnc, vals, round_args=c(-log10(C$appx$a))))
  
}



#' Construct derivative estimator theta'_n
#' 
#' 
#' @param theta_n An estimator of theta_0 (usually theta_n or gcomp_n)
#' @param type One of c("gcomp", "linear", "spline")
construct_deriv_theta_n <- function(theta_n, type) {
  
  if (type=="linear") {
    
    # Estimate entire function on grid
    grid <- seq(0,1,0.01)
    theta_ns <- theta_n(grid)
    
    grid_width <- grid[2] - grid[1]
    points_x <- c(grid[1])
    points_y <- c(theta_ns[1])
    for (i in 2:length(grid)) {
      if (theta_ns[i]-theta_ns[i-1]!=0) {
        points_x <- c(points_x, grid[i]-(grid_width/2))
        points_y <- c(points_y, mean(c(theta_ns[i],theta_ns[i-1])))
      }
    }
    points_x <- c(points_x, grid[length(grid)])
    points_y <- c(points_y, theta_ns[length(grid)])
    points_sl <- c()
    for (i in 2:length(points_x)) {
      slope <- (points_y[i]-points_y[i-1]) /
        (points_x[i]-points_x[i-1])
      points_sl <- c(points_sl, slope)
    }
    
    fnc <- function(x) {
      if (x==0) {
        index <- 1
      } else {
        index <- which(x<=points_x)[1]-1
      }
      return(max(points_sl[index],0))
    }
    
  } else if (type=="line") {
    
    theta_n_left <- theta_n(0.1)
    theta_n_right <- theta_n(0.9)
    
    fnc <- function(x) {
      (theta_n_right-theta_n_left)/0.8
    }
    
  } else if (type=="spline") {
    
    # Estimate entire function on grid
    grid <- seq(0,1,0.01)
    theta_ns <- theta_n(grid)
    
    # Identify jump points of step function
    # jump_points <- c(0)
    jump_points <- c()
    for (i in 2:length(grid)) {
      if (theta_ns[i]!=theta_ns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    # jump_points <- c(jump_points,grid[length(grid)])
    
    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)
    
    # Fit cubic smoothing spline
    theta_n_smoothed <- smooth.spline(x=midpoints, y=theta_n(midpoints))
    
    # Construct derivative function
    fnc <- function(x) {
      
      width <- 0.3
      x1 <- x - width/2
      x2 <- x + width/2
      
      if (x1<0) {
        x2 <- x2 - x1
        x1 <- 0
      }
      if (x2>1) {
        x1 <- x1 - x2 + 1
        x2 <- 1
      }
      
      y1 <- predict(theta_n_smoothed, x=x1)$y
      y2 <- predict(theta_n_smoothed, x=x2)$y
      
      return(max((y2-y1)/width,0))
      
    }
    
  } else if (type=="m-spline") {
    
    # Estimate entire function on grid
    grid <- seq(0,1,0.01)
    theta_ns <- theta_n(grid)
    
    # Identify jump points of step function
    # jump_points <- c(0)
    jump_points <- c()
    for (i in 2:length(grid)) {
      if (theta_ns[i]!=theta_ns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    # jump_points <- c(jump_points,grid[length(grid)])
    
    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)
    
    # Fit monotone cubic smoothing spline
    theta_n_smoothed <- splinefun(
      x = midpoints,
      y = theta_n(midpoints),
      method = "monoH.FC"
    )
    
    # Construct derivative function
    fnc <- function(x) {
      
      width <- 0.3
      x1 <- x - width/2
      x2 <- x + width/2
      
      if (x1<0) {
        x2 <- x2 - x1
        x1 <- 0
      }
      if (x2>1) {
        x1 <- x1 - x2 + 1
        x2 <- 1
      }
      
      y1 <- theta_n_smoothed(x1)
      y2 <- theta_n_smoothed(x2)
      # y1 <- predict(theta_n_smoothed, x=x1)$y
      # y2 <- predict(theta_n_smoothed, x=x2)$y
      
      return(max((y2-y1)/width,0))
      
    }
    
  } else if (type=="gcomp") {
    
    fnc <- function(x) {
      
      # Set derivative appx x-coordinates
      width <- 0.2
      p1 <- x - width/2
      p2 <- x + width/2
      if (p1<0) {
        p2 <- p2 - p1
        p1 <- 0
      }
      if (p2>1) {
        p1 <- p1 - p2 + 1
        p2 <- 1
      }
      c(p1,p2)
      
      return(max((theta_n(p2)-theta_n(p1))/width,0))
      
    }
    
  }
  
  return(memoise(Vectorize(fnc)))

}



#' Construct tau_n Chernoff scale factor function
#' 
#' @param deriv_theta_n A derivative estimator returned by
#'     deriv_theta_n()
#' @param gamma_n Nuisance function estimator returned by construct_gamma_n()
#' @param f_a_n Density estimator returned by construct_f_a_n()
#' @return Chernoff scale factor estimator function
construct_tau_n <- function(deriv_theta_n, gamma_n, f_a_n) {
  
  return(memoise(Vectorize(function(x){
    (4*deriv_theta_n(x)*f_a_n(x)*gamma_n(x))^(1/3)
  })))
  
}



#' Construct gamma_n nuisance estimator function
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param type Type of regression; one of c("linear", "cubic")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, vals, type, omega_n, f_aIw_n, f_a_n,
                              f_a_delta1_n) {
  
  # Estimate marginal delta probability
  delta_prob <- mean(dat_orig$delta)
  
  # Construct weights
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  # Construct pseudo-outcomes
  dat %<>% mutate(
    po = (
      (weights*omega_n(w1,w2,a,y_star,delta_star)) / f_aIw_n(a,w1,w2)
    )^2
  )
  
  # Remove outliers to prevent errors (revisit this)
  cutoff <- as.numeric(quantile(dat$po, 0.75) + 10000*iqr(dat$po))
  dat %<>% filter(is.finite(po))
  dat %<>% filter(po<cutoff)
  
  # Run regression
  if (type=="cubic") {

    model <- lm(po~a+I(a^2)+I(a^3), data=dat)
    coeff <- as.numeric(model$coefficients)
    
    fnc <- function(x) {
      delta_prob*(f_a_delta1_n(x)/f_a_n(x)) * (
        coeff[1] + coeff[2]*x + coeff[3]*(x^2) + coeff[4]*(x^3)
      )
    }
    
  } else if (type=="kernel") {
    
    ks <- ksmooth(
      x = dat$a,
      y = dat$po,
      kernel = "normal",
      bandwidth = 0.2,
      range.x = c(0,1),
      x.points = vals$a
    )
    
    fnc <- function(x) {
      index <- which.min(abs(x-ks$x))
      return(ks$y[index])
    }
    
  }
  
  return (create_htab(fnc, vals, round_args=-log10(C$appx$a)))
  
}



#' Construct estimator of conditional density f_{A|W}
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param type One of c("parametric", "binning")
#' @param k Number of bins for the binning estimator (if k=0, then the number of
#'     bins will be selected via cross-validation); ignored for the parametric
#'     estimator
#' @param delta1 Compute the density conditional on Delta=1
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of A is [0,1]
construct_f_aIw_n <- function(dat_orig, vals, type, k=0, delta1=FALSE) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  if (delta1) {
    dat$weights <- rep(1, nrow(dat))
  }
  
  if (type=="true") {
    
    # Note: this is not accurate if edge!="none"
    
    if (L$distr_A=="Unif(0,1)") {
      fnc <- function(a, w1, w2) { 1 }
    } else if (L$distr_A=="Beta(1.5+w1,1.5+w2)") {
      fnc <- function(a, w1, w2) {
        dbeta(a, shape1=1.5+w1, shape2=1.5+w2)
      }
    } else if (L$distr_A=="N(0.5,0.01)") {
      fnc <- function(a, w1, w2) {
        dnorm(a, mean=0.5, sd=0.1)
      }
    } else if (L$distr_A=="N(0.5,0.04)") {
      fnc <- function(a, w1, w2) {
        dnorm(a, mean=0.5, sd=0.2)
      }
    } else if (L$distr_A=="N(0.4+0.2w1+0.1w2,0.01)") {
      fnc <- function(a, w1, w2) {
        dnorm(a, mean=0.4+0.2*w1+0.1*w2, sd=0.1)
      }
    }
    
    
  } else if (type=="parametric") {
    
    # Set up weighted likelihood
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:nrow(dat)), function(i) {
        shape1 <- par[1] + par[2]*dat$w1[i]
        shape2 <- par[3] + par[4]*dat$w2[i]
        loglik <- dbeta(ifelse(dat$a[i]==0,1e-4,dat$a[i]),
                        shape1=shape1, shape2=shape2, log=TRUE)
        return(loglik*dat$weights[i])
      }))
      
      return(-1*sum_loglik)
      
    }
    
    # Run optimizer
    opt <- optim(par=c(a1=1, a2=0.1, a3=1, a4=0.1), fn=wlik)
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: optim() did not converge")
    }
    
    fnc <- function(a, w1, w2){
      shape1 <- opt$par[1] + opt$par[2]*w1
      shape2 <- opt$par[3] + opt$par[4]*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    }
    
  } else if (type=="binning") {
    
    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # par[1] through par[k-1] are the hazard components for the bins 1 to k-1
    # par[k] and par[k+1] are the coefficients for W1 and W2
    create_dens <- function(k, dat) {
      
      # Cut points
      alphas <- seq(0, 1, length.out=k+1)
      
      # Density for a single observation
      dens <- Vectorize(function(a, w1, w2, par) {
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        hz <- expit(
          par[c(1:(ifelse(bin==k,k-1,bin)))] + par[k]*w1 + par[k+1]*w2
        )
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        dens <- k*p1*p2
        return(dens)
      }, vectorize.args=c("a","w1","w2"))
      
      # Set up weighted likelihood
      wlik <- function(par) {
        -1 * sum(dat$weights *
                   log(pmax(dens(a=dat$a, w1=dat$w1, w2=dat$w2, par), 1e-8)))
      }
      
      # Run optimizer
      opt <- solnp(
        pars = rep(0.001,k+1),
        fun = wlik
      )
      if (opt$convergence!=0) {
        warning("construct_f_aIw_n: solnp() did not converge")
      }
      
      fnc <- Vectorize(function(a, w1, w2){
        
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        par <- opt$pars
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(par[j] + par[k]*w1 + par[k+1]*w2)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        
        return(k*p1*p2)
        
      })
      
    }
    
    # Select k via cross-validation
    if (k==0) {
      
      # Prep
      n_folds <- 5
      folds <- sample(cut(c(1:length(dat)), breaks=n_folds, labels=FALSE))
      ks <- c(5,10,15,20,25)
      best <- list(k=999, max_log_lik=999)
      
      # Cross-validation
      for (k in ks) {
        
        sum_log_lik <- 0
        for (i in c(1:n_folds)) {
          dat_train <- dat[-which(folds==i),]
          dat_test <- dat[which(folds==i),]
          dens <- create_dens(k, dat_train)
          sum_log_lik <- sum_log_lik + sum(log(
            dens(dat_test$a, dat_test$w1, dat_test$w2)
          ))
        }
        
        if (sum_log_lik>best$max_log_lik || best$max_log_lik==999) {
          best$k <- k
          best$max_log_lik <- sum_log_lik
        }
        
      }
      
      k <- best$k
      
    }
    
    fnc <- create_dens(k, dat)
    
  }
  
  round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return (create_htab(fnc, vals, round_args=round_args))
  
}



#' Construct estimator of marginal density f_A
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return Marginal density estimator function
construct_f_a_n <- function(dat_orig, vals, f_aIw_n) {
  
  fnc <- function(a) {
    mean(f_aIw_n(a,dat_orig$w1,dat_orig$w2))
  }
  
  return (create_htab(fnc, vals, round_args=-log10(C$appx$a)))
  
}



#' Construct density ratio estimator g_n
#' 
#' @param vals Dataframe of values to run function on
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
construct_g_n <- function(vals, f_aIw_n, f_a_n) {
  
  fnc <- function(a,w1,w2) {
    f_aIw_n(a,w1,w2) / f_a_n(a)
  }
  
  round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return (create_htab(fnc, vals, round_args=round_args))
  
}



#' Construct estimator of nuisance influence function omega_n
#' 
#' @param vals Dataframe of values to run function on
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param Sc_n Conditional censoring survival function estimator returned by
#'     construct_S_n
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(vals, S_n, Sc_n) {
  
  # First, construct cumulative hazard estimator
  H_n <- function(t, w1, w2, a) {
    -1 * log(S_n(t, w1, w2, a))
  }
  
  fnc <- function(w1,w2,a,y_star,delta_star) {
    
    k <- round(min(y_star,C$t_e))
    if (k==0) {
      integral <- 0
    } else {
      i <- c(1:k)
      # i <- c(1:m)
      # k <- min(y_star,C$t_e)
      
      integral <- 0.5 * sum(
        ( H_n(i,w1,w2,a) - H_n(i-1,w1,w2,a) ) * (
          ( S_n(i,w1,w2,a) * Sc_n(i,w1,w2,a) )^-1 +
            ( S_n(i-1,w1,w2,a) * Sc_n(i-1,w1,w2,a))^-1
        )
      )
      
      integral_righthand <- sum(
        ( H_n(i,w1,w2,a) - H_n(i-1,w1,w2,a) ) *
          ( (S_n(i,w1,w2,a)) * Sc_n(i,w1,w2,a) )^-1
      )
      
      # integral_diffgrid <- 0.5 * sum(
      #   ( H_n((i*k)/m,w1,w2,a) - H_n(((i-1)*k)/m,w1,w2,a) ) * (
      #     (S_n((i*k)/m,w1,w2,a))^-1 * (Sc_n((i*k)/m,w1,w2,a))^-1 +
      #       (S_n(((i-1)*k)/m,w1,w2,a))^-1 * (Sc_n(((i-1)*k)/m,w1,w2,a))^-1
      #   )
      # )
      
    }
    
    return(S_n(C$t_e,w1,w2,a) * (
        (delta_star * as.integer(y_star<=C$t_e)) /
          (S_n(k,w1,w2,a) * Sc_n(k,w1,w2,a)) -
        integral
    ))
    
  }
  
  return (create_htab(fnc, vals))
  
}



#' Construct nuisance estimator eta_n
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @return Estimator function of nuisance eta_0
construct_eta_n <- function(dat_orig, vals, S_n) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  fnc <- function(x,w1,w2) {
    (1/n_orig) * sum(
      weights * as.integer(dat$a<=x) *
        (1-S_n(C$t_e,w1,w2,dat$a))
    )
  }
  
  round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return (create_htab(fnc, vals, round_args=round_args))
  
}



#' lambda estimator
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @return Value of lambda
lambda <- function(dat_orig, k, G) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  lambda <- (1/n_orig) * sum(
    weights * (G(dat$a))^k
  )
  return(lambda)
  
}



#' Construct Gamma_n primitive one-step estimator
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param S_n A conditional survival function returned by construct_S_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @param type One of c("one-step", "plug-in")
#' @return Gamma_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Gamma_n <- function(dat_orig, vals, omega_n, S_n, g_n,
                              type="one-step") {
  
  dat_orig$weights <- wts(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- filter(dat_orig,!is.na(a))
  n <- nrow(dat)
  
  i_long <- rep(c(1:n), each=n)
  j_long <- rep(c(1:n), times=n)
  a_i_long <- dat$a[i_long]
  w1_i_long <- dat$w1[i_long]
  w1_j_long <- dat$w1[j_long]
  w2_i_long <- dat$w2[i_long]
  w2_j_long <- dat$w2[j_long]
  delta_star_i_long <- dat$delta_star[i_long]
  delta_star_j_long <- dat$delta_star[j_long]
  weights_i_long <- dat$weights[i_long]
  weights_j_long <- dat$weights[j_long]
  
  if (type=="one-step") {
    
    subpiece_1a <- dat$weights * ( 1 + (
      omega_n(dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star) /
        g_n(dat$a,dat$w1,dat$w2)
    ) )
    subpiece_2a <- (weights_i_long*weights_j_long) *
      S_n(C$t_e,w1_j_long,w2_j_long,a_i_long)
    
    fnc <- function(x) {
      
      subpiece_1b <- as.integer(dat$a<=x)
      piece_1 <- (1/n_orig) * sum(subpiece_1a*subpiece_1b)
      
      subpiece_2b <- as.integer(a_i_long<=x)
      piece_2 <- (1/(n_orig^2)) * sum(subpiece_2a*subpiece_2b)
      
      return(piece_1-piece_2)
      
    }
    
  }
  
  if (type=="plug-in") {
    
    fnc <- function(x) {
      return(
        (1/n_orig^2) * sum(
          (weights_i_long*weights_j_long) *
            as.integer(a_i_long<=x) *
            (1 - S_n(C$t_e,w1_j_long,w2_j_long,a_i_long))
        )
      )
    }
    
  }
  
  return (create_htab(fnc, vals, round_args=-log10(C$appx$a)))
  
}






#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param type One of c("ecdf", "inverse").
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
#'   - This accesses wts() globally
construct_Phi_n <- function (dat_orig, type="ecdf") {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  dat %<>% arrange(a)
  vals_x <- unique(dat$a)
  vals_y <- c()
  
  for (j in 1:length(vals_x)) {
    indices <- which(dat$a==vals_x[j])
    weights_j <- dat$weights[indices]
    new_y_val <- (1/n_orig) * sum(weights_j)
    vals_y <- c(vals_y, new_y_val)
  }
  vals_y <- cumsum(vals_y)
  
  if (type=="ecdf") {
    rval <- approxfun(vals_x, vals_y, method="constant", yleft=0,
                      yright=1, f=0, ties="ordered")
  } else if (type=="inverse") {
    rval <- approxfun(vals_y, vals_x, method="constant", yleft=min(vals_x),
                      yright=max(vals_x), f=1, ties="ordered")
  }
  return(rval)
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_rho_n <- function(dat_orig, Phi_n) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- dat$weights
  a_j <- dat$a
  
  return(memoise(Vectorize(function(a) {
    
    return((1/n_orig) * sum(
      weights_j * (Phi_n(a_j)^3) * (as.integer(a<=a_j) - Phi_n(a_j))
    ))
    
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_xi_n <- function(Phi_n, lambda_2, lambda_3) {
  
  return(memoise(Vectorize(function(a_i,a_j) {
    return(
      (2*as.integer(a_i<=a_j) - 3*Phi_n(a_j))*Phi_n(a_j)*lambda_2 +
      (2*Phi_n(a_j) - as.integer(a_i<=a_j))*lambda_3
    )
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_1 <- function(dat_orig, Gamma_n, Phi_n, xi_n, rho_n,
                                lambda_2, lambda_3) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- dat$weights
  a_j <- dat$a
  
  return(memoise(Vectorize(function(w1,w2,y_star,delta_star,a) {
    
    piece_1 <- (1/n_orig) * sum(
      weights_j * (xi_n(a,a_j)) * Gamma_n(a_j)
    )
    
    piece_2 <- (lambda_2*(Phi_n(a)^2) - lambda_3*Phi_n(a)) * Gamma_n(a)
    
    return(piece_1+piece_2)
    
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Gamma <- function(dat_orig, omega_n, g_n, gcomp_n, eta_n,
                                    Gamma_n) {
  
  return(memoise(Vectorize(function(x,w1,w2,y_star,delta_star,a) {
    as.integer(a<=x)*(
      (omega_n(w1,w2,a,y_star,delta_star)/g_n(a,w1,w2)) +
        gcomp_n(a)
    ) +
      eta_n(x,w1,w2) - 
      2*Gamma_n(x)
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_2 <- function(dat_orig, Phi_n, infl_fn_Gamma, lambda_2,
                                lambda_3) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- dat$weights
  a_j <- dat$a
  
  return(memoise(Vectorize(function(w1,w2,y_star,delta_star,a) {
    (1/n_orig) * sum(weights_j * (
      ( lambda_2*(Phi_n(a_j)^2) - lambda_3*Phi_n(a_j) ) *
        infl_fn_Gamma(a_j,w1,w2,y_star,delta_star,a)
    ))
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
beta_n_var_hat <- function(dat_orig, infl_fn_1, infl_fn_2) {
  
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  
  b_sum <- 0
  for (i in c(1:nrow(dat))) {
    b_sum <- b_sum + (dat$weights[i] * (
      infl_fn_1(dat$w1[i], dat$w2[i], dat$y_star[i],
                dat$delta_star[i], dat$a[i]) +
        infl_fn_2(dat$w1[i], dat$w2[i], dat$y_star[i],
                  dat$delta_star[i], dat$a[i])
    ))^2
  }
  
  return( (1/n_orig)*sum(b_sum) )
  
}



#' Construct dataframes of values to pre-compute functions on
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param appx Approximation spec used for grids (stored in C$appx)
#' @return A list of dataframes
create_val_list <- function(dat_orig, appx) {
  
  dat <- dat_orig %>% filter(!is.na(a))
  
  omega <- subset(dat, select=c(w1,w2,a,y_star,delta_star))
  omega_copy <- omega
  omega_copy$a <- 0
  omega <- rbind(omega, omega_copy)
  
  return(list(
    A = data.frame(a=dat$a),
    AW = data.frame(a=dat$a, w1=dat$w1, w2=dat$w2),
    A_grid = data.frame(a=seq(0,1,appx$a)),
    W_grid = expand.grid(w1=seq(0,1,appx$w1), w2=c(0,1)),
    AW_grid = expand.grid(a=seq(0,1,appx$a), w1=seq(0,1,appx$w1),
                          w2=c(0,1)),
    S_n = expand.grid(t=seq(0,C$t_e,appx$t_e), w1=seq(0,1,appx$w1b),
                      w2=c(0,1), a=seq(0,1,appx$a)),
    omega = omega
  ))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_Gamma_cf_k <- function(dat_train, dat_test, vals, omega_n, g_n,
                                 gcomp_n, eta_n) {
  
  n_test <- nrow(dat_test)
  dat_test$weights <- wts(dat_test)
  d1 <- dat_test %>% filter(!is.na(a))
  weights_1 <- d1$weights
  
  n_train <- nrow(dat_train)
  dat_train$weights <- wts(dat_train)
  d2 <- dat_train %>% filter(!is.na(a))
  weights_2 <- d2$weights
  
  fnc <- function(x) {
    
    piece_1 <- (1/n_test) * sum(weights_1 * (
      as.integer(d1$a<=x) *
        (
          (omega_n(d1$w1,d1$w2,d1$a,d1$y_star,d1$delta_star) /
             g_n(d1$a,d1$w1,d1$w2)
          ) +
            gcomp_n(d1$a)
        ) +
        eta_n(x,d1$w1,d1$w2)
    ))
    
    piece_2 <- (1/n_train) * sum(weights_2 * as.integer(d2$a<=x)*gcomp_n(d2$a))
    
    return(piece_1-piece_2)
    
  }
  
  return (create_htab(fnc, vals, round_args=-log10(C$appx$a)))
  
}



#' Construct cross-fitted Gamma_0 estimator
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param params The same params object passed to est_curve
#' @param vlist A list of dataframes returned by create_val_list()
#' @return A cross-fitted Gamma_0 estimator
construct_Gamma_cf <- function(dat_orig, params, vlist) {
  
  # Prep for cross-fitting
  Gamma_cf_k <- list()
  n_orig <- nrow(dat_orig)
  rows <- c(1:n_orig)
  folds <- sample(cut(rows, breaks=params$cf_folds, labels=FALSE))
  
  # Loop through folds
  for (k in 1:params$cf_folds) {
    
    # Split data
    dat_train <- dat_orig[-which(folds==k),]
    dat_test <- dat_orig[which(folds==k),]
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat_train)
    Phi_n_inv <- construct_Phi_n(dat_train, type="inverse")
    S_n <- construct_S_n(dat_train, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_train, vlist$S_n, type=params$S_n_type, csf=TRUE)
    gcomp_n <- construct_gcomp_n(dat_train, vlist$A_grid, S_n)
    eta_n <- construct_eta_n(dat_train, vlist$AW_grid, S_n)
    f_aIw_n <- construct_f_aIw_n(dat_train, vlist$AW_grid, type=params$g_n_type,
                                 k=15)
    f_a_n <- construct_f_a_n(dat_train, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    
    # Construct K functions
    Gamma_cf_k[[k]] <- construct_Gamma_cf_k(
      dat_train, dat_test, vlist$A_grid, omega_n, g_n, gcomp_n, eta_n
    )
    
    # Remove objects
    rm(Phi_n,Phi_n_inv,S_n,Sc_n,gcomp_n,eta_n,f_aIw_n,f_a_n,g_n,omega_n)
    
  }
  
  # Construct cross-fitted Gamma_n
  return(Vectorize(function(x) {
    mean(sapply(c(1:params$cf_folds), function(k) {
      Gamma_cf_k[[k]](x)
    }))
  }))
  
}



#' Construct propensity score estimator of pi_0 = P(A=0|W=w)
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param type One of c("logistic","")
#' @return Propensity score estimator of pi_0
construct_pi_n <- function(dat_orig, vals, type) {
  
  # Construct weights
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  # Construct indicator I{A=0}
  dat$ind_A0 <- as.integer(dat$a==0)
  
  if (type=="true") {
    
    if (L$edge=="expit") {
      fnc <- function(w1,w2) { expit(w1+w2-3.3) }
    } else if (L$edge=="complex") {
      fnc <- function(w1,w2) { 0.84*w2*pmax(0,1-4*abs(w1-0.5)) }
    } else if (L$edge=="none") {
      fnc <- function(w1,w2) { 0 }
    }
    
  }
  
  if (type=="logistic") {
    
    suppressWarnings({
      model <- glm(
        ind_A0~w1+w2,
        data = dat,
        family = "binomial",
        weights = weights
      )
    })
    coeffs <- model$coefficients
    
    fnc <- function(w1,w2) {
      expit( coeffs[[1]] + coeffs[[2]]*w1 + coeffs[[3]]*w2 )
    }
    
  }
  
  if (type=="SL") {
    
    sl <- SuperLearner(
      Y = dat$ind_A0,
      X = subset(dat, select=c(w1,w2)),
      newX = vals,
      family = binomial(),
      SL.library = "SL.earth", # SL.glm SL.gbm SL.ranger SL.earth
      # SL.library = c("SL.earth", "SL.gam", "SL.ranger"), # SL.glm SL.gbm SL.ranger SL.earth
      obsWeights = weights,
      control = list(saveFitLibrary=FALSE)
    )
    assign("sl", sl, envir=.GlobalEnv)
    
    fnc <- function(w1, w2) {
      r1 <- which(abs(w1-vals$w1)<1e-10)
      r2 <- which(abs(w2-vals$w2)<1e-10)
      index <- intersect(r1,r2)
      return(sl$SL.predict[index])
    }
    
  }
  
  return (create_htab(fnc, vals, round_args=c(-log10(C$appx$w1), 0)))
  
}



#' COmpute one-step estimator of counterfactual survival at A=0
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @return Value of one-step estiamtor
theta_os_n <- function(dat_orig, pi_n, S_n, omega_n) {
  
  # Construct weights
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  # Return estimate
  return(
    1 - (1/n_orig) * sum(weights * (
      S_n(C$t_e,dat$w1,dat$w2,a=0) - (
        (as.integer(dat$a==0)/pi_n(dat$w1,dat$w2)) *
          omega_n(dat$w1,dat$w2,a=0,dat$y_star,dat$delta_star)
      )
    ))
  )
  
}



#' Compute asymptotic variance of one-step estimator theta_os_n
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals Dataframe of values to run function on
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param theta_os_n_est Estimate returned by one-step estimator theta_os_n()
#' @return Asymptotic variance estimate
sigma2_os_n <- function(dat_orig, pi_n, S_n, omega_n, theta_os_n_est) {
  
  # Construct weights
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  # Return estimate
  return(
    (1/n_orig) * sum((weights * (
      S_n(C$t_e,dat$w1,dat$w2,a=0) - (
        (as.integer(dat$a==0)/pi_n(dat$w1,dat$w2)) *
          omega_n(dat$w1,dat$w2,a=0,dat$y_star,dat$delta_star)
      ) -
        theta_os_n_est
    ))^2)
  )
  
}


