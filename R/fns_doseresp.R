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



#' Given a function, return a "superfunction" that is:
#'     (1) memoised,
#'     (2) vectorized for both scalar and vector base function inputs,
#'     (3) allows for automatic argument rounding (for approximations), and
#'     (4) allows for pre-computation on a given set of values
#' 
#' @param fnc The base function
#' @param aux Any auxillary data that the function needs to access in its
#'     environment (e.g. a model object)
#' @param vec Vectorization pattern. If the function inputs are all scalars, use
#'     vec=TRUE. If the function inputs are a combination of scalars and
#'     vectors, use a vector of 0s, 1s, and 2s corresponding to whether each
#'     input is non-vectorized (0), a scalar (1) or a vector (2). For example,
#'     for a function f(A,B,C,D) that takes scalars A and B, a vector C, and
#'     non-vectorized argument D, use vec=c(1,1,2,0). If the entire function
#'     should not be vectorized, use vec=FALSE.
#' @param vals Values to precompute the function on. This should be a list
#'     with keys corresponding to the function arguments and values supplied
#'     in the same way they would be supplied to the function
#' @return A function with the above modifications. This can be called in the
#'     same way as the original function or in vectorized format. For example,
#'     the function f(A,B,C,D) described above can be called via
#'     f(2,4,c(1,2,3),9) (same as original) or via
#'     f(c(2,3),c(4,5),data.frame(x=c(1,5),y=c(2,6),z=c(3,7)),9) (vectorized).
#'     The number of rows of each data frame should be the same as the number of
#'     rows of each vector. The data frame column names are ignored.
#' @notes
#'   - !!!!! vec=c(...,0,...) option not yet fully implemented; workaround is to
#'     wrap the argument in a list
#'   - The order of inputs matters and argument names are ignored.
#'   - This function only works on numeric data and will coerce all data to be
#'     numeric.
construct_superfunc <- function(fnc, aux=NA, vec=TRUE, vals=NA, rnd=NA) {
  
  htab <- new.env()
  ..new_fnc <- function() {
    
    # !!!!! add .. in front of all local variables
    
    ..e <- parent.env(environment())
    args <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    arg_names <- names(args)
    args1 <- args[[1]]
    if (identical(..e$vec[1],T)) {
      ..e$vec <- rep(1, length(args))
    }
    if (!is.numeric(..e$vec) || any(!(..e$vec %in% c(0,1,2)))) {
      stop("vec must be a vector of 0s, 1s, and 2s")
    }
    if (max(..e$vec)!=1) {
      index <- min(which(..e$vec==2))
      args2 <- args[[index]]
    }
    if (max(..e$vec)==1 && length(args1)==1) {
      vec_mode <- F
    } else if (max(..e$vec)==1 && length(args1)>1) {
      vec_mode <- T
    } else if (max(..e$vec)!=1 && class(args2)!="data.frame") {
      vec_mode <- F
    } else if (max(..e$vec)!=1 && class(args2)=="data.frame") {
      vec_mode <- T
    } else {
      stop("Function inputs incorrect")
    }
    
    # First, we create a key and check to see if a key/value pair is already
    #     stored. If it is, we retrieve the value. If it is not, we run the
    #     function and store the value. This is done separately based on
    #     whether the current function call is in "vectorized mode".
    # Note: memoising will not work if the function itself returns NULL
    if (vec_mode) {
      
      # Vectorized operation
      for (j in 1:length(..e$vec)) {
        if (..e$vec[j]==2) {
          args[[j]] <- as.list(as.data.frame(t(eval(args[[j]])),
                                                row.names=NA))
        }
      }
      FUN <- function(...) {
        keylist <- sapply(arg_names, function(arg_name) {
          as.numeric(list(...)[[arg_name]])
        })
        key <- paste(keylist, collapse=";")
        val <- ..e$htab[[key]]
        if (is.null(val)) {
          val <- do.call(..e$fnc, list(...))
          ..e$htab[[key]] <- val
        }
        return(val)
      }
      dovec <- (..e$vec!=0)
      res <- do.call(
        mapply,
        c(FUN=FUN, args[dovec], MoreArgs=list(args[!dovec]))
      )
      
    } else {
      
      # Non-vectorized operation
      keylist <- sapply(arg_names, function(arg_name) {
        as.numeric(args[[arg_name]])
      })
      key <- paste(keylist, collapse=";")
      val <- ..e$htab[[key]]
      if (is.null(val)) {
        val <- do.call(..e$fnc, keylist)
        ..e$htab[[key]] <- val
      }
      res <- val
      
    }
    
    # Return value
    return(res)
    
  }
  
  # Set formals and set up environment
  formals(..new_fnc) <- formals(fnc)
  f_env <- new.env(parent=environment(fnc))
  f_env$vec <- vec
  f_env$rnd <- rnd
  f_env$htab <- htab
  f_env$aux <- aux
  f_env$fnc <- fnc
  environment(..new_fnc) <- f_env

  # Run function on vals list
  if (is.list(vals)) {
    do.call(..new_fnc, vals)
  }
  
  return(..new_fnc)
  
}



#' Probability of sampling
#' 
#' @param sampling One of c("iid", "two-phase (6%)", "two-phase (72%)")
#' @param delta_star Component of dataset returned by generate_data()
#' @param y_star Component of dataset returned by generate_data()
#' @param w Component of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
Pi <- function(sampling, delta_star, y_star, w) {
  
  # if (!(sampling %in% c("iid", "two-phase (6%)", "two-phase (72%)"))) {
  #   stop("sampling incorrectly specified")
  # }
  
  if (sampling=="two-phase (10% compl random)") {
    return(rep(0.1, length(y_star)))
  }
  
  if (sampling=="iid") {
    return(rep(1, length(y_star)))
  } else {
    if (sampling=="two-phase (6%)") {
      pi_w <- function(w) { expit(w$w1+w$w2-3.85) }
    } else if (sampling=="two-phase (6% random)") {
      pi_w <- function(w) { 0.06 }
    } else if (sampling=="two-phase (72%)") {
      pi_w <- function(w) { expit(w$w1+w$w2-0.1) }
    }
    ev <- as.integer(delta_star==1 & y_star<=200)
    return(ev + (1-ev)*pi_w(w))
  }
  
}



#' Return IP weights
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param scale One of c("none", "stabilized")
#' @return A sum-to-one vector of weights
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
wts <- function(dat_orig, scale="stabilized") {
  
  sampling <- attr(dat_orig,"sampling")
  
  weights <- dat_orig$delta /
    Pi(sampling, dat_orig$delta_star, dat_orig$y_star, dat_orig$w)
  
  if (scale=="stabilized") {
    if (sampling=="iid") {
      s <- 1
    } else {
      s <- sum(weights) / length(dat_orig$delta)
    }
    weights <- weights / s
  }
  
  return(weights)
  
}



#' Filter dat_orig according to indices
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param indices Indices to filter dataset by
#' @return Filtered subsample of dataset
ss <- function(dat_orig, indices) {
  
  i <- indices
  
  return(list(
    w = dat_orig$w[i,],
    a = dat_orig$a[i],
    delta = dat_orig$delta[i],
    y_star = dat_orig$y_star[i],
    delta_star = dat_orig$delta_star[i],
    weights = dat_orig$weights[i]
  ))
  
}



#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param which One of c("ecdf", "inverse")
#' @param type Defaults to "estimated". Override with "true" for debugging.
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
construct_Phi_n <- function (dat, which="ecdf", type="estimated") {
  
  if (type!="true") {
    
    n_orig <- sum(dat$weights)
    a <- sort(dat$a)
    vals_x <- unique(a)
    vals_y <- c()
    
    for (j in 1:length(vals_x)) {
      indices <- which(a==vals_x[j])
      weights_j <- dat$weights[indices]
      new_y_val <- (1/n_orig) * sum(weights_j)
      vals_y <- c(vals_y, new_y_val)
    }
    vals_y <- cumsum(vals_y)
    
    if (type=="step") {
      method <- "constant"
    } else if (type=="linear (top)") {
      method <- "linear"
    } else if (type=="linear (mid)") {
      vals_x <- c(vals_x[1], vals_x[1:(length(vals_x)-1)]+diff(vals_x)/2,
                  vals_x[length(vals_x)])
      vals_y <- c(0, vals_y[1:(length(vals_y)-1)], 1)
      method <- "linear"
    }
    
    if (which=="ecdf") {
      rval <- approxfun(vals_x, vals_y, method=method, yleft=0,
                        yright=1, f=0, ties="ordered")
    } else if (which=="inverse") {
      rval <- approxfun(vals_y, vals_x, method=method, yleft=min(vals_x),
                        yright=max(vals_x), f=1, ties="ordered")
    }
    return(rval)
    
  } else if (type=="true") {
    
    if (L$distr_A=="Unif(0,1)") {
      return(function(x) {x})
    } else if (L$distr_A=="N(0.5,0.01)") {
      if (which=="ecdf") {
        return(function(x) { ptruncnorm(x, a=0, b=1, mean=0.5, sd=0.1) })
      } else {
        return(function(x) { qtruncnorm(x, a=0, b=1, mean=0.5, sd=0.1) })
      }
    } else if (L$distr_A=="N(0.5,0.04)") {
      if (which=="ecdf") {
        return(function(x) { ptruncnorm(x, a=0, b=1, mean=0.5, sd=0.2) })
      } else {
        return(function(x) { qtruncnorm(x, a=0, b=1, mean=0.5, sd=0.2) })
      }
    }
    
  }
  
}



#' Construct conditional survival estimator S_n
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param type One of c("true", "Cox PH", "KM", "Random Forest")
#' @param csf Logical; if TRUE, estimate the conditional survival
#'     function of the censoring distribution instead
#' @return Conditional density estimator function
construct_S_n <- function(dat, vals, type, csf=FALSE) {
  
  if (csf) { dat$delta_star <- 1 - dat$delta_star }
  
  if (type %in% c("Cox PH", "Random Forest")) {
    fml <- "Surv(y_star,delta_star)~a"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind(dat$y_star, dat$delta_star, dat$a, dat$w, dat$weights)
  }
  
  if (type=="Cox PH") {
    
    weights_m1 <- dat$weights * (length(dat$weights)/sum(dat$weights))
    
    # Fit Cox model
    model <- coxph(fml, data=df, weights=weights_m1)
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
    
  } else if (type=="Random Forest") {
    
    model <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=100,
                   splitrule="logrank", nsplit=0, case.wt=weights,
                   samptype="swor")
    
    newX <- subset(filter(vals, t==0), select=-c(t))
    pred <- predict(model, newdata=newX)
    
    fnc <- function(t, w, a) {
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-pred$time.interest))
      if (length(row)!=1) { stop("Error in construct_S_n (B)") }
      if (length(col)!=1) { stop("Error in construct_S_n (C)") }
      return(pred$survival[row,col])
    }
    
  } else if (type=="Random Forest Ted") {
    
    # !!!!! Adapt to new structure
    
    # method <- "survSL.rfsrc"
    # 
    # newX <- subset(filter(vals, t==0), select=-c(t))
    # new.times <- unique(vals$t)
    # 
    # srv <- survSuperLearner(
    #   time = dat$y_star,
    #   event = dat$delta_star,
    #   X = subset(dat, select=c(w1,w2,a)),
    #   newX = newX,
    #   new.times = new.times,
    #   event.SL.library = c(method),
    #   cens.SL.library = c(method),
    #   obsWeights = dat$weights,
    #   control = list(
    #     initWeightAlg = method,
    #     max.SL.iter = 10
    #   )
    # )
    # 
    # fnc <- function(t, w1, w2, a) {
    #   r1 <- which(abs(w1-newX$w1)<1e-8)
    #   r2 <- which(abs(w2-newX$w2)<1e-8)
    #   r3 <- which(abs(a-newX$a)<1e-8)
    #   row <- intersect(r1,intersect(r2,r3))
    #   col <- which(t==new.times)
    #   return(srv$event.SL.predict[row,col])
    # }
    
  } else if (type=="true") {
    
    surv_true <- L$surv_true
    alpha_3 <- L$alpha_3
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    lmbd2 <- L$sc_params$lmbd2
    v2 <- L$sc_params$v2
    if (csf) {
      fnc <- function(t, w, a) {
        if (L$surv_true=="Cox PH") {
          lin <- C$alpha_1*w[1] + C$alpha_2*w[2] - 1
        } else if (L$surv_true=="complex") {
          lin <- C$alpha_1*pmax(0,2-8*abs(w[1]-0.5)) - 0.35
        }
        return(exp(-1*lmbd2*(t^v2)*exp(lin)))
      }
    } else {
      fnc <- function(t, w, a) {
        if (L$surv_true=="Cox PH") {
          lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*a - 1.7
        } else if (L$surv_true=="complex") {
          lin <- C$alpha_1*pmax(0,2-8*abs(w[1]-0.5)) + 1.2*alpha_3*w[2]*a - 1
        }
        return(exp(-1*lmbd*(t^v)*exp(lin)))
      }
    }
    
  }  
  
  # round_args <- c(-log10(C$appx$t_e), -log10(C$appx$w1b), 0, -log10(C$appx$a))
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2,1), vals=vals))
  
}



#' Construct g-computation estimator function of theta_0
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @return G-computation estimator of theta_0
construct_gcomp_n <- function(dat_orig, vals, S_n) {
  
  fnc <- function(a) {
    1 - mean(S_n(C$t_e,dat_orig$w,a))
  }
  
  # round_args <- c(-log10(C$appx$a)))
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
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
      
      width <- 0.1
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
  
  return(Vectorize(function(x){
    (4*deriv_theta_n(x)*f_a_n(x)*gamma_n(x))^(1/3)
  }))
  
}



#' Construct estimator of nuisance influence function omega_n
#' 
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param Sc_n Conditional censoring survival function estimator returned by
#'     construct_S_n
#' @param type Defaults to "estimated". Override with "true" for debugging. Note
#'     that type="true" only works for surv_true="Cox PH" and assumes that S_0
#'     and Sc_0 (i.e. the true functions) are passed in.
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(vals, S_n, Sc_n, type="estimated") {
  
  if (type=="estimated") {
    
    # First, construct cumulative hazard estimator
    H_n <- function(t,w,a) { -1 * log(S_n(t,w,a)) }
    
    fnc <- function(w,a,y_star,delta_star) {
      
      k <- round(min(y_star,C$t_e))
      if (k==0) {
        integral <- 0
      } else {
        i <- c(1:k)
        # i <- c(1:m)
        # k <- min(y_star,C$t_e)
        
        integral <- 0.5 * sum(
          ( H_n(i,w,a) - H_n(i-1,w,a) ) * (
            ( S_n(i,w,a) * Sc_n(i,w,a) )^-1 +
              ( S_n(i-1,w,a) * Sc_n(i-1,w,a))^-1
          )
        )
        
        # integral_righthand <- sum(
        #   ( H_n(i,w,a) - H_n(i-1,w,a) ) *
        #     ( (S_n(i,w,a)) * Sc_n(i,w,a) )^-1
        # )
        
        # integral_diffgrid <- 0.5 * sum(
        #   ( H_n((i*k)/m,w,a) - H_n(((i-1)*k)/m,w,a) ) * (
        #     (S_n((i*k)/m,w,a))^-1 * (Sc_n((i*k)/m,w,a))^-1 +
        #       (S_n(((i-1)*k)/m,w,a))^-1 * (Sc_n(((i-1)*k)/m,w,a))^-1
        #   )
        # )
        
      }
      
      return(S_n(C$t_e,w,a) * (
        (delta_star * as.integer(y_star<=C$t_e)) /
          (S_n(k,w,a) * Sc_n(k,w,a)) -
          integral
      ))
      
    }
    
  } else if (type=="true") {
    
    fnc <- function(w,a,y_star,delta_star) {
      
      # Shorten parameter variable names
      lmbd <- L$sc_params$lmbd
      v <- L$sc_params$v
      lmbd2 <- L$sc_params$lmbd2
      v2 <- L$sc_params$v2
      alpha_3 <- L$alpha_3
      
      # Construct linear predictors
      lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*a - 1.7
      lin2 <- C$alpha_1*w[1] + C$alpha_2*w[2] - 1
      
      # Compute omega_0
      piece_1 <- exp(-1*lmbd*(C$t_e^v)*exp(lin))
      piece_2 <- (delta_star*as.integer(y_star<=C$t_e)) /
        exp(-1*(lmbd*y_star^v*exp(lin)+lmbd2*y_star^v2*exp(lin2)))
      integral <- integrate(
        function(t) {
          t^(v-1) * exp(lin+lmbd*t^v*exp(lin)+lmbd2*t^v2*exp(lin2))
        },
        lower = 0,
        upper = min(C$t_e,y_star)
      )$value
      piece_3 <- lmbd*v*integral
      
      return(piece_1*(piece_2-piece_3))
      
    }
    
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(2,1,1,1), vals=vals))
  
}



#' Construct gamma_n nuisance estimator function
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param type Type of regression; one of c("linear", "cubic")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n() among the observations for which delta==1
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, dat, vals, type, omega_n, f_aIw_n,
                              f_a_n, f_a_delta1_n) {
  
  # Estimate marginal delta probability
  delta_prob <- mean(dat_orig$delta)
  
  # Construct I_star variable and pseudo-outcomes
  I_star <- dat$delta_star*as.integer(dat$y_star<=C$t_e)
  po <- ((dat$weights*omega_n(dat$w,dat$a,dat$y_star,dat$delta_star)) /
               f_aIw_n(dat$a,dat$w))^2
  
  # Create dataframe for regression
  df <- data.frame(
    a=dat$a,
    po=po,
    I_star=I_star
  )
  
  # Remove outliers to prevent errors (revisit this)
  df %<>% filter(is.finite(po))
  # if (type %in% c("cubic", "kernel")) {
  #   cutoff <- as.numeric(quantile(df$po, 0.75) + 10000*iqr(df$po))
  #   df %<>% filter(po<cutoff)
  # }
  
  # Run regression
  if (type=="cubic") {
    
    model <- lm(po~a+I(a^2)+I(a^3), data=df)
    coeff <- as.numeric(model$coefficients)
    
    reg <- function(x) {
      coeff[1] + coeff[2]*x + coeff[3]*(x^2) + coeff[4]*(x^3)
    }

  } else if (type=="kernel") {
    
    ks <- ksmooth(
      x = df$a,
      y = df$po,
      kernel = "normal",
      bandwidth = 0.2,
      x.points = vals$a
    )
    
    reg <- function(x) {
      index <- which.min(abs(x-ks$x))
      return(ks$y[index])
    }
    
  } else if (type=="kernel2") {
    
    df_0 <- filter(df, I_star==0)
    df_1 <- filter(df, I_star==1)
    
    ks_0 <- ksmooth(
      x = df_0$a,
      y = df_0$po,
      kernel = "normal",
      bandwidth = 0.2,
      x.points = vals$a
    )
    
    ks_1 <- ksmooth(
      x = df_1$a,
      y = df_1$po,
      kernel = "normal",
      bandwidth = 0.2,
      x.points = vals$a
    )
    
    reg_0 <- function(x) {
      index <- which.min(abs(x-ks_0$x))
      return(ks_0$y[index])
    }
    
    reg_1 <- function(x) {
      index <- which.min(abs(x-ks_1$x))
      return(ks_1$y[index])
    }
    
    # Run logistic regression
    model <- glm(I_star~a, data=df, family="binomial")
    coeff <- as.numeric(model$coefficients)
    
    prob_1 <- function(x) { expit(coeff[1]+coeff[2]*x) }
    prob_0 <- function(x) { 1 - prob_1(x) }
    
    reg <- function(x) {
      reg_0(x)*prob_0(x) + reg_1(x)*prob_1(x)
    }
    
  }
  
  fnc <- function(x) {
    delta_prob*(f_a_delta1_n(x)/f_a_n(x)) * reg(x)
  }
  
  # round_args <- -log10(C$appx$a))
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct estimator of conditional density f_{A|W}
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param type One of c("parametric", "binning")
#' @param k Number of bins for the binning estimator (if k=0, then the number of
#'     bins will be selected via cross-validation); ignored for the parametric
#'     estimator
#' @param delta1 Compute the density conditional on Delta=1
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of A is [0,1]
construct_f_aIw_n <- function(dat, vals, type, k=0, delta1=FALSE) {
  
  if (delta1) { dat$weights <- rep(1, length(dat$weights)) }
  
  if (type=="true") {
    
    # Note: this is not accurate if edge!="none"
    
    if (L$distr_A=="Unif(0,1)") {
      fnc <- function(a, w) { 1 }
    } else if (L$distr_A=="N(0.5,0.01)") {
      fnc <- function(a, w) {
        dnorm(a, mean=0.5, sd=0.1)
      }
    } else if (L$distr_A=="N(0.5,0.04)") {
      fnc <- function(a, w) {
        dnorm(a, mean=0.5, sd=0.2)
      }
    } else if (L$distr_A=="N(0.4+0.2w1+0.1w2,0.01)") {
      fnc <- function(a, w) {
        dnorm(a, mean=0.4+0.2*w[1]+0.1*w[2], sd=0.1)
      }
    }
    
    
  } else if (type=="parametric") {
    
    # !!!!! Update this to estimate Normals instead of Beta
    
    # # Set up weighted likelihood
    # wlik <- function(prm) {
    #   
    #   sum_loglik <- sum(sapply(c(1:length(dat$a)), function(i) {
    #     shape1 <- prm[1] + prm[2]*dat$w[i,1]
    #     shape2 <- prm[3] + prm[4]*dat$w[i,2]
    #     loglik <- dbeta(ifelse(dat$a[i]==0,1e-4,dat$a[i]),
    #                     shape1=shape1, shape2=shape2, log=TRUE)
    #     return(loglik*dat$weights[i])
    #   }))
    #   
    #   return(-1*sum_loglik)
    #   
    # }
    # 
    # # Run optimizer
    # opt <- optim(prm=c(a1=1, a2=0.1, a3=1, a4=0.1), fn=wlik)
    # if (opt$convergence!=0) {
    #   warning("construct_f_aIw_n: optim() did not converge")
    # }
    # 
    # fnc <- function(a, w1, w2){
    #   shape1 <- opt$prm[1] + opt$prm[2]*w1
    #   shape2 <- opt$prm[3] + opt$prm[4]*w2
    #   return(dbeta(a, shape1=shape1, shape2=shape2))
    # }
    
  } else if (type=="binning") {
    
    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # prm[1] through prm[k-1] are the hazard components for the bins 1 to k-1
    # prm[k] and prm[k+1] are the coefficients for W1 and W2
    create_dens <- function(k, dat) {
      
      # Cut points
      alphas <- seq(0, 1, length.out=k+1)
      
      # Density for a single observation
      dens_s <- construct_superfunc(function(a, w, prm) {
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        if (x==0) {
          print("hey")
          print("prm")
          print(prm)
          print("c(k:(k+length(w)-1))")
          print(c(k:(k+length(w)-1)))
          print("prm[c(k:(k+length(w)-1))]")
          print(prm[c(k:(k+length(w)-1))])
          print("w")
          print(w)
          print("as.numeric(prm[c(k:(k+length(w)-1))] %*% w)")
          print(as.numeric(prm[c(k:(k+length(w)-1))] %*% w))
          x<<-1
        }
        hz <- expit(
          prm[c(1:(ifelse(bin==k,k-1,bin)))] +
            as.numeric(prm[c(k:(k+length(w)-1))] %*% w)
            # (par[k]*w[1] + par[k+1]*w[2])
        )
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        return(k*p1*p2)
      }, vec=c(1,2,0))
      
      # Set up weighted likelihood
      wlik <- function(prm) {
        # if (x==0) {
        #   print("dat$a")
        #   print(dat$a)
        #   print("dat$w")
        #   print(dat$w)
        #   x<<-1
        # }
        -1 * sum(dat$weights * log(pmax(dens_s(a=dat$a, w=dat$w, prm),1e-8)))
      }
      
      # Run optimizer
      opt <- solnp(
        pars = rep(0.001,k+length(dat$w)-1),
        fun = wlik
      )
      if (opt$convergence!=0) {
        warning("construct_f_aIw_n: solnp() did not converge")
      }
      
      fnc <- function(a, w) {
        
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        prm <- opt$pars
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(prm[j] + prm[c(k:(k+length(w)-1))] %*% w)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        
        return(k*p1*p2)
        
      }
      
    }
    
    # Select k via cross-validation
    if (k==0) {
      
      # Prep
      n_folds <- 5
      folds <- sample(cut(c(1:length(dat$a)), breaks=n_folds, labels=FALSE))
      ks <- c(5,10,15,20,25)
      best <- list(k=999, max_log_lik=999)
      
      # Cross-validation
      for (k in ks) {
        
        sum_log_lik <- 0
        for (i in c(1:n_folds)) {
          dat_train <- list(
            a = dat$a[-which(folds==i)],
            w = dat$w[-which(folds==i),]
          )
          dat_test <- list(
            a = dat$a[which(folds==i)],
            w = dat$w[which(folds==i),]
          )
          dens <- create_dens(k, dat_train)
          sum_log_lik <- sum_log_lik +
            sum(log(dens(dat_test$a, dat_test$w)))
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
  
  # round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' Construct estimator of marginal density f_A
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return Marginal density estimator function
construct_f_a_n <- function(dat_orig, vals, f_aIw_n) {
  
  fnc <- function(a) {
    mean(f_aIw_n(a,dat_orig$w))
  }
  
  # round_args <- -log10(C$appx$a))
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct density ratio estimator g_n
#' 
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
construct_g_n <- function(vals, f_aIw_n, f_a_n) {
  
  fnc <- function(a,w) {
    f_aIw_n(a,w) / f_a_n(a)
  }
  
  # round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' Construct nuisance estimator eta_n
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @return Estimator function of nuisance eta_0
construct_eta_n <- function(dat, vals, S_n) {
  
  n_orig <- sum(dat$weights)
  
  fnc <- function(x,w) {
    (1/n_orig) * sum(
      dat$weights * as.integer(dat$a<=x) *
        (1-S_n(C$t_e,w,dat$a))
    )
  }
  
  # round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' lambda estimator
#' 
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @return Value of lambda
lambda <- function(n_orig, dat, k, G) {
  
  n_orig <- sum(dat$weights)
  
  lambda <- (1/n_orig) * sum(
    dat$weights * (G(dat$a))^k
  )
  return(lambda)
  
}



#' Construct Gamma_n primitive one-step estimator
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param S_n A conditional survival function returned by construct_S_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @param type One of c("one-step", "plug-in")
#' @return Gamma_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Gamma_n <- function(dat, vals, omega_n, S_n, g_n, type="one-step") {
  
  n_orig <- sum(dat$weights)
  n <- length(dat$a)
  i_long <- rep(c(1:n), each=n)
  j_long <- rep(c(1:n), times=n)
  a_i_long <- dat$a[i_long]
  w_i_long <- dat$w[i_long,]
  w_j_long <- dat$w[j_long,]
  delta_star_i_long <- dat$delta_star[i_long]
  delta_star_j_long <- dat$delta_star[j_long]
  weights_i_long <- dat$weights[i_long]
  weights_j_long <- dat$weights[j_long]
  
  if (type=="one-step") {
    
    subpiece_1a <- dat$weights * ( 1 + (
      omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
        g_n(dat$a,dat$w)
    ) )
    subpiece_2a <- (weights_i_long*weights_j_long) *
      S_n(C$t_e,w_j_long,a_i_long)
    
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
            (1 - S_n(C$t_e,w_j_long,a_i_long))
        )
      )
    }
    
  }
  
  # round_args <- -log10(C$appx$a)
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_rho_n <- function(dat, Phi_n) {
  
  n_orig <- sum(dat$weights)
  
  return(memoise(Vectorize(function(a) {
    (1/n_orig) * sum(
      dat$weights * (Phi_n(dat$a)^3) * (as.integer(a<=dat$a) - Phi_n(dat$a))
    )
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_xi_n <- function(Phi_n, lambda_2, lambda_3) {
  
  return(memoise(Vectorize(function(a_i,a_j) {
    (2*as.integer(a_i<=a_j) - 3*Phi_n(a_j))*Phi_n(a_j)*lambda_2 +
    (2*Phi_n(a_j) - as.integer(a_i<=a_j))*lambda_3
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_1 <- function(dat, Gamma_n, Phi_n, xi_n, rho_n,
                                lambda_2, lambda_3) {
  
  n_orig <- sum(dat$weights)
  weights_j <- dat$weights
  a_j <- dat$a
  
  fnc <- function(a) {
    
    piece_1 <- (1/n_orig) * sum(
      weights_j * (xi_n(a,a_j)) * Gamma_n(a_j)
    )
    
    piece_2 <- (lambda_2*(Phi_n(a)^2) - lambda_3*Phi_n(a)) * Gamma_n(a)
    
    return(piece_1+piece_2)
    
  }
  
  return(construct_superfunc(fnc))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Gamma <- function(omega_n, g_n, gcomp_n, eta_n,
                                    Gamma_n) {
  
  fnc <- function(x,w,y_star,delta_star,a) {
    as.integer(a<=x)*(
      (omega_n(w,a,y_star,delta_star)/g_n(a,w)) +
        gcomp_n(a)
    ) +
      eta_n(x,w) - 
      2*Gamma_n(x)
  }
  
  return(construct_superfunc(fnc, vec=c(1,2,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_2 <- function(dat, Phi_n, infl_fn_Gamma, lambda_2, lambda_3) {
  
  n_orig <- sum(dat$weights)
  weights_j <- dat$weights
  a_j <- dat$a
  
  fnc <- function(w,y_star,delta_star,a) {
    (1/n_orig) * sum(weights_j * (
      ( lambda_2*(Phi_n(a_j)^2) - lambda_3*Phi_n(a_j) ) *
        infl_fn_Gamma(a_j,w,y_star,delta_star,a)
    ))
  }
  
  return(construct_superfunc(fnc, vec=c(2,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
beta_n_var_hat <- function(dat, infl_fn_1, infl_fn_2) {
  
  n_orig <- sum(dat$weights)
  
  b_sum <- 0
  for (i in c(1:length(dat$a))) {
    b_sum <- b_sum + (dat$weights[i] * (
      infl_fn_1(dat$a[i]) +
      # infl_fn_1(dat$w1[i], dat$w2[i], dat$y_star[i],
      #           dat$delta_star[i], dat$a[i]) +
        infl_fn_2(dat$w[i,], dat$y_star[i],
                  dat$delta_star[i], dat$a[i])
    ))^2
  }
  
  return( (1/n_orig)*sum(b_sum) )
  
}



#' Construct lists of values to pre-compute functions on
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param appx Approximation spec used for grids (stored in C$appx)
#' @return A list of dataframes
create_val_list <- function(dat, appx) {
  
  return(list(
    A = list(a=dat$a),
    # AW = data.frame(a=dat$a, w1=dat$w1, w2=dat$w2),
    A_grid = list(a=seq(0,1,appx$a)),
    # W_grid = expand.grid(w1=seq(0,1,appx$w1), w2=c(0,1)),
    # AW_grid = expand.grid(a=seq(0,1,appx$a), w1=seq(0,1,appx$w1),
                          # w2=c(0,1)),
    # S_n = expand.grid(t=seq(0,C$t_e,appx$t_e), w1=seq(0,1,appx$w1b),
                      # w2=c(0,1), a=seq(0,1,appx$a)),
    omega = list(
      w = rbind(dat$w,dat$w),
      a = c(dat$a, rep(0, length(dat$a))),
      y_star = rep(dat$y_star,2),
      delta_star = rep(dat$delta_star,2)
    )
  ))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_Gamma_cf_k <- function(dat_train, dat_test, vals, omega_n, g_n,
                                 gcomp_n, eta_n) {
  
  # !!!!! Needs to be updated
  
  n_test <- length(dat_test$a)
  dat_test$weights <- wts(dat_test) # !!!!! Weights need to be re-stabilized here
  d1 <- ss(dat_test, which(dat_test$delta==1))
  weights_1 <- d1$weights
  
  n_train <- length(dat_train$a)
  dat_train$weights <- wts(dat_train) # !!!!! Weights need to be re-stabilized here
  d2 <- ss(dat_train, which(dat_train$delta==1))
  weights_2 <- d2$weights
  
  fnc <- function(x) {
    
    piece_1 <- (1/n_test) * sum(weights_1 * (
      as.integer(d1$a<=x) *
        ((omega_n(d1$w,d1$a,d1$y_star,d1$delta_star) / g_n(d1$a,d1$w)) +
           gcomp_n(d1$a)
        ) +
        eta_n(x,d1$w)
    ))
    
    piece_2 <- (1/n_train) * sum(weights_2 * as.integer(d2$a<=x)*gcomp_n(d2$a))
    
    return(piece_1-piece_2)
    
  }
  
  # round_args <- -log10(C$appx$a)
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
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
  n_orig <- length(dat_orig$delta)
  rows <- c(1:n_orig)
  folds <- sample(cut(rows, breaks=params$cf_folds, labels=FALSE))
  
  # Loop through folds
  for (k in 1:params$cf_folds) {
    
    # Split data
    dat_orig_train <- ss(dat_orig, -which(folds==k))
    dat_orig_test <- ss(dat_orig, which(folds==k))
    dat_train <- ss(dat_orig_train, which(dat_orig_train$delta==1))
    dat_test <- ss(dat_orig_test, which(dat_orig_test$delta==1))
        
    # Construct component functions
    Phi_n <- construct_Phi_n(dat_train, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat_train, which="inverse",
                                 type=params$ecdf_type)
    S_n <- construct_S_n(dat_train, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_train, vlist$S_n, type=params$S_n_type, csf=TRUE)
    gcomp_n <- construct_gcomp_n(dat_orig_train, vlist$A_grid, S_n)
    eta_n <- construct_eta_n(dat_train, vlist$AW_grid, S_n)
    f_aIw_n <- construct_f_aIw_n(dat_train, vlist$AW_grid, type=params$g_n_type,
                                 k=15)
    f_a_n <- construct_f_a_n(dat_orig_train, vlist$A_grid, f_aIw_n)
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
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param type One of c("logistic","")
#' @return Propensity score estimator of pi_0
construct_pi_n <- function(dat, vals, type) {
  
  # Construct indicator I{A=0}
  ind_A0 <- as.integer(dat$a==0)
  
  if (type=="true") {
    
    if (L$edge=="expit") {
      fnc <- function(w) { expit(w[1]+w[2]-3.3) }
    } else if (L$edge=="expit2") {
      fnc <- function(w) { expit(w[1]+w[2]-1) }
    } else if (L$edge=="complex") {
      fnc <- function(w) { 0.84*w[2]*pmax(0,1-4*abs(w[1]-0.5)) }
    } else if (L$edge=="none") {
      fnc <- function(w) { 0 }
    }
    
  } else if (type=="logistic") {
    
    fml <- "ind_A0~1"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind(ind_A0, dat$w, dat$weights)
    suppressWarnings({
      model <- glm(
        fml,
        data = df,
        family = "binomial",
        weights = weights
      )
    })
    coeffs <- model$coefficients
    
    fnc <- function(w) {
      # expit( coeffs[[1]] + coeffs[[2]]*w1 + coeffs[[3]]*w2 )
      expit( as.numeric(coeffs) %*% c(1,w) )
    }
    
  } else if (type=="SL") {
    
    sl <- SuperLearner(
      Y = ind_A0,
      X = dat$w,
      newX = vals,
      family = binomial(),
      SL.library = "SL.earth", # SL.glm SL.gbm SL.ranger SL.earth
      # SL.library = c("SL.earth", "SL.gam", "SL.ranger"), # SL.glm SL.gbm SL.ranger SL.earth
      obsWeights = dat$weights,
      control = list(saveFitLibrary=FALSE)
    )
    assign("sl", sl, envir=.GlobalEnv) # ?????
    
    fnc <- function(w) {
      
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      index <- Reduce(intersect, r)
      return(sl$SL.predict[index])
    }
    
  }
  
  # round_args <- c(-log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(2), vals=vals))
  
}



#' COmpute one-step estimator of counterfactual survival at A=0
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @return Value of one-step estiamtor
theta_os_n <- function(dat, pi_n, S_n, omega_n) {
  
  n_orig <- sum(dat$weights)
  
  return(
    1 - (1/n_orig) * sum(dat$weights * (
      S_n(C$t_e,dat$w,a=0) - (
        (as.integer(dat$a==0)/pi_n(dat$w)) *
          omega_n(dat$w,a=0,dat$y_star,dat$delta_star)
      )
    ))
  )
  
}



#' Compute asymptotic variance of one-step estimator theta_os_n
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param theta_os_n_est Estimate returned by one-step estimator theta_os_n()
#' @return Asymptotic variance estimate
sigma2_os_n <- function(dat, pi_n, S_n, omega_n, theta_os_n_est) {
  
  n_orig <- sum(dat$weights)
  
  return(
    (1/n_orig) * sum((dat$weights * (
      S_n(C$t_e,dat$w,a=0) - (
        (as.integer(dat$a==0)/pi_n(dat$w)) *
          omega_n(dat$w,a=0,dat$y_star,dat$delta_star)
      ) -
        (1-theta_os_n_est)
    ))^2)
  )
  
}
