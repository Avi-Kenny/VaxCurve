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
expit <- function(x) { 1 / (1+exp(-x)) }



#' Derivative of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_expit <- function(x) { exp(x) / ((1+exp(x))^2) }



#' Alias for indicator (as.integer) function
#' 
#' @param x Logical input
#' @return BInary (integer) output
In <- as.integer



#' Derivative of logit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_logit <- function(x) { 1 / (x-x^2) }



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
#' @param rnd If provided, a list of length equal to the number of arguments in
#'     fnc; each value should be an integer or vector of integers that will be
#'     passed to round()
#' @return A function with the above modifications. This can be called in the
#'     same way as the original function or in vectorized format. For example,
#'     the function f(A,B,C,D) described above can be called via
#'     f(2,4,c(1,2,3),9) (same as original) or via
#'     f(c(2,3),c(4,5),data.frame(x=c(1,5),y=c(2,6),z=c(3,7)),9) (vectorized).
#'     The number of rows of each data frame should be the same as the number of
#'     rows of each vector. The data frame column names are ignored.
#' @notes
#'   - The order of inputs matters and argument names are ignored.
#'   - This function only works on numeric data and will coerce all data to be
#'     numeric.
construct_superfunc <- function(fnc, aux=NA, vec=TRUE, vals=NA, rnd=NA) {
  
  htab <- new.env()
  ..new_fnc <- function() {
    
    ..e <- parent.env(environment())
    ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
    ..mc1 <- ..mc[[1]]
    if (max(..e$vec)!=1) {
      index <- min(which(..e$vec==2))
      ..mc2 <- ..mc[[..e$arg_names[index]]]
    }
    if (max(..e$vec)==1 && length(..mc1)==1) {
      vec_mode <- F
    } else if (max(..e$vec)==1 && length(..mc1)>1) {
      vec_mode <- T
    } else if (max(..e$vec)!=1 && class(..mc2)!="data.frame") {
      vec_mode <- F
    } else if (max(..e$vec)!=1 && class(..mc2)=="data.frame") {
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
      for (j in 1:length(..e$arg_names)) {
        if (..e$vec[j]==2) {
          ..mc[[..e$arg_names[j]]] <- as.list(
            as.data.frame(t(..mc[[..e$arg_names[j]]]), row.names=NA)
          )
        }
      }
      FUN <- function(...) {
        keylist <- lapply(..e$arg_names, function(arg_name) {
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
      res <- do.call(mapply, c(
        FUN = FUN,
        ..mc[..e$arg_names[dovec]],
        MoreArgs = list(..mc[..e$arg_names[!dovec]]),
        USE.NAMES = F
      ))
      
    } else {
      
      # # !!!!!
      # if (is.na(rnd[1])) {
      #   l <- list(...)
      # } else {
      #   l <- lapply(..e$arg_names, function(arg_name) {
      #     as.numeric(list(...)[[arg_name]])
      #   })
      # }
      
      # Non-vectorized operation
      keylist <- lapply(..e$arg_names, function(arg_name) {
        as.numeric(..mc[[arg_name]])
      })
      key <- paste(keylist, collapse=";")
      val <- ..e$htab[[key]]
      if (is.null(val)) {
        val <- do.call(..e$fnc, ..mc)
        ..e$htab[[key]] <- val
      }
      res <- val
      
    }
    
    # Return value
    return(res)
    
  }
  
  # Transform `vec` (if needed) and validate
  if (identical(vec[1],T)) { vec <- rep(1, length(names(formals(fnc)))) }
  if (!is.numeric(vec) || any(!(vec %in% c(0,1,2)))) {
    stop("`vec` must be a vector of 0s, 1s, and 2s")
  }
  
  # Set formals and set up environment
  formals(..new_fnc) <- formals(fnc)
  f_env <- new.env(parent=environment(fnc))
  f_env$arg_names <- names(formals(fnc))
  f_env$vec <- vec
  f_env$rnd <- rnd
  f_env$htab <- htab
  f_env$aux <- aux
  f_env$fnc <- fnc
  environment(..new_fnc) <- f_env
  
  # Run function on vals list
  if (is.list(vals)) { do.call(..new_fnc, vals) }
  
  return(..new_fnc)
  
}



#' Probability of sampling
#' 
#' @param sampling One of c("iid", "two-phase (6%)", "two-phase (72%)",
#'     "two-phase (70% random)", "two-phase (6% random)")
#' @param delta_star Component of dataset returned by generate_data()
#' @param y_star Component of dataset returned by generate_data()
#' @param w Component of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
Pi <- function(sampling, delta_star, y_star, w) {
  
  if (sampling=="iid") {
    probs <- rep(1, length(y_star))
  } else if (sampling=="two-phase (70% random)") {
    probs <- rep(0.7, length(y_star))
  } else if (sampling=="two-phase (50% random)") {
    probs <- rep(0.5, length(y_star))
  } else if (sampling=="two-phase (25% random)") {
    probs <- rep(0.25, length(y_star))
  } else if (sampling=="two-phase (6% random)") {
    probs <- rep(0.06, length(y_star))
  } else if (sampling=="cycle") {
    probs <- rep(c(0.8,0.6,0.4,0.2), length.out=length(y_star))
  } else {
    ev <- In(delta_star==1 & y_star<=C$t_0)
    if (sampling=="two-phase (6%)") {
      probs <- ev + (1-ev)*expit(w$w1+w$w2-3.85)
    } else if (sampling=="two-phase (72%)") {
      probs <- ev + (1-ev)*expit(w$w1+w$w2-0.1)
    } else if (sampling=="two-phase (50%)") {
      probs <- ev + (1-ev)*expit(w$w1+w$w2-1)
    } else if (sampling=="two-phase (25%)") {
      probs <- ev + (1-ev)*expit(w$w1+w$w2-2.2)
    } else if (sampling=="w1") {
      probs <- (0.2 + 0.6*w$w1)
    } else if (sampling=="w2") {
      probs <- (0.2 + 0.6*w$w2)
    }
    
  }
  
  return(probs)
  
}



#' Return IP weights
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param scale One of c("none", "stabilized")
#' @param type One of c("true", "estimated")
#' @param return_strata Whether the discrete two-phase sampling strata
#'     membership variable should be returned
#' @return A sum-to-one vector of weights
#' @notes
#'   - Only used for simulation; for the real analysis, the weights are
#'     calculated separately
wts <- function(dat_orig, scale="stabilized", type="true", return_strata=F) {
  
  sampling <- attr(dat_orig,"sampling")
  Pi_0 <- Pi(sampling, dat_orig$delta_star, dat_orig$y_star, dat_orig$w)
  strata1 <- In(factor(Pi_0))
  
  if (type=="true") {
    
    weights <- dat_orig$delta / Pi_0
    
  } else if (type=="estimated") {
    
    Pi_n_vals <- c()
    for (i in c(1:max(strata1))) {
      Pi_n_vals[i] <- sum(In(strata1==i)*dat_orig$delta) / sum(In(strata1==i))
      if (Pi_n_vals[i]==0) {
        # Hack to avoid NA values in small sample sizes
        warning(paste0("stratum ", i, " had no one sampled."))
        Pi_n_vals[i] <- 1
      }
    }
    weights <- dat_orig$delta / Pi_n_vals[strata1]
    
  }
  
  if (scale=="stabilized") {
    s <- sum(weights) / length(dat_orig$delta)
    weights <- weights / s
  }
  
  if (!return_strata) {
    return(weights)
  } else {
    return(list(weights=weights, strata=strata1))
  }
  
}



#' Subset dat_orig according to indices
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param indices Indices to filter dataset by
#' @return Filtered subsample of dataset
ss <- function(dat_orig, indices) {
  
  i <- indices
  
  return(list(
    id = dat_orig$id[i],
    w = dat_orig$w[i,, drop=F],
    a = dat_orig$a[i],
    delta = dat_orig$delta[i],
    y_star = dat_orig$y_star[i],
    delta_star = dat_orig$delta_star[i],
    weights = dat_orig$weights[i],
    strata = dat_orig$strata[i]
  ))
  
}



#' Round values in dat_orig
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @return Dataset with values rounded
round_dat <- function(dat_orig) {
  
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$y_star <- round(dat_orig$y_star, -log10(C$appx$t_0))
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
  
  return(dat_orig)
  
}



#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param which One of c("ecdf", "inverse")
#' @param type One of c("true", "step", "linear (top)", "linear (mid)")
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
construct_Phi_n <- function (dat, which="ecdf", type="step") {
  
  if (type!="true") {
    
    n_orig <- sum(dat$weights)
    df <- data.frame(a=dat$a, weights=dat$weights)
    df %<>% arrange(a)
    vals_x <- unique(df$a)
    vals_y <- c()
    
    for (j in 1:length(vals_x)) {
      indices <- which(df$a==vals_x[j])
      weights_j <- df$weights[indices]
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
    
    if (L$distr_S=="Unif(0,1)") {
      return(function(a) {a})
    } else if (L$distr_S=="N(0.5,0.01)") {
      if (which=="ecdf") {
        return(function(a) { ptruncnorm(a, a=0, b=1, mean=0.5, sd=0.1) })
      } else {
        return(function(a) { qtruncnorm(a, a=0, b=1, mean=0.5, sd=0.1) })
      }
    } else if (L$distr_S=="N(0.5,0.04)") {
      if (which=="ecdf") {
        return(function(a) { ptruncnorm(a, a=0, b=1, mean=0.5, sd=0.2) })
      } else {
        return(function(a) { qtruncnorm(a, a=0, b=1, mean=0.5, sd=0.2) })
      }
    } else if (L$distr_S=="N(0.3+0.4x2,0.04)") {
      stop("ecdf not programmed for N(0.3+0.4x2,0.04)")
    }
    
  }
  
}



#' Construct conditional survival estimator Q_n
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc(); REQUIRED FOR SUPERLEARNER
#' @param type One of c("true", "Cox PH", "Random Forest", "Super Learner")
#' @param return_model Logical; if TRUE, return the model object instead of the
#'     function
#' @param print_coeffs Logical; if TRUE, print the algorithm coefficients for
#'     the conditional survival and censoring functions (applicable only if
#'     type=="Super Learner")
#' @return Conditional density estimator function
construct_Q_n <- function(dat, vals, type, return_model=F, print_coeffs=F) {
  
  if (type %in% c("Cox PH", "Super Learner")) {
    if (type=="Cox PH") {
      methods <- c("survSL.coxph")
    } else if (type=="Super Learner") {
      # Excluding "survSL.rfsrc" for now. survSL.pchSL gives errors.
      methods <- c("survSL.coxph", "survSL.expreg", "survSL.km",
                   "survSL.loglogreg", "survSL.pchreg", "survSL.weibreg")
      # methods <- c("survSL.km", "survSL.pchreg", "survSL.rfsrc") # !!!!!
    }
    
    newX <- cbind(vals$w, a=vals$a)[which(vals$t==0),]
    new.times <- unique(vals$t)
    srv <- survSuperLearner(
      time = dat$y_star,
      event = dat$delta_star,
      X = cbind(dat$w, a=dat$a),
      newX = newX,
      new.times = new.times,
      event.SL.library = methods,
      cens.SL.library = methods,
      obsWeights = dat$weights,
      control = list(
        initWeightAlg = methods[1],
        max.SL.iter = 10
      )
    )
    
    if (print_coeffs && type=="Super Learner") {
      cat("\n------------------------------\n")
      cat("SuperLearner algorithm weights\n")
      cat("------------------------------\n\n")
      cat("event.coef\n")
      cat("----------\n")
      print(sort(srv$event.coef, decr=T))
      cat("\ncens.coef\n")
      cat("---------\n")
      print(sort(srv$cens.coef, decr=T))
      cat("\n------------------------------\n")
    }
    
    srv_pred <- srv$event.SL.predict
    cens_pred <- srv$cens.SL.predict
    rm(srv)
    
    # !!!!! Later consolidate these via a wrapper/constructor function
    fnc_srv <- function(t, w, a) {
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      if (class(newX[["a"]][1])=="factor") {
        r[[length(w)+1]] <- which(a==newX[["a"]])
      } else {
        r[[length(w)+1]] <- which(abs(a-newX[["a"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",w=(",
                    paste(w,collapse=","),"),a=",a,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",w=(",
                    paste(w,collapse=","),"),a=",a,""))
      }
      return(srv_pred[row,col])
    }
    
    fnc_cens <- function(t, w, a) {
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      if (class(newX[["a"]][1])=="factor") {
        r[[length(w)+1]] <- which(a==newX[["a"]])
      } else {
        r[[length(w)+1]] <- which(abs(a-newX[["a"]])<1e-8)
      }
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-new.times))
      if (length(row)!=1) {
        stop(paste0("Error in construct_Q_n (B); ", "t=",t,",w=(",
                    paste(w,collapse=","),"),a=",a,""))
      }
      if (length(col)!=1) {
        stop(paste0("Error in construct_Q_n (C); ", "t=",t,",w=(",
                    paste(w,collapse=","),"),a=",a,""))
      }
      return(cens_pred[row,col])
    }
    
  }
  
  if (type=="true") {
    
    surv_true <- L$surv_true
    alpha_3 <- L$alpha_3
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    lmbd2 <- L$sc_params$lmbd2
    v2 <- L$sc_params$v2
    if (L$surv_true=="exp") {
      v <- 1   # Overwrite
      v2 <- 1 # Overwrite
    }
    
    fnc_srv <- function(t, w, a) {
      if (L$surv_true=="Non PH") {
        if (L$dir=="decr") {
          # lin <- C$alpha_1*log(t)*w[1]/3 + C$alpha_2*log(t)*w[2]/3 +
          #   log(t)*alpha_3*a/2 - 2
        } else {
          # lin <- C$alpha_1*log(t)*w[1]/3 + C$alpha_2*log(t)*w[2]/3 +
          #   log(t)*alpha_3*(1-a)/2 - 2
        }
      } else {
        if (L$surv_true=="Cox PH") {
          if (L$dir=="decr") {
            lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*a - 1.7
          } else {
            lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*(1-a) - 1.7
          }
        } else if (L$surv_true=="Complex") {
          if (L$dir=="decr") {
            lin <- alpha_3*expit(10*(2*a-1)) - 0.5
          } else {
            lin <- alpha_3*expit(10*(2*(1-a)-1)) - 0.5
          }
        } else if (L$surv_true=="exp") {
          lin <- 0
        }
      }
      return(exp(-1*lmbd*(t^v)*exp(lin)))
    }
    
    fnc_cens <- function(t, w, a) {
      if (L$surv_true %in% c("Cox PH", "Complex", "Non PH")) {
        lin <- C$alpha_1*w[1] + C$alpha_2*w[2] - 1
      } else if (L$surv_true=="exp") {
        lin <- 0
      }
      return(exp(-1*lmbd2*(t^v2)*exp(lin)))
    }
    
  }
  
  if (type=="constant") {
    fnc_srv <- function(t, w, a) { 0.5 }
    fnc_cens <- function(t, w, a) { 0.5 }
  }
  
  if (type=="Random Forest") {
    
    # Note: not running this through Super Learner for now
    
    # Setup
    fml <- "Surv(y_star,delta_star)~a"
    for (i in 1:length(dat$w)) { fml <- paste0(fml, "+w",i) }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star, "a"=dat$a,
                dat$w, "weights"=dat$weights)
    
    # Survival function
    model_srv <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=20,
                   splitrule="logrank", nsplit=0, case.wt=df$weights,
                   samptype="swor")
    newX <- cbind(a=vals$a, vals$w)[which(vals$t==0),]
    pred_srv <- predict(model_srv, newdata=newX)
    fnc_srv <- function(t, w, a) {
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      r[[length(w)+1]] <- which(abs(a-newX[["a"]])<1e-8)
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-pred_srv$time.interest))
      if (length(row)!=1) { stop("Error in construct_Q_n (B)") }
      if (length(col)!=1) { stop("Error in construct_Q_n (C)") }
      return(pred_srv$survival[row,col])
    }
    
    # Censoring function
    df$delta_star <- round(1-df$delta_star)
    model_cens <- rfsrc(fml, data=df, ntree=500, mtry=2, nodesize=20,
                   splitrule="logrank", nsplit=0, case.wt=df$weights,
                   samptype="swor")
    pred_cens <- predict(model_cens, newdata=newX)
    fnc_cens <- function(t, w, a) {
      r <- list()
      for (i in 1:length(w)) {
        r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
      }
      r[[length(w)+1]] <- which(abs(a-newX[["a"]])<1e-8)
      row <- Reduce(intersect, r)
      col <- which.min(abs(t-pred_cens$time.interest))
      if (length(row)!=1) { stop("Error in construct_Q_n (B)") }
      if (length(col)!=1) { stop("Error in construct_Q_n (C)") }
      return(pred_cens$survival[row,col])
    }
    
  }
  
  sfnc_srv <- construct_superfunc(fnc_srv, aux=NA, vec=c(1,2,1), vals=vals)
  sfnc_cens <- construct_superfunc(fnc_cens, aux=NA, vec=c(1,2,1), vals=vals)
  rm("vals", envir=environment(get("fnc_srv",envir=environment(sfnc_srv))))
  
  return(list(srv=sfnc_srv, cens=sfnc_cens))
  
}



#' Construct g-computation estimator function of theta_0
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @return G-computation estimator of theta_0
construct_gcomp_n <- function(dat_orig, vals=NA, Q_n) {
  
  fnc <- function(a) {
    1 - mean(Q_n(
      rep(C$t_0, nrow(dat_orig$w)),
      dat_orig$w,
      rep(a, nrow(dat_orig$w))
    ))
  }
  
  # round_args <- c(-log10(C$appx$a)))
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct derivative estimator theta'_n
#' 
#' 
#' @param r_Mn An estimator of r_M0
#' @param type One of c("gcomp", "linear", "spline")
#' @param dir Direction of monotonicity; one of c("incr", "decr")
construct_deriv_r_Mn <- function(r_Mn, type, dir="incr") {
  
  # Estimate entire function on grid
  grid <- round(seq(0,1,0.01),2)
  r_Mns <- r_Mn(grid)
  
  if (type=="linear") {
    
    grid_width <- grid[2] - grid[1]
    points_x <- c(grid[1])
    points_y <- c(r_Mns[1])
    for (i in 2:length(grid)) {
      if (r_Mns[i]-r_Mns[i-1]!=0) {
        points_x <- c(points_x, grid[i]-(grid_width/2))
        points_y <- c(points_y, mean(c(r_Mns[i],r_Mns[i-1])))
      }
    }
    points_x <- c(points_x, grid[length(grid)])
    points_y <- c(points_y, r_Mns[length(grid)])
    points_sl <- c()
    for (i in 2:length(points_x)) {
      slope <- (points_y[i]-points_y[i-1]) /
        (points_x[i]-points_x[i-1])
      points_sl <- c(points_sl, slope)
    }
    
    fnc <- function(a) {
      if (a==0) {
        index <- 1
      } else {
        index <- which(round(a,6)<=round(points_x,6))[1]-1
      }
      if (dir=="incr") {
        return(max(points_sl[index],0))
      } else {
        return(min(points_sl[index],0))
      }
    }
    
  } else if (type=="line") {
    
    r_Mn_left <- r_Mn(0) # 0.1
    r_Mn_right <- r_Mn(1) # 0.9
    fnc <- function(a) { (r_Mn_right-r_Mn_left)/1 } # 0.8
    
  } else if (type=="spline") {
    
    # Identify jump points of step function
    jump_points <- c(0)
    for (i in 2:length(grid)) {
      if (r_Mns[i]!=r_Mns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    jump_points <- c(jump_points,grid[length(grid)])
    
    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)
    
    # Fit cubic smoothing spline
    r_Mn_smoothed <- smooth.spline(x=midpoints, y=r_Mn(midpoints))
    
    # Construct derivative function
    fnc <- function(a) {
      
      width <- 0.3
      x1 <- a - width/2
      x2 <- a + width/2
      if (x1<0) { x2<-width; x1<-0 }
      if (x2>1) { x1<-1-width; x2<-1; }
      
      y1 <- predict(r_Mn_smoothed, x=x1)$y
      y2 <- predict(r_Mn_smoothed, x=x2)$y
      
      if (dir=="incr") {
        return(max((y2-y1)/width,0))
      } else {
        return(min((y2-y1)/width,0))
      }
      
    }
    
  } else if (type=="m-spline") {
    
    # Identify jump points of step function
    jump_points <- c(0)
    for (i in 2:length(grid)) {
      if (r_Mns[i]!=r_Mns[i-1]) {
        jump_points <- c(jump_points, mean(c(grid[i],grid[i-1])))
      }
    }
    jump_points <- c(jump_points,grid[length(grid)])
    
    # Identify midpoints of jump points
    midpoints <- jump_points[1:(length(jump_points)-1)]+(diff(jump_points)/2)
    
    if (length(midpoints)>=2) {
      # Fit monotone cubic smoothing spline
      r_Mn_smoothed <- splinefun(
        x = midpoints,
        y = r_Mn(midpoints),
        method = "monoH.FC"
      )
    } else {
      # Fit a straight line instead if there are <2 midpoints
      r_Mn_smoothed <- function(x) {
        slope <- r_Mn(1) - r_Mn(0)
        intercept <- r_Mn(0)
        return(intercept + slope*x)
      }
    }
    
    # Construct derivative function
    fnc <- function(a) {
      
      width <- 0.2
      x1 <- a - width/2
      x2 <- a + width/2
      
      if (x1<0) { x2<-width; x1<-0; }
      if (x2>1) { x1<-1-width; x2<-1; }
      
      y1 <- r_Mn_smoothed(x1)
      y2 <- r_Mn_smoothed(x2)
      
      if (dir=="incr") {
        return(max((y2-y1)/width,0))
      } else {
        return(min((y2-y1)/width,0))
      }
      
    }
    
  } else if (type=="gcomp") {
    
    fnc <- function(a) {
      
      # Set derivative appx x-coordinates
      width <- 0.2
      p1 <- a - width/2
      p2 <- a + width/2
      if (p1<0) { p2 <- width; p1 <- 0; }
      if (p2>1) { p1<-1-width; p2<-1; }
      c(p1,p2)
      
      if (dir=="incr") {
        return(max((r_Mn(p2)-r_Mn(p1))/width,0))
      } else {
        return(min((r_Mn(p2)-r_Mn(p1))/width,0))
      }
      
    }
    
  }
  
  return(Vectorize(fnc))

}



#' Construct tau_n Chernoff scale factor function
#' 
#' @param deriv_r_Mn A derivative estimator returned by construct_deriv_r_Mn()
#' @param gamma_n Nuisance function estimator returned by construct_gamma_n()
#' @param f_a_n Density estimator returned by construct_f_a_n()
#' @return Chernoff scale factor estimator function
construct_tau_n <- function(deriv_r_Mn, gamma_n, f_a_n,
                            pi_star_n=NA, g_n=NA, dat_orig=NA) {
  
  n_orig <- length(dat_orig$a)
  w <- dat_orig$w
  return(Vectorize(function(x) {
    abs(
      ((4*deriv_r_Mn(x))/(n_orig*f_a_n(x))) *
        sum((gamma_n(w,x)*pi_star_n(w,x))/g_n(x,w))
    )^(1/3)
  }))
  
}



#' Construct estimator of nuisance influence function omega_n
#' 
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param Qc_n Conditional censoring survival function estimator returned by
#'     construct_Q_n
#' @param type Defaults to "estimated". Override with "true" for debugging. Note
#'     that type="true" only works for surv_true="Cox PH" and assumes that S_0
#'     and Sc_0 (i.e. the true functions) are passed in.
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(vals=NA, Q_n, Qc_n, type="estimated") {
  
  if (type=="estimated") {
    
    # Construct cumulative hazard estimator
    H_n <- function(t,w,a) { -1 * log(Q_n(t,w,a)) }
    
    # Construct memoised integral function
    omega_integral <- function(k,w,a) {
      if (k==0) {
        integral <- 0
      } else {
        i <- round(seq(C$appx$t_0,k,C$appx$t_0), -log10(C$appx$t_0))
        incr <- C$appx$t_0
        w_long <- as.data.frame(
          matrix(rep(w,length(i)), ncol=length(w), byrow=T)
        )
        a_long <- rep(a,length(i))
        integral <- 0.5 * sum(
          ( H_n(i,w_long,a_long) - H_n(i-incr,w_long,a_long) ) * (
            ( Q_n(i,w_long,a_long) * Qc_n(i,w_long,a_long) )^-1 +
              ( Q_n(i-incr,w_long,a_long) * Qc_n(i-incr,w_long,a_long))^-1
          )
        )
      }
    }
    omega_integral <- construct_superfunc(omega_integral, aux=NA,
                                          vec=c(1,2,1), vals=NA)
    
    fnc <- function(w,a,y_star,delta_star) {
      
      k <- round(min(y_star,C$t_0), -log10(C$appx$t_0))
      
      if (F) {
        return (0)
      } # DEBUG
      
      return(Q_n(C$t_0,w,a) * (
        (delta_star * In(y_star<=C$t_0)) / (Q_n(k,w,a) * Qc_n(k,w,a)) -
          omega_integral(k,w,a)
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
      if (L$dir=="decr") {
        lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*a - 1.7
      } else {
        lin <- C$alpha_1*w[1] + C$alpha_2*w[2] + alpha_3*(1-a) - 1.7
      }
      lin2 <- C$alpha_1*w[1] + C$alpha_2*w[2] - 1
      
      # Compute omega_0
      piece_1 <- exp(-1*lmbd*(C$t_0^v)*exp(lin))
      piece_2 <- (delta_star*In(y_star<=C$t_0)) /
        exp(-1*(lmbd*y_star^v*exp(lin)+lmbd2*y_star^v2*exp(lin2)))
      integral <- integrate(
        function(t) {
          t^(v-1) * exp(lin+lmbd*t^v*exp(lin)+lmbd2*t^v2*exp(lin2))
        },
        lower = 0,
        upper = min(C$t_0,y_star)
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
#' @param type Type of regression; one of c("cubic", "kernel", "kernel2")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n() among the observations for which delta==1
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, dat, type="Super Learner", vals=NA,
                              omega_n=NA, f_aIw_n=NA, f_a_n=NA,
                              f_a_delta1_n=NA) {
  
  # Construct pseudo-outcomes
  po <- (dat$weights*omega_n(dat$w,dat$a,dat$y_star,dat$delta_star))^2
  
  # Create dataframe for regression
  df <- data.frame(a=dat$a, po=po)
  if (sum(!is.finite(df$po))!=0) {
    df %<>% filter(is.finite(po))
    warning(paste("construct_gamma_n:", sum(!is.finite(df$po)),
                  "non-finite po values"))
  }
  
  # Setup
  x_grid <- round(seq(0,1,C$appx$a), -log10(C$appx$a))
  X <- cbind(dat$w, a=dat$a)
  W_reduced <- distinct(dat_orig$w)
  W_reduced <- cbind("w_index"=c(1:nrow(W_reduced)), W_reduced)
  newX <- expand.grid(
    w_index = W_reduced$w_index,
    a = x_grid
  )
  newX <- inner_join(W_reduced, newX, by="w_index")
  newX$w_index <- NULL
  
  # Run regression
  if (type=="Super Learner") {
    
    # Fit SuperLearner regression
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.loess",
                    "SL.nnet", "SL.ksvm", "SL.caret", "SL.rpartPrune",
                    "SL.svm")
    
    model_sl <- SuperLearner(Y=df$po, X=X, newX=newX, family="gaussian",
                             SL.library=SL.library, verbose=F)
    pred <- as.numeric(model_sl$SL.predict)
    if (sum(pred<0)!=0) {
      warning(paste("construct_gamma_n:", sum(pred<0),
                    "negative predicted values"))
    }
    if (F) {
      print(coef(model_sl))
    } # DEBUG: Print algorithm coeffs
    
    rm(model_sl)
    
    newX$index <- c(1:nrow(newX))
    reg <- function(w,a) {
      # Dynamically filter to select index
      cond <- paste0("a==",a,"")
      for (i in c(1:length(w))) {
        cond <- paste0(cond," & w",i,"==",w[i])
      }
      index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
      if (length(index)!=1) {
        stop(paste0("Error in construct_gamma_n; ", "w=(",
                    paste(w,collapse=","), "), a=",a))
      }
      
      # Return prediction
      return(pred[index])
    }
    
  }
  
  # Remove large intermediate objects
  rm(dat_orig,dat,omega_n,f_aIw_n)
  
  fnc <- reg
  return(construct_superfunc(fnc, aux=NA, vec=c(2,1), vals=vals))
  
}



#' Construct q_n nuisance estimator function
#' 
#' @param type One of c("new", "zero")
#' @return q_n nuisance estimator function
#' !!!!! Eventually phase out type="Super Learner"
construct_q_n <- function(type="Super Learner", dat, dat_orig,
                          omega_n=NA, g_n=NA, z_n=NA, gcomp_n=NA, # g_n_star=NA
                          alpha_star_n=NA, f_aIw_n=NA, Q_n=NA, Qc_n=NA,
                          vals=NA) {
  
  if (type=="new") {
    
    surv_deriv <- function(Q_n) {

      fnc <- function(t, w, a) {

        width <- 40
        t1 <- t - width/2
        t2 <- t + width/2
        if (t1<0) { t2<-width; t1<-0; }
        if (t2>C$t_0) { t1<-C$t_0-width; t2<-C$t_0; }
        t1 <- round(t1,-log10(C$appx$t_0))
        t2 <- round(t2,-log10(C$appx$t_0))
        ch_y <- Q_n(t2,w,a) - Q_n(t1,w,a)

        return(ch_y/width)

      }

      return(construct_superfunc(fnc, aux=NA, vec=c(1,2,1), vals=NA))

    }
    
    Q_n_deriv <- surv_deriv(Q_n)
    Qc_n_deriv <- surv_deriv(Qc_n)
    
    q_n_star_inner <- function(y_star, delta_star, w, a) {
      In(a!=0)*omega_n(w,a,y_star,delta_star)/g_n(a,w) + gcomp_n(a)
    }
    q_n_star_inner <- construct_superfunc(q_n_star_inner, aux=NA, vec=c(1,1,2,1), vals=NA)
    
    q_n_star <- function(y_star, delta_star, w, a, x) {
      In(a<=x)*q_n_star_inner(y_star, delta_star, w, a) -
        In(a!=0)*alpha_star_n(x)
    }
    q_n_star <- construct_superfunc(q_n_star, aux=NA, vec=c(1,1,2,1,1), vals=NA)
    
    f_n_srv <- function(y_star, delta_star, w, a) {
      if (delta_star==1) {
        return(-1*Qc_n(y_star,w,a)*Q_n_deriv(y_star,w,a))
      } else {
        return(-1*Q_n(y_star,w,a)*Qc_n_deriv(y_star,w,a))
      }
    }
    f_n_srv <- construct_superfunc(f_n_srv, aux=NA, vec=c(1,1,2,1), vals=NA)
    
    n <- length(dat$a)
    
    fnc <- function(w, y_star, delta_star, x) {
      
      y_star_ <- rep(y_star,n)
      delta_star_ <- rep(delta_star,n)
      w_ <- as.data.frame(matrix(rep(w,n), ncol=length(w), byrow=T))
      x_ <- rep(x,n)
      f_n_srv_ <- f_n_srv(y_star_, delta_star_, w_, dat$a)
      q_n_star_ <- q_n_star(y_star, delta_star_, w_, dat$a, x_)
      g_n_ <- g_n(dat$a,w_)
      
      denom <- sum(dat$weights*f_n_srv_*g_n_)
      if (denom==0) {
        return (0)
      } else {
        num <- sum(dat$weights*q_n_star_*f_n_srv_*g_n_)
        return(num/denom)
      }
      
    }
    
  }
  
  if (type=="zero") {
    
    fnc <- function(w, y_star, delta_star, x) { 0 }
    
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(2,1,1,0), vals=vals))
  
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
#' @param edge_corr Edge correction (see est_curve)
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of A is [0,1]
construct_f_aIw_n <- function(dat, vals=NA, type, k=0, delta1=F,
                              edge_corr="none") {
  
  if (delta1) { dat$weights <- rep(1, length(dat$weights)) }
  
  # Create normalization factor (to multiply density by if there is a point mass
  # at the edge)
  if (edge_corr=="none") {
    norm_factor <- 1
  } else {
    norm_factor <- sum(dat$a!=0)/length(dat$a)
    dat <- ss(dat, which(dat$a!=0))
  }
  
  if (type=="true") {
    
    # Note: this is not accurate if edge!="none"
    
    if (L$distr_S=="Unif(0,1)") {
      fnc <- function(a, w) { 1 }
    } else if (L$distr_S=="Unif(0.3,0.7)") {
      fnc <- function(a, w) { 2.5 * In(a>=0.3 & a<=0.7) }
    } else if (L$distr_S=="N(0.5,0.01)") {
      fnc <- function(a, w) { dtruncnorm(a, a=0, b=1, mean=0.5, sd=0.1) }
    } else if (L$distr_S=="N(0.5,0.04)") {
      fnc <- function(a, w) { dtruncnorm(a, a=0, b=1, mean=0.5, sd=0.2) }
    } else if (L$distr_S=="N(0.3+0.4x2,0.04)") {
      fnc <- function(a, w) {
        dtruncnorm(a, a=0, b=1, mean=0.3+0.4*w[2], sd=0.2)
      }
    } else if (L$distr_S=="Tri UP") {
      fnc <- function(a, w) { 2*a }
    } else if (L$distr_S=="Tri DN") {
      fnc <- function(a, w) { 2-2*a }
    }
    
    
  } else if (type=="parametric") {
    
    # Density for a single observation
    dens_s <- function(a, w, prm) {
      mu <- sum( c(1,w) * c(expit(prm[1]),prm[2:(length(w)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( dtruncnorm(a, a=0, b=1, mean=mu, sd=sigma) )
      # return( dnorm(a, mean=mu, sd=sigma) )
    }
    
    # Set up weighted likelihood
    wlik <- function(prm) {
      -1 * sum(sapply(c(1:length(dat$a)), function(i) {
        dat$weights[i] * log(pmax(dens_s(a=dat$a[i], w=as.numeric(dat$w[i,]), prm),1e-8))
      }))
    }
    
    # Run optimizer
    opt <- solnp(
      pars = c(rep(0, length(dat$w)+1), 0.15),
      fun = wlik
    )
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: solnp() did not converge")
    }
    prm <- opt$pars
    
    # Remove large intermediate objects
    rm(dat,dens_s,opt)
    
    fnc <- function(a, w) {
      mu <- sum( c(1,w) * c(expit(prm[1]),prm[2:(length(w)+1)]) )
      sigma <- 10*expit(prm[length(prm)])
      return( dtruncnorm(a, a=0, b=1, mean=mu, sd=sigma) )
    }
    
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
        hz <- expit(
          prm[c(1:(ifelse(bin==k,k-1,bin)))] +
            as.numeric(prm[c(k:(k+length(w)-1))] %*% w)
        )
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        return(k*p1*p2)
      }, vec=c(1,2,0))
      
      # Set up weighted likelihood
      wlik <- function(prm) {
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
      prm <- opt$pars
      
      # Remove large intermediate objects
      rm(dat,dens_s,opt)
      
      fnc <- function(a, w) {
        
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(prm[j] + prm[c(k:(k+length(w)-1))] %*% w)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        
        return(norm_factor*k*p1*p2)
        
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
construct_f_a_n <- function(dat_orig, vals=NA, f_aIw_n) {
  
  fnc <- function(a) {
    mean(f_aIw_n(rep(a,nrow(dat_orig$w)),dat_orig$w))
  }
  
  # round_args <- -log10(C$appx$a))
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct density ratio estimator g_n
#' 
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
construct_g_n <- function(f_aIw_n, f_a_n, vals=NA) {
  
  fnc <- function(a,w) {
    f_aIw_n(a,w) / f_a_n(a)
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' Construct nuisance estimator eta_n
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @return Estimator function of nuisance eta_0
construct_eta_n <- function(dat, vals=NA, Q_n) {
  
  n_orig <- sum(dat$weights)
  
  fnc <- function(x,w) {
    
    w_long <- as.data.frame(
      matrix(rep(w,length(dat$a)), ncol=length(w), byrow=T)
    )

    return(
      (1/n_orig) * sum(
        dat$weights * In(dat$a<=x) *
          (1-Q_n(rep(C$t_0,length(dat$a)),w_long,dat$a))
      )
    )
  }
  
  # round_args <- c(-log10(C$appx$a), -log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' Construct nuisance estimator eta*_n
#' 
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @return Estimator function of nuisance eta*_0
construct_etastar_n <- function(Q_n, vals=NA) {
  
  fnc <- function(x,w) {
    x <- round(x,-log10(C$appx$a))
    if (x==0) {
      integral <- 0
    } else {
      integral <- sum(sapply(seq(C$appx$a,x,C$appx$a), function(a) {
        C$appx$a * Q_n(C$t_0, w, round(a,-log10(C$appx$a)))
      }))
    }
    return(x-integral)
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
lambda <- function(dat, k, G) {
  
  n_orig <- sum(dat$weights)
  lambda <- (1/n_orig) * sum( dat$weights * (G(dat$a))^k )
  return(lambda)
  
}



#' Construct Theta_os_n primitive one-step estimator
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param etastar_n A nuisance estimator returned by construct_etastar_n()
#' @return Gamma_os_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Theta_os_n <- function(dat, vals=NA, omega_n, f_aIw_n, etastar_n) {
  
  weights_i <- dat$weights
  n_orig <- sum(weights_i)
  a_i <- dat$a
  w_i <- dat$w
  piece_1 <- omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
    f_aIw_n(dat$a,dat$w)
  
  # Remove large intermediate objects
  rm(dat,omega_n,f_aIw_n)
  
  fnc <- function(x) {
    (1/n_orig) * sum(weights_i * (
      In(a_i<=x) * piece_1 + etastar_n(rep(x,nrow(w_i)),w_i)
    ))
  }
  
  # round_args <- -log10(C$appx$a)
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct Gamma_os_n primitive one-step estimator
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc()
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param Q_n A conditional survival function returned by construct_Q_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @param type One of c("one-step", "plug-in")
#' @return Gamma_os_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Gamma_os_n <- function(dat, vals=NA, omega_n, Q_n, g_n,
                                 type="one-step") {
  
  n_orig <- sum(dat$weights)
  n <- length(dat$a)
  i_long <- rep(c(1:n), each=n)
  j_long <- rep(c(1:n), times=n)
  a_i_short <- dat$a
  a_i_long <- dat$a[i_long]
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
      Q_n(rep(C$t_0, length(a_i_long)),w_j_long,a_i_long)
    
    # Remove large intermediate objects
    rm(dat,delta_star_i_long,delta_star_j_long,i_long,j_long,omega_n,Q_n,
       w_i_long,w_j_long,weights_i_long,weights_j_long)
    
    fnc <- function(a) {
      
      subpiece_1b <- In(round(a_i_short,-log10(C$appx$a))<=
                          round(a,-log10(C$appx$a)))
      piece_1 <- (1/n_orig) * sum(subpiece_1a*subpiece_1b)
      
      subpiece_2b <- In(round(a_i_long,-log10(C$appx$a))<=
                                  round(a,-log10(C$appx$a)))
      piece_2 <- (1/(n_orig^2)) * sum(subpiece_2a*subpiece_2b)
      
      return(piece_1-piece_2)
      
    }
    
  }
  
  if (type=="plug-in") {
    
    piece <- weights_i_long*weights_j_long *
      (1 - Q_n(rep(C$t_0, length(a_i_long)),w_j_long,a_i_long))
    
    # Remove large intermediate objects
    rm(dat,delta_star_i_long,delta_star_j_long,i_long,j_long,omega_n,Q_n,
       w_i_long,w_j_long,weights_i_long,weights_j_long)
    
    fnc <- function(a) {
      return((1/n_orig^2) * sum(piece * In(a_i_long<=a)))
    }
    
  }
  
  # !!!!! Use aux; also dat should not be in the environment of returned function but it is
  
  # round_args <- -log10(C$appx$a)
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_rho_n <- function(dat, Phi_n, vals=NA) {
  
  n_orig <- sum(dat$weights)
  
  fnc <- function(a) {
    (1/n_orig) * sum(
      dat$weights * (Phi_n(dat$a)^3) * (In(a<=dat$a) - Phi_n(dat$a))
    )
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_xi_n <- function(Phi_n, lambda_2, lambda_3, vals=NA) {
  
  fnc <- function(a_i,a_j) {
    (2*In(a_i<=a_j) - 3*Phi_n(a_j))*Phi_n(a_j)*lambda_2 +
    (2*Phi_n(a_j) - In(a_i<=a_j))*lambda_3
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(1,1), vals=vals))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_1 <- function(dat, Gamma_os_n, Phi_n, lambda_2,
                                lambda_3, vals=NA) {
  
  n_orig <- sum(dat$weights)
  weights_j <- dat$weights
  a_j <- dat$a
  
  rho_n <- function(a,x) {
    (1/n_orig) * sum(
      weights_j * (Phi_n(a))^x * Gamma_os_n(round(a,-log10(C$appx$a)))
    )
  }
  
  rho_1 <- rho_n(a_j,1)
  rho_2 <- rho_n(a_j,2)
  piece_01 <- Gamma_os_n(round(a_j,-log10(C$appx$a)))
  piece_20 <- (Phi_n(a_j))^2
  piece_10 <- Phi_n(a_j)
  piece_11 <- Phi_n(a_j) * Gamma_os_n(round(a_j,-log10(C$appx$a)))
  
  fnc <- function(a_i) {
    
    piece_1 <- (2*(1/n_orig)*sum(weights_j*In(a_i<=a_j)*piece_10)+
                  (Phi_n(a_i))^2-6*lambda_2)*rho_2
    piece_2 <- lambda_2*(
      2*(1/n_orig)*sum(weights_j*In(a_i<=a_j)*piece_11)+
        (Phi_n(a_i))^2*Gamma_os_n(round(a_i,-log10(C$appx$a)))
    )
    piece_3 <- (3*(1/n_orig)*sum(weights_j*In(a_i<=a_j)*piece_20)+
                  (Phi_n(a_i))^3-6*lambda_3)*rho_1
    piece_4 <- lambda_3*((1/n_orig)*sum(weights_j*In(a_i<=a_j)*
                                          piece_01)+Phi_n(a_i)*
                           Gamma_os_n(round(a_i,-log10(C$appx$a))))
    
    return(piece_1+piece_2-piece_3-piece_4)
    
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=vals))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Gamma <- function(omega_n, g_n, gcomp_n, eta_n,
                                    Gamma_os_n) {
  
  fnc <- function(x,w,y_star,delta_star,a) {
    In(a<=x)*((omega_n(w,a,y_star,delta_star)/g_n(a,w)) + gcomp_n(a)) +
      eta_n(x,w) -
      2*Gamma_os_n(round(x,-log10(C$appx$a)))
  }
  
  return(construct_superfunc(fnc, vec=c(1,2,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Gamma2 <- function(omega_n, g_n_star, gcomp_n, z_n,
                                     alpha_star_n, q_n, eta_ss_n,
                                     Gamma_os_n_star) {
  
  fnc <- function(x,w,y_star,delta_star,a,wt) {
    if (wt==0) {
      piece_1 <- 0
      piece_2 <- 0
      piece_3 <- 0
    } else {
      piece_1 <- In(a!=0 & a<=x)
      piece_2 <- omega_n(w,a,y_star,delta_star)/g_n_star(a,w) + gcomp_n(a)/z_n
      piece_3 <- In(a!=0)*alpha_star_n(x)
    }
    wt*(piece_1*piece_2-piece_3/z_n) +
      (1-wt)*q_n(w,y_star,delta_star,x) +
      eta_ss_n(x,w) -
      Gamma_os_n_star(round(x,-log10(C$appx$a)))
  }
  
  return(construct_superfunc(fnc, vec=c(1,2,1,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Theta <- function(omega_n, f_aIw_n, q_star_n, etastar_n,
                                    Theta_os_n) {
  
  fnc <- function(x,w,y_star,delta_star,a,wt) {
    if (wt==0) {
      piece_1 <- 0
      piece_2 <- 0
    } else {
      piece_1 <- In(a<=x)
      piece_2 <- omega_n(w,a,y_star,delta_star)/f_aIw_n(a,w)
    }
    wt*piece_1*piece_2 +
      (1-wt) * q_star_n(w,y_star,delta_star,x) +
      etastar_n(x,w) -
      Theta_os_n(round(x,-log10(C$appx$a)))
  }
  
  return(construct_superfunc(fnc, vec=c(1,2,1,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_2 <- function(dat, Phi_n, infl_fn_Gamma, lambda_2, lambda_3) {
  
  n_orig <- sum(dat$weights)
  weights_j <- dat$weights
  a_j <- dat$a
  len <- length(a_j)
  
  # fnc <- function(w,y_star,delta_star,a) {
  fnc <- function(w,y_star,delta_star,a,wt) {
    w_long <- as.data.frame(
      matrix(rep(w,len), ncol=length(w), byrow=T)
    )
    (1/n_orig) * sum(weights_j * (
      ( lambda_2*(Phi_n(a_j)^2) - lambda_3*Phi_n(a_j) ) *
        # infl_fn_Gamma(a_j,w_long,rep(y_star,len),rep(delta_star,len),rep(a,len))
        infl_fn_Gamma(a_j, w_long, rep(y_star,len), rep(delta_star,len),
                      rep(a,len), rep(wt,len))
    ))
  }
  
  # return(construct_superfunc(fnc, vec=c(2,1,1,1)))
  return(construct_superfunc(fnc, vec=c(2,1,1,1,1)))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
beta_n_var_hat <- function(dat, dat_orig, infl_fn_1, infl_fn_2) {
  
  n_orig <- sum(dat$weights)
  
  b_sum <- 0
  for (i in c(1:length(dat_orig$a))) {
    if (dat_orig$delta[i]==0) {
      b_sum <- b_sum + (infl_fn_2(dat_orig$w[i,], dat_orig$y_star[i],
                                  dat_orig$delta_star[i], dat_orig$a[i],
                                  dat_orig$weights[i]))^2
    } else {
      b_sum <- b_sum + (
        dat_orig$weights[i]*infl_fn_1(dat_orig$a[i]) +
          infl_fn_2(dat_orig$w[i,], dat_orig$y_star[i], dat_orig$delta_star[i],
                    dat_orig$a[i], dat_orig$weights[i])
      )^2
    }
  }
  
  return( (1/n_orig)*b_sum )
  
}



#' Construct lists of values to pre-compute functions on
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param factor_A If true, factor_A is used instead of seq(0,1,C$appx$a) for the
#'     Q_n component only
#' @return A list of lists
create_val_list <- function(dat_orig, factor_A=NA) {
  
  names(dat_orig$w) <- paste0("w", c(1:length(dat_orig$w)))
  W_reduced <- distinct(dat_orig$w)
  W_reduced <- cbind("w_index"=c(1:nrow(W_reduced)), W_reduced)
  if (is.na(factor_A[1])) {
    a <- round(seq(0,1,C$appx$a),-log10(C$appx$a))
  } else {
    a <- factor_A
  }
  Q_n_pre <- expand.grid(
    t = round(seq(0,C$t_0,C$appx$t_0),-log10(C$appx$t_0)),
    w_index = W_reduced$w_index,
    a = a
  )
  Q_n_pre <- inner_join(Q_n_pre, W_reduced, by="w_index")
  
  return(list(
    A = NA, # list(a=dat$a),
    AW = NA, # data.frame(a=dat$a, w1=dat$w1, w2=dat$w2),
    A_grid = NA, # list(a=round(seq(0,1,C$appx$a),-log10(C$appx$a)),
    W_grid = NA,
    AW_grid = NA,
    Q_n = list(
      t = Q_n_pre$t,
      w = subset(Q_n_pre, select=-c(t,w_index,a)),
      a = Q_n_pre$a
    ),
    omega = NA # list(w=rbind(dat$w,dat$w), a=c(dat$a, rep(0, length(dat$a))), y_star=rep(dat$y_star,2), delta_star=rep(dat$delta_star,2))
  ))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_Gamma_cf_k <- function(dat_train, dat_test, vals=NA, omega_n, g_n,
                                 gcomp_n, eta_n) {
  
  # !!!!! Needs to be updated
  
  n_test <- length(dat_test$a)
  dat_test$weights <- wts(dat_test) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d1 <- ss(dat_test, which(dat_test$delta==1))
  weights_1 <- d1$weights
  
  n_train <- length(dat_train$a)
  dat_train$weights <- wts(dat_train) # !!!!! Weights need to be re-calculated and/or re-stabilized here
  d2 <- ss(dat_train, which(dat_train$delta==1))
  weights_2 <- d2$weights
  
  fnc <- function(a) {
    
    piece_1 <- (1/n_test) * sum(weights_1 * (
      In(d1$a<=a) *
        ((omega_n(d1$w,d1$a,d1$y_star,d1$delta_star) / g_n(d1$a,d1$w)) +
           gcomp_n(d1$a)
        ) +
        eta_n(a,d1$w)
    ))
    
    piece_2 <- (1/n_train) * sum(weights_2 * In(d2$a<=a)*gcomp_n(d2$a))
    
    return(piece_1-piece_2)
    
  }
  
  # round_args <- -log10(C$appx$a)
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct cross-fitted Gamma_0 estimator
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param params The same params object passed to est_curve
#' @param vlist A list of dataframes returned by create_val_list(); Q_n REQUIRED
#'     but others can be NA
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
    srvSL <- construct_Q_n(dat_train, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    gcomp_n <- construct_gcomp_n(dat_orig_train, vlist$A_grid, Q_n)
    eta_n <- construct_eta_n(dat_train, vlist$AW_grid, Q_n)
    f_aIw_n <- construct_f_aIw_n(dat_train, vlist$AW_grid, type=params$g_n_type,
                                 k=15)
    f_a_n <- construct_f_a_n(dat_orig_train, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n,
                                 type=params$omega_n_type)
    
    # Construct K functions
    Gamma_cf_k[[k]] <- construct_Gamma_cf_k(
      dat_train, dat_test, vlist$A_grid, omega_n, g_n, gcomp_n, eta_n
    )
    
    # Remove objects
    rm(Phi_n,Phi_n_inv,Q_n,Qc_n,gcomp_n,eta_n,f_aIw_n,f_a_n,g_n,omega_n)
    
  }
  
  # Construct cross-fitted Gamma_os_n
  return(Vectorize(function(a) {
    mean(sapply(c(1:params$cf_folds), function(k) {
      Gamma_cf_k[[k]](a)
    }))
  }))
  
}



#' Construct propensity score estimator of pi_0 = P(A=0|W=w)
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param vals List of values to pre-compute function on; passed to
#'     construct_superfunc(); REQUIRED FOR SUPERLEARNER
#' @param type One of c("true", "logistic", "Super Learner", "generalized"). If
#'     type=="true", the only valid value is zero. If type=="generalized", the
#'     arguments `f_aIw_n` and `cutoffs` must also be supplied.
#' @return Propensity score estimator of pi_0
#' @notes For all types except for "generalized", this function constructs the
#'     probability P(A=0|W=w). The type "generalized" constructs the probability
#'     P(A=a|W=w) for a generic value a that has positive mass.
construct_pi_n <- function(dat, vals=NA, type, f_aIw_n=NA, cutoffs=NA) {
  
  # Construct indicator I{A=0}
  ind_A0 <- In(dat$a==0)
  
  if (type=="true") {
    
    # Note: These are only valid for val==0
    if (L$edge=="expit") {
      fnc <- function(w, val) {
        if (val==0) { expit(w[1]+w[2]-3.3) }
      }
    } else if (L$edge=="expit2") {
      fnc <- function(w, val) {
        if (val==0) { expit(w[1]+w[2]-1) }
      }
    } else if (L$edge=="Complex") {
      fnc <- function(w, val) {
        if (val==0) { 0.84*w[2]*pmax(0,1-4*abs(w[1]-0.5)) }
      }
    } else if (L$edge=="none") {
      fnc <- function(w, val) {
        if (val==0) { 0 }
      }
    }
    
  } else if (type=="logistic") {
    
    fml <- "ind_A0~1"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind(ind_A0, dat$w)
    suppressWarnings({
      model <- glm(
        fml,
        data = df,
        family = "binomial",
        weights = dat$weights
      )
    })
    coeffs <- model$coefficients
    
    # !!!!! val is unused
    fnc <- function(w, val) { as.numeric(expit(coeffs %*% c(1,w))) }
    
  } else if (type=="SL") {
    
    # sl <- SuperLearner(
    #   Y = ind_A0,
    #   X = dat$w,
    #   newX = vals,
    #   family = binomial(),
    #   SL.library = "SL.earth", # SL.glm SL.gbm SL.ranger SL.earth
    #   # SL.library = c("SL.earth", "SL.gam", "SL.ranger"), # SL.glm SL.gbm SL.ranger SL.earth
    #   obsWeights = dat$weights,
    #   control = list(saveFitLibrary=FALSE)
    # )
    # assign("sl", sl, envir=.GlobalEnv) # ?????
    # 
    # fnc <- function(w, val) {
    #   
    #   r <- list()
    #   for (i in 1:length(w)) {
    #     r[[i]] <- which(abs(w[i]-newX[[paste0("w",i)]])<1e-8)
    #   }
    #   index <- Reduce(intersect, r)
    #   return(sl$SL.predict[index])
    # }
    
  } else if (type=="generalized") {
    
    fnc <- function(w, val) {
      
      # val should be an index of the bin
      bin_start <- cutoffs[as.numeric(val)]
      bin_stop <- cutoffs[as.numeric(val)+1]
      grid <- seq(bin_start, bin_stop, length.out=11)[1:10]
      avg_height <- mean(sapply(grid, function(a) { f_aIw_n(a, w=w) }))
      
      return((bin_stop-bin_start)*avg_height)
      
    }
    
  }
  
  # round_args <- c(-log10(C$appx$w1), 0)
  return(construct_superfunc(fnc, aux=NA, vec=c(2,1), vals=vals))
  
}



#' Construct estimator of pi*_0 = P(Delta=0|W=w,A=a)
#' 
#' @param x !!!!!
construct_pi_star_n <- function(dat_orig, vals=NA, type="Super Learner",
                                f_aIw_n=NA, f_aIw_delta1_n=NA) {
  
  # Set library
  if (type=="Super Learner") {
    SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
                    "SL.svm", "SL.glmnet")
  } else if (type=="logistic") {
    SL.library <- c("SL.glm")
  }
  
  # Create data objects
  X <- dat_orig$w
  newX <- distinct(dat_orig$w)
  
  # Fit SuperLearner regression
  model_sl <- SuperLearner(Y=dat_orig$delta, X=X, newX=newX, family="binomial",
                           SL.library=SL.library, verbose=F)
  
  pred <- as.numeric(model_sl$SL.predict)
  rm(model_sl)
  
  # Construct function
  newX$index <- c(1:nrow(newX))
  reg <- function(w) {
    
    # Dynamically filter to select index
    cond <- "1==1"
    for (i in c(1:length(w))) {
      cond <- paste0(cond," & round(w",i,",8)==",round(w[i],8))
    }
    index <- (dplyr::filter(newX, eval(parse(text=cond))))$index
    if (length(index)!=1) {
      stop(paste0("Error in construct_pi_star_n; ",
                  "w=(", paste(w,collapse=","), ")"))
    }
    
    # Return prediction
    return(pred[index])
    
  }
  
  # Check results
  if (F) {
    d0 <- ss(dat_orig, which(dat_orig$w$w2==0))
    d1 <- ss(dat_orig, which(dat_orig$w$w2==0))
    grid <- round(seq(0,1,0.1),1)
    preds0 <- sapply(grid, function(w1) { reg(c(w1,0)) })
    preds1 <- sapply(grid, function(w1) { reg(c(w1,1)) })
    ggplot(
      data.frame(x=d0$w$w1, y=d0$delta),
      aes(x=x, y=y)
    ) +
      geom_jitter(alpha=0.3, width=0.04, height=0.04) +
      geom_line(data=data.frame(x=grid, y=preds0), color="forestgreen")
    ggplot(
      data.frame(x=d1$w$w1, y=d1$delta),
      aes(x=x, y=y)
    ) +
      geom_jitter(alpha=0.3, width=0.04, height=0.04) +
      geom_line(data=data.frame(x=grid, y=preds1), color="forestgreen")
  }
  
  fnc <- function(w,a) {
    (f_aIw_delta1_n(a,w) / f_aIw_n(a,w)) * reg(w)
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(2,1), vals=vals))
  
}



#' Compute one-step estimator of counterfactual survival at A=0
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param val Value of A
#' @return Value of one-step estiamtor
r_Mn_edge <- function(dat, pi_n, Q_n, omega_n, val=0) {
  
  n_orig <- sum(dat$weights)
  n_dat <- nrow(dat$w)
  
  return(
    1 - (1/n_orig) * sum(dat$weights * (
      Q_n(rep(C$t_0,n_dat),dat$w,a=rep(val,n_dat)) - (
        (In(dat$a==val)/pi_n(dat$w, rep(val,n_dat))) *
          omega_n(dat$w,a=rep(val,n_dat),dat$y_star,dat$delta_star)
      )
    ))
  )
  
}



#' Compute asymptotic variance of one-step estimator r_Mn_edge
#' 
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param pi_n Propensity score estimator returned by construct_pi_n()
#' @param Q_n Conditional survival function estimator returned by construct_Q_n
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param r_Mn_edge_est Estimate returned by one-step estimator r_Mn_edge()
#' @param val Value of A
#' @return Asymptotic variance estimate
sigma2_edge <- function(dat, pi_n, Q_n, omega_n, r_Mn_edge_est, val=0) {
  
  n_orig <- sum(dat$weights)
  n_dat <- nrow(dat$w)
  
  return(
    (1/n_orig) * sum((dat$weights * (
      Q_n(rep(C$t_0,n_dat),dat$w,a=rep(val,n_dat)) - (
        (In(dat$a==val)/pi_n(dat$w, rep(val,n_dat))) *
          omega_n(dat$w,a=rep(val,n_dat),dat$y_star,dat$delta_star)
      ) -
        (1-r_Mn_edge_est)
    ))^2)
  )
  
}



#' Construct density ratio estimator g*_n
#' 
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
#' @notes The indicator functions are for computational conveniene; the actual
#'     function is technically undefined there
construct_g_n_star <- function(f_aIw_n, f_a_n, z_n) {
  
  function(a,w) {
    In(a==0) + In(a!=0) * (z_n*(f_aIw_n(a,w)/f_a_n(a)))
  }
  
}



#' Construct nuisance estimator alpha*_n
#' 
#' @param dat !!!!!
#' @param gcomp !!!!!
#' @param z_n !!!!!
#' @param vals !!!!!
#' @return Nuisance estimator function
construct_alpha_star_n <- function(dat, gcomp_n, z_n, vals=NA) {
  
  n_orig <- sum(dat$weights)
  piece_1 <- dat$weights * In(dat$a!=0) * gcomp_n(dat$a)
  
  fnc <- function(x) {
    (1/(n_orig*z_n)) * sum(In(dat$a<=x)*piece_1)
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct nuisance estimator eta_ss_n
#' 
#' @param dat !!!!!
#' @param Q_n !!!!!
#' @param z_n !!!!!
#' @param vals !!!!!
#' @return Nuisance estimator function
construct_eta_ss_n <- function(dat, Q_n, z_n, vals=NA) {
  
  n_orig <- sum(dat$weights)
  piece_1 <- dat$weights * In(dat$a!=0)
  
  fnc <- function(x,w) {
    w_long <- as.data.frame(
      matrix(rep(w,length(dat$a)), ncol=length(w), byrow=T)
    )
    (1/(n_orig*z_n)) * sum(
      piece_1 * In(dat$a<=x) *
        (1-Q_n(rep(C$t_0,length(dat$a)),w_long,dat$a))
    )
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  
}



#' Construct Gamma_os_n_star primitive one-step estimator
#'
#' @param x !!!!!
construct_Gamma_os_n_star <- function(dat, omega_n, g_n_star, eta_ss_n, z_n,
                                      gcomp_n, alpha_star_n, vals=NA) {

  weights_i <- dat$weights
  n_orig <- sum(weights_i)
  a_i <- dat$a
  w_i <- dat$w
  piece_1 <- omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
    g_n_star(dat$a,dat$w)
  piece_2 <- In(a_i!=0)
  piece_3 <- gcomp_n(a_i)

  # Remove large intermediate objects
  rm(dat,omega_n,g_n_star,gcomp_n)
  
  fnc <- function(x) {
    (1/n_orig) * sum(weights_i * (
      piece_2*In(a_i<=x)*piece_1 +
        eta_ss_n(rep(x,nrow(w_i)),w_i) +
        (piece_2/z_n)*(In(a_i<=x)*piece_3-alpha_star_n(x))
    ))
  }

  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))

}



#' Construct Gamma_os_n_star primitive one-step estimator (based on EIF)
#' 
#' @param x !!!!!
construct_Gamma_os_n_star2 <- function(dat, dat_orig, omega_n, g_n, # g_n_star
                                       eta_ss_n, z_n, q_n, gcomp_n,
                                       alpha_star_n, vals=NA) {
  
  n_orig <- round(sum(dat$weights))
  piece_1 <- In(dat$a!=0)
  piece_2 <- (omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
                g_n(dat$a,dat$w)) + gcomp_n(dat$a)
  piece_3 <- (1-dat_orig$weights)
  
  # Remove large intermediate objects
  rm(omega_n,g_n,gcomp_n)
  
  fnc <- function(x) {
    
    if (F) {
      eta_ss_n(x,c(0,0))
      return(
        (1/n_orig) * sum(
          eta_ss_n(rep(x,n_orig),dat_orig$w)
        )
      )
    } # DEBUG
    
    (1/(n_orig*z_n)) * sum(dat$weights * (
      piece_1*In(dat$a<=x)*piece_2 - piece_1*alpha_star_n(x)
    )) +
      (1/n_orig) * sum(
        piece_3 * (
          q_n(dat_orig$w,dat_orig$y_star,dat_orig$delta_star,x)/z_n
        ) +
          eta_ss_n(rep(x,n_orig),dat_orig$w)
      )
    
  }
  
  return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}



#' Construct Theta_os_n primitive one-step estimator (based on EIF)
#' 
#' @param x !!!!!
construct_Theta_os_n2 <- function(dat, dat_orig, omega_n, f_aIw_n, q_star_n,
                                  etastar_n, vals=NA) { # eta_ss_n, z_n, q_n, gcomp_n, alpha_star_n
  
  # n_orig <- round(sum(dat$weights))
  # piece_1 <- In(dat$a!=0)
  # piece_2 <- (omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
                # g_n(dat$a,dat$w)) + gcomp_n(dat$a)
  # piece_3 <- (1-dat_orig$weights)
  
  # Remove large intermediate objects
  # rm(omega_n,g_n,gcomp_n)
  
  # fnc <- function(x) {
  #   (1/(n_orig*z_n)) * sum(dat$weights * (
  #     piece_1*In(dat$a<=x)*piece_2 - piece_1*alpha_star_n(x)
  #   )) +
  #     (1/n_orig) * sum(
  #       piece_3 * (
  #         q_n(dat_orig$w,dat_orig$y_star,dat_orig$delta_star,x)/z_n
  #       ) +
  #         eta_ss_n(rep(x,n_orig),dat_orig$w)
  #     )
  # }
  
  # n_orig <- round(sum(dat$weights))
  # piece_2 <- omega_n(dat$w,dat$a,dat$y_star,dat$delta_star)/f_aIw_n(dat$a,dat$w)
  # piece_3 <- (1-dat_orig$weights)
  
  # Remove large intermediate objects
  # rm(omega_n,f_aIw_n)
  
  # fnc <- function(x) {
  #   (1/n_orig) * sum(dat$weights * In(dat$a<=x) * piece_2) +
  #     (1/n_orig) * sum(
  #       piece_3 * q_star_n(dat_orig$w,dat_orig$y_star,dat_orig$delta_star,x) +
  #         etastar_n(rep(x,n_orig),dat_orig$w)
  #     )
  # }
  
  # return(construct_superfunc(fnc, aux=NA, vec=T, vals=vals))
  
}





# Estimate the variance of the Cox model based marginalized survival estimator
#' 
#' @param dat_orig Dataset returned by generate_data()
#' @param dat Subsample of dataset returned by ss() for which delta==1
#' @param t The end time of interest
#' @param points The A-values of interest
#' @param se_beta Estimates related to parameter vector beta_n
#' @param se_bshz Estimates related to Breslow baseline cuml. hazard estimator
#' @param se_surv Estimates related to survival at a point
#' @param se_marg Estimates related to marginalized survival
#' @param return_extras For debugging
#' @param parallel Parallelize internal lapply calls
#' @param verbose Print status updates and progress bars
#' @return A list containing the selected return values
cox_var <- function(dat_orig, dat, t, points, se_beta=F, se_bshz=F,
                    se_surv=F, se_marg=F, return_extras=F, parallel=F,
                    verbose=F) {
  # Setup
  if (verbose) {
    print(paste("Check 0 (start):", Sys.time()))
    pbapply::pboptions(type="txt", char="#", txt.width=40, style=3)
  } else {
    pbapply::pboptions(type="none")
  }
  if (parallel) {
    n_cores <- parallel::detectCores() - 1
    cl <- parallel::makeCluster(n_cores)
    if (verbose) {
      print(cl)
    }
  } else {
    cl <- NULL
  }
  
  # Fit a Cox model
  # Note: scaling the weights affects the SEs but not the estimates; thus, this
  #       is only needed for debugging
  model <- coxph(
    formula = formula(paste0("Surv(y_star,delta_star)~",
                             paste(names(dat$w),collapse="+"),"+a")),
    data = cbind(y_star=dat$y_star, delta_star=dat$delta_star,
                 dat$w, a=dat$a),
    weights = dat$weights * (length(dat$weights)/sum(dat$weights))
  )
  beta_n <- as.numeric(model$coefficients)
  
  if (verbose) { print(paste("Check 1 (Cox model fit):", Sys.time())) }
  
  # !!!!! Stabilize weights
  
  # Alias random variables
  WT <- dat$weights
  ST <- dat$strata
  N <- length(dat_orig$a)
  n <- length(dat$a)
  Z_ <- t(as.matrix(cbind(dat$w,a=dat$a)))
  T_ <- dat$y_star
  Ds_ <- dat$delta_star
  lin <- as.numeric(t(beta_n)%*%Z_)
  d <- dim(Z_)[1]
  
  # Intermediate functions
  {
    S_2n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            res <- matrix(NA, nrow=d, ncol=d)
            for (i in c(1:d)) {
              for (j in c(1:d)) {
                if (!is.na(res[j,i])) {
                  res[i,j] <- res[j,i]
                } else {
                  res[i,j] <- (1/N)*sum(WT*In(T_>=x)*Z_[i,]*Z_[j,]*exp(lin))
                }
              }
            }
            return(res)
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()
    
    S_0n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
              (1/N) * sum(WT*In(T_>=x)*exp(lin))
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()
    
    S_1n <- (function() {
      .cache <- new.env()
      function(x) {
        val <- .cache[[as.character(x)]]
        if (is.null(val)) {
          val <- (function(x) {
            (1/N) * as.numeric(Z_ %*% (WT*In(T_>=x)*exp(lin)))
          })(x)
          .cache[[as.character(x)]] <- val
        }
        return(val)
      }
    })()
    
    m_n <- function(x) {
      S_1n(x) / S_0n(x)
    }
    
    h <- function(x) {
      (S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x))
    }
    
  }
  
  # Create set of event times
  i_ev <- which(Ds_==1)
  t_ev <- T_[which(Ds_==1)]
  
  # Estimated information matrix (for an individual)
  I_tilde <- Reduce("+", lapply(i_ev, function(i) {
    WT[i] * h(T_[i])
  }))
  I_tilde <- (1/N)*I_tilde
  I_tilde_inv <- solve(I_tilde)
  
  # Score function (Cox model)
  l_star <- function(z_i,ds_i,t_i) {
    ds_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*exp(sum(z_i*beta_n))*In(T_[j]<=t_i)*(z_i-m_n(T_[j]))) /
        S_0n(T_[j])
    }))
  }
  
  # Influence function: beta_hat (regular Cox model)
  l_tilde <- (function() {
    .cache <- new.env()
    function(z_i,ds_i,t_i) {
      key <- paste(c(z_i,ds_i,t_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,ds_i,t_i) {
          I_tilde_inv %*% l_star(z_i,ds_i,t_i)
        })(z_i,ds_i,t_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Create p_n and p1_n vectors
  n_strata <- max(dat_orig$strata)
  p_n <- c()
  p1_n <- c()
  for (c in c(1:n_strata)) {
    p_n[c] <- mean(c==dat_orig$strata)
    p1_n[c] <- mean(c==dat_orig$strata & dat_orig$delta==1)
  }
  
  # Influence function: beta_hat (est. weights)
  lstar_tilde <- (function() {
    .cache <- new.env()
    
    exp_terms <- list()
    for (i in c(1:max(dat_orig$strata))) {
      js <- which(ST==i)
      if (length(js)>0) {
        exp_terms[[i]] <- (1/length(js)) * Reduce("+", lapply(js, function(j) {
          l_tilde(Z_[,j],Ds_[j],T_[j])
        }))
      } else {
        exp_terms[[i]] <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
      }
    }
    
    function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          if (d_i) {
            piece_1 <- wt_i*l_tilde(z_i,ds_i,t_i)
          } else {
            piece_1 <- as.matrix(rep(0,d))
          }
          
          if (p1_n[st_i]!=0) {
            piece_2 <- (1 - (d_i*p_n[st_i])/p1_n[st_i]) * exp_terms[[st_i]]
          } else {
            piece_2 <- as.matrix(rep(0,d)) # Hack to avoid NAs in small samples
          }
          
          return(as.numeric(piece_1+piece_2))
        })(z_i,d_i,ds_i,t_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Nuisance constant: mu_n
  mu_n <- (1/N) * as.numeric(Reduce("+", lapply(t_ev, function(t_j) {
    (In(t_j<=t) * S_1n(t_j)) / (S_0n(t_j))^2
  })))
  
  # Nuisance function: v_n
  v_n <- (function() {
    .cache <- new.env()
    function(st_i,d_i,t_j) {
      key <- paste(c(st_i,d_i,t_j), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(st_i,d_i,t_j) {
          k_set <- which(ST==st_i)
          if (length(k_set)>0) {
            return(
              (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
                unlist(lapply(k_set, function(k) {
                  In(T_[k]>=t_j) * exp(sum(beta_n*Z_[,k]))
                }))
              )
            )
          } else { return(0) }
        })(st_i,d_i,t_j)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Nuisance function: vstar_n
  vstar_n <- function(st_i,d_i,t) {
    k_set <- which(ST==st_i & Ds_==1)
    if (length(k_set)>0) {
      return(
        (1/N) * (p1_n[st_i]-p_n[st_i]*d_i)/(p1_n[st_i])^2 * sum(
          unlist(lapply(k_set, function(k) {
            In(T_[k]<=t) / S_0n(T_[k])
          }))
        )
      )
    } else { return(0) }
  }
  
  # Breslow estimator
  Lambda_n <- Vectorize((function() {
    .cache <- new.env()
    function(x) {
      key <- as.character(x)
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(x) {
          (1/N) * sum(unlist(lapply(i_ev, function(i) {
            (WT[i] * In(T_[i]<=x)) / S_0n(T_[i])
          })))
        })(x)
        .cache[[key]] <- val
      }
      return(val)
    }
  })())
  
  # Survival estimator (at a point)
  Q_n <- (function() {
    .cache <- new.env()
    function(z) {
      key <- paste(z, collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z) {
          exp(-exp(sum(z*beta_n))*Lambda_n(t))
        })(z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Influence function: Breslow estimator (est. weights)
  infl_fn_Bres <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i) {
          pc_2 <- vstar_n(st_i,d_i,t)
          pc_4 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
            ( In(t_j<=t) * v_n(st_i,d_i,t_j) ) / (S_0n(t_j))^2
          })))
          pc_5 <- sum(mu_n*lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))
          
          if (d_i==1) {
            pc_1 <- ( wt_i * ds_i * In(t_i<=t) ) / S_0n(t_i)
            pc_3 <- (1/N) * sum(unlist(lapply(t_ev, function(t_j) {
              (In(t_j<=t)*wt_i*In(t_i>=t_j)*exp(sum(beta_n*z_i))) /
                (S_0n(t_j))^2
            })))
            return(pc_1+pc_2-pc_3-pc_4-pc_5)
          } else {
            return(pc_2-pc_4-pc_5)
          }
        })(z_i,d_i,ds_i,t_i,wt_i,st_i)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Influence function: survival at a point (est. weights)
  omega_n <- (function() {
    .cache <- new.env()
    function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
      key <- paste(c(z_i,d_i,ds_i,t_i,wt_i,st_i,z), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z_i,d_i,ds_i,t_i,wt_i,st_i,z) {
          explin <- exp(sum(z*beta_n))
          piece_1 <- Lambda_n(t) * explin *
            (t(z)%*%lstar_tilde(z_i,d_i,ds_i,t_i,wt_i,st_i))[1]
          piece_2 <- explin * infl_fn_Bres(z_i,d_i,ds_i,t_i,wt_i,st_i)
          return(-1*Q_n(z)*(piece_1+piece_2))
        })(z_i,d_i,ds_i,t_i,wt_i,st_i,z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Influence function: marginalized survival (est. weights)
  infl_fn_marg <- (function() {
    .cache <- new.env()
    w_ <- as.list(as.data.frame(t(dat_orig$w)))
    function(w_i,a_i,d_i,ds_i,t_i,wt_i,st_i,a) {
      key <- paste(c(w_i,a_i,d_i,ds_i,t_i,wt_i,st_i,a), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(w_i,a_i,d_i,ds_i,t_i,wt_i,st_i,a) {
          piece_1 <- Q_n(c(w_i,a))
          piece_2 <- (1/N) * sum(unlist(lapply(c(1:N), function(j) {
            w_j <- w_[[j]]
            omgn <- omega_n(c(w_i,a_i),d_i,ds_i,t_i,wt_i,st_i,c(w_j,a))
            return(omgn-Q_n(c(w_j,a)))
          })))
          return(piece_1+piece_2)
        })(w_i,a_i,d_i,ds_i,t_i,wt_i,st_i,a)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  # Construct results object
  res <- list(model=model, beta_n=beta_n)
  
  if (verbose) { print(paste("Check 2 (functions declared):", Sys.time())) }
  
  # Variance estimate: beta_hat
  if (se_beta) {
    
    res$est_beta <- beta_n
    
    if (verbose) { print(paste("Check 3a (var est START: beta):", Sys.time())) }
    res$var_est_beta <- (1/N^2) * Reduce("+", pblapply(c(1:N), function(i) {
      (lstar_tilde(
        z_i = as.numeric(c(dat_orig$w[i,], dat_orig$a[i])),
        d_i = dat_orig$delta[i],
        ds_i = dat_orig$delta_star[i],
        t_i = dat_orig$y_star[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl))
    if (verbose) { print(paste("Check 4a (var est END: beta):", Sys.time())) }
    
  }
  
  # Variance estimate: Breslow estimator
  if (se_bshz) {
    
    res$est_bshz <- Lambda_n(t)
    
    if (verbose) { print(paste("Check 3b (var est START: bshz):", Sys.time())) }
    res$var_est_bshz <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (infl_fn_Bres(
        z_i = as.numeric(c(dat_orig$w[i,], dat_orig$a[i])),
        d_i = dat_orig$delta[i],
        ds_i = dat_orig$delta_star[i],
        t_i = dat_orig$y_star[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i]
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4b (var est END: bshz):", Sys.time())) }
    
  }
  
  # Variance estimate: survival at a point
  if (se_surv) {
    
    z_0 <- c(0.3,1,0.5)
    
    res$est_surv <- Q_n(z_0)
    
    if (verbose) { print(paste("Check 3c (var est START: surv):", Sys.time())) }
    res$var_est_surv <- (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
      (omega_n(
        z_i = as.numeric(c(dat_orig$w[i,], dat_orig$a[i])),
        d_i = dat_orig$delta[i],
        ds_i = dat_orig$delta_star[i],
        t_i = dat_orig$y_star[i],
        wt_i = dat_orig$weights[i],
        st_i = dat_orig$strata[i],
        z = z_0
      ))^2
    }, cl=cl)))
    if (verbose) { print(paste("Check 4c (var est END: surv):", Sys.time())) }
    
  }
  
  # Variance estimate: marginalized survival
  if (se_marg) {
    
    # !!!! basehaz vs. Lambda_hat
    bh <- basehaz(model, centered=FALSE)
    index <- max(which((bh$time<C$t_0)==T))
    est_bshz <- bh$hazard[index]
    N <- sum(dat$weights)
    res$est_marg <- unlist(lapply(points, function(a) {
      (1/N) * sum(unlist(lapply(c(1:N), function(i) {
        exp(-1*exp(sum(beta_n*c(as.numeric(dat_orig$w[i,]),a)))*est_bshz)
      })))
    }))
    
    # !!!!! Pre-calculate omega_n values
    
    if (verbose) { print(paste("Check 3d (var est START: marg):", Sys.time())) }
    res$var_est_marg <- unlist(lapply(points, function(a) {
      if (verbose) { print(paste0("Check 4d (point=",a,"): ", Sys.time())) }
      (1/N^2) * sum(unlist(pblapply(c(1:N), function(i) {
        (infl_fn_marg(
          w_i = as.numeric(dat_orig$w[i,]),
          a_i = dat_orig$a[i],
          d_i = dat_orig$delta[i],
          ds_i = dat_orig$delta_star[i],
          t_i = dat_orig$y_star[i],
          wt_i = dat_orig$weights[i],
          st_i = dat_orig$strata[i],
          a = a
        ))^2
      }, cl=cl)))
    }))
    if (verbose) { print(paste("Check 5d (var est END: marg):", Sys.time())) }
    
  }
  
  if (return_extras) {
    res$S_0n <- S_0n
    res$S_1n <- S_1n
    res$S_2n <- S_2n
    res$I_tilde <- I_tilde
    res$I_tilde_inv <- I_tilde_inv
    res$l_tilde <- l_tilde
    res$omega_n <- omega_n
    res$Lambda_n <- Lambda_n
  }
  
  return(res)
  
}
