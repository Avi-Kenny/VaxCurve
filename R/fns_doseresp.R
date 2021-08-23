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



#' Helper function to create hash table keys
#' 
#' @param ... A vector of function arguments
#' @return Hashed environment key
z <- function(...) {
  paste(c(...), collapse=";")
}



#' Helper function to allow for vectorized access of hash table
#' 
#' @param fn The hash table (function) to access
#' @param ... Function arguments, to construct hash table keys
#' @return Function values
# v <- function(fn, ...) {
#   apply(
#     X = cbind(...),
#     MARGIN = 1,
#     FUN = function(x) {
#       key <- paste(c(x), collapse=";")
#       val <- fn[[key]]
#       return(val)
#     }
#   )
# }
v <- memoise(function(fn, ...) {
  do.call("mapply", c(
    FUN = function(...) {
      key <- paste(c(...), collapse=";")
      return(get(fn, envir=parent.frame(n=3))[[key]])
    },
    list(...)
  ))
})



#' Helper function to create hash table from a function
#' 
#' @param fn The function to convert to a hash table
#' @param vals Dataframe of values to run function on
#' @return Hash table
create_htab <- function(fn, vals) {
  
  # Create, populate, and return hash table
  htab <- new.env()
  for (i in 1:nrow(vals)) {
    row <- vals[i,]
    key <- paste(row, collapse=";")
    htab[[key]] <- do.call(fn, as.list(as.numeric(row)))
    # v <- do.call(fn, as.list(as.numeric(row)))
    # if (is.null(v)) {
    #   stop(paste0("key: ",key))
    # }
    # htab[[key]] <- v
  }
  return (htab)
  
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
    
    pi <- function(w1,w2) { expit(w1+w2) }
    return(delta_star + (1-delta_star)*pi(w1,w2))
    
  }
  
}



#' Return IP weights
#' 
#' @param dat Dataset returned by generate_data()
#' @param scale One of c("none", "sum 1", "mean 1")
#' @return A sum-to-one vector of weights
#' @notes
#'   - Pi is accessed globally
wts <- function(dat, scale) {
  
  weights <- 1 / Pi(attr(dat,"sampling"), dat$delta_star, dat$w1, dat$w2)
  
  if (scale=="sum 1") {
    weights <- weights / sum(weights)
  }
  
  if (scale=="mean 1") {
    weights <- weights / mean(weights)
  }
  
  return(weights)
  
}



#' Return IP weight stabilization factor
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return IP weight stabilization factor "s"
#' @notes
#'   - wts is accessed globally (and therefore so is Pi)
stab <- function(dat_orig) {
  
  if (attr(dat_orig,"sampling")=="iid") {
    return(1)
  } else if (attr(dat_orig,"sampling")=="two-phase") {
    n_orig <- nrow(dat_orig)
    dat_orig %<>% filter(!is.na(a))
    weights <- wts(dat_orig, scale="none")
    return(sum(weights)/n_orig)
  }
  
}



#' Construct conditional survival estimator S_n
#' 
#' @param dat Dataset returned by generate_data(); accepts either full data or
#'     truncated data
#' @param vals Dataframe of values to run function on
#' @param type Currently only "Cox PH" is implemented
#' @param csf Logical; if TRUE, estimate the conditional survival
#'     function of the censoring distribution instead
#' @return Conditional density estimator function
construct_S_n <- function(dat, vals, type, csf=FALSE) {
  
  if (type=="Cox PH") {
    
    weights <- wts(dat, scale="mean 1")
    
    if (csf) dat$delta_star <- 1-dat$delta_star
    
    # Fit Cox model
    model <- coxph(
      Surv(y_star, delta_star)~w1+w2+a,
      data = dat,
      weights = weights
    )
    c_1 <- model$coefficients[[1]]
    c_2 <- model$coefficients[[2]]
    c_3 <- model$coefficients[[3]]
    
    # Get cumulative hazard estimate
    bh <- basehaz(model, centered=FALSE)
    
    # Pre-calculate H_0 vector
    H_0 <- c()
    for (t in 0:C$t_e) {
      index <- which.min(abs(bh$time-t))
      H_0[t+1] <- bh$hazard[index]
    }
    
    # return(memoise(Vectorize(function(t, w1, w2, a){
    #   return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
    # })))
    fn <- function(t, w1, w2, a){
      return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
    }
    
    # Run function on vals and return hash table
    return (create_htab(fn, vals))
    
  }
  
}



#' Construct g-computation estimator function of theta_0
#' 
#' @param dat_orig Dataset returned by generate_data(); must be full data
#' @param vals Dataframe of values to run function on
#' @param mu_n A regression estimator returned by construct_mu_n()
#' @return G-computation estimator of theta_0
construct_gcomp_n <- function(dat_orig, vals, S_n) {
  
  # Declare function
  fn <- function(a) {
    1 - mean(v("S_n",C$t_e,round(dat_orig$w1,1),dat_orig$w2,round(a,1)))
  }
  
  # Run function on vals and return hash table
  return (create_htab(fn, vals))
  
}



#' Construct derivative estimator theta'_n
#' 
#' @param gcomp_n G-comp estimator of theta_0 returned by construct_gcomp_n()
construct_deriv_theta_n <- function(gcomp_n) {
  
  deriv_theta_n <- function(a) {
    
    # Set derivative appx x-coordinates
    width <- 0.05
    p1 <- a - width/2
    p2 <- a + width/2
    if (p1<0) {
      p2 <- p2 - p1
      p1 <- 0
    }
    if (p2>1) {
      p1 <- p1 - p2 + 1
      p2 <- 1
    }
    
    return( (gcomp_n[[z(round(p2,3))]]-gcomp_n[[z(round(p1,3))]])/width )
    
  }
  
  return(memoise(Vectorize(deriv_theta_n)))

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
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param type Type of regression; one of c("linear", "cubic")
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return gamma_n nuisance estimator function
construct_gamma_n <- function(dat_orig, type, omega_n, f_aIw_n) {
  
  # Construct weights
  s <- stab(dat_orig)
  dat_orig$weights <- wts(dat_orig, scale="none")
  
  # Construct pseudo-outcomes
  dat_orig %<>% mutate(
    po = ifelse(
      is.na(a),
      0,
      ( (omega_n(w1,w2,y_star,delta_star,a))*(weights/s)) / f_aIw_n(a,w1,w2) )^2
  )
  
  # Run regression
  if (type=="linear") {
    model <- lm(po~a, data=filter(dat_orig,!is.na(a)), weights=weights)
    coeff <- as.numeric(model$coefficients)
    
    return(memoise(Vectorize(function(x){
      coeff[1] + coeff[2]*x
    })))
  }
  
  # Run regression
  if (type=="cubic") {
    model <- lm(po~a+I(a^2)+I(a^3), data=filter(dat_orig,!is.na(a)), weights=weights)
    coeff <- as.numeric(model$coefficients)
    
    return(memoise(Vectorize(function(x){
      coeff[1] + coeff[2]*x + coeff[3]*(x^2) + coeff[4]*(x^3)
    })))
  }
  
}



#' Construct estimator of marginal density f_A
#' 
#' @param dat_orig Dataset returned by generate_data(); must be full data
#' @param f_aIw_n A conditional density estimator returned by
#'     construct_f_aIw_n()
#' @return Marginal density estimator function
construct_f_a_n <- function(dat_orig, f_aIw_n) {
  
  return(memoise(Vectorize(function(a) {
    mean(f_aIw_n(a,dat_orig$w1,dat_orig$w2))
  })))

}



#' Construct estimator of conditional density f_{A|W}
#' 
#' @param dat Dataset returned by generate_data(); accepts either full data or
#'     truncated data
#' @param type One of c("parametric", "binning")
#' @return Conditional density estimator function
#' @notes
#'   - Assumes support of A is [0,1]
construct_f_aIw_n <- function(dat, type) {
  
  dat_trunc <- filter(dat, !is.na(a))
  n_trunc <- nrow(dat_trunc)
  weights <- wts(dat_trunc, scale="none")
  
  if (type=="parametric") {
    
    # Set up weighted likelihood
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        shape1 <- par[1] + par[2]*dat_trunc$w1[i]
        shape2 <- par[3] + par[4]*dat_trunc$w2[i]
        loglik <- dbeta(dat_trunc$a[i], shape1=shape1, shape2=shape2, log=TRUE)
        return(loglik*weights[i])
      }))
      
      return(-1*sum_loglik)
      
    }
    
    # Run optimizer
    opt <- optim(par=c(a1=0.5, a2=0.1, a3=0.5, a4=0.1), fn=wlik)
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: optim() did not converge")
    }
    
    f_aIw_n <- function(a, w1, w2){
      shape1 <- opt$par[1] + opt$par[2]*w1
      shape2 <- opt$par[3] + opt$par[4]*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    }
    
  }
  
  if (type=="binning") {
    
    # k is fixed for now; later choose via cross-validation
    k <- 10
    alphas <- seq(0, 1, length.out=k+1)
    
    # Set up binning density (based on Diaz and Van Der Laan 2011)
    # par[1] through par[k-1] are the hazard components for the bins 1 to k-1
    # par[k] and par[k+1] are the coefficients for W1 and W2
    dens <- Vectorize(function(a, w1, w2, par) {
      bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
      hz <- sapply(c(1:(ifelse(bin==k,k-1,bin))), function(j) {
        expit(par[j] + par[k]*w1 + par[k+1]*w2)
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      dens <- k*p1*p2
      return(dens)
    }, vectorize.args=c("a","w1","w2"))
    
    # Set up weighted likelihood
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        lik <- dens(a=dat_trunc$a[i], w1=dat_trunc$w1[i], w2=dat_trunc$w2[i],
                    par)
        return(weights[i]*log(lik))
      }))
      
      return(-1*sum_loglik)
      
    }
    
    # Run optimizer
    opt <- optim(par=rep(0,k+1), fn=wlik, method="CG")
    if (opt$convergence!=0) {
      warning("construct_f_aIw_n: optim() did not converge")
    }
    
    f_aIw_n <- function(a, w1, w2){
      
      bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
      par <- opt$par
      hz <- sapply(c(1:(k-1)), function(j) {
        expit(par[j] + par[k]*w1 + par[k+1]*w2)
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      
      return(k*p1*p2)
      
    }
    
  }
  
  return(memoise(Vectorize(f_aIw_n)))
  
}



#' Construct density ratio estimator g_n
#' 
#' @param vals Dataframe of values to run function on
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Marginal density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
construct_g_n <- function(vals, f_aIw_n, f_a_n) {
  
  # Declare function
  fn <- function(a,w1,w2) {
    f_aIw_n(a,w1,w2) / f_a_n(a)
  }
  
  # Run function on vals and return hash table
  return (create_htab(fn, vals))

}



#' Construct estimator of nuisance influence function omega_n
#' 
#' @param vals Dataframe of values to run function on
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @param Sc_n Conditional censoring survival function estimator returned by
#'     construct_S_n
#' @param m Number of partitions used for the integral approximation
#' @return Estimator function of nuisance omega_0
construct_omega_n <- function(vals, S_n, Sc_n, m=400) {
  
  # Declare function
  fn <- function(w1,w2,a,y_star,delta_star) {
    
    i <- c(1:m)
    k <- min(y_star,C$t_e)
    w1 <- round(w1,1)
    a <- round(a,1)
    
    integral <- sum(
      (v("S_n",round((i*k)/m),w1,w2,a) - v("S_n",round(((i-1)*k)/m),w1,w2,a)) *
      ((v("S_n",round((i*k)/m),w1,w2,a))^2 * v("Sc_n",round((i*k)/m),w1,w2,a))^-1
    )
    S_n[[z(C$t_e,w1,w2,a)]] * (
      (delta_star * as.integer(y_star<=C$t_e)) /
        (S_n[[z(round(k),w1,w2,a)]] * Sc_n[[z(round(k),w1,w2,a)]]) +
        integral
    )

  }
  
  # Run function on vals and return hash table
  return (create_htab(fn, vals))
  
}



#' Construct nuisance estimator eta_n
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param vals Dataframe of values to run function on
#' @param S_n Conditional survival function estimator returned by construct_S_n
#' @return Estimator function of nuisance eta_0
construct_eta_n <- function(dat_orig, vals, S_n) {
  
  # Prep
  s <- stab(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- wts(dat, scale="none")
  
  # Declare temporary memoised version of v(S_n,...)
  # v_mem <- memoise(function(w1,w2) {
  #   v(S_n,C$t_e,round(w1,1),w2,round(dat$a,1))
  # })
  
  # Declare function
  fn <- function(x,w1,w2) {
    sum(
      weights * as.integer(dat$a<=x) *
        (1-v("S_n",C$t_e,round(w1,1),w2,round(dat$a,1)))
        # (1-v_mem(w1,w2))
        # (1-S_n(C$t_e,round(w1,1),w2,round(dat$a,1)))
    ) /
      (n_orig*s)
  }
  
  # Run function on vals and return hash table
  return (create_htab(fn, vals))
  
}



#' lambda estimator
#' 
#' @param k Power k
#' @param G Transformation function G; usually returned by a function
#'     constructed by construct_Phi_n()
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @return Value of lambda
lambda <- function(k, G, dat_orig) {
  
  n_orig <- nrow(dat_orig)
  s <- stab(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- wts(dat, scale="none")
  lambda <- (1/n_orig) * sum(
    (weights/s) * (G(dat$a))^k
  )
  return(lambda)
  
}



#' Construct Gamma_n primitive one-step estimator
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param omega_n A nuisance influence function returned by construct_omega_n()
#' @param S_n A conditional survival function returned by construct_S_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @return Gamma_n estimator
#' @notes This is a generalization of the one-step estimator from Westling &
#'     Carone 2020
construct_Gamma_n <- function(dat_orig, vals, omega_n, S_n, g_n) {
  
  # Prep
  dat <- dat_orig
  s <- stab(dat_orig)
  n_orig <- nrow(dat)
  dat %<>% filter(!is.na(a))
  n <- nrow(dat)
  weights <- wts(dat, scale="none")
  i_long <- rep(c(1:n), each=n)
  j_long <- rep(c(1:n), times=n)
  a_i_long <- dat$a[i_long]
  w1_i_long <- dat$w1[i_long]
  w1_j_long <- dat$w1[j_long]
  w2_i_long <- dat$w2[i_long]
  w2_j_long <- dat$w2[j_long]
  delta_star_i_long <- dat$delta_star[i_long]
  delta_star_j_long <- dat$delta_star[j_long]
  sampling <- attr(dat,"sampling")
  subpiece_1a <- ( 1 + (
    v("omega_n",dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star) /
    v("g_n",dat$a,dat$w1,dat$w2)
  ) ) * ( weights / s )
  subpiece_2a <- v("S_n",C$t_e,round(w1_j_long,1),w2_j_long,round(a_i_long,1)) /
  (
    s^2 * Pi(sampling,delta_star_i_long,w1_i_long,w2_i_long) *
    Pi(sampling,delta_star_j_long,w1_j_long,w2_j_long)
  )
  
  # Declare function
  fn <- function(x) {
    
    subpiece_1b <- as.integer(dat$a<=x)
    piece_1 <- sum(subpiece_1a*subpiece_1b) * (1/n_orig)
    
    subpiece_2b <- as.integer(a_i_long<=x)
    piece_2 <- sum(subpiece_2a*subpiece_2b) * (1/(n_orig^2))
    
    return(piece_1-piece_2)
    
  }
  
  # Run function on vals and return hash table
  return (create_htab(fn, vals))
  
}



#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat_orig Dataset returned by generate_data(); this must be the FULL
#'     DATA, including the "missing" observations
#' @param type One of c("ecdf", "inverse").
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - Adaptation of stats::ecdf() source code
#'   - This accesses wts() globally
construct_Phi_n <- function (dat_orig, type="ecdf") {
  
  dat <- dat_orig
  s <- stab(dat_orig)
  dat <- cbind(dat, wts=wts(dat, scale="none"))
  dat %<>% arrange(a)
  n_orig <- nrow(dat)
  dat %<>% filter(!is.na(a))
  vals_x <- unique(dat$a)
  
  vals_y <- c()
  for (j in 1:length(vals_x)) {
    indices <- which(dat$a==vals_x[j])
    weights_j <- dat$wts[indices]
    new_y_val <- sum(weights_j) / (n_orig*s)
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



#' Hypothesis test based on logistic regression
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_wald <- function(dat, alt_type="incr", params) {
  
  model <- glm(y~w1+w2+a, data=dat, family="binomial")
  one_sided_p <- pnorm(summary(model)$coefficients["a",3], lower.tail=F)
  reject <- as.integer(one_sided_p<0.05)
  
  return (reject)
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_rho_n <- function(dat_orig, Phi_n) {
  
  s <- stab(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- wts(dat, scale="none")
  a_j <- dat$a
  
  return(memoise(Vectorize(function(a) {
    
    return((1/n_orig) * sum(
      (weights_j/s) * (Phi_n(a_j)^3) * (as.integer(a<=a_j) - Phi_n(a_j))
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
  
  s <- stab(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- wts(dat, scale="none")
  a_j <- dat$a
  sampling <- attr(dat,"sampling")
  
  return(memoise(Vectorize(function(w1,w2,y_star,delta_star,delta,a) {
    
    piece_1 <- (1/n_orig) * sum(
      # (weights_j/s) * (xi_n(a,a_j)-rho_n(a)) * Gamma_n(a_j)
      (weights_j/s) * (xi_n(a,a_j)) * v("Gamma_n",a_j)
    )
    
    piece_2 <- (lambda_2*(Phi_n(a)^2) - lambda_3*Phi_n(a)) * Gamma_n[[z(a)]]
    
    return(
      (delta/(s*Pi(sampling, delta_star,w1,w2))) * (piece_1+piece_2)
    )
    
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_Gamma <- function(dat_orig, omega_n, g_n, gcomp_n, eta_n,
                                    Gamma_n) {
  
  s <- stab(dat_orig)
  sampling <- attr(dat_orig,"sampling")
  
  return(memoise(Vectorize(function(x,w1,w2,y_star,delta_star,delta,a) {
    (delta/(s*Pi(sampling, delta_star,w1,w2))) * (
      as.integer(a<=x)*(
        (omega_n[[z(w1,w2,a,y_star,delta_star)]]/g_n[[z(a,w1,w2)]]) +
        gcomp_n[[z(round(a,3))]]
      ) +
        eta_n[[z(round(x,2),round(w1,2),w2)]] - 
        2*Gamma_n[[z(x)]]
    )
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
construct_infl_fn_2 <- function(dat_orig, Phi_n, infl_fn_Gamma, lambda_2,
                                lambda_3) {
  
  s <- stab(dat_orig)
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights_j <- wts(dat, scale="none")
  a_j <- dat$a
  
  return(memoise(Vectorize(function(w1,w2,y_star,delta_star,delta,a) {
    (1/n_orig) * sum((weights_j/s) * (
      ( lambda_2*(Phi_n(a_j)^2) - lambda_3*Phi_n(a_j) ) *
        infl_fn_Gamma(a_j,w1,w2,y_star,delta_star,delta,a)
    ))
  })))
  
}



#' !!!!! document
#' 
#' @param x x
#' @return x
beta_n_var_hat <- function(dat_orig, infl_fn_1, infl_fn_2) {
  
  n_orig <- nrow(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  
  b_sum <- 0
  for (i in c(1:nrow(dat))) {
    s <- (
      infl_fn_1(dat$w1[i], dat$w2[i], dat$y_star[i], dat$delta_star[i],
                dat$delta[i], dat$a[i]) +
      infl_fn_2(dat$w1[i], dat$w2[i], dat$y_star[i], dat$delta_star[i],
                dat$delta[i], dat$a[i])
    )^2
    b_sum <- b_sum + s
  }
  
  return(
    (1/n_orig) * sum(b_sum)
  )
  
}
