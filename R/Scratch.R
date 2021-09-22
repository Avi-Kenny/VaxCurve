
# Histogram for edge values
if (F) {
  
  r_cox <- filter(sim$results, surv_true=="Cox PH")
  r_com <- filter(sim$results, surv_true=="complex")
  cox <- r_cox$est_0.0
  com <- r_com$est_0.0
  tr_cox <- r_cox[1,"theta_0.0"]
  tr_com <- r_com[1,"theta_0.0"]
  
  (mean(cox)-tr_cox)/mean(cox)
  (mean(com)-tr_com)/mean(com)
  
  df <- data.frame(
    x = c(cox,com),
    grp = c(rep("Cox PH", length(cox)),
            rep("complex", length(com)))
  )
  
  ggplot(df, aes(x=x, group=grp, fill=factor(grp))) +
    geom_histogram(color="white") +
    facet_wrap(~grp, ncol=2) +
    geom_vline(
      aes(xintercept=tr),
      data = data.frame(tr=c(tr_cox,tr_com), grp=c("Cox PH","complex"))
    )
  
}

# Monotone spline for derivative
if (F) {
  
  x <- sort(runif(5))
  y <- sort(runif(5))
  
  theta_n_smoothed <- splinefun(x=x, y=y, method="monoH.FC")
  
  ggplot(
    data.frame(x=seq(0,1,0.01), y=theta_n_smoothed(seq(0,1,0.01))),
    aes(x=x,y=y)) +
    geom_line() +
    labs(title=mmm) +
    geom_point(data=data.frame(x=x, y=y))
  
}

# Use softer cutoffs instead of step functions
if (F) {
  
  ff <- function(w1) { as.integer(abs(w1-0.5)<0.11) }
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=ff)
  
  ff2 <- function(w1) { pmax(0,1-4*abs(w1-0.5)) }
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=ff2)
  
}

#
if (F) {
  
  # Generate data
  dat_orig <- generate_data(
    n = 2000,
    alpha_3 = 0.7,
    distr_A = "Beta(1.5+w1,1.5+w2)", # "Unif(0,1)" "Beta(1.5+w1,1.5+w2)"
    edge = "none",
    surv_true = "complex",
    sampling = "two-phase"
  )
  
  # Prep
  n_orig <- nrow(dat_orig)
  dat_orig$weights <- wts(dat_orig)
  dat <- dat_orig %>% filter(!is.na(a))
  weights <- dat$weights
  
  # Construct dataframes of values to pre-compute functions on
  vlist <- create_val_list(dat_orig, C$appx)
  
  # Construct component functions
  S_n <- construct_S_n(dat_orig, vlist$S_n, type="true")
  Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="true",
                        csf=TRUE)
  f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                               type="true", k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  
  offset_25 <- (dat$weights *
     as.integer(dat$a<=0.25) *
     omega_n(dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star)
  ) / g_n(dat$a,dat$w1,dat$w2)
  offset_75 <- (dat$weights *
                  as.integer(dat$a<=0.75) *
                  omega_n(dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star)
  ) / g_n(dat$a,dat$w1,dat$w2)
  print(mean(offset_25))
  print(mean(offset_75))
  # print(head(sort(offset)))
  # print(tail(sort(offset)))
  # print(mean(sort(offset)[2:length(offset)]))
    
}

# Variance of one-step edge estimator
if (F) {
  
  ##### Edge mass #####
  
  n_reps <- 30
  edge <- "expit"
  ests <- c()
  sigma2s <- c()
  
  for (i in 1:n_reps) {
    
    # Generate data
    dat_orig <- generate_data(
      n = 1000,
      alpha_3 = 0.7,
      distr_A = "Beta(1.5+w1,1.5+w2)",
      edge = edge,
      surv_true = "Cox PH",
      sampling = "two-phase"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig, C$appx)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    pi_n <- construct_pi_n(dat_orig, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat_orig, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat_orig, pi_n, S_n, omega_n,
                                   theta_os_n_est)
    
    ests <- c(ests, theta_os_n_est)
    sigma2s <- c(sigma2s, sigma2_os_n_est/n_orig)
    
  }
  
  v1 <- var(ests)
  m1 <- mean(sigma2s)
  
  ##### Edge mass #####
  
  n_reps <- 30
  edge <- "none"
  ests <- c()
  sigma2s <- c()
  
  for (i in 1:n_reps) {
    
    # Generate data
    dat_orig <- generate_data(
      n = 1000,
      alpha_3 = 0.7,
      distr_A = "Beta(1.5+w1,1.5+w2)",
      edge = edge,
      surv_true = "Cox PH",
      sampling = "two-phase"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig, C$appx)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    pi_n <- construct_pi_n(dat_orig, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat_orig, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat_orig, pi_n, S_n, omega_n,
                                   theta_os_n_est)
    
    ests <- c(ests, theta_os_n_est)
    sigma2s <- c(sigma2s, sigma2_os_n_est/n_orig)
    
  }
  
  v2 <- var(ests)
  m2 <- mean(sigma2s)
  
  print("Edge mass")
  print(v1)
  print(m1)
  print("No edge mass")
  print(v2)
  print(m2)

}

# Testing cross-fitting
if (F) {
  
  system.time({
    ests_reg <- est_curve(
      dat_orig = dat_tp,
      estimator = "Grenander",
      params = list(S_n_type=params$S_n_type, g_n_type=params$g_n_type, ci_type="logit",
                    cf_folds=1),
      points = C$points
    )
  })
  
  system.time({
    ests_cf <- est_curve(
      dat_orig = dat_tp,
      estimator = "Grenander",
      params = list(S_n_type=params$S_n_type, g_n_type=params$g_n_type, ci_type="logit",
                    cf_folds=3),
      points = C$points
    )
  })
  
  system.time({
    reject_reg <- test_2(
      dat_orig = dat_tp,
      alt_type = "incr",
      params = list(
        var = "asymptotic",
        S_n_type = params$S_n_type,
        g_n_type = params$g_n_type,
        est_known_nuis = FALSE,
        cf_folds = 1
      )
    )
  })
  
  system.time({
    reject_cf <- test_2(
      dat_orig = dat_tp,
      alt_type = "incr",
      params = list(
        var = "asymptotic",
        S_n_type = params$S_n_type,
        g_n_type = params$g_n_type,
        est_known_nuis = FALSE,
        cf_folds = 3
      )
    )
  })
  
  
}

# Does memoising work if function is passed?
if (F) {
  
  f <- function() {
    print(get("hey", envir=parent.frame()))
  }
  
  (function(hey) {
    f()
  })("there")
  
  fn_mem <- memoise(function(x) {
    Sys.sleep(2)
    return(9)
  })
  fn_mem(1)
  
  fn2 <- function(fn) {
    fn(1)
  }
  fn3 <- function(fn) {
    fn(1)
  }
  fn2(fn_mem)
  fn3(fn_mem)
  
}

# Fastest way of getting a vectorized/memoised function (part 2)
if (F) {
  
  construct_S_n_mem <- function(dat, vals) {
      model <- coxph(Surv(y_star, delta_star)~w1+w2+a, data=dat)
      c_1 <- model$coefficients[[1]]
      c_2 <- model$coefficients[[2]]
      c_3 <- model$coefficients[[3]]
      bh <- basehaz(model, centered=FALSE)
      H_0 <- c()
      for (t in 0:C$t_e) {
        index <- which.min(abs(bh$time-t))
        H_0[t+1] <- bh$hazard[index]
      }
      # return(Vectorize(memoise(function(t, w1, w2, a){
      #   return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
      # })))
      return(memoise(Vectorize(function(t, w1, w2, a){
        return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
      })))
  }
  
  construct_S_n_htab <- function(dat, vals) {
    model <- coxph(Surv(y_star, delta_star)~w1+w2+a, data=dat)
    c_1 <- model$coefficients[[1]]
    c_2 <- model$coefficients[[2]]
    c_3 <- model$coefficients[[3]]
    bh <- basehaz(model, centered=FALSE)
    H_0 <- c()
    for (t in 0:C$t_e) {
      index <- which.min(abs(bh$time-t))
      H_0[t+1] <- bh$hazard[index]
    }
    fn <- function(t, w1, w2, a){
      return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
    }
    return (create_htab(fn, vals))
  }
  
  vals_S_n <- expand.grid(t=c(0:C$t_e),w1=seq(0,1,0.1),w2=c(0,1),a=seq(0,1,0.1))
  S_n_mem <- construct_S_n_mem(dat, vals_S_n)
  for (i in 1:nrow(vals_S_n)) {
    row <- vals_S_n[i,]
    x <- do.call(S_n_mem, as.list(as.numeric(row)))
  }
  S_n_htab <- construct_S_n_htab(dat, vals_S_n)
  
  v1 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        # key <- "200;0;0;0.1"
        key <- paste(x, collapse=";")
        val <- fn[[key]]
        return(val)
      }
    )
  }
  v2 <- function(fn, ...) {
    do.call("mapply", c(
      FUN = function(...) {
        # key <- "200;0;0;0.1"
        key <- paste(c(...), collapse=";")
        return(fn[[key]])
      },
      list(...)
    ))
  }
  v3 <- memoise(function(...) {
    do.call("mapply", c(
      FUN = function(...) {
        key <- paste(c(...), collapse=";")
        return(S_n_htab[[key]])
      },
      list(...)
    ))
  })
  
  S_n_mem(C$t_e,0,0,round(dat$a,1))[1:10]
  v1(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  v2(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  v3(C$t_e,0,0,round(dat$a,1))[1:10]
  v4("S_n_htab",C$t_e,0,0,round(dat$a,1))[1:10]
  
  microbenchmark({
    # S_n_mem(C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    S_n_mem(C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v1(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v1(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v2(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v3(C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v4("S_n_htab",C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  
}

# Fastest way of getting a vectorized/memoised function (part 1)
if (F) {
  
  n <- 4*10^4
  fn <- function(x,y) { x^2+y }
  
  # Method 1: memoise/Vectorize
  # mean: 838 microsec
  fn_1 <- Vectorize(memoise(fn))
  for(i in 1:n) {
    x <- fn_1(i,10)
  }
  microbenchmark({
    fn_1(c(100:200),10)
  }, times=100L)
  
  # Method 2: hash table
  # mean: 205 microsec
  htab <- new.env()
  for(i in 1:n) {
    htab[[paste(c(i,10), collapse=";")]] <- fn_1(i,10)
  }
  fn_2 <- function(...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        key <- paste(c(x), collapse=";")
        val <- htab[[key]]
        return(val)
      }
    )
  }
  microbenchmark({
    fn_2(c(100:200),10)
  }, times=100L)
  
  # Create, populate, and return hash table
  for (i in 1:n) {
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

# !!!! TEMP
if (F) {
  
  # fn_mem <- memoise(Vectorize(function(x) {x^2}))
  # htab <- new.env()
  # 
  # for(i in 1:10000) {
  #   x <- fn_mem(i)
  #   row <- vals[i,]
  #   key <- paste(row, collapse=";")
  #   htab[[key]] <- do.call(fn, as.list(as.numeric(row)))
  # }
  
  S_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type)
  
  
  
  
  S_n_htab_vc <- Vectorize(function(w,x,y,z) {
    S_n_htab[["200;0;0;0.1"]]
    # S_n_htab[[paste(c(w,x,y,z), collapse=";")]]
  })
  
  
  microbenchmark({
    S_n_htab_vc(200,0,0,rep(0.1,100))
    # S_n_htab[[paste0("200;0;0;",0.1)]]
    # v(S_n_htab,C$t_e,0,0,round(dat$a,1))
  }, times=100L)
  
  microbenchmark({
    S_n_mem(200,0,0,rep(0.1,100))
    # S_n_mem(C$t_e,0,0,round(dat$a,1))
  }, times=100L)
  
  (function(a,b,c) {
    mc <- match.call()
    called_args <- as.list(mc)[-1]
    print(called_args)
  })(1,2,3)

  function (t, w1, w2, a) 
  {
    mc <- match.call()
    encl <- parent.env(environment())
    called_args <- as.list(mc)[-1]
    default_args <- encl$`_default_args`
    default_args <- default_args[setdiff(names(default_args), 
                                         names(called_args))]
    called_args[encl$`_omit_args`] <- NULL
    args <- c(lapply(called_args, eval, parent.frame()), lapply(default_args, 
                                                                eval, envir = environment()))
    key <- encl$`_hash`(c(encl$`_f_hash`, args, lapply(encl$`_additional`, 
                                                       function(x) eval(x[[2L]], environment(x)))))
    res <- encl$`_cache`$get(key)
    if (inherits(res, "key_missing")) {
      mc[[1L]] <- encl$`_f`
      res <- withVisible(eval(mc, parent.frame()))
      encl$`_cache`$set(key, res)
    }
    if (res$visible) {
      res$value
    }
    else {
      invisible(res$value)
    }
  }
  
  
  
  
  
  library(stringi)
  
  v1 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[paste(c(x), collapse=";")]]
      }
    )
  }
  
  v2 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[paste(x, collapse=";")]]
      }
    )
  }
  
  v3 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[stri_join(x,collapse=";")]]
      }
    )
  }
  
  microbenchmark({
    v1(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
  microbenchmark({
    v2(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
  microbenchmark({
    v3(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
}

# New hashing/memoising structure for creators
if (F) {
  
  # Creator option #1
  creator_1 <- function(vals) {
    
    # Declare function
    fn <- function(x,y) { x*y }
    
    # Set environment
    htab <- new.env()
    
    # Run function on vals
    for (i in 1:nrow(vals)) {
      row <- vals[i,]
      key <- paste(row, collapse=";")
      htab[[key]] <- do.call("fn", as.list(as.numeric(row)))
    }
    
    return (htab)
    
  }
  
  # Usage
  vals <- expand.grid(a=c(1:n), w=c(1:n))
  fn_1 <- creator_1(vals)
  z <- function(...) {
    paste(c(...), collapse=";")
  }
  v <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {fn[[z(x)]]}
    )
  }
  
  fn_1[[z(1,2)]]
  fn_1[[z(4,5)]]
  fn_1[[z(7,8)]]
  v(fn_1,c(1,4,7),c(2,5,8))

  fn_1[[z(2,2)]]
  fn_1[[z(2,5)]]
  fn_1[[z(2,8)]]
  v(fn_1,3,c(2,5,8))
  
  n <- 10000
  vals <- data.frame(a=c(1:n), w=c(2:(n+1)))
  
  microbenchmark({
    fn_1[["77;98"]]
  }, times=1000L)
  microbenchmark({
    fn_1[[paste(c(77,98), collapse=";")]]
  }, times=1000L)
  microbenchmark({
    fn_1[[z(77,98)]]
  }, times=1000L)

  # memoise() vs. new.env()
  mem_fn <- memoise(fn)
  htab_old <- new.env()
  
  # Initialize hash structures
  # n= 200:  0.4 sec
  # n= 400:  1.8 sec
  # n= 800: 10.4 sec
  # n=1600: 10.4 sec
  n <- 1600
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        htab_old[[paste0(i,";",j)]] <- fn(i,j)
      }
    }
  })
  
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        x <- mem_fn(i,j)
      }
    }
  })
  
  # Run benchmarks
  microbenchmark({x <- mem_fn(5,6)}, times=1000L)
  microbenchmark({x <- htab[[paste0(5,";",6)]]}, times=1000L)
  
}

# Checking lambda factors with mixture distribution
if (F) {
  
  n <- 10^6
  prob <- 0.9 # 0.9
  mix <- rbinom(n, size=1, prob=prob)
  unif <- runif(n)
  a <- mix*unif
  F_a <- ecdf(a)
  mean((F_a(a))^2)
  mean((F_a(a))^3)
  
}

# Probability integral transform
if (F) {
  
  x <- rbeta(5,0.1,0.1)
  F_n <- ecdf(x)
  print(x)
  print(F_n(x))
  
  n <- 100
  i <- c(1:n)
  sum(i^2)
  n^3/3 + n^2/2 + n/6
  sum(i^3)
  n^4/4 + n^3/2 + n^2/4
  
}

# Benchmarking hashing/memoising structures
if (F) {
  
  library(memoise)
  library(microbenchmark)
  
  # Unmemoised function
  fn <- function(x,y) { x*y*1000 }
  
  # memoise() vs. new.env()
  mem_fn <- memoise(fn)
  htab <- new.env()
  htab_old <- new.env()
  
  # Initialize hash structures
  # n= 200:  0.4 sec
  # n= 400:  1.8 sec
  # n= 800: 10.4 sec
  # n=1600: 10.4 sec
  n <- 1600
  htab_old <- new.env()
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        htab_old[[paste0(i,";",j)]] <- fn(i,j)
      }
    }
  })
  
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        x <- mem_fn(i,j)
      }
    }
  })
  
  # Run benchmarks
  microbenchmark({x <- mem_fn(5,6)}, times=1000L)
  microbenchmark({x <- htab[[paste0(5,";",6)]]}, times=1000L)
  
}

# TEMP
if (F) {
  
  int <- (3/5)*3^5 + (2/3)*3^3
  print(int)
  
  G <- function(x) { x^3 + 2*x }
  h <- function(x) { x^2 }
  m <- 1000000
  i <- c(1:m)
  int_appx <- sum(
    (G((3*i)/m)-G((3*i-3)/m)) * (h((3*i)/m))
  )
  print(int_appx)
  
}

# TEMP
if (F) {
  
  dat <- c(0,0,0,0,1,1,1,2,2,3)+0.5
  n <- length(dat)
  I <- as.integer
  
  # Version 1
  {
    dens <- Vectorize(function(x, phi1, phi2, phi3) {
      phi4 <- 4-(phi1+phi2+phi3)
      bin <- which.min(x>c(0,1,2,3,4))-1
      lik <- phi1*I(bin==1)+phi2*I(bin==2)+phi3*I(bin==3)+phi4*I(bin==4)
      return(lik)
    }, vectorize.args=c("x"))
    
    wlik <- function(par) {
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        lik <- dens(dat[i], par[1], par[2], par[3])
        return(log(lik))
      }))
      return(-1*sum_loglik)
    }
    
    opt <- optim(par=c(1,1,1), fn=wlik, method="CG") # BFGS CG SANN
    if (opt$convergence!=0) { warning("nonconvergence of optim()") }
    dens(c(0,1,2,3)+0.5, phi1=opt$par[1], phi2=opt$par[2], phi3=opt$par[3])
    # sum(d)
  }
  
  # Version 2 (only dens)
  {
    dens <- Vectorize(function(x, par) {
      bin <- which.min(x>c(0,1,2,3,4))-1
      k <- 4 # num bins
      hz <- sapply(c(1:(ifelse(bin==k,k-1,bin))), function(j) {
        expit(par[j])
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      lik <- k*p1*p2
      return(lik)
    }, vectorize.args=c("x"))

    wlik <- function(par) {
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        lik <- dens(dat[i], par)
        return(log(lik))
      }))
      return(-1*sum_loglik)
    }
    
    opt <- optim(par=c(0.5,0.5,0.5), fn=wlik, method="CG") # BFGS CG SANN
    if (opt$convergence!=0) { warning("nonconvergence of optim()") }
    dens(c(0,1,2,3)+0.5, par=c(opt$par[1],opt$par[2],opt$par[3]))
    sum(d)
    
  }
  
}

# Testing Ted's density estimator function
if (F) {
  
  library(ctsCausal)
  
  # Sample data
  n <- 1000
  WW <- data.frame(W1 = runif(n))
  Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-WW$W1)))
  AA <- (1-Z) * rnorm(n, mean = WW$W1, sd = abs(1 + WW$W1))
  fit <- cmdSuperLearner(
    A = AA,
    W = WW,
    control = list(
      SL.library = c("SL.mean","SL.glm","SL.gam","SL.earth"),
      verbose = FALSE,
      n.bins = c(2:10)
    )
  )  
  
}


# Testing G-comp variance estimator
if (F) {
  
  # Generate data
  ests <- c()
  for (i in 1:100) {
    n <- 100
    alpha_0 <- -1.5
    alpha_1 <- 0.3
    alpha_2 <- 0.7
    alpha_3 <- 0.7
    w1 <- runif(n)
    w2 <- rbinom(n, size=1, prob=0.5)
    a <- rbeta(n, shape1=0.9, shape2=(1.1 + 0.4*w2))
    probs <- expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*a)
    y <- rbinom(n, size=1, prob=probs)
    dat <- data.frame(w1=w1, w2=w2, a=a, y=y)
    
    # Fit logistic regression model
    model <- glm(y~w1+w2+a, data=dat, family="binomial")
    coefs <- as.numeric(summary(model)$coefficients[,1])
    
    # Run gcomp
    x <- 0.5
    gcomp_i <- apply(
      X = dat,
      MARGIN = 1,
      FUN = function(r) {
        M_i <- c(1,r[["w1"]],r[["w2"]],x)
        return(expit(as.numeric(M_i %*% coefs)))
      }
    )
    
    ests <- c(ests, mean(gcomp_i))
  }
  sd(ests)
  
  # Run histogram
  (function(x){ggplot(data.frame(x=x),aes(x=x))+geom_histogram()})(ests)
  
  # Calculate variance estimate
  I_inv <- n*vcov(model)
  mu_n <- function(w1,w2,a) {
    expit(coefs[1] + coefs[2]*w1 + coefs[3]*w2 + coefs[4]*a)
  }
  exp_lin <- function(w1,w2,a) {
    exp(coefs[1] + coefs[2]*w1 + coefs[3]*w2 + coefs[4]*a)
  }
  v1 <- function(w1,w2,a) {
    t(as.matrix(c(1,w1,w2,a)))
  }
  v2 <- function(y_i,w1_i,w2_i,a_i,w1_j,w2_j,a) {
    y_mins_mu <- y_i - mu_n(w1_j,w2_j,a)
    as.matrix(c(
      y_mins_mu,
      w1_i*y_mins_mu,
      w2_i*y_mins_mu,
      a_i*y_mins_mu
    ))
  }
  
  w1 <- dat$w1
  w2 <- dat$w2
  y <- dat$y
  a <- dat$a
  var_est <- (1/(n^4)) * sum(sapply(c(1:n), function (i) {
    (sum(sapply(c(1:n), function (j) {
      
      mu_i <- mu_n(w1[i],w2[i],x)
      mu_j <- mu_n(w1[j],w2[j],x)
      mu_i - mu_j + (
        ( mu_j / (1+exp_lin(w1[j],w2[j],x)) ) *
        ( v1(w1[j],w2[j],x) %*%
          I_inv %*%
          v2(y[i],w1[i],w2[i],a[i],w1[j],w2[j],x) )
      )
      
    })))^2
  }))
  print(var_est)
  
  # var_est <- (1/(n^4)) * sum(
  #   
  #   mu_n(w1_i_long,w2_i_long,x) - mu_n(w1_j_long,w2_j_long,x) + (
  #     ( mu_n(w1_j_long,w2_j_long,x) / (1 + exp_lin(w1_j_long,w2_j_long,x)) ) *
  #     ( v1(w1_j_long,w2_j_long,x) %*%
  #       I_inv %*%
  #       v2(y_i_long,w1_i_long,w2_i_long,a_i_long,w1_j_long,w2_j_long,x) )
  #   )
  #   
  # )
  
}

# Testing conditional distribution estimators
if (F) {
  
  w <- 10
  gamma <- c(0,-3,-3)
  hz1 <- expit(gamma[1]+0.2*w)
  hz2 <- expit(gamma[2]+0.2*w)
  hz3 <- expit(gamma[3]+0.2*w)
  hz4 <- expit(gamma[4]+0.2*w)
  p1 <- hz1
  p2 <- hz2 * (1-hz1)
  p3 <- hz3 * (1-hz2)*(1-hz1)
  p4 <- (1-hz3)*(1-hz2)*(1-hz1) # Modification of Diaz and VDL
  ggplot(data.frame(
    x = c(1:4),
    y = c(p1,p2,p3,p4)
  ), aes(x=x, y=y)) + geom_bar(stat="identity") + ylim(c(0,1))
  #
  
  
  n <- 10000
  w1 <- runif(n)
  w2 <- runif(n)
  a <- rbeta(n, shape1=0.3+w1, shape2=0.3+w2)
  ggplot(data.frame(a=a), aes(x=a)) + geom_histogram(bins=50)
  
}

# Experimenting with the profiler
if (F) {
  
  y <- function() {
    rnorm(10^7)
  }
  z <- function() {
    y()
  }
  x1 <- rnorm(10^7)
  x2 <- 999
  x3 <- y()
  x4 <- 999
  x5 <- z()
  x6 <- 999
  
}

# Testing whether memoising works with constructors
if (F) {
  
  construct_test <- function() {
    
    return(memoise(Vectorize(function(a) {
      mean(mu_n(a,dat$w1,dat$w2))
    })))
    
  }
  t1 <- construct_test()
  t2 <- construct_test()
  environment(t1)
  environment(t2)
  
  construct_test2 <- function(test_func) {
    
    return(function(a) {
      print(environment(test_func))
    })
    
  }
  
  tt1 <- construct_test2(t1)
  tt2 <- construct_test2(t2)
  tt1()
  tt2()
  
}

# Adapting ECDF and inverse ECDF to two-phase sampling
if (F) {
  
  dat <- data.frame(
    a = c(0.5,0.75,0.75,NA,NA),
    wts = c(6,1,1,NA,NA)
  )
  
  # Old functions
  Phi_n <- ecdf(dat$a)
  Phi_n_inv <- function(x) {
    dat %<>% filter(!is.na(a))
    qemp(p=x, obs=dat$a, discrete=T)
  }
  
  # New functions
  construct_Phi_n2 <- function(dat) {
    dat <- cbind(dat, wts=wts(dat, scale="sum 1"))
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    s <- sum(dat$wts) / n
    Vectorize(function(x) { (1/(n*s))*sum(dat$wts*as.integer(dat$a<=x)) })
  }
  construct_Phi_n3 <- function (dat, type="ecdf") {
    # Adaptation of ecdf() source code
    dat <- cbind(dat, wts=wts(dat, scale="sum 1"))
    dat %<>% arrange(a)
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    vals_x <- unique(dat$a)
    vals_y <- c()
    s <- sum(dat$wts) / n
    for (j in 1:length(vals_x)) {
      indices <- which(dat$a==vals_x[j])
      weights_j <- dat$wts[indices]
      new_y_val <- sum(weights_j) / (n*s)
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
  Phi_n2 <- construct_Phi_n2(dat)
  Phi_n3 <- construct_Phi_n3(dat)
  Phi_n_inv3 <- construct_Phi_n3(dat, type="inverse")
  
  
  # # Run benchmarks
  # microbenchmark(Phi_n(0.5), times=100)
  # microbenchmark(Phi_n2(0.5), times=100)
  # microbenchmark(Phi_n3(0.5), times=100)
  
  # ECDF
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n2)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n3)
  
  # Inverse ECDF
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n_inv)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n_inv3)
  
  # Correct behavior at knots
  Phi_n(0.5)
  Phi_n2(0.5)
  Phi_n3(0.5)
  Phi_n_inv(1/3)
  Phi_n_inv3(1/3)
  
}

# Beta density estimator (weighted)
if (F) {
  
  library(kdensity)
  library(ks)
  
  n <- 500
  shape <- c(0.6,0.6) # c(0.3,0.3)
  sample <- rbeta(n, shape1=shape[1], shape2=shape[2]) %>% sort()
  grid <- seq(0,1,0.01)
  
  wts <- c(
    # rep(1,length(sample))
    rep(0.5,length(sample)/2),
    rep(1.5,length(sample)/2)
  )
  wts <- wts / sum(wts)
  kd_gauss2_ests <- density(
    x = sample,
    kernel = "gaussian",
    from = 0,
    to = 1,
    weights = wts
  )
  kd_gauss2 <- function(x) {
    index <- which.min(abs(kd_gauss2_ests$x - x))
    area <- mean(kd_gauss2_ests$y)
    return(kd_gauss2_ests$y[index]/area)
  }
  kd_gauss <- kdensity(
    x = sample,
    start = "gumbel",
    kernel = "gaussian",
    support = c(0,1)
  )
  # kd_beta <- kdensity( # Weights unavailable
  #   x = sample,
  #   start = "gumbel",
  #   kernel = "beta",
  #   support = c(0,1)
  # )
  kd_beta2_ests <- kde.boundary(
    x = sample,
    boundary.kernel = "beta",
    w = wts,
    xmin = 0,
    xmax = 1
  )
  kd_beta2 <- function(x) {
    len <- length(kd_beta2_ests$eval.points)
    k_x <- kd_beta2_ests$eval.points[2:(len-1)]
    k_dens <- kd_beta2_ests$estimate[2:(len-1)]
    index <- which.min(abs(k_x - x))
    return(k_dens[index])
  }
  
  # Graph of true Beta distribution against estimates
  dens_true <- sapply(grid, function(x) {
    dbeta(x, shape1=shape[1], shape2=shape[2])
  })
  dens_est_gauss <- sapply(grid, kd_gauss)
  dens_est_gauss2 <- sapply(grid, kd_gauss2)
  dens_est_beta <- sapply(grid, kd_beta)
  dens_est_beta2 <- sapply(grid, kd_beta2)
  dat_plot <- data.frame(
    x = rep(grid, 5),
    density = c(dens_true, dens_est_gauss, dens_est_gauss2,
                dens_est_beta, dens_est_beta2),
    which = rep(c("True", "Est: Gaussian", "Est: Gaussian 2",
                  "Est: Beta", "Est: Beta 2"), each=length(grid))
  )
  ggplot(dat_plot, aes(x=x, y=density, color=factor(which))) +
    geom_line() +
    labs(color="Which")
  
  # # Histogram of sample
  # ggplot(data.frame(x=sample), aes(x=x)) + geom_histogram()
  
  # # Check to see if densities are properly normalized (should appx equal 1)
  # mean(dens_est_beta[1:100])
  # # mean(dens_est_beta2[1:100])
  # mean(dens_est_3[1:100])
  # mean(dens_est_gauss[1:100])
  

}
