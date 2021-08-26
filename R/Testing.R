
# !!!!! Right now, I'm collecting previous ad-hoc tests here. Need to make this
#     more structured

#################.
##### Setup #####
#################.

{
  
  # 1. Run CONFIG and SETUP within MAIN.R
  
  # 2. Set global constants
  params <- list(g_n_type="parametric", S_n_type="Cox PH")
  C <- list(lambda=10^-4, v=1.5, lambda2=0.5*10^-4, v2=1.5,
            points=seq(0,1,0.1), alpha_1=0.3, alpha_2=0.7, t_e=200,
            appx=list(t_e=10,w1=0.01,w1b=0.1,a=0.01))
  
}



##########################.
##### Code profiling #####
##########################.

{
  
  # Generate datasets
  dat_iid <- generate_data(n=400, alpha_3=0.7, distr_A="Unif(0,1)",
                      surv_true="Cox PH", sampling="iid")
  dat_tp <- generate_data(n=400, alpha_3=0.7, distr_A="Unif(0,1)",
                      surv_true="Cox PH", sampling="two-phase")
  
  # Estimation
  ests <- est_curve(
    dat_orig = dat_tp,
    estimator = "Grenander",
    params = list(
      S_n_type = params$S_n_type,
      g_n_type = params$g_n_type,
      ci_type = "logit",
      cf_folds = 1
    ),
    points = C$points
  )
  
  # Testing
  reject <- test_2(
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
  
}



#############################.
##### Empirical cdf G_n #####
#############################.

{
  
  # Generate datasets
  d1 <- generate_data(n=1000, alpha_3=0.7, distr_A="Unif(0,1)",
                      surv_true="Cox PH", sampling="iid")
  d2 <- generate_data(n=1000, alpha_3=0.7, distr_A="Unif(0,1)",
                      surv_true="Cox PH", sampling="two-phase")
  
  for (dat_orig in list(d1,d2)) {
    
    # Curve 1: Phi_n
    Phi_n <- construct_Phi_n(dat_orig)
    Phi_0 <- function(x) {x}
    
    # Curve 2: Phi_n_inv
    Phi_n_inv <- construct_Phi_n(dat_orig, type="inverse")
    Phi_0_inv <- function(x) {x}
    
    # Plot true curves against estimated curve
    grid <- seq(0,1,0.01)
    curves <- c("Phi","Phi_inv")
    estimators <- c("Estimate","True")
    len <- length(grid)
    n_curves <- length(curves)
    n_estimators <- length(estimators)
    df <- data.frame(
      x = rep(grid, n_curves*n_estimators),
      y = c(
        Phi_n(grid),
        Phi_0(grid),
        Phi_n_inv(grid),
        Phi_0_inv(grid)
      ),
      curve = rep(rep(curves, each=n_estimators*len)),
      which = rep(rep(estimators, each=len), n_curves)
    )
    print(ggplot(df, aes(x=x, y=y, color=which)) +
      geom_line() +
      facet_wrap(~curve, ncol=2) +
      labs(title="Estimation of empirical CDFs (and inverse ECDFs)",
           color="Which"))

  }
  
}



######################################################.
##### Checking the influence function of Gamma_n #####
######################################################.

{
  
  # Set up data
  dat <- generate_data(
    n = 400, # 5000
    alpha_3 = 0,
    distr_A = "Unif(0,1)",
    surv_true = "Cox PH",
    sampling = "iid" # iid two-phase
  )
  
  # Define the statistic to bootstrap
  bootstat <- function(dat_orig, indices) {
    
    dat_orig <- dat_orig[indices,]
    
    # Construct component functions
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat_orig)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    Gamma_n <- construct_Gamma_n(dat_orig, omega_n, S_n, g_n)
    
    # Return the value of Gamma_n(0.5)
    return (Gamma_n(0.5))
    
  }
  
  # Run bootstrap
  boot_obj <- boot(data=dat, statistic=bootstat, R=50)
  
  # See actual distribution of test statistic
  Gammas <- c()
  for (i in 1:50) {
    
    dat <- generate_data(
      n = 400, # 5000
      alpha_3 = 0,
      distr_A = "Unif(0,1)",
      surv_true = "Cox PH",
      sampling = "iid" # iid two-phase
    )
    
    Gammas <- c(Gammas, bootstat(dat, c(1:nrow(dat))))
    
  }
  
  # Estimate variance using influence function
  {
    
    # Construct component functions
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat_orig)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    gcomp_n <- construct_gcomp_n(dat_orig, S_n)
    eta_n <- construct_eta_n(dat_orig, S_n)
    Gamma_n <- construct_Gamma_n(dat_orig, omega_n, S_n, g_n)
    
    # !!!!!
    if (F) {
      
      plot1 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point() +
        geom_smooth(formula=y~x, method=lm)
      
      plot2 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point()
      model <- lm(po~a, data=filter(dat_orig,!is.na(a)))
      coeff1 <- as.numeric(model$coefficients)
      plot2$layers <- c({
        stat_function(fun = function(x) {
          coeff1[1] + coeff1[2]*x
        }, color="turquoise")
      }, plot2$layers)
      
      plot3 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point()
      model <- lm(po~a+I(a^2)+I(a^3), data=filter(dat_orig,!is.na(a)))
      coeff <- as.numeric(model$coefficients)
      plot3$layers <- c({
        stat_function(fun = function(x) {
          coeff[1] + coeff[2]*x + coeff[3]*(x^2) + coeff[4]*(x^3)
        }, color="turquoise")
      }, plot3$layers)
      
      
      plot1
      plot2
      plot3
      #
      
      
      
      
    }
    
    # Construct influence function
    infl_fn_Gamma <- construct_infl_fn_Gamma(dat_orig, omega_n, g_n, gcomp_n,
                                             eta_n, Gamma_n)
    
    # Estimate variance and SD
    var_hat <- mean((
      infl_fn_Gamma(x=0.5,dat$w1,dat$w2,dat$y_star,dat$delta_star,dat$delta,dat$a)
    )^2)
    sd_hat <- sqrt(var_hat/nrow(dat))
    
  }
  
  # Bootstrap SE
  # n=400: 0.0315789
  #        0.03146832
  print(sd(boot_obj$t))
  
  # Empirical SE
  # n=400: 0.03128305
  #        0.03187881
  print(sd(Gammas))
  
  # IF-based variance estimator
  # n=400: 0.02280982
  #        0.02301398
  print(sd_hat)
  
}



#################################################.
##### New test statistic variance estimator #####
#################################################.

{
  
  # Generate data
  dat_orig <- generate_data(n=400, alpha_3=0, distr_A="Unif(0,1)",
                            surv_true="Cox PH",sampling="iid")
  
  # Define the test statistic
  test_stat <- function(dat_orig, indices) {
    
    dat_orig <- dat_orig[indices,]
    
    # Prep
    n_orig <- nrow(dat_orig)
    dat_orig$weights <- wts(dat_orig)
    dat <- dat_orig %>% filter(!is.na(a))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vals_A <- data.frame(a=dat$a)
    vals_AW <- data.frame(a=dat$a, w1=dat$w1, w2=dat$w2)
    vals_A_grid <- data.frame(a=seq(0,1,C$appx$a))
    vals_AW_grid <- expand.grid(a=seq(0,1,C$appx$a), w1=seq(0,1,C$appx$w1),
                                w2=c(0,1))
    vals_S_n <- expand.grid(t=seq(0,C$t_e,C$appx$t_e), w1=seq(0,1,C$appx$w1b),
                            w2=c(0,1), a=seq(0,1,C$appx$a))
    vals_omega <- subset(dat, select=-c(delta,weights))
    
    # Construct component functions
    G_n <- construct_Phi_n(dat_orig)
    f_aIw_n <- construct_f_aIw_n(dat, vals_AW_grid, type=params$g_n_type)
    f_a_n <- construct_f_a_n(dat_orig, vals_A_grid, f_aIw_n)
    g_n <- construct_g_n(vals_AW_grid, f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vals_omega, S_n, Sc_n)
    Gamma_n <- construct_Gamma_n(dat_orig, vals_A_grid, omega_n, S_n, g_n)
    gcomp_n <- construct_gcomp_n(dat_orig, vals_A_grid, S_n)
    eta_n <- construct_eta_n(dat_orig, vals_AW_grid, S_n)
    xi_n <- construct_xi_n(Phi_n=G_n, lambda_2, lambda_3)
    rho_n <- construct_rho_n(dat_orig, Phi_n=G_n)
    lambda_2 <- lambda(dat_orig,2,G_n)
    lambda_3 <- lambda(dat_orig,3,G_n)
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * (
        lambda_2*(G_n(dat$a))^2 -
          lambda_3*G_n(dat$a)
      ) *
        (Gamma_n(dat$a))
    )
    
    return (beta_n)
    
  }
  
  
  # See actual distribution of test statistic
  betas <- c()
  for (i in 1:50) {
    
    dat_orig <- generate_data(n=400, alpha_3=0, distr_A="Unif(0,1)",
                              surv_true="Cox PH",sampling="iid")
    
    beta_n <- test_stat(dat_orig, c(1:nrow(dat_orig)))
    betas <- c(betas, beta_n)
    
    print(paste("i:",i))
    
  }
  
  # Histogram of betas
  ggplot(data.frame(x=betas), aes(x=x)) + geom_histogram()
  
  # New variance estimator
  sd_hat <- test_2(
    dat_orig = dat_orig,
    alt_type = "incr",
    params = list(
      var = "asymptotic",
      S_n_type = params$S_n_type,
      g_n_type = params$g_n_type,
      est_known_nuis = FALSE,
      cf_folds = 1
    ),
    return_sd = TRUE
  )
  
  # Empirical SE
  # n=400: 0.0002369075
  print(sd(betas))
  
  # New variance estimator
  # n=400: 0.0002077375
  # n=400: 0.000221076
  print(sd_hat)
  
  run_bootstrap <- FALSE
  if (run_bootstrap) {
    
    # # Bootstrap SE
    # # OLD: n=400: 0.0002719989
    # # OLD: n=800: 0.0001763408
    # print(sd(boot_obj$t))
    
    # Run bootstrap
    boot_obj <- boot(data=dat, statistic=test_stat, R=50)
    
  }
  
}



####################################################.
##### Conditional survival function estimators #####
####################################################.

{
  
  # Generate data
  alpha_3 <- 0.7
  surv_true <- "Complex"
  dat <- generate_data(
    n = 2000,
    alpha_3 = alpha_3,
    distr_A = "Unif(0,1)",
    surv_true = surv_true,
    sampling = "two-phase"
  )
  
  S_n <- construct_S_n(dat, type="Cox PH")
  
  S_0 <- Vectorize(function(t, w1, w2, a) {
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
      return(exp(-1*C$lambda*(t^C$v)*exp(lin)))
    } else if (surv_true=="Complex") {
      # lin <- C$alpha_2*w2*as.numeric(abs(w1-0.5)<0.2) + alpha_3*w1*a
      lin <- w2*w1*a
      return(exp(-1*C$lambda*(t^C$v)*exp(lin)))
    }
  })
  
  # !!!!!
  {
    # survlistWrappers()
    dat2 <- filter(dat, !is.na(a))
    weights <- wts(dat2)
    newX <- expand.grid(w1=c(0,1), w2=c(1), a=c(0,1))
    # sl_library <- c("survSL.coxph", "survSL.weibreg", "survSL.rfsrc", "survSL.km")
    sl_library <- c("survSL.coxph", "survSL.km", "survSL.rfsrc")
    
    # [1] "survSL.coxph"     "survSL.expreg"    ""      
    # [4] "survSL.km"        "survSL.loglogreg" "survSL.pchreg"   
    # [7] "survSL.pchSL"     "survSL.require"   ""    
    # [10] "survSL.template"  "survSL.weibreg"    
    
    srv <- survSuperLearner(
      time = dat2$y_star,
      event = dat2$delta_star,
      X = subset(dat2,select=c(w1,w2,a)),
      newX = newX,
      new.times = c(1:200),
      # event.SL.library = c("survSL.coxph"),
      # cens.SL.library = c("survSL.coxph"),
      event.SL.library = sl_library,
      cens.SL.library = sl_library,
      control = list(max.SL.iter=20),
      obsWeights = weights
    )
    
  }
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      S_n(t=times, w1=0, w2=1, a=0),
      S_n(t=times, w1=1, w2=1, a=0),
      S_n(t=times, w1=0, w2=1, a=1),
      S_n(t=times, w1=1, w2=1, a=1),
      srv$event.SL.predict[1,],
      srv$event.SL.predict[2,],
      srv$event.SL.predict[3,],
      srv$event.SL.predict[4,],
      S_0(t=times, w1=0, w2=1, a=0),
      S_0(t=times, w1=1, w2=1, a=0),
      S_0(t=times, w1=0, w2=1, a=1),
      S_0(t=times, w1=1, w2=1, a=1)
    ),
    which = rep(c("Cox PH","SL","True S_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0,a=0","w1=1,a=0","w1=0,a=1","w1=1,a=1"),3),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: S_0[t|W,A]",
         color="Estimator")
  
}



#####################################################.
##### Conditional censoring function estimators #####
#####################################################.

{
  
  # Generate data
  alpha_3 <- 0.7
  dat <- generate_data(
    n = 5000,
    alpha_3 = alpha_3,
    distr_A = "Unif(0,1)",
    surv_true = "Cox PH",
    sampling = "two-phase"
  )
  
  Sc_n <- construct_S_n(dat, type="Cox PH", csf=TRUE)
  
  Sc_0 <- Vectorize(function(t, w1, w2, a) {
    lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a
    return(exp(-1*C$lambda2*(t^C$v2)*exp(lin)))
  })
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 8),
    survival = c(
      Sc_n(t=times, w1=0, w2=1, a=0),
      Sc_n(t=times, w1=1, w2=1, a=0),
      Sc_n(t=times, w1=0, w2=1, a=1),
      Sc_n(t=times, w1=1, w2=1, a=1),
      Sc_0(t=times, w1=0, w2=1, a=0),
      Sc_0(t=times, w1=1, w2=1, a=0),
      Sc_0(t=times, w1=0, w2=1, a=1),
      Sc_0(t=times, w1=1, w2=1, a=1)
    ),
    which = rep(c("Cox PH","True S^C_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0,a=0","w1=1,a=0","w1=0,a=1","w1=1,a=1"),2),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: S_0[t|W,A]",
         color="Estimator")
  
}



#############################.
##### theta_n estimator #####
#############################.

{
  
  # Generate data
  set.seed(1)
  dat <- generate_data(
    n = 1000,
    alpha_3 = 0.7,
    distr_A = "Unif(0,1)",
    surv_true = "Cox PH",
    sampling = "two-phase"
  )
  
  # Obtain estimates
  ests <- est_curve(
    dat = dat,
    estimator = "Grenander",
    params = list(
      S_n_type = params$S_n_type,
      g_n_type = params$g_n_type,
      ci_type = "logit", # none
      cf_folds = 1
    ),
    points = C$points
  )
  
  # Return results
  theta_true <- attr(dat, "theta_true")
  theta_ests <- c()
  ci_lo <- c()
  ci_hi <- c()
  len <- length(C$points)
  for (i in 1:len) {
    theta_ests <- c(theta_ests, ests[[i]]$est)
    ci_lo <- c(ci_lo, ests[[i]]$ci_lo)
    ci_hi <- c(ci_hi, ests[[i]]$ci_hi)
  }
  
  plot_data <- data.frame(
    x = rep(C$points, 2),
    theta = c(theta_ests, theta_true),
    which = rep(c("Estimate","Truth"), each=len),
    ci_lo = c(ci_lo, theta_true),
    ci_hi = c(ci_hi, theta_true)
  )
  ggplot(plot_data, aes(x=x, y=theta, color=factor(which))) +
    geom_line() +
    ylim(c(0,1)) +
    labs(color="Which", fill="Which") +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_hi, fill=factor(which)),
      alpha = 0.2,
      linetype = "dotted"
    )
  
}



##########################################.
##### Conditional density estimators #####
##########################################.

{
  
  # Set levels here
  n <- 5000
  distr_A <- "Beta(0.8+0.9*w1,0.8+0.4*w2)"  # Unif(0,1) Beta(0.9,1.1+0.4*w2) Beta(0.8+0.9*w1,0.8+0.4*w2)
  sampling <- "two-phase"                     # iid two-phase
  
  # Generate data
  dat <- generate_data(
    n = n,
    alpha_3 = 0.7,
    distr_A = distr_A,
    surv_true = "Cox PH",
    sampling = sampling
  )
  
  # True conditional density function
  f_aIw_0 <- function(a,w1,w2) {
    if (distr_A=="Unif(0,1)") {
      return(1)
    } else if (distr_A=="Beta(0.9,1.1+0.4*w2)") {
      shape1 <- 0.9
      shape2 <- 1.1 + 0.4*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    } else if (distr_A=="Beta(0.8+0.9*w1,0.8+0.4*w2)") {
      shape1 <- 0.8 + 0.9*w1
      shape2 <- 0.8 + 0.4*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    }
  }
  
  # Parametric estimate
  {
    dat_trunc <- filter(dat, !is.na(a))
    n_trunc <- nrow(dat_trunc)
    weights <- wts(dat_trunc)
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        shape1 <- par[1] + par[2]*dat_trunc$w1[i]
        shape2 <- par[3] + par[4]*dat_trunc$w2[i]
        loglik <- dbeta(dat_trunc$a[i], shape1=shape1, shape2=shape2, log=TRUE)
        return(loglik*weights[i])
      }))
      
      return(-1*sum_loglik)
      
    }
    opt <- optim(par=c(a1=0.5, a2=0.1, a3=0.5, a4=0.1), fn=wlik)
    f_aIw_n_para <- Vectorize(function(a, w1, w2){
      shape1 <- opt$par[1] + opt$par[2]*w1
      shape2 <- opt$par[3] + opt$par[4]*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    })
  }
  
  # Nonparametric estimate (adapted from Diaz & VDL)
  {
    
    # k is the number of bins
    construct_f_aIw_n_nonpar <- function(dat, k) {
      
      alphas <- seq(0, 1, length.out=k+1)
      dat_trunc <- filter(dat, !is.na(a))
      n_trunc <- nrow(dat_trunc)
      weights <- wts(dat_trunc)
      
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
      
      wlik <- function(par) {
        
        # par[1] through par[k-1] are the hazard components for bins 1 to k-1
        # par[k] and par[k+1] correspond to W1 and W2
        sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
          lik <- dens(a=dat_trunc$a[i], w1=dat_trunc$w1[i], w2=dat_trunc$w2[i],
                      par)
          return(weights[i]*log(lik))
        }))
        
        return(-1*sum_loglik)
        
      }
      
      opt <- optim(par=rep(0,k+1), fn=wlik, method="CG")
      if (opt$convergence!=0) {
        warning("Nonpar conditional density: optim() did not converge")
      }
      
      return(Vectorize(memoise(function(a, w1, w2){
        
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        par <- opt$par
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(par[j] + par[k]*w1 + par[k+1]*w2)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        
        return(k*p1*p2)
        
      })))
      
    }
    
    f_aIw_n_nonpar <- construct_f_aIw_n_nonpar(dat, k=10)
    
  }
  
  # Generate plot data
  grid <- seq(0.01,0.99,0.01)
  f_aIw_models <- c("Truth", "Parametric", "Nonparametric")
  n_models <- length(f_aIw_models)
  len <- length(grid)
  plot_data <- data.frame(
    a = rep(grid, 4*n_models),
    density = c(
      sapply(grid, function(a) { f_aIw_0(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.8, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.8, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.8, w2=1) })
    ),
    which = rep(f_aIw_models, each=len*4),
    covariates = rep(c(
      rep("W1=0.2, W2=0",len),
      rep("W1=0.8, W2=0",len),
      rep("W1=0.2, W2=1",len),
      rep("W1=0.8, W2=1",len)
    ), n_models)
  )
  ggplot(plot_data, aes(x=a, y=density, color=factor(which))) +
    geom_line() +
    facet_wrap(~covariates, ncol=4) +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of conditional density: f(A|W)") +
    ylim(c(0,NA))
  
}
