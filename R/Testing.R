
#################.
##### Setup #####
#################.

{
  
  # 1. Run beginning of MAIN.R
  
  # 2. Set global constants
  params <- list(S_n_type="true", g_n_type="true", ci_type="regular",
                 cf_folds=1, edge_corr="none", deriv_type="spline")
  C <- sim$constants
  L <- list(n=1000, alpha_3=0.8, lambda2=0.0001, distr_A="N(0.5,0.01)",
            edge="none", surv_true="complex", sampling="two-phase")
  L$estimator <- list(est="Grenander", params=params)
  
}



#############################.
##### Empirical cdf G_n #####
#############################.

{
  
  # Generate datasets
  distr_A <- "Unif(0,1)"
  d1 <- generate_data(n=1000, alpha_3=0.7, distr_A=distr_A, edge="none",
                      surv_true="Cox PH", sampling="iid")
  d2 <- generate_data(n=1000, alpha_3=0.7, distr_A=distr_A, edge="none",
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



#########################.
##### deriv_theta_n #####
#########################.

{
  
  # Generate true deriv_thetas
  {
    library(numDeriv)
    
    surv_true <- "Cox PH"
    m <- 10^5
    w1 <- runif(m)
    w2 <- rbinom(m, size=1, prob=0.5)
    alpha_3 <- 0.75
    
    theta_true_f <- Vectorize(function(a) {
      
      lin <- function(w1,w2,a) {
        if (surv_true=="Cox PH") {
          C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1
        } else if (surv_true=="complex") {
          C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + alpha_3*w2*a - 1
        }
      }
      
      S_0 <- function(t, w1, w2, a) {
        exp( -1 * C$lambda * (t^C$v) * exp(lin(w1,w2,a)) )
      }
      
      return(1 - mean(S_0(C$t_e, w1, w2, a)))
      
    })
    
    deriv_theta_0 <- function (x) { grad(func=theta_true_f, x=x) }
  }
  
  # Esimate deriv_theta_n n_samples times
  C$appx=list(t_e=10,w1=0.01,w1b=0.1,a=0.01)
  n_samples <- 20
  deriv_ests <- c()
  grid <- seq(0,1,0.1)
  for (i in 1:n_samples) {
    
    distr_A <- "Beta(1.5+w1,1.5+w2)" # "Unif(0,1)"
    dat_orig <- generate_data(n=1000, alpha_3=0.75, distr_A=distr_A,
                              edge="none", surv_true="Cox PH",
                              sampling="two-phase")
    vlist <- create_val_list(dat_orig, C$appx)
    
    Phi_n <- construct_Phi_n(dat_orig)
    Phi_n_inv <- construct_Phi_n(dat_orig, type="inverse")
    S_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest")
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest",
                          csf=TRUE)
    # S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    # Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
    #                       csf=TRUE)
    f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                 type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    Gamma_n <- construct_Gamma_n(dat_orig, vlist$A_grid, omega_n, S_n, g_n)
    
    Psi_n <- Vectorize(function(x) { return(Gamma_n(Phi_n_inv(x))) })
    gcm <- gcmlcm(x=seq(0,1,0.01), y=Psi_n(seq(0,1,0.01)), type="gcm")
    dGCM <- Vectorize(function(x) {
      if (x==0) {
        index <- 1
      } else {
        index <- which(x<=gcm$x.knots)[1]-1
      }
      return(gcm$slope.knots[index])
    })
    theta_n <- function(x) {
      x_trans <- Phi_n(x)
      return(dGCM(x_trans))
    }
    deriv_theta_n <- construct_deriv_theta_n(theta_n, type="gcomp")
    deriv_ests <- c(deriv_ests, deriv_theta_n(grid))
    
  }
  
  # Plot true curves against estimated curve
  len <- length(grid)
  df <- data.frame(
    x = rep(grid, n_samples+1),
    y = c(deriv_ests, deriv_theta_0(grid)),
    which = c(rep("Estimate", len*n_samples), rep("True", len))
  )
  ggplot(df, aes(x=x, y=y, color=which)) +
    geom_point() +
    ylim(c(0,1)) +
    labs(title="Estimation of derivative of theta (RF)",
         color="Which")
  
}



###################.
##### omega_n #####
###################.

{
  
  # Note: this check only works for surv_true=="Cox PH"
  
  # Generate data
  dat_orig <- generate_data(
    n = 3000,
    alpha_3 = L$alpha_3,
    distr_A = "Beta(1.5+w1,1.5+w2)", # "Unif(0,1)" "Beta(1.5+w1,1.5+w2)"
    edge = "none",
    surv_true = L$surv_true,
    sampling = "two-phase" # two-phase iid
  )
  vlist <- create_val_list(dat_orig, C$appx)
  
  # True omega_0
  S_0 <- construct_S_n(dat_orig, vlist$S_n, type="true")
  Sc_0 <- construct_S_n(dat_orig, vlist$S_n, type="true", csf=TRUE)
  omega_0 <- Vectorize(function(w1,w2,a,y_star,delta_star) {
    
    if (L$surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + L$alpha_3*a - 1
    } else if (L$surv_true=="complex") {
      lin <- C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + L$alpha_3*w2*a - 1
    }
    
    piece_1 <- exp(-1*C$lambda*(C$t_e^C$v)*exp(lin))
    
    piece_2 <- (delta_star*as.integer(y_star<=C$t_e)) /
      exp(-exp(lin)*(C$lambda*y_star^C$v+L$lambda2*y_star^C$v2))
    
    integral <- integrate(
      function(t) {
        t^(C$v-1) * exp(lin+exp(lin)*(C$lambda*t^C$v+L$lambda2*t^C$v2))
      },
      lower = 0,
      upper = min(C$t_e,y_star)
    )$value
    
    piece_3 <- C$lambda*C$v*integral
    
    return(piece_1*(piece_2-piece_3))
    
  })
  
  # Construct omega_n
  S_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest")
  Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest", csf=TRUE)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  
  omegas_est <- c()
  omegas_true <- c()
  for (i in 1:nrow(dat_orig)) {
    omegas_est <- c(omegas_est, do.call(omega_0, as.list(vlist$omega[i,])))
    omegas_true <- c(omegas_true, do.call(omega_n, as.list(vlist$omega[i,])))
  }
  
  for (i in 1:500) {
    est <- do.call(omega_0, as.list(vlist$omega[i,]))
    true <- do.call(omega_n, as.list(vlist$omega[i,]))
    if (abs(est-true)>0.4) {
      print(paste("i:",i))
      print(vlist$omega[i,])
      print(paste("Est:",est))
      print(paste("True:",true))
    }
  }
  
  # Scatterplot of estimated vs true
  ggplot(data.frame(x=omegas_true, y=omegas_est), aes(x=x, y=y)) +
    geom_abline(slope=1, intercept=0, color="green") +
    geom_point(alpha=0.1) +
    labs(title="omega_n vs omega_0", x="omega_0", y="omega_n")
  
}



###################.
##### gamma_n #####
###################.

{
  
  # Generate datasets
  alpha_3 <- 0.7
  distr_A <- "Unif(0,1)" # "Unif(0,1)" "Beta(1.5+w1,1.5+w2)"
  surv_true <- "Cox PH"
  dat_orig <- generate_data(
    n = 2000,
    alpha_3 = alpha_3,
    distr_A = distr_A,
    edge = "none",
    surv_true = surv_true,
    sampling = "two-phase"
  )
  
  # True omega_0
  omega_0 <- Vectorize(function(w1,w2,a,y_star,delta_star) {
    
    if (surv_true=="Cox PH") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1
    } else if (surv_true=="complex") {
      lin <- C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + alpha_3*w2*a - 1
    }
    
    piece_1 <- exp(-1*C$lambda*(C$t_e^C$v)*exp(lin))
    
    piece_2 <- (delta_star*as.integer(y_star<=C$t_e)) /
      exp(-exp(lin)*(C$lambda*y_star^C$v+L$lambda2*y_star^C$v2))
    
    integral <- integrate(
      function(t) {
        t^(C$v-1) * exp(lin+exp(lin)*(C$lambda*t^C$v+L$lambda2*t^C$v2))
      },
      lower = 0,
      upper = min(C$t_e,y_star)
    )$value
    
    piece_3 <- C$lambda*C$v*integral
    
    return(piece_1*(piece_2-piece_3))
    
  })
  
  # Approximate true gamma_0
  construct_gamma_0 <- function() {
    
    epsilon <- 0.02
    
    dat_mc <- generate_data(
      n = 50000,
      alpha_3 = alpha_3,
      distr_A = distr_A,
      edge = "none",
      surv_true = surv_true,
      sampling = "two-phase"
    )
    
    f_aIw_0 <- function(a,w1,w2) {
      if (distr_A=="Unif(0,1)") {
        return(1)
      } else if (distr_A=="Beta(1.5+w1,1.5+w2)") {
        return(dbeta(a, shape1=1.5+w1, shape2=1.5+w2))
      }
    }
    
    return(Vectorize(function(x) {
      dat_mc %<>% filter(abs(a-x)<=epsilon)
      d <- dat_mc
      return(mean(
        (omega_0(d$w1,d$w2,d$a,d$y_star,d$delta_star) /
           f_aIw_0(d$a,d$w1,d$w2))^2
      ))
    }))
    
    
  }
  gamma_0 <- construct_gamma_0()
  
  # Construct gamma_n
  vlist <- create_val_list(dat_orig, C$appx)
  S_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest")
  Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest", csf=TRUE)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                               type=params$g_n_type, k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  f_aIw_delta1_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                      type=params$g_n_type, delta1=TRUE, k=15)
  f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                  f_aIw_delta1_n)
  gamma_n <- construct_gamma_n(dat_orig, vlist$A_grid, type="kernel",
                               omega_n, f_aIw_n, f_a_n, f_a_delta1_n)
  
  # Plot true curves against estimated curve
  grid <- seq(0,1,0.05)
  df <- data.frame(
    x = rep(grid, 2),
    y = c(gamma_n(grid), gamma_0(grid)),
    which = rep(c("Estimate","Truth"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) +
    geom_line() +
    # ylim(c(0,10)) +
    labs(title=paste("Estimation of gamma_0; distr_A:",distr_A),
         color="Which")
  
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
    edge = "none",
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
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type, k=15)
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
      edge = "none",
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
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    gcomp_n <- construct_gcomp_n(dat_orig, S_n)
    eta_n <- construct_eta_n(dat_orig, S_n)
    Gamma_n <- construct_Gamma_n(dat_orig, omega_n, S_n, g_n)
    
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
  dat_orig <- generate_data(n=400, alpha_3=0, distr_A="Unif(0,1)", edge="none",
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
    f_aIw_n <- construct_f_aIw_n(dat, vals_AW_grid, type=params$g_n_type, k=15)
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
                              edge="none", surv_true="Cox PH", sampling="iid")
    
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
  dat_orig <- generate_data(
    n = 3000,
    alpha_3 = L$alpha_3,
    distr_A = L$distr_A,
    edge = L$edge,
    surv_true = L$surv_true,
    sampling = L$sampling
  )
  vlist <- create_val_list(dat_orig, C$appx)
  
  # S_n_CoxPH <- construct_S_n(dat_orig, vlist$S_n, type="Cox PH")
  # S_n_RF <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest")
  S_0 <- construct_S_n(dat_orig, vlist$S_n, type="true")
  S_n_CoxPH <- S_0
  S_n_RF <- S_0
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      S_n_CoxPH(t=times, w1=0.2, w2=1, a=0.2),
      S_n_CoxPH(t=times, w1=0.5, w2=1, a=0.2),
      S_n_CoxPH(t=times, w1=0.2, w2=1, a=0.8),
      S_n_CoxPH(t=times, w1=0.5, w2=1, a=0.8),
      S_n_RF(t=times, w1=0.2, w2=1, a=0.2),
      S_n_RF(t=times, w1=0.5, w2=1, a=0.2),
      S_n_RF(t=times, w1=0.2, w2=1, a=0.8),
      S_n_RF(t=times, w1=0.5, w2=1, a=0.8),
      S_0(t=times, w1=0.2, w2=1, a=0.2),
      S_0(t=times, w1=0.5, w2=1, a=0.2),
      S_0(t=times, w1=0.2, w2=1, a=0.8),
      S_0(t=times, w1=0.5, w2=1, a=0.8)
    ),
    which = rep(c("Cox PH","RF","True S_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0.2,w2=1,a=0.2","w1=0.5,w2=1,a=0.2",
                     "w1=0.2,w2=1,a=0.8","w1=0.5,w2=1,a=0.8"),3), each=length(times))
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
  dat_orig <- generate_data(
    n = 3000,
    alpha_3 = L$alpha_3,
    distr_A = L$distr_A,
    edge = L$edge,
    surv_true = L$surv_true,
    sampling = L$sampling
  )
  vlist <- create_val_list(dat_orig, C$appx)
  
  # Sc_n_CoxPH <- construct_S_n(dat_orig, vlist$S_n, type="Cox PH", csf=TRUE)
  # Sc_n_RF <- construct_S_n(dat_orig, vlist$S_n, type="Random Forest", csf=TRUE)
  Sc_0 <- construct_S_n(dat_orig, vlist$S_n, type="true", csf=TRUE)
  Sc_n_CoxPH <- Sc_0
  Sc_n_RF <- Sc_0
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      Sc_n_CoxPH(t=times, w1=0.2, w2=1, a=0.2),
      Sc_n_CoxPH(t=times, w1=0.1, w2=1, a=0.2),
      Sc_n_CoxPH(t=times, w1=0.2, w2=1, a=0.8),
      Sc_n_CoxPH(t=times, w1=0.1, w2=1, a=0.8),
      Sc_n_RF(t=times, w1=0.2, w2=1, a=0.2),
      Sc_n_RF(t=times, w1=0.1, w2=1, a=0.2),
      Sc_n_RF(t=times, w1=0.2, w2=1, a=0.8),
      Sc_n_RF(t=times, w1=0.1, w2=1, a=0.8),
      Sc_0(t=times, w1=0.2, w2=1, a=0.2),
      Sc_0(t=times, w1=0.1, w2=1, a=0.2),
      Sc_0(t=times, w1=0.2, w2=1, a=0.8),
      Sc_0(t=times, w1=0.1, w2=1, a=0.8)
    ),
    which = rep(c("Cox PH","RF","True Sc_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0.2,w2=1,a=0.2","w1=0.5,w2=1,a=0.2",
                     "w1=0.2,w2=1,a=0.8","w1=0.5,w2=1,a=0.8"),3), each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional (censoring) survival: Sc_0[t|W,A]",
         color="Estimator")
  
}



#############################.
##### theta_n estimator #####
#############################.

{
  
  # Generate data
  set.seed(10)
  dat_orig <- generate_data(
    n = 2000,
    alpha_3 = 0.7,
    distr_A = "Beta(1.5+w1,1.5+w2)", # "Unif(0,1)" "Beta(1.5+w1,1.5+w2)"
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase"
  )
  
  # Obtain estimates
  ests <- est_curve(
    dat_orig = dat_orig,
    estimator = "Grenander",
    params = list(S_n_type="true", g_n_type="true",
                  deriv_type="spline", ci_type="regular", cf_folds=1,
                  edge_corr="none"),
    points = C$points
  )
  
  # Return results
  theta_true <- attr(dat_orig, "theta_true")
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
    which = rep(c("Est (RF)","Truth"), each=len),
    ci_lo = c(ci_lo, theta_true),
    ci_hi = c(ci_hi, theta_true)
  )
  ggplot(plot_data, aes(x=x, y=theta, color=factor(which))) +
    geom_line() +
    labs(color="Which") +
    # ylim(c(0.3,0.7)) +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_hi),
      alpha = 0.2,
      fill = NA,
      linetype = "dotted"
    )
  
}



##########################################.
##### Conditional density estimators #####
##########################################.

{
  
  # Set levels here
  n <- 2000
  distr_A <- "Unif(0,1)"
  # distr_A <- "Beta(1.5+w1,1.5+w2)"
  edge <- "none"
  sampling <- "two-phase"
  
  # Generate data
  dat_orig <- generate_data(
    n = n,
    alpha_3 = 0.7,
    distr_A = distr_A,
    edge = edge,
    surv_true = "Cox PH",
    sampling = sampling
  )
  
  # True conditional density function
  f_aIw_0 <- function(a,w1,w2) {
    if (distr_A=="Unif(0,1)") {
      return(1)
    } else if (distr_A=="Beta(1.5+w1,1.5+w2)") {
      return(dbeta(a, shape1=1.5+w1, shape2=1.5+w2))
    }
  }
  
  # Parametric estimate
  vlist <- create_val_list(dat_orig, C$appx)
  f_aIw_n_para <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                    type="parametric", k=0)
  f_aIw_n_semi <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                                    # type="binning", k=15)
                                    type="binning", k=0)
  
  # Generate plot data
  grid <- seq(0,1,0.01)
  f_aIw_models <- c("Truth", "Parametric", "Semiparametric")
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
      sapply(grid, function(a) { f_aIw_n_semi(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w1=0.8, w2=1) })
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
  
  # Check marginal density
  f_a_0 <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_0)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n_semi)
  
  # Generate plot data
  grid <- seq(0,1,0.01)
  f_a_models <- c("Truth", "Semiparametric")
  len <- length(grid)
  plot_data <- data.frame(
    a = rep(grid, 2),
    density = c(f_a_0(grid),f_a_n(grid)),
    which = rep(f_a_models, each=len)
  )
  ggplot(plot_data, aes(x=a, y=density, color=factor(which))) +
    geom_line() +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of marginal density: f(A)") +
    ylim(c(0,NA))
  
}



#######################################.
##### Propensity score estimators #####
#######################################.

{
  
  # Set levels here
  n <- 1000
  edge <- "expit"
  
  # Generate data
  dat_orig <- generate_data(
    n = n,
    alpha_3 = 0.7,
    distr_A = "Beta(1.5+w1,1.5+w2)",
    edge = edge,
    surv_true = "Cox PH",
    sampling = "two-phase"
  )
  
  # True propensity score function
  pi_0 <- function(w1,w2) {
    if (edge=="expit") {
      return(expit(w1+w2-3.3))
    } else if (edge=="complex") {
      return(0.84*w2*pmax(0,1-4*abs(w1-0.5)))
    }
  }
  
  # Construct estimators
  vlist <- create_val_list(dat_orig, C$appx)
  pi_n_logistic <- construct_pi_n(dat_orig, vlist$W_grid, type="logistic")
  pi_n_SL <- construct_pi_n(dat_orig, vlist$W_grid, type="SL")
  
  # Curve 1: W2=0
  pi_n_log_0 <- function(w1) { pi_n_logistic(w1,w2=0) }
  pi_n_SL_0 <- function(w1) { pi_n_SL(w1,w2=0) }
  pi_0_0 <- function(w1) { pi_0(w1,w2=0) }
  
  # Curve 2: W2=1
  pi_n_log_1 <- function(w1) { pi_n_logistic(w1,w2=1) }
  pi_n_SL_1 <- function(w1) { pi_n_SL(w1,w2=1) }
  pi_0_1 <- function(w1) { pi_0(w1,w2=1) }
  
  # Plot true curves against estimated curve
  grid <- seq(0,1,0.01)
  curves <- c("W2=0","W2=1")
  estimators <- c("Logistic","SL","True")
  len <- length(grid)
  n_curves <- length(curves)
  n_estimators <- length(estimators)
  df <- data.frame(
    x = rep(grid, n_curves*n_estimators),
    y = c(
      pi_n_log_0(grid),
      pi_n_SL_0(grid),
      pi_0_0(grid),
      pi_n_log_1(grid),
      pi_n_SL_1(grid),
      pi_0_1(grid)
    ),
    curve = rep(rep(curves, each=n_estimators*len)),
    which = rep(rep(estimators, each=len), n_curves)
  )
  print(ggplot(df, aes(x=x, y=y, color=which)) +
          geom_line() +
          facet_wrap(~curve, ncol=2) +
          labs(title="Estimation of propensity score function",
               color="Which"))
  
}
