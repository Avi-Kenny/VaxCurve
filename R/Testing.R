
#################.
##### Setup #####
#################.

{
  
  # 1. Run beginning of MAIN.R
  
  # 2. Set global constants
  params <- list(
    S_n_type="Cox PH", g_n_type="binning", ci_type="regular", cf_folds=1,
    edge_corr="none", ecdf_type="linear (mid)", deriv_type="m-spline",
    gamma_type="kernel", gamma_which="new", omega_n_type="estimated",
    marg="Gamma_star"
  )
  C <- list(points=seq(0,1,0.02), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10, w_tol=25, y_star=1, a=0.01))
  L <- list(
    n=500, alpha_3=-2, dir="decr",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
    # sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-7, v2=1.5), # Uncomment this for (almost) no censoring
    distr_A="Unif(0,1)", edge="none", surv_true="Cox PH",
    sampling="two-phase (50%)", estimator=list(est="Grenander",params=params)
  )
  
  # 3. Generate dataset
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir)
  
}



##########################################.
##### Dataset descriptive statistics #####
##########################################.

{
  
  res <- function(title, s, rnd) {
    print(title)
    print(paste0(round(mean(s),rnd), " (", round(min(s),rnd),
                 " -- ", round(max(s),rnd), ")"))
  }
  
  # Set up counters
  pct_rc <- num_inf <- pct_inf <- num_tp <- c()
  pct_inf_a0 <- pct_inf_a3 <- pct_inf_a5 <- pct_inf_a7 <- pct_inf_a1 <- c()
  
  n_reps <- 10
  for (i in 1:n_reps) {
    
    # # !!!!!
    # dat_orig <- sim$creators$generate_data(
    #   n = 5000,
    #   alpha_3 = 3,
    #   distr_A = "N(0.5,0.01)",
    #   edge = "none",
    #   surv_true = "Cox PH",
    #   sc_params = sim$levels$sc_params$sc_params,
    #   sampling = "two-phase (72%)",
    #   dir = "decr"
    # )
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Calculate summary stats
    pct_rc <- c(pct_rc, mean(as.integer(dat_orig$delta_star==0)*
                                     as.integer(dat_orig$y_star<=200)))
    num_inf <- c(num_inf, sum(as.integer(dat_orig$delta_star==1)*
                                as.integer(dat_orig$y_star<=200)))
    pct_inf <- c(pct_inf, mean(as.integer(dat_orig$delta_star==1)*
                                 as.integer(dat_orig$y_star<=200)))
    num_tp <- c(num_tp, sum(dat_orig$delta))
    pct_inf_a0 <- c(pct_inf_a0, attr(dat_orig,"theta_true")[1])
    pct_inf_a3 <- c(pct_inf_a3, attr(dat_orig,"theta_true")[16])
    pct_inf_a5 <- c(pct_inf_a5, attr(dat_orig,"theta_true")[26])
    pct_inf_a7 <- c(pct_inf_a7, attr(dat_orig,"theta_true")[36])
    pct_inf_a1 <- c(pct_inf_a1, attr(dat_orig,"theta_true")[51])
    
  }
  
  res("% right-censored by t_e", pct_rc, 3)
  res("# of infections by t_e", num_inf, 0)
  res("% of infections by t_e", pct_inf, 4)
  res("# in phase-2 sample", num_tp, 0)
  res("True inf rate (A=0)", pct_inf_a0, 5)
  res("True inf rate (A=0.3)", pct_inf_a3, 5)
  res("True inf rate (A=0.5)", pct_inf_a5, 5)
  res("True inf rate (A=0.7)", pct_inf_a7, 5)
  res("True inf rate (A=1.0)", pct_inf_a1, 5)
  
}



#############################.
##### Empirical cdf G_n #####
#############################.

{
  
  # !!!!!
  p1a <- c(); p1b <- c(); p2a <- c(); p2b <- c();
  n_reps <- 100
  for (i in c(1:n_reps)) {
    print(paste("rep:",i))
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    Phi_n1 <- construct_Phi_n(dat, type="step")
    Phi_n2 <- construct_Phi_n(dat, type="linear (mid)")
    p1a <- c(p1a, Phi_n1(0.5))
    p1b <- c(p1b, Phi_n1(0.8))
    p2a <- c(p2a, Phi_n2(0.5))
    p2b <- c(p2b, Phi_n2(0.8))
  }
  p1a <- p1a - 0.5
  p1b <- p1b - 0.8
  p2a <- p2a - 0.5
  p2b <- p2b - 0.8
  ggplot(
    data.frame(
      x = c(p1a,p1b,p2a,p2b),
      grp = rep(c("step (0.5)","step (0.8)","linear (0.5)","linear (0.8)"),
                each=n_reps)
    ),
    aes(x=x, group=grp, fill=factor(grp))) +
    geom_histogram(color="white") +
    geom_vline(xintercept=0, color="grey") +
    facet_wrap(~grp, ncol=2)
  
  # !!!!!
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, "iid", L$dir) # L$sampling
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  vlist <- create_val_list(dat, C$appx)
  Phi_n1 <- construct_Phi_n(dat, type="step")
  Phi_n2 <- ecdf(dat$a)
  round(Phi_n1(seq(0,1,0.1))-Phi_n2(seq(0,1,0.1)),10) # !!!!!
  
  # Curve 1: Phi_n
  Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
  Phi_0 <- construct_Phi_n(dat, type="true")
  
  # Curve 2: Phi_n_inv
  Phi_n_inv <- construct_Phi_n(dat, which="inverse",
                               type=params$ecdf_type)
  Phi_0_inv <- construct_Phi_n(dat, which="inverse",
                               type="true")
  
  # # Plot true curves against estimated curve
  # grid <- seq(0,1,0.01)
  # curves <- c("Phi","Phi_inv")
  # estimators <- c("True","Estimate")
  # len <- length(grid)
  # n_curves <- length(curves)
  # n_estimators <- length(estimators)
  # df <- data.frame(
  #   x = rep(grid, n_curves*n_estimators),
  #   y = c(
  #     Phi_0(grid),
  #     Phi_n(grid),
  #     Phi_0_inv(grid),
  #     Phi_n_inv(grid)
  #   ),
  #   curve = rep(rep(curves, each=n_estimators*len)),
  #   which = rep(rep(estimators, each=len), n_curves)
  # )
  # ggplot(df, aes(x=x, y=y, color=which)) +
  #   geom_line() +
  #   facet_wrap(~curve, ncol=2) +
  #   # xlim(c(0.25,0.35)) + # !!!!!
  #   # ylim(c(0.25,0.35)) + # !!!!!
  #   labs(title="Estimation of empirical CDFs (and inverse ECDFs)",
  #        color="Which")
  
  # !!!!!
  grid <- seq(0,1,0.01)
  curves <- c("Phi","Phi_inv")
  len <- length(grid)
  n_curves <- length(curves)
  df <- data.frame(
    x = rep(grid, n_curves),
    y = c(
      Phi_n(grid),
      Phi_n_inv(grid)
    ),
    curve = rep(rep(curves, each=len))
  )
  ggplot(df, aes(x=x, y=y, color=curve)) +
    geom_line() +
    # xlim(c(0.25,0.35)) + # !!!!!
    # ylim(c(0.25,0.35)) + # !!!!!
    labs(title="Estimation of empirical CDFs (and inverse ECDFs)",
         color="Which") +
    geom_abline(slope=1, color="grey")
  
}



#########################.
##### deriv_theta_n #####
#########################.

{
  
  # Generate true deriv_thetas
  {
    library(numDeriv)
    
    m <- 10^5
    w1 <- sample(seq(0,1,0.1), size=m, replace=T)
    w2 <- rbinom(m, size=1, prob=0.5)
    alpha_3 <- L$alpha_3
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    
    theta_true_f <- Vectorize(function(a) {
      
      lin <- function(w1,w2,a) {
        if (L$surv_true=="Cox PH") {
          if (L$dir=="decr") {
            C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1.7
          } else {
            C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*(1-a) - 1.7
          }
        } else if (L$surv_true=="complex") {
          if (L$dir=="decr") {
            C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 2.5*alpha_3*w2*a +
              0.7*alpha_3*(1-w2)*a - 1.3
          } else {
            C$alpha_1*pmax(0,2-8*abs(w1-0.5)) + 2.5*alpha_3*w2*(1-a) +
              0.7*alpha_3*(1-w2)*(1-a) - 1.3
          }
        }
      }
      
      S_0 <- function(t, w1, w2, a) {
        exp( -1 * lmbd * (t^v) * exp(lin(w1,w2,a)) )
      }
      
      return(1 - mean(S_0(C$t_e, w1, w2, a)))
      
    })
    
    deriv_theta_0 <- function (x) { grad(func=theta_true_f, x=x) }
  }
  
  # Esimate deriv_theta_n n_samples times
  deriv_ests <- c()
  grid <- seq(0,1,0.1)
  n_samples <- 3 # !!!!!
  for (i in 1:n_samples) {
    
    # Construct dat_orig and vlist
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    vlist <- create_val_list(dat, C$appx)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    
    # Construct additional component functions
    Psi_n <- Vectorize(function(x) {
      Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
    })
    gcm <- gcmlcm(x=seq(0,1,0.0001), y=Psi_n(seq(0,1,0.0001)), type="lcm")
    dGCM <- Vectorize(function(x) {
      if (x==0) {
        index <- 1
      } else {
        index <- which(round(x,6)<=gcm$x.knots)[1]-1
      }
      return(gcm$slope.knots[index])
    })
    theta_n_Gr <- function(x) { dGCM(Phi_n(x)) }
    f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                        type=params$g_n_type, k=15, delta1=TRUE)
    f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                    f_aIw_delta1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, type=params$gamma_type,
                                 vals=vlist$A_grid, omega_n=omega_n,
                                 f_aIw_n=f_aIw_n, f_a_n=f_a_n,
                                 f_a_delta1_n=f_a_delta1_n)
    
    theta_n <- theta_n_Gr
    deriv_theta_n <- construct_deriv_theta_n(theta_n, type=params$deriv_type,
                                             L$dir="decr")
    
    # Compute estimates
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
    ylim(c(0,0.1)) +
    labs(title="Estimation of derivative of theta",
         color="Which")
  
}



###############################.
##### omega_n and gamma_n #####
###############################.

{
  
  # Note: this check only works for surv_true=="Cox PH"
  
  # Create vlist
  vlist <- create_val_list(dat, C$appx)
  
  # True omega_0
  S_0 <- construct_S_n(dat, vlist$S_n, type="true")
  Sc_0 <- construct_S_n(dat, vlist$S_n, type="true", csf=TRUE)
  omega_0 <- Vectorize(function(w1,w2,a,y_star,delta_star) {
    
    # Shorten parameter variable names
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    lmbd2 <- L$sc_params$lmbd2
    v2 <- L$sc_params$v2
    alpha_3 <- L$alpha_3
    
    # Construct linear predictors
    if (L$dir=="decr") {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a - 1.7
    } else {
      lin <- C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*(1-a) - 1.7
    }
    lin2 <- C$alpha_1*w1 + C$alpha_2*w2 - 1
    
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
    
  })
  
  # Construct omega_n
  # S_n <- S_0
  # Sc_n <- Sc_0
  S_n <- construct_S_n(dat, vlist$S_n, type="Random Forest")
  Sc_n <- construct_S_n(dat, vlist$S_n, type="Random Forest", csf=TRUE)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  
  omegas_est <- c()
  omegas_true <- c()
  for (i in 1:nrow(vlist$omega)) {
    omegas_est <- c(omegas_est, do.call(omega_0, as.list(vlist$omega[i,])))
    omegas_true <- c(omegas_true, do.call(omega_n, as.list(vlist$omega[i,])))
  }
  
  highs <- c()
  for (i in c(1:nrow(vlist$omega))) {
    omg0 <- do.call(omega_0, as.list(vlist$omega[i,]))
    if (omg0>0.5) {
      highs <- c(highs, i)
    }
  }
  vlist$omega[highs,]
  # highs have delta_star=1, y_star<C$t_e, and A=0 or A>0.19
  
  # Scatterplot of estimated vs true
  ggplot(data.frame(x=omegas_true, y=omegas_est), aes(x=x, y=y)) +
    geom_abline(slope=1, intercept=0, color="green") +
    geom_point(alpha=0.1) +
    # xlim(-0.03,0.01)+ylim(-0.03,0.01)+
    # xlim(0.95,1.15)+ylim(0.95,1.15)+
    labs(title="omega_n vs omega_0", x="omega_0", y="omega_n")
  
  # Approximate true gamma_0
  construct_gamma_0 <- function() {
    epsilon <- 0.02
    dat_mc <- generate_data(n=50000, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, sampling="iid", L$dir)
    
    f_aIw_0 <- construct_f_aIw_n(dat, vlist$AW_grid, type="true", k=15)
    return(Vectorize(function(x) {
      i <- which(abs(dat_mc$a-x)<=epsilon)
      d <- list(
        w = dat_mc$w[i,],
        a = dat_mc$a[i],
        delta = dat_mc$delta[i],
        y_star = dat_mc$y_star[i],
        delta_star = dat_mc$delta_star[i]
      )
      return(mean(
        (omega_0(d$w,d$a,d$y_star,d$delta_star) /
           f_aIw_0(d$a,d$w))^2
      ))
    }))
  }
  gamma_0 <- construct_gamma_0()
  
  # Construct gamma_n
  f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid,
                                      type=params$g_n_type, delta1=TRUE, k=15)
  f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_delta1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, type="kernel",
                               vals=vlist$A_grid, omega_n=omega_n,
                               f_aIw_n=f_aIw_n, f_a_n=f_a_n,
                               f_a_delta1_n=f_a_delta1_n)
  gamma_n2 <- construct_gamma_n(dat_orig, dat, type="kernel2",
                               vals=vlist$A_grid, omega_n=omega_n,
                               f_aIw_n=f_aIw_n, f_a_n=f_a_n,
                               f_a_delta1_n=f_a_delta1_n)
  
  
  # Plot true curves against estimated curve
  grid <- seq(0,1,0.05)
  df <- data.frame(
    x = rep(grid, 3),
    y = c(gamma_n(grid), gamma_n2(grid), gamma_0(grid)),
    which = rep(c("Est (kernel)", "Est (kernel2)", "Truth"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) +
    geom_line() +
    # ylim(c(0,0.3)) +
    labs(title="Estimation of gamma_0", color="Which")
  
}



#########################################################.
##### Checking the influence function of Gamma_os_n #####
#########################################################.

{
  
  # Set up data
  dat <- generate_data(
    n = 400, # 5000
    alpha_3 = 0,
    distr_A = "Unif(0,1)",
    edge = "none",
    surv_true = "Cox PH",
    sc_params = L$sc_params,
    sampling = "iid",
    dir = L$dir
  )
  
  # Define the statistic to bootstrap
  bootstat <- function(dat_orig, indices) {
    
    dat_orig <- dat_orig[indices,]
    
    # Construct component functions
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat, type=params$ecdf_type)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, omega_n, S_n, g_n)
    
    # Return the value of Gamma_os_n(0.5)
    return (Gamma_os_n(0.5))
    
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
      sc_params = L$sc_params,
      sampling = "iid",
      dir = L$dir
    )
    
    Gammas <- c(Gammas, bootstat(dat, c(1:length(dat$a))))
    
  }
  
  # Estimate variance using influence function
  {
    
    # Construct component functions
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat, type=params$ecdf_type)
    f_aIw_n <- construct_f_aIw_n(dat, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(S_n, Sc_n)
    gcomp_n <- construct_gcomp_n(dat_orig, vals_A_grid, S_n)
    eta_n <- construct_eta_n(dat, vals_AW_grid, S_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, omega_n, S_n, g_n)
    
    if (F) {
      
      plot1 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point() +
        geom_smooth(formula=y~x, method=lm)
      
      plot2 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point()
      model <- lm(po~a, data=filter(dat_orig,!is.na(a))) # !!!!!
      coeff1 <- as.numeric(model$coefficients)
      plot2$layers <- c({
        stat_function(fun = function(x) {
          coeff1[1] + coeff1[2]*x
        }, color="turquoise")
      }, plot2$layers)
      
      plot3 <- ggplot(dat_orig, aes(x=a, y=po)) + geom_point()
      model <- lm(po~a+I(a^2)+I(a^3), data=filter(dat_orig,!is.na(a))) # !!!!!
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
    infl_fn_Gamma <- construct_infl_fn_Gamma(omega_n, g_n, gcomp_n,
                                             eta_n, Gamma_os_n)
    
    # Estimate variance and SD
    var_hat <- mean((
      infl_fn_Gamma(x=0.5,dat$w,dat$y_star,dat$delta_star,dat$delta,dat$a)
    )^2)
    sd_hat <- sqrt(var_hat/length(dat$a))
    
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



#############################################.
##### Test statistic variance estimator #####
#############################################.

{
  
  # Generate data
  dat_orig <- generate_data(n=400, alpha_3=0, distr_A="Unif(0,1)", edge="none",
                            surv_true="Cox PH", sc_params=L$sc_params,
                            sampling="two-phase (72%)", dir=L$dir)
  
  # Define the test statistic
  test_stat <- function(dat_orig, indices) {
    
    # # Uncomment this if running bootstrap
    # dat_orig <- ss(dat_orig, indices)
    
    # Prep
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat, C$appx)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    # lambda_2 <- lambda(dat,2,Phi_n)
    # lambda_3 <- lambda(dat,3,Phi_n)
    lambda_2 <- 1/3
    lambda_3 <- 1/4
    f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
    f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
    g_n <- construct_g_n(f_aIw_n, f_a_n)
    S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * (
        lambda_2*(Phi_n(dat$a))^2 -
          lambda_3*Phi_n(dat$a)
      ) *
        Gamma_os_n(round(dat$a, -log10(C$appx$a)))
    )
    
    return (beta_n)
    
  }
  
  
  # See actual distribution of test statistic
  betas <- c()
  for (i in 1:50) {
    
    dat_orig <- generate_data(n=400, alpha_3=0, distr_A="Unif(0,1)", edge="none",
                              surv_true="Cox PH", sc_params=L$sc_params,
                              sampling="two-phase (72%)", dir=L$dir) # "iid" "two-phase (72%)"
    
    beta_n <- test_stat(dat_orig, c(1:length(dat_orig$delta)))
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
    )
  )$sd_n
  
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
  
  # Set C$appx$t_e to 10 in "Setup" section
  
  S_0 <- construct_S_n(dat, vlist$S_n, type="true")
  S_n_CoxPH <- construct_S_n(dat, vlist$S_n, type="Cox PH")
  S_n_SL <- construct_S_n(dat, vlist$S_n, type="Super Learner")
  # S_n_RF <- construct_S_n(dat, vlist$S_n, type="Random Forest")
  
  # Plot true curve against estimated curve
  times <- seq(0,200,10)
  n <- length(times)
  w_a <- as.data.frame(cbind(w1=rep(0.2,n), w2=rep(1,n)))
  w_b <- as.data.frame(cbind(w1=rep(0.5,n), w2=rep(1,n)))
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      S_n_CoxPH(t=times, w=w_a, a=rep(0.2,n)),
      S_n_CoxPH(t=times, w=w_b, a=rep(0.2,n)),
      S_n_CoxPH(t=times, w=w_a, a=rep(0.8,n)),
      S_n_CoxPH(t=times, w=w_b, a=rep(0.8,n)),
      S_n_SL(t=times, w=w_a, a=rep(0.2,n)),
      S_n_SL(t=times, w=w_b, a=rep(0.2,n)),
      S_n_SL(t=times, w=w_a, a=rep(0.8,n)),
      S_n_SL(t=times, w=w_b, a=rep(0.8,n)),
      S_0(t=times, w=w_a, a=rep(0.2,n)),
      S_0(t=times, w=w_b, a=rep(0.2,n)),
      S_0(t=times, w=w_a, a=rep(0.8,n)),
      S_0(t=times, w=w_b, a=rep(0.8,n))
    ),
    which = rep(c("Cox PH","SL","True S_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0.2,w2=1,a=0.2","w1=0.5,w2=1,a=0.2",
                     "w1=0.2,w2=1,a=0.8","w1=0.5,w2=1,a=0.8"),3),
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
  dat_orig <- generate_data(
    n = 3000,
    alpha_3 = L$alpha_3,
    distr_A = L$distr_A,
    edge = L$edge,
    surv_true = L$surv_true,
    sampling = L$sampling,
    dir = L$dir
  )
  vlist <- create_val_list(dat, C$appx)
  
  # Sc_n_CoxPH <- construct_S_n(dat, vlist$S_n, type="Cox PH", csf=TRUE)
  # Sc_n_RF <- construct_S_n(dat, vlist$S_n, type="Random Forest", csf=TRUE)
  Sc_0 <- construct_S_n(dat, vlist$S_n, type="true", csf=TRUE)
  Sc_n_CoxPH <- Sc_0
  Sc_n_RF <- Sc_0
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      Sc_n_CoxPH(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_n_CoxPH(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_n_CoxPH(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.8,n)),
      Sc_n_CoxPH(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.8,n)),
      Sc_n_RF(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_n_RF(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_n_RF(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.8,n)),
      Sc_n_RF(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.8,n)),
      Sc_0(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_0(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.2,n)),
      Sc_0(t=times, cbind(w1=rep(0.2,n), w2=rep(1,n)), a=rep(0.8,n)),
      Sc_0(t=times, cbind(w1=rep(0.5,n), w2=rep(1,n)), a=rep(0.8,n))
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
  
  # Obtain estimates
  for (i in c(1:10)) {
    
    # !!!!! START !!!!!
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true, # !!!!!
                              L$sc_params, L$sampling, L$dir) # !!!!!
    dat <- ss(dat_orig, which(dat_orig$delta==1)) # !!!!!
    vlist <- create_val_list(dat, C$appx) # !!!!!
    # !!!!! END !!!!!
    
    # dat_orig1 <- dat_orig
    ests <- est_curve(
      dat_orig = dat_orig,
      estimator = "Grenander",
      params = params,
      points = C$points,
      dir = L$dir,
      return_extra = c("Phi_n_inv", "deriv_theta_n", "f_a_n", "gamma_n", "Psi_n",
                       "omega_n", "f_aIw_n", "S_n", "gcm", "dGCM", "etastar_n")
    )
    
    assign(paste0("ests",i), ests)
    assign(paste0("dat_orig",i), dat_orig)
    
  }
  
  # Return results
  # Good: 1,2
  # Bad: 5,7
  {
    i=7
    ests<-eval(as.name(paste0("ests",i)))
    theta_true <- attr(dat_orig, "theta_true")
    theta_ests <- ests$est
    ci_lo <- ests$ci_lo
    ci_hi <- ests$ci_hi
    len <- length(C$points)
    
    # Plot theta_n (estimate vs. truth)
    plot_data <- data.frame(
      x = rep(C$points, 2),
      theta = c(theta_ests, theta_true),
      which = rep(c("Est","Truth"), each=len),
      ci_lo = c(ci_lo, theta_true),
      ci_hi = c(ci_hi, theta_true)
    )
    ggplot(plot_data, aes(x=x, y=theta, color=factor(which))) +
      geom_line() +
      # labs(color="Which", title=paste("ests",i)) +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.2,
        fill = NA,
        linetype = "dotted"
      )
    
    # Analyze intermediate objects
    grid <- round(seq(0,1,0.01),2)
    gcm <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots, ties="ordered")
    plot_data <- data.frame(
      x = rep(grid,3),
      y = c(ests$Psi_n(grid), gcm(grid), ests$dGCM(grid)),
      which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*theta_n)"), each=101)
    )
    ggplot(plot_data, aes(x=x, y=y, color=which)) +
      geom_line() +
      labs(title=paste("ests",i)) +
      theme(legend.position="bottom") + ylim(c(-0.55,0.1))
    
    # # omega_n, etastar_n, and S_n all look fine
    # ggplot(data.frame(x=grid, y=omega_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="omega_n (ests 7)")
    # ggplot(data.frame(x=grid, y=etastar_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="etastar_n (ests 7)")
    # ggplot(data.frame(x=grid, y=S_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="S_n (ests 7)")
    
  }
  
  # Analysis of ests 7
  
  # > quantile(dat$a)
  # 0%  25%  50%  75% 100% 
  # 0.15 0.43 0.49 0.56 0.83  
  
  # > grid[c(15,16)]
  # [1] 0.14 0.15
  # > ests$Psi_n(grid)[c(15,16)]
  # [1] -0.08357442  0.01014071
  
  # What is causing the shape of the curve on [0,0.14] ?????
  
  # piece_1 <- omega_n(dat$w,dat$a,dat$y_star,dat$delta_star) /
  #   f_aIw_n(dat$a,dat$w)
  # fnc <- function(x) {
  #   (1/n_orig) * sum(weights_i * (
  #     as.integer(a_i<=x) * piece_1 + etastar_n(rep(x,nrow(w_i)),w_i)
  #   ))
  # }
  # etastar_n(x,c(0,1))
  
}



##########################################.
##### Conditional density estimators #####
##########################################.

{
  
  # Construct vlist
  vlist <- create_val_list(dat, C$appx)
  
  # # !!!!!
  # R.utils::withTimeout({
  #   f_aIw_n_semi <- construct_f_aIw_n(dat, NA, type="binning", k=15)
  # }, timeout=5)
  
  # True conditional density function
  f_aIw_0 <- construct_f_aIw_n(dat, NA, type="true", k=15)
  f_aIw_n_semi <- construct_f_aIw_n(dat, NA, type="binning", k=15)
  # f_aIw_0 <- construct_f_aIw_n(dat, vlist$AW_grid, type="true", k=15)
  # f_aIw_n_semi <- construct_f_aIw_n(dat, vlist$AW_grid, type="binning", k=15)
  f_aIw_n_para <- f_aIw_0 # !!!!!
  
  # Generate plot data
  grid <- seq(0,1,0.01)
  f_aIw_models <- c("Truth", "Parametric", "Semiparametric")
  n_models <- length(f_aIw_models)
  len <- length(grid)
  plot_data <- data.frame(
    a = rep(grid, 4*n_models),
    density = c(
      sapply(grid, function(a) { f_aIw_0(a, w=c(0.2,0)) }),
      sapply(grid, function(a) { f_aIw_0(a, w=c(0.8,0)) }),
      sapply(grid, function(a) { f_aIw_0(a, w=c(0.2,1)) }),
      sapply(grid, function(a) { f_aIw_0(a, w=c(0.8,1)) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w=c(0.2,0)) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w=c(0.8,0)) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w=c(0.2,1)) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w=c(0.8,1)) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w=c(0.2,0)) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w=c(0.8,0)) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w=c(0.2,1)) }),
      sapply(grid, function(a) { f_aIw_n_semi(a, w=c(0.8,1)) })
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
  # f_a_0 <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_0)
  # f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n_semi)
  f_a_0 <- construct_f_a_n(dat_orig, NA, f_aIw_0)
  f_a_n <- construct_f_a_n(dat_orig, NA, f_aIw_n_semi)
  
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
    alpha_3 = -4,
    distr_A = "N(0.5,0.04)",
    edge = edge,
    surv_true = "Cox PH",
    sampling = L$sampling,
    dir = L$dir
  )
  
  # True propensity score function
  pi_0 <- construct_superfunc(
    function(w) {
      if (edge=="expit") {
        return(expit(w[1]+w[2]-3.3))
      } else if (edge=="complex") {
        return(0.84*w[2]*pmax(0,1-4*abs(w[1]-0.5)))
      }
    },
    vec=c(2)
  )
  
  # Construct estimators
  vlist <- create_val_list(dat, C$appx)
  pi_n_logistic <- construct_pi_n(dat, vlist$W_grid, type="logistic")
  pi_n_SL <- construct_pi_n(dat, vlist$W_grid, type="SL")
  
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
