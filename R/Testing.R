
#################.
##### Setup #####
#################.

{
  
  # 1. Run beginning of MAIN.R
  
  # 2. Set global constants
  params <- list(
    Q_n_type="Cox PH", g_n_type="parametric", ci_type="regular", cf_folds=1,
    edge_corr="none", ecdf_type="linear (mid)", deriv_type="linear",
    gamma_type="Super Learner", omega_n_type="estimated",
    marg="Gamma_star2" # Gamma_star
  )
  C <- list(points=round(seq(0,1,0.02),2), alpha_1=0.5, alpha_2=0.7, t_0=200,
            appx=list(t_0=10, x_tol=25, s=0.01))
  L <- list(
    n=500, alpha_3=-2, dir="decr",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5), # sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-7, v2=1.5), # Uncomment this for (almost) no censoring
    distr_S="Unif(0,1)", edge="none", surv_true="Cox PH", # "N(0.3+0.4x2,0.04)"
    sampling="two-phase (50%)", estimator=list(est="Grenander",params=params) # two-phase (72%)
    # n=15000, alpha_3=-4, dir="decr",
    # sc_params=list(lmbd=3e-5, v=1.5, lmbd2=3e-5, v2=1.5), # sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-7, v2=1.5), # Uncomment this for (almost) no censoring
    # distr_S="N(0.5,0.04)", edge="none", surv_true="Cox PH",
    # sampling="two-phase (6%)", estimator=list(est="Grenander",params=params)
  )
  
  # 3. Generate dataset
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir) # wts_type="estimated"
  dat <- ss(dat_orig, which(dat_orig$z==1))
  vlist <- create_val_list(dat_orig)
  
}



##########################################.
##### Dataset descriptive statistics #####
##########################################.

{
  
  res <- function(title, ss, rnd) {
    print(title)
    print(paste0(round(mean(ss),rnd), " (", round(min(ss),rnd),
                 " -- ", round(max(ss),rnd), ")"))
  }
  
  # Set up counters
  pct_rc <- num_inf <- pct_inf <- num_tp <- c()
  pct_inf_s0 <- pct_inf_s3 <- pct_inf_s5 <- pct_inf_s7 <- pct_inf_s1 <- c()
  
  n_reps <- 10
  for (i in 1:n_reps) {
    
    # Generate dataset
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    
    # Calculate summary stats
    pct_rc <- c(pct_rc, mean(as.integer(dat_orig$delta==0)*
                                     as.integer(dat_orig$y<=200)))
    num_inf <- c(num_inf, sum(as.integer(dat_orig$delta==1)*
                                as.integer(dat_orig$y<=200)))
    pct_inf <- c(pct_inf, mean(as.integer(dat_orig$delta==1)*
                                 as.integer(dat_orig$y<=200)))
    num_tp <- c(num_tp, sum(dat_orig$z))
    pct_inf_s0 <- c(pct_inf_s0, attr(dat_orig,"r_M0")[1])
    pct_inf_s3 <- c(pct_inf_s3, attr(dat_orig,"r_M0")[16])
    pct_inf_s5 <- c(pct_inf_s5, attr(dat_orig,"r_M0")[26])
    pct_inf_s7 <- c(pct_inf_s7, attr(dat_orig,"r_M0")[36])
    pct_inf_s1 <- c(pct_inf_s1, attr(dat_orig,"r_M0")[51])
    
  }
  
  res("% right-censored by t_0", pct_rc, 3)
  res("# of infections by t_0", num_inf, 0)
  res("% of infections by t_0", pct_inf, 4)
  res("# in phase-2 sample", num_tp, 0)
  res("True inf rate (S=0)", pct_inf_s0, 5)
  res("True inf rate (S=0.3)", pct_inf_s3, 5)
  res("True inf rate (S=0.5)", pct_inf_s5, 5)
  res("True inf rate (S=0.7)", pct_inf_s7, 5)
  res("True inf rate (S=1.0)", pct_inf_s1, 5)
  
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
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    dat <- ss(dat_orig, which(dat_orig$z==1))
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
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                            L$sc_params, "iid", L$dir) # L$sampling
  dat <- ss(dat_orig, which(dat_orig$z==1))
  vlist <- create_val_list(dat_orig)
  Phi_n1 <- construct_Phi_n(dat, type="step")
  Phi_n2 <- ecdf(dat$s)
  round(Phi_n1(round(seq(0,1,0.1),1))-Phi_n2(round(seq(0,1,0.1),1)),10) # !!!!!
  
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
  grid <- round(seq(0,1,0.01),2)
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



######################.
##### deriv_r_Mn #####
######################.

{
  
  # Generate true deriv_thetas
  {
    library(numDeriv)
    
    m <- 10^5
    x1 <- sample(round(seq(0,1,0.1),1), size=m, replace=T)
    x2 <- rbinom(m, size=1, prob=0.5)
    alpha_3 <- L$alpha_3
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    
    r_M0_f <- Vectorize(function(s) {
      
      lin <- function(x1,x2,s) {
        if (L$surv_true=="Cox PH") {
          if (L$dir=="decr") {
            C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*s - 1.7
          } else {
            C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*(1-s) - 1.7
          }
        } else if (L$surv_true=="Complex") {
          if (L$dir=="decr") {
            C$alpha_1*pmax(0,2-8*abs(x1-0.5)) + 2.5*alpha_3*x2*s +
              0.7*alpha_3*(1-x2)*s - 1.3
          } else {
            C$alpha_1*pmax(0,2-8*abs(x1-0.5)) + 2.5*alpha_3*x2*(1-s) +
              0.7*alpha_3*(1-x2)*(1-s) - 1.3
          }
        }
      }
      
      Q_0 <- function(t, x1, x2, s) {
        exp( -1 * lmbd * (t^v) * exp(lin(x1,x2,s)) )
      }
      
      return(1 - mean(Q_0(C$t_0, x1, x2, s)))
      
    })
    
    deriv_theta_0 <- function (u) { grad(func=r_M0_f, u=u) }
    
  }
  
  # Esimate deriv_r_Mn n_samples times
  deriv_ests <- c()
  grid <- round(seq(0,1,0.1),1)
  n_samples <- 3
  for (i in 1:n_samples) {
    
    # Construct dat_orig and vlist
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=params$g_n_type, k=15)
    f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$S_grid, omega_n, Q_n, g_n)
    
    # Construct additional component functions
    Psi_n <- Vectorize(function(u) {
      Gamma_os_n(round(Phi_n_inv(u), -log10(C$appx$s)))
    })
    gcm <- gcmlcm(x=seq(0,1,0.0001), y=Psi_n(seq(0,1,0.0001)), type="lcm")
    dGCM <- Vectorize(function(u) {
      if (u==0) {
        index <- 1
      } else {
        index <- which(round(u,6)<=gcm$x.knots)[1]-1
      }
      return(gcm$slope.knots[index])
    })
    r_Mn_Gr <- function(u) { dGCM(Phi_n(u)) }
    f_sIx_z1_n <- construct_f_sIx_n(dat, vlist$SX_grid,
                                        type=params$g_n_type, k=15, z1=TRUE)
    f_s_z1_n <- construct_f_s_n(dat_orig, vlist$S_grid,
                                    f_sIx_z1_n)
    gamma_n <- construct_gamma_n(dat_orig, dat, type=params$gamma_type,
                                 vals=vlist$S_grid, omega_n=omega_n,
                                 f_sIx_n=f_sIx_n, f_s_n=f_s_n,
                                 f_s_z1_n=f_s_z1_n)
    
    r_Mn <- r_Mn_Gr
    deriv_r_Mn <- construct_deriv_r_Mn(r_Mn, type=params$deriv_type,
                                             L$dir="decr")
    
    # Compute estimates
    deriv_ests <- c(deriv_ests, deriv_r_Mn(grid))
    
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
  vlist <- create_val_list(dat_orig)
  
  # True omega_0
  srvSL <- construct_Q_n(dat, vlist$Q_n, type=params$Q_n_type)
  Q_0 <- srvSL$srv
  Qc_0 <- srvSL$cens
  omega_0 <- Vectorize(function(x1,x2,s,y,delta) {
    
    # Shorten parameter variable names
    lmbd <- L$sc_params$lmbd
    v <- L$sc_params$v
    lmbd2 <- L$sc_params$lmbd2
    v2 <- L$sc_params$v2
    alpha_3 <- L$alpha_3
    
    # Construct linear predictors
    if (L$dir=="decr") {
      lin <- C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*s - 1.7
    } else {
      lin <- C$alpha_1*x1 + C$alpha_2*x2 + alpha_3*(1-s) - 1.7
    }
    lin2 <- C$alpha_1*x1 + C$alpha_2*x2 - 1
    
    # Compute omega_0
    piece_1 <- exp(-1*lmbd*(C$t_0^v)*exp(lin))
    piece_2 <- (delta*as.integer(y<=C$t_0)) /
      exp(-1*(lmbd*y^v*exp(lin)+lmbd2*y^v2*exp(lin2)))
    integral <- integrate(
      function(t) {
        t^(v-1) * exp(lin+lmbd*t^v*exp(lin)+lmbd2*t^v2*exp(lin2))
      },
      lower = 0,
      upper = min(C$t_0,y)
    )$value
    piece_3 <- lmbd*v*integral
    return(piece_1*(piece_2-piece_3))
    
  })
  
  # Construct omega_n
  # Q_n <- Q_0
  # Qc_n <- Qc_0
  # Q_n <- construct_Q_n(dat, vlist$Q_n, type="Random Forest")
  # Qc_n <- construct_Q_n(dat, vlist$Q_n, type="Random Forest", csf=TRUE)
  omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n)
  
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
  # highs have delta=1, y<C$t_0, and S=0 or S>0.19
  
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
    dat_mc <- generate_data(n=50000, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                            L$sc_params, sampling="iid", L$dir)
    
    f_sIx_0 <- construct_f_sIx_n(dat, vlist$SX_grid, type="true", k=15)
    return(Vectorize(function(u) {
      i <- which(abs(dat_mc$s-u)<=epsilon)
      d <- list(
        x = dat_mc$x[i,],
        s = dat_mc$s[i],
        z = dat_mc$z[i],
        y = dat_mc$y[i],
        delta = dat_mc$delta[i]
      )
      return(mean(
        (omega_0(d$x,d$s,d$y,d$delta) /
           f_sIx_0(d$s,d$x))^2
      ))
    }))
  }
  gamma_0 <- construct_gamma_0()
  
  # Construct gamma_n
  f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=params$g_n_type, k=15)
  f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
  f_sIx_z1_n <- construct_f_sIx_n(dat, vlist$SX_grid,
                                      type=params$g_n_type, z1=TRUE, k=15)
  f_s_z1_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_z1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, type="kernel",
                               vals=vlist$S_grid, omega_n=omega_n,
                               f_sIx_n=f_sIx_n, f_s_n=f_s_n,
                               f_s_z1_n=f_s_z1_n)
  gamma_n2 <- construct_gamma_n(dat_orig, dat, type="kernel2",
                               vals=vlist$S_grid, omega_n=omega_n,
                               f_sIx_n=f_sIx_n, f_s_n=f_s_n,
                               f_s_z1_n=f_s_z1_n)
  
  # Plot true curves against estimated curve
  grid <- round(seq(0,1,0.05),2)
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
    distr_S = "Unif(0,1)",
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
    n_orig <- length(dat_orig$z)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat, type=params$ecdf_type)
    f_sIx_n <- construct_f_sIx_n(dat, type=params$g_n_type, k=15)
    f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(Q_n, Qc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, omega_n, Q_n, g_n)
    
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
      distr_S = "Unif(0,1)",
      edge = "none",
      surv_true = "Cox PH",
      sc_params = L$sc_params,
      sampling = "iid",
      dir = L$dir
    )
    
    Gammas <- c(Gammas, bootstat(dat, c(1:length(dat$s))))
    
  }
  
  # Estimate variance using influence function
  {
    
    # Construct component functions
    n_orig <- length(dat_orig$z)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    weights <- dat$weights
    G_n <- construct_Phi_n(dat, type=params$ecdf_type)
    f_sIx_n <- construct_f_sIx_n(dat, type=params$g_n_type, k=15)
    f_s_n <- construct_f_s_n(dat_orig, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(Q_n, Qc_n)
    gcomp_n <- construct_gcomp_n(dat_orig, vals_S_grid, Q_n)
    # eta_n <- construct_eta_n(dat, vals_SX_grid, Q_n) # ARCHIVED
    Gamma_os_n <- construct_Gamma_os_n(dat, omega_n, Q_n, g_n)
    
    if (F) {
      
      plot1 <- ggplot(dat_orig, aes(x=s, y=po)) + geom_point() +
        geom_smooth(formula=y~x, method=lm)
      
      plot2 <- ggplot(dat_orig, aes(x=s, y=po)) + geom_point()
      model <- lm(po~s, data=filter(dat_orig,!is.na(s))) # !!!!!
      coeff1 <- as.numeric(model$coefficients)
      plot2$layers <- c({
        stat_function(fun = function(u) {
          coeff1[1] + coeff1[2]*u
        }, color="turquoise")
      }, plot2$layers)
      
      plot3 <- ggplot(dat_orig, aes(x=s, y=po)) + geom_point()
      model <- lm(po~s+I(s^2)+I(s^3), data=filter(dat_orig,!is.na(s))) # !!!!!
      coeff <- as.numeric(model$coefficients)
      plot3$layers <- c({
        stat_function(fun = function(u) {
          coeff[1] + coeff[2]*u + coeff[3]*(u^2) + coeff[4]*(u^3)
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
      infl_fn_Gamma(u=0.5,dat$x,dat$y,dat$delta,dat$z,dat$s)
    )^2)
    sd_hat <- sqrt(var_hat/length(dat$s))
    
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
  dat_orig <- generate_data(n=400, alpha_3=0, distr_S="Unif(0,1)", edge="none",
                            surv_true="Cox PH", sc_params=L$sc_params,
                            sampling="two-phase (72%)", dir=L$dir)
  
  # Define the test statistic
  test_stat <- function(dat_orig, indices) {
    
    # # Uncomment this if running bootstrap
    # dat_orig <- ss(dat_orig, indices)
    
    # Prep
    n_orig <- length(dat_orig$z)
    dat <- ss(dat_orig, which(dat_orig$z==1))
    weights <- dat$weights
    
    # Construct dataframes of values to pre-compute functions on
    vlist <- create_val_list(dat_orig)
    
    # Construct component functions
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    # lambda_2 <- lambda(dat,2,Phi_n)
    # lambda_3 <- lambda(dat,3,Phi_n)
    lambda_2 <- 1/3
    lambda_3 <- 1/4
    f_sIx_n <- construct_f_sIx_n(dat, vlist$SX_grid, type=params$g_n_type, k=15)
    f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n)
    g_n <- construct_g_n(f_sIx_n, f_s_n)
    srvSL <- construct_Q_n(dat, vlist$Q_n, type=params$Q_n_type)
    Q_n <- srvSL$srv
    Qc_n <- srvSL$cens
    omega_n <- construct_omega_n(vlist$omega, Q_n, Qc_n)
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$S_grid, omega_n, Q_n, g_n)
    
    # Compute the test statistic
    beta_n <- (1/n_orig) * sum(
      weights * (
        lambda_2*(Phi_n(dat$s))^2 -
          lambda_3*Phi_n(dat$s)
      ) *
        Gamma_os_n(round(dat$s, -log10(C$appx$s)))
    )
    
    return (beta_n)
    
  }
  
  
  # See actual distribution of test statistic
  betas <- c()
  for (i in 1:50) {
    
    dat_orig <- generate_data(n=400, alpha_3=0, distr_S="Unif(0,1)", edge="none",
                              surv_true="Cox PH", sc_params=L$sc_params,
                              sampling="two-phase (72%)", dir=L$dir) # "iid" "two-phase (72%)"
    
    beta_n <- test_stat(dat_orig, c(1:length(dat_orig$z)))
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
      Q_n_type = params$Q_n_type,
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
  # Set C$appx$t_0 to 10 in "Setup" section
  Q_0 <- (construct_Q_n(dat, vlist$Q_n, type="true"))$srv
  Q_CoxPH <- (construct_Q_n(dat, vlist$Q_n, type="Cox PH"))$srv
  # Q_SL <- (construct_Q_n(dat, vlist$Q_n, type="Super Learner"))$srv
  Q_SL <- (construct_Q_n(dat, vlist$Q_n, type="Random Forest"))$srv
  
  # Plot true curve against estimated curve (as function of T)
  times <- round(seq(0,200,10))
  n <- length(times)
  x_a <- as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
  x_b <- as.data.frame(cbind(x1=rep(0.5,n), x2=rep(1,n)))
  df <- data.frame(
    time = rep(times, 12),
    survival = c(Q_CoxPH(t=times, x=x_a, s=rep(0.2,n)),
                 Q_CoxPH(t=times, x=x_b, s=rep(0.2,n)),
                 Q_CoxPH(t=times, x=x_a, s=rep(0.8,n)),
                 Q_CoxPH(t=times, x=x_b, s=rep(0.8,n)),
                 Q_SL(t=times, x=x_a, s=rep(0.2,n)),
                 Q_SL(t=times, x=x_b, s=rep(0.2,n)),
                 Q_SL(t=times, x=x_a, s=rep(0.8,n)),
                 Q_SL(t=times, x=x_b, s=rep(0.8,n)),
                 Q_0(t=times, x=x_a, s=rep(0.2,n)),
                 Q_0(t=times, x=x_b, s=rep(0.2,n)),
                 Q_0(t=times, x=x_a, s=rep(0.8,n)),
                 Q_0(t=times, x=x_b, s=rep(0.8,n))),
    which = rep(c("Cox PH","SL","True Q_0"), each=4*length(times)),
    covs = rep(rep(c("x1=0.2,x2=1,s=0.2","x1=0.5,x2=1,s=0.2",
                     "x1=0.2,x2=1,s=0.8","x1=0.5,x2=1,s=0.8"),3),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: Q_0[t|X,S]",
         color="Estimator")
  
  # Plot true curve against estimated curve (as function of S)
  a_grid <- round(seq(0,1,0.02),2)
  n <- length(a_grid)
  x_a <- as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
  df <- data.frame(
    s = rep(s_grid, 3),
    survival = c(Q_CoxPH(t=rep(C$t_0,n), x=x_a, s=s_grid),
                 Q_SL(t=rep(C$t_0,n), x=x_a, s=s_grid),
                 Q_0(t=rep(C$t_0,n), x=x_a, s=s_grid)),
    which = rep(c("Cox PH","SL","True Q_0"), each=n)
  )
  ggplot(df, aes(x=s, y=survival, color=which)) +
    geom_line() +
    labs(title="Estimation of conditional survival: Q_0[t|X,S]",
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
    distr_S = L$distr_S,
    edge = L$edge,
    surv_true = L$surv_true,
    sc_params = L$sc_params,
    sampling = L$sampling,
    dir = L$dir
  )
  vlist <- create_val_list(dat_orig)
  
  Qc_0 <- (construct_Q_n(dat, vlist$Q_n, type="true"))$cens
  Qc_n_Cox <- (construct_Q_n(dat, vlist$Q_n, type="Cox PH"))$cens
  Qc_n_SL <- (construct_Q_n(dat, vlist$Q_n, type="Super Learner"))$cens
  
  # Plot true curve against estimated curve
  times <- round(seq(0,200,10))
  n <- length(times)
  x_02 <- as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
  x_05 <- as.data.frame(cbind(x1=rep(0.5,n), x2=rep(1,n)))
  df <- data.frame(
    time = rep(times, 12),
    survival = c(
      Qc_n_Cox(t=times, x_02, s=rep(0.2,n)),
      Qc_n_Cox(t=times, x_05, s=rep(0.2,n)),
      Qc_n_Cox(t=times, x_02, s=rep(0.8,n)),
      Qc_n_Cox(t=times, x_05, s=rep(0.8,n)),
      Qc_n_SL(t=times, x_02, s=rep(0.2,n)),
      Qc_n_SL(t=times, x_05, s=rep(0.2,n)),
      Qc_n_SL(t=times, x_02, s=rep(0.8,n)),
      Qc_n_SL(t=times, x_05, s=rep(0.8,n)),
      Qc_0(t=times, x_02, s=rep(0.2,n)),
      Qc_0(t=times, x_05, s=rep(0.2,n)),
      Qc_0(t=times, x_02, s=rep(0.8,n)),
      Qc_0(t=times, x_05, s=rep(0.8,n))
    ),
    which = rep(c("Old","New","True Qc_0"), each=4*length(times)),
    covs = rep(rep(c("x1=0.2,x2=1,s=0.2","x1=0.5,x2=1,s=0.2",
                     "x1=0.2,x2=1,s=0.8","x1=0.5,x2=1,s=0.8"),3), each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional (censoring) survival: Qc_0[t|X,S]",
         color="Estimator")
  
}



#############################.
##### r_Mn estimator #####
#############################.

{
  
  # Obtain estimates
  for (i in c(1:10)) {
    
    # !!!!! START !!!!!
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true, # !!!!!
                              L$sc_params, L$sampling, L$dir) # !!!!!
    dat <- ss(dat_orig, which(dat_orig$z==1)) # !!!!!
    vlist <- create_val_list(dat_orig) # !!!!!
    # !!!!! END !!!!!
    
    # dat_orig1 <- dat_orig
    ests <- est_curve(
      dat_orig = dat_orig,
      estimator = "Grenander",
      params = params,
      points = C$points,
      dir = L$dir,
      return_extra = c("Phi_n_inv", "deriv_r_Mn", "f_s_n", "gamma_n", "Psi_n",
                       "omega_n", "f_sIx_n", "Q_n", "gcm", "dGCM", "etastar_n")
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
    r_M0 <- attr(dat_orig, "r_M0")
    theta_ests <- ests$est
    ci_lo <- ests$ci_lo
    ci_hi <- ests$ci_hi
    len <- length(C$points)
    
    # Plot r_Mn (estimate vs. truth)
    plot_data <- data.frame(
      x = rep(C$points, 2),
      theta = c(theta_ests, r_M0),
      which = rep(c("Est","Truth"), each=len),
      ci_lo = c(ci_lo, r_M0),
      ci_hi = c(ci_hi, r_M0)
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
      which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*r_Mn)"), each=101)
    )
    ggplot(plot_data, aes(x=x, y=y, color=which)) +
      geom_line() +
      labs(title=paste("ests",i)) +
      theme(legend.position="bottom") + ylim(c(-0.55,0.1))
    
    # # omega_n, etastar_n, and Q_n all look fine
    # ggplot(data.frame(x=grid, y=omega_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="omega_n (ests 7)")
    # ggplot(data.frame(x=grid, y=etastar_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="etastar_n (ests 7)")
    # ggplot(data.frame(x=grid, y=Q_n(grid)), aes(x=x, y=y)) +
    #   geom_line()+theme(legend.position="bottom")+labs(title="Q_n (ests 7)")
    
  }
  
  # Analysis of ests 7
  
  # > quantile(dat$s)
  # 0%  25%  50%  75% 100% 
  # 0.15 0.43 0.49 0.56 0.83  
  
  # > grid[c(15,16)]
  # [1] 0.14 0.15
  # > ests$Psi_n(grid)[c(15,16)]
  # [1] -0.08357442  0.01014071
  
  # What is causing the shape of the curve on [0,0.14] ?????
  
  # piece_1 <- omega_n(dat$x,dat$s,dat$y,dat$delta) /
  #   f_sIx_n(dat$s,dat$x)
  # fnc <- function(u) {
  #   (1/n_orig) * sum(weights_i * (
  #     as.integer(s_i<=u) * piece_1 + etastar_n(rep(u,nrow(x_i)),x_i)
  #   ))
  # }
  # etastar_n(u,c(0,1))
  
}



##########################################.
##### Conditional density estimators #####
##########################################.

{
  # True conditional density function
  f_sIx_0 <- construct_f_sIx_n(dat, NA, type="true", k=15)
  f_sIx_n_semi <- construct_f_sIx_n(dat, NA, type="binning", k=15)
  f_sIx_n_para <- construct_f_sIx_n(dat, NA, type="parametric")
  
  # Generate plot data
  grid <- round(seq(0,1,0.01),2)
  f_sIx_models <- c("Truth", "Parametric", "Semiparametric")
  n_models <- length(f_sIx_models)
  len <- length(grid)
  plot_data <- data.frame(
    s = rep(grid, 4*n_models),
    density = c(
      sapply(grid, function(s) { f_sIx_0(s, x=c(0.2,0)) }),
      sapply(grid, function(s) { f_sIx_0(s, x=c(0.8,0)) }),
      sapply(grid, function(s) { f_sIx_0(s, x=c(0.2,1)) }),
      sapply(grid, function(s) { f_sIx_0(s, x=c(0.8,1)) }),
      sapply(grid, function(s) { f_sIx_n_para(s, x=c(0.2,0)) }),
      sapply(grid, function(s) { f_sIx_n_para(s, x=c(0.8,0)) }),
      sapply(grid, function(s) { f_sIx_n_para(s, x=c(0.2,1)) }),
      sapply(grid, function(s) { f_sIx_n_para(s, x=c(0.8,1)) }),
      sapply(grid, function(s) { f_sIx_n_semi(s, x=c(0.2,0)) }),
      sapply(grid, function(s) { f_sIx_n_semi(s, x=c(0.8,0)) }),
      sapply(grid, function(s) { f_sIx_n_semi(s, x=c(0.2,1)) }),
      sapply(grid, function(s) { f_sIx_n_semi(s, x=c(0.8,1)) })
    ),
    which = rep(f_sIx_models, each=len*4),
    covariates = rep(c(
      rep("X1=0.2, X2=0",len),
      rep("X1=0.8, X2=0",len),
      rep("X1=0.2, X2=1",len),
      rep("X1=0.8, X2=1",len)
    ), n_models)
  )
  ggplot(plot_data, aes(x=s, y=density, color=factor(which))) +
    geom_line() +
    facet_wrap(~covariates, ncol=4) +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of conditional density: f(S|X)") +
    ylim(c(0,NA))
  
  # Check marginal density
  # f_s_0 <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_0)
  # f_s_n <- construct_f_s_n(dat_orig, vlist$S_grid, f_sIx_n_semi)
  f_s_0 <- construct_f_s_n(dat_orig, NA, f_sIx_0)
  f_s_n <- construct_f_s_n(dat_orig, NA, f_sIx_n_semi)
  
  # Generate plot data
  grid <- round(seq(0,1,0.01),2)
  f_s_models <- c("Truth", "Semiparametric")
  len <- length(grid)
  plot_data <- data.frame(
    s = rep(grid, 2),
    density = c(f_s_0(grid),f_s_n(grid)),
    which = rep(f_s_models, each=len)
  )
  ggplot(plot_data, aes(x=s, y=density, color=factor(which))) +
    geom_line() +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of marginal density: f(S)") +
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
    distr_S = "N(0.5,0.04)",
    edge = edge,
    surv_true = "Cox PH",
    sampling = L$sampling,
    dir = L$dir
  )
  
  # True propensity score function
  pi_0 <- construct_superfunc(
    function(x) {
      if (edge=="expit") {
        return(expit(x[1]+x[2]-3.3))
      } else if (edge=="Complex") {
        return(0.84*x[2]*pmax(0,1-4*abs(x[1]-0.5)))
      }
    },
    vec=c(2)
  )
  
  # Construct estimators
  vlist <- create_val_list(dat_orig)
  pi_n_logistic <- construct_pi_n(dat, vlist$X_grid, type="logistic")
  pi_n_SL <- construct_pi_n(dat, vlist$X_grid, type="SL")
  
  # Curve 1: X2=0
  pi_n_log_0 <- function(x1) { pi_n_logistic(x1,x2=0) }
  pi_n_SL_0 <- function(x1) { pi_n_SL(x1,x2=0) }
  pi_0_0 <- function(x1) { pi_0(x1,x2=0) }
  
  # Curve 2: X2=1
  pi_n_log_1 <- function(x1) { pi_n_logistic(x1,x2=1) }
  pi_n_SL_1 <- function(x1) { pi_n_SL(x1,x2=1) }
  pi_0_1 <- function(x1) { pi_0(x1,x2=1) }
  
  # Plot true curves against estimated curve
  grid <- round(seq(0,1,0.01),2)
  curves <- c("X2=0","X2=1")
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
