
# Cross-validation wrapper
if (F) {
  
  # CV prep
  n_folds <- 5
  folds <- sample(cut(c(1:nrow(df)), breaks=n_folds, labels=FALSE))
  bws <- seq(2*C$appx$a, 0.3, length.out=30)
  
  best <- list(bw=999, sum_sse=999)
  for (bw in bws) {
    
    sum_sse <- 0
    for (i in c(1:n_folds)) {
      
      df_train <- df[-which(folds==i),]
      df_test <- df[which(folds==i),]
      
      ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal", bandwidth=bw)
      reg <- Vectorize(function(a) {
        index <- which.min(abs(a-ks$x))
        return(ks$y[index])
      })
      
      sum_sse <- sum_sse + sum((reg(df_test$a)-df_test$po)^2, na.rm=T)
      
    }
    
    # print(paste0("bw: ",bw,", sum_sse: ",sum_sse))
    if (sum_sse<best$sum_sse || best$sum_sse==999) {
      best$bw <- bw
      best$sum_sse <- sum_sse
    }
    
    
  }
  
  # Construct optimal function from true data
  ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal",
                bandwidth=best$bw)
  reg <- Vectorize(function(a) {
    index <- which.min(abs(a-ks$x))
    return(ks$y[index])
  })
  
  # !!!!!
  # x <- runif(100)
  # y <- sin(x*10)+rnorm(100,sd=0.1)+1
  # df=data.frame(a=x,po=y)
  print(best)
  grid <- seq(0,1,0.01)
  y_range <- c(0,1)
  ggplot(df, aes(x=a, y=po)) +
    geom_point() +
    geom_line(
      aes(x=x, y=y),
      data.frame(x=grid, y=reg(grid))
    ) +
    ylim(y_range) +
    labs(title=paste0("Regression (CV; ",y_range[1],"--",y_range[2],")"))
  
}

# Checking variance components
if (F) {
  
  sim %>% summarize(
    mean = list(
      list(x="psi1psi2_var", name="var_est_psi1psi2"),
      list(x="Psi_1_var_est", name="var_est_psi1"),
      list(x="Psi_2_var_est", name="var_est_psi2")
    ),
    var = list(
      list(x="psi1psi2", name="var_emp_psi1psi2"),
      list(x="Psi_1_est", name="var_emp_psi1"),
      list(x="Psi_2_est", name="var_emp_psi2")
    )
  )
  
}

# Debugging infl_fn_1
if (F) {
  
  ######################################################.
  ##### Old plot of estimated Var against true Var #####
  ######################################################.
  
  if (F) {
    # Var estimates (True Var is 0.00025)
    ggplot(data.frame(x=Psi_6_var_est), aes(x=x)) +
      geom_histogram(alpha=0.7) +
      geom_vline(xintercept=var(Psi_6_est), color="forestgreen", linetype="dashed")
    
    # Var estimates (True Var is ???)
    print(paste("Observed variance:", round(var(Psi_4_est),7)))
    print(paste("Estimated variance:", round(mean(Psi_4_var_est),7)))
  }
  
  
  
  #################.
  ##### Psi_6 #####
  #################.
  
  # IF constructor
  construct_infl_fn_6 <- function(Phi_n) {
    fnc <- function(a_i) {
      as.integer(a_i<=0.5) - Phi_n(0.5)
    }
    return(fnc) # return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_6_est <- rep(NA, n_reps)
  Psi_6_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF
    n_orig <- 1000
    a <- round(runif(n_orig), 3)
    Phi_n <- ecdf(a)
    
    # Construct estimator, IF, and variance estimator
    Psi_6_est[i] <- Phi_n(0.5)
    infl_fn_Psi_6 <- construct_infl_fn_6(Phi_n)
    Psi_6_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_6(a))^2)
    
    # Hypothesis test (Null: Phi_0(0.5)=0.5)
    Phi_0 <- 0.5
    test_stat <- (Psi_6_est[i]-Phi_0)^2/Psi_6_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_4 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_4 <- function(a,Phi_n,Gamma_0) {
    a_j <- a
    piece_1 <- Gamma_0(a_j)
    piece_2 <- 2*mean(Phi_n(a_j)*Gamma_0(a_j))
    fnc <- function(a_i) {
      Phi_n(a_i)*Gamma_0(a_i) + mean(as.integer(a_i<=a_j)*piece_1) - piece_2
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_4_est <- rep(NA, n_reps)
  Psi_4_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF
    n_orig <- 1000
    a <- round(runif(n_orig), 3)
    Phi_n <- ecdf(a)
    
    # Construct estimator, IF, and variance estimator
    Psi_4_est[i] <- mean(Phi_n(a)*Gamma_0(a))
    infl_fn_Psi_4 <- construct_infl_fn_4(a,Phi_n,Gamma_0)
    Psi_4_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_4(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.13335
    test_stat <- (Psi_4_est[i]-Phi_0)^2/Psi_4_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_7 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_7 <- function(a,Phi_n,Gamma_0,lambda_3n) {
    a_j <- a
    piece_1 <- (Phi_n(a_j))^2
    piece_2 <- mean(Phi_n(a_j)*Gamma_0(a_j))
    piece_3 <- Gamma_0(a_j)
    fnc <- function(a_i) {
      (3*mean(as.integer(a_i<=a_j)*piece_1) + (Phi_n(a_i))^3 - 6*lambda_3n) *
        piece_2 +
        lambda_3n * (
          mean(as.integer(a_i<=a_j)*piece_3) + Phi_n(a_i)*Gamma_0(a_i)
        )
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_7_est <- rep(NA, n_reps)
  Psi_7_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_3n <- mean((Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_7_est[i] <- mean(lambda_3n*Phi_n(a)*Gamma_0(a))
    infl_fn_Psi_7 <- construct_infl_fn_7(a,Phi_n,Gamma_0,lambda_3n)
    Psi_7_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_7(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.03333
    test_stat <- (Psi_7_est[i]-Phi_0)^2/Psi_7_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_8 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_8 <- function(a,Phi_n,Gamma_0,lambda_2n) {
    a_j <- a
    piece_1 <- Phi_n(a_j)
    piece_2 <- mean((Phi_n(a_j))^2*Gamma_0(a_j))
    piece_3 <- Phi_n(a_j)*Gamma_0(a_j)
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_1) + (Phi_n(a_i))^2 - 6*lambda_2n) * 
        piece_2 +
        lambda_2n *
        (2*mean(as.integer(a_i<=a_j)*piece_3) + (Phi_n(a_i))^2*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 100
  Psi_8_est <- rep(NA, n_reps)
  Psi_8_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_2n <- mean((Phi_n(a))^2)
    
    # Construct estimator, IF, and variance estimator
    Psi_8_est[i] <- mean(lambda_2n*(Phi_n(a))^2*Gamma_0(a))
    infl_fn_Psi_8 <- construct_infl_fn_8(a,Phi_n,Gamma_0,lambda_2n)
    Psi_8_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_8(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.031746
    test_stat <- (Psi_8_est[i]-Phi_0)^2/Psi_8_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_1 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  rho_n <- function(a,Phi_n,Gamma_0,x) {
    mean( (Phi_n(a))^x * Gamma_0(a) )
  }
  construct_infl_fn_1 <- function(a,Phi_n,Gamma_0,lambda_2n,lambda_3n) {
    a_j <- a
    rho_1 <- rho_n(a,Phi_n,Gamma_0,1)
    rho_2 <- rho_n(a,Phi_n,Gamma_0,2)
    # xi_01 <- construct_xi_n(a,Phi_n,Gamma_0,0,1)
    # xi_20 <- construct_xi_n(a,Phi_n,Gamma_0,2,0)
    # xi_10 <- construct_xi_n(a,Phi_n,Gamma_0,1,0)
    # xi_11 <- construct_xi_n(a,Phi_n,Gamma_0,1,1)
    piece_01 <- Gamma_0(a_j)
    piece_20 <- (Phi_n(a_j))^2
    piece_10 <- Phi_n(a_j)
    piece_11 <- Phi_n(a_j) * Gamma_0(a_j)
    
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_10)+(Phi_n(a_i))^2-6*lambda_2n)*rho_2 +
        lambda_2n*(2*mean(as.integer(a_i<=a_j)*piece_11)+(Phi_n(a_i))^2*Gamma_0(a_i)) -
        (3*mean(as.integer(a_i<=a_j)*piece_20)+(Phi_n(a_i))^3-6*lambda_3n)*rho_1 -
        lambda_3n*(mean(as.integer(a_i<=a_j)*piece_01)+Phi_n(a_i)*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_1_est <- rep(NA, n_reps)
  Psi_1_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_2n <- mean((Phi_n(a))^2)
    lambda_3n <- mean((Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_1_est[i] <- mean(
      lambda_2n*(Phi_n(a))^2*Gamma_0(a) - lambda_3n*Phi_n(a)*Gamma_0(a)
    )
    infl_fn_Psi_1 <- construct_infl_fn_1(a,Phi_n,Gamma_0,lambda_2n,lambda_3n)
    Psi_1_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_1(a))^2)
    
    # Hypothesis test
    Phi_0 <- -0.0015878
    test_stat <- (Psi_1_est[i]-Phi_0)^2/Psi_1_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  ################################.
  ##### Psi_1 (with weights) #####
  ################################.
  
  # Setup
  
  C <- list(points=seq(0,1,0.02), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10,w1=0.1,a=0.01))
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  rho_n <- function(a,Phi_n,Gamma_0,x) {
    mean( (Phi_n(a))^x * Gamma_0(a) )
  }
  construct_infl_fn_1 <- function(a,Phi_n,Gamma_0,lambda_2n,lambda_3n) {
    a_j <- a
    rho_1 <- rho_n(a,Phi_n,Gamma_0,1)
    rho_2 <- rho_n(a,Phi_n,Gamma_0,2)
    # xi_01 <- construct_xi_n(a,Phi_n,Gamma_0,0,1)
    # xi_20 <- construct_xi_n(a,Phi_n,Gamma_0,2,0)
    # xi_10 <- construct_xi_n(a,Phi_n,Gamma_0,1,0)
    # xi_11 <- construct_xi_n(a,Phi_n,Gamma_0,1,1)
    piece_01 <- Gamma_0(a_j)
    piece_20 <- (Phi_n(a_j))^2
    piece_10 <- Phi_n(a_j)
    piece_11 <- Phi_n(a_j) * Gamma_0(a_j)
    
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_10)+(Phi_n(a_i))^2-6*lambda_2n)*rho_2 +
        lambda_2n*(2*mean(as.integer(a_i<=a_j)*piece_11)+(Phi_n(a_i))^2*Gamma_0(a_i)) -
        (3*mean(as.integer(a_i<=a_j)*piece_20)+(Phi_n(a_i))^3-6*lambda_3n)*rho_1 -
        lambda_3n*(mean(as.integer(a_i<=a_j)*piece_01)+Phi_n(a_i)*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 10
  Psi_1_est <- rep(NA, n_reps)
  Psi_1_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    dat_orig <- generate_data(1000, 0, "Unif(0,1)", "none", "Cox PH",
                              list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
                              "iid", "decr") # "two-phase (72%)"
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    a <- dat$a
    weights <- dat$weights
    
    Phi_n <- construct_Phi_n(dat, type="step")
    lambda_2n <- (1/n_orig) * sum(weights*(Phi_n(a))^2)
    lambda_3n <- (1/n_orig) * sum(weights*(Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_1_est[i] <- (1/n_orig) * sum(weights*(
      lambda_2n*(Phi_n(a))^2*Gamma_0(a) - lambda_3n*Phi_n(a)*Gamma_0(a)
    ))
    infl_fn_Psi_1 <- construct_infl_fn_1(a,Phi_n,Gamma_0,lambda_2n,lambda_3n)
    Psi_1_var_est[i] <- (1/n_orig^2) * sum(weights*(infl_fn_Psi_1(a))^2)
    
    # Hypothesis test
    Phi_0 <- -0.0015878
    test_stat <- (Psi_1_est[i]-Phi_0)^2/Psi_1_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
}



# Testing infl_fn_1
if (F) {
  
  L <- list(n=1000, alpha_3=-2, dir="decr",
            sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
            distr_A="Unif(0,1)", edge="none", surv_true="Cox PH", # Unif(0,1) N(0.5,0.04)
            ecdf_type="true", sampling="iid", # two-phase (72%)
            estimator=list(est="Grenander",params=params)
  )
  
  n_reps <- 50
  ests <- rep(NA,n_reps)
  vars <- rep(NA,n_reps)
  for (i in c(1:n_reps)) {
    
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    vlist <- create_val_list(dat, C$appx)
    
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    lambda_2 <- lambda(dat,2,Phi_n)
    lambda_3 <- lambda(dat,3,Phi_n)
    xi_n <- construct_xi_n(Phi_n, lambda_2, lambda_3)
    rho_n <- function(x) { 0 }
    
    # Partial stat
    Gamma_0 <- Vectorize(function(a) {
      Theta_true <- attr(dat_orig,"Theta_true")
      grid <- seq(0,1,0.02)
      index <- which.min(abs(a-seq(0,1,0.02)))
      return(Theta_true[index])
    })
    partial_est <- (1/n_orig) * sum( weights * (
      (lambda_2*(Phi_n(dat$a))^2 - lambda_3*Phi_n(dat$a)) * Gamma_0(dat$a)
    ))
    
    # Variance estimate of partial stat
    infl_fn_1b <- construct_infl_fn_1(dat, Gamma_0, Phi_n, xi_n, rho_n,
                                      lambda_2, lambda_3)
    partial_var <- (1/n_orig^2) * sum((weights * (infl_fn_1b(dat$a)))^2)
    
    ests[i] <- partial_est
    vars[i] <- partial_var
    
    print(paste0("Rep ", i, " of ", n_reps))
    
  }
  
  # Process results
  print(paste0("sd(ests), iid:", sd(ests)))
  print(paste0("mean(sqrt(vars)), iid:", mean(sqrt(vars))))
  ggplot(data.frame(x=ests), aes(x=x)) + geom_histogram() + labs(title="ests")
  ggplot(data.frame(x=vars), aes(x=x)) + geom_histogram() + labs(title="var")
  
}

# Figure out why returned functions are so large
if (F) {
  
  S_n <- readRDS("705 (SL, marker 8)/S_n.rds")
  
  objs <- ls(environment(get("fnc",envir=environment(gamma_n))))
  for (obj in objs) {
    print(obj)
    print(object.size(get(
      obj,
      envir = environment(get("fnc",envir=environment(gamma_n)))
    )))
  }

}

# Debugging
if (F) {
  
  # Set up containers
  df <- data.frame(
    "rep" = integer(),
    "point" = double(),
    "theta_n" = double(),
    "Gamma_n" = double(),
    "Psi_n" = double()
  )
  df_list <- list()
  
  # Run for-loop
  for (i in 1:2) {
    
    print(paste0("Rep ",i,": ",Sys.time()))
    set.seed(i)
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat, C$appx)
    vlist$AW_grid <- NA; vlist$omega <- NA; vlist$W_grid <- NA;
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    Psi_n <- Vectorize(function(x) {
      -1 * Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
    })
    gcm <- gcmlcm(
      x = seq(0,1,C$appx$a),
      y = Psi_n(seq(0,1,C$appx$a)), # rev(Psi_n(seq(0,1,C$appx$a)))
      type = "gcm"
    )
    dGCM <- Vectorize(function(x) {
      # The round deals with a floating point issue
      index <- which(round(x,5)<=gcm$x.knots)[1]-1
      if (index==0) { index <- 1 }
      return(gcm$slope.knots[index])
    })
    theta_n_Gr <- Vectorize(function(x) { min(max(-1 * dGCM(Phi_n(x)),0),1) })
    
    for (j in 1:51) {
      df[nrow(df)+1,] <- c(
        i, points[j], theta_n_Gr(points[j]), Gamma_os_n(points[j]),
        Psi_n(points[j])
      )
    }
    
    df_list[[i]] <- list(
      dat_orig = dat_orig,
      Phi_n = Phi_n,
      Phi_n_inv = Phi_n_inv,
      S_n = S_n,
      Sc_n = Sc_n,
      g_n = g_n,
      omega_n = omega_n,
      Gamma_os_n = Gamma_os_n
    )
    
  }
  
  # Plot results
  df$Psi_n <- -1 * df$Psi_n
  ggplot(df, aes(x=point, y=Psi_n, group=rep)) + # Gamma_n Psi_n
    geom_line(alpha=0.4) +
    ylim(c(0,0.4)) # 0.4
  
  # !!!!!
  ggplot(
    data.frame(
      x = c(gcm$x.knots, seq(0,1,0.01)),
      y = c(gcm$y.knots, rev(Psi_n(seq(0,1,0.01)))),
      grp = c(rep("gcm",length(gcm$x.knots)),rep("Psi_n",101))
    ),
    aes(x=x, y=y, color=factor(grp))) +
    geom_line()
  
  # !!!!!
  identical(
    sim$results_complex$sim_uid_1$dat_orig,
    df_list[[1]]$dat_orig
  )
  sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.1))
  df_list[[1]]$Phi_n(seq(0,1,0.1))
  
  # !!!!!
  sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.1))
  ggplot(
    data.frame(
      x = seq(0,1,0.01),
      y = sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.01))
    ),
    aes(x=x, y=y)
  ) + geom_line()
  df_list[[1]]$Phi_n(seq(0,1,0.1))
  ggplot(
    data.frame(
      x = seq(0,1,0.01),
      y = df_list[[1]]$Phi_n(seq(0,1,0.01))
    ),
    aes(x=x, y=y)
  ) + geom_line()
  
  # return(function(a) { ptruncnorm(a, a=0, b=1, mean=0.5, sd=0.2) })
  
  
}

# Unit tests for superfunc
if (F) {
  
  f <- function(a,b,c) {
    if (T) { Sys.sleep(1) }
    (a-b)+sum(c)/5
  }
  fs <- construct_superfunc(f, vec=c(0,1,2))
  # f(9,2,c(10,40))
  # f(9,3,c(20,50))
  # f(9,4,c(30,60))
  # system.time({x<-fs(9,2,c(10,40)); print(x);})
  # system.time({x<-fs(9,3,c(20,50)); print(x);})
  # system.time({x<-fs(9,4,c(30,60)); print(x);})
  # system.time({x<-fs(9,c(2,3,4),data.frame(c(10,20,30),c(40,50,60))); print(x)})
  
  # Test rounding
  fs(9,2,c(10,40))
  fs(9,c(2,3,4),data.frame(c(10,20,30),c(40,50,60)))
  
}

# Deconstruct vectorize
if (F) {
  
  # # !!!!!
  # (function(){
  #   (function(){
  #     parent.frame()
  #   })()
  #   # parent.frame()
  # })()
  
  Vectorize2 <- function (FUN, vectorize.args = arg.names, SIMPLIFY = TRUE, 
            USE.NAMES = TRUE) 
  {
    arg.names <- as.list(formals(FUN))
    arg.names[["..."]] <- NULL
    arg.names <- names(arg.names)
    vectorize.args <- as.character(vectorize.args)
    if (!length(vectorize.args)) 
      return(FUN)
    if (!all(vectorize.args %in% arg.names)) 
      stop("must specify names of formal arguments for 'vectorize'")
    collisions <- arg.names %in% c("FUN", "SIMPLIFY", "USE.NAMES", 
                                   "vectorize.args")
    if (any(collisions)) 
      stop(sQuote("FUN"), " may not have argument(s) named ", 
           paste(sQuote(arg.names[collisions]), collapse = ", "))
    rm(arg.names, collisions)
    (function() {
      FUNV <- function() {
        print("parent.frame()")
        print(parent.frame())
        args <- lapply(as.list(match.call())[-1L], eval, 
                       parent.frame())
        print("args")
        print(args)
        names <- if (is.null(names(args))) 
          character(length(args))
        else names(args)
        print("names")
        print(names)
        dovec <- names %in% vectorize.args
        print("dovec")
        print(dovec)
        print("args[dovec]")
        print(args[dovec])
        print("args[!dovec]")
        print(args[!dovec])
        do.call("mapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
                            SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
      }
      formals(FUNV) <- formals(FUN)
      environment(FUNV) <- parent.env(environment())
      FUNV
    })()
  }
  
  x <- Vectorize2(function(a,b,c,d) {a+b+c+d}, vectorize.args=c("a","b"))
  x(1,2,3,4)
  
}

# Deconstruct memoise
if (F) {
  
  mem2 <- function (f) {
    memo_f <- function(...) {
      encl <- parent.env(environment())
      v <- encl$`_val`
    }
    memo_f_env <- new.env(parent=environment(f))
    memo_f_env$`_val` <- 3
    environment(memo_f) <- memo_f_env
    memo_f
  }
  
  f1 <- function() { print("there") }
  f2 <- mem2(f1)
  f2()
  
}

# construct_superfunc: fixing environments issue
if (F) {
  
  x <- 1
  f0 <- function() { return(x) }
  f1 <- memoise(f)
  f2 <- construct_superfunc(f)
  f3 <- Vectorize(f)
  
  f2 <- function() {
    x <- 99
    f()
    # f1()
  }
  f2()
  
  f2(x=4) # This should give same answer as f2b
  f2b <- function(x) {
    f1b <- function() { return(x) }
    f1b(a=3)
  }
  f2b(x=4)
  
  # !!!!! Try both passing function anonymously and named
  # !!!!! Test use of aux
  
}

# Smoothed ecdf
if (F) {
  
  # Create dataset
  n_orig <- 100
  dat <- data.frame(weights=rep(1,n_orig), a=runif(n_orig, min=0, max=1))
  
  # Calculate estimators
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
  
  vals_x_ <- c(vals_x[1], vals_x[1:(length(vals_x)-1)]+diff(vals_x)/2,
               vals_x[length(vals_x)])
  vals_y_ <- c(0, vals_y[1:(length(vals_y)-1)], 1)
  
  rval <- approxfun(vals_x_, vals_y_, method="linear", yleft=0,
                    yright=1, ties="ordered")
  
  rval_inv <- approxfun(vals_y_, vals_x_, method="linear", yleft=min(vals_x),
                        yright=max(vals_x), ties="ordered")
  
  grid <- seq(0,1,0.01)
  jitter <- 0
  # jitter <- rnorm(2*length(grid), sd=0.003)
  df <- data.frame(
    x = rep(grid,2),
    y = c(rval(grid), rval_inv(grid)) + jitter,
    which = rep(c("ecdf", "smoothed"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line() +
    geom_abline(slope=1, color="grey")
  
  rval(rval_inv(seq(0,1,0.1)))
  rval_inv(rval(seq(0,1,0.1)))
  
}

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
    alpha_3 = -3,
    distr_A = "N(1.5+w1,1.5+w2)", # "Unif(0,1)"
    edge = "none",
    surv_true = "complex",
    sc_params = L$sc_params,
    sampling = "two-phase (72%)",
    dir = "decr"
  )
  
  # Prep
  n_orig <- nrow(dat_orig)

  # Construct dataframes of values to pre-compute functions on
  vlist <- create_val_list(dat_orig, C$appx)
  
  # Construct component functions
  S_n <- construct_S_n(dat_orig, vlist$S_n, type="true")
  Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="true",
                        csf=TRUE)
  f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                               type="true", k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  g_n <- construct_g_n(f_aIw_n, f_a_n)
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
      alpha_3 = -3,
      distr_A = "N(0.5,0.04)",
      edge = edge,
      surv_true = "Cox PH",
      sc_params = L$sc_params,
      sampling = "two-phase (72%)",
      dir = "decr"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig, C$appx)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat, pi_n, S_n, omega_n, theta_os_n_est)
    
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
      alpha_3 = -3,
      distr_A = "N(0.5,0.04)",
      edge = edge,
      surv_true = "Cox PH",
      sc_params = L$sc_params,
      sampling = "two-phase (72%)",
      dir = "decr"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig, C$appx)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    pi_n <- construct_pi_n(dat, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat, pi_n, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat, pi_n, S_n, omega_n, theta_os_n_est)
    
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
    alpha_3 <- -3
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
  construct_Phi_n2 <- function(dat, type=params$ecdf_type) {
    dat <- cbind(dat, wts=wts(dat, scale="sum 1"))
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    s <- sum(dat$wts) / n
    Vectorize(function(x) { (1/(n*s))*sum(dat$wts*as.integer(dat$a<=x)) })
  }
  construct_Phi_n3 <- function (dat, which="ecdf", type=params$ecdf_type) {
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
  Phi_n2 <- construct_Phi_n2(dat, type=params$ecdf_type)
  Phi_n3 <- construct_Phi_n3(dat, type=params$ecdf_type)
  Phi_n_inv3 <- construct_Phi_n3(dat, which="inverse", type=params$ecdf_type)
  
  
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
