############################.
##### Define functions #####
############################.

# Estimates the integral \int_a^b h(x)dG(x) based on a Riemann sum approximation
# m is the number of rectangles in the Riemann sum
int_appx <- function(h, G, a, b, m=500) {
  i <- c(1:m)
  x1 <- (i*(b-a))/m
  x2 <- ((i-1)*(b-a))/m
  sum((G(a+x1)-G(a+x2))*h(a+x1))
}

# Estimates the integral \int_a^b h(x)dG(x)
# G should be a dataframe with columns `x` and `y`, where x is a vector of step
#     function jump points and y=G(x)
int_step <- function(h, G, a, b) {
  x1 <- G$x[-1]
  y1 <- G$y[-1]
  y2 <- G$y[-length(G$y)]
  sum(as.integer(x1>=a & x1<=b)*(y1-y2)*h(x1))
}

# Estimates the variance of the cumulative hazard at time t for a subject with
#     covariate vector z_0
cox_var <- function(dat, theta_hat, bh, z_0, t, tau=NA,
                    return_extra=NA, breslow=NA) {
  
  # theta_hat=as.numeric(coeffs);z_0=c(0.5,1,0.5);t=100;
  
  # Alias random variables
  WT <- dat$weights
  N <- sum(WT)
  Z_ <- t(matrix(c(dat$w$w1, dat$w$w2, dat$a), ncol=3))
  T_ <- dat$y_star
  D_ <- dat$delta_star
  lin <- as.numeric(t(theta_hat)%*%Z_)
  d <- dim(Z_)[1]
  
  S_2n <- memoise(function(x) {
    Y_ <- as.integer(T_>=x)
    res <- matrix(NA, nrow=d, ncol=d)
    for (i in c(1:d)) {
      for (j in c(1:d)) {
        if (!is.na(res[j,i])) {
          res[i,j] <- res[j,i]
        } else {
          res[i,j] <- (1/N)*sum(Y_*Z_[i,]*Z_[j,]*exp(lin))
        }
      }
    }
    return(res)
  })
  
  S_0n <- Vectorize(memoise(function(x) {
    Y_ <- as.integer(T_>=x)
    return((1/N)*sum(Y_*exp(lin)))
  }))
  
  S_1n <- memoise(function(x) {
    Y_ <- as.integer(T_>=x)
    return((1/N)*as.numeric(Z_ %*% (Y_*exp(lin))))
  })
  
  m_n <- function(x) {
    S_1n(x) / S_0n(x)
  }
  
  h <- function(x) {
    return((S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x)))
  }
  
  # Create set of event times
  t_ev <- T_[D_==1]
  
  # Create estimated information matrix (for an individual)
  I_tilde <- matrix(0, nrow=d, ncol=d)
  for (t_i in t_ev) {
    I_tilde <- I_tilde + h(t_i)
  }
  I_tilde <- (1/N)*I_tilde
  I_tilde_inv <- solve(I_tilde)
  
  # Create score function
  l_star <- function(z_i,delta_i,t_i) {
    vec_list <- lapply(t_ev, function(t_j) {
      (exp(sum(z_i*theta_hat))*as.integer(t_j<=t_i)*(z_i-m_n(t_j)))/S_0n(t_j)
    })
    delta_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+",vec_list)
  }
  l_tilde <- function(z_i,delta_i,t_i) {
    I_tilde_inv %*% l_star(z_i,delta_i,t_i)
  }
  
  # Create omega function
  omega_n <- (function(){
    
    cmp1 <- as.integer(t_ev<=t)/S_0n(t_ev)
    cmp2 <- matrix(NA, nrow=d, ncol=length(t_ev))
    for (j in c(1:length(t_ev))) {
      cmp2[,j] <- z_0 - m_n(t_ev[j])
    }
    piece_3a <- t((1/N)*(cmp2%*%cmp1))
    piece_4 <- exp(sum(z_0*theta_hat))
    
    return(memoise(function(z_i,delta_i,t_i) {
      piece_1 <- (delta_i*as.integer(t_i<=t)) / S_0n(t_i)
      piece_2 <- (1/N) * exp(sum(z_i*theta_hat)) * sum(
        sapply(t_ev, function(t_j) {
          as.integer(t_j<=min(t,t_i))/((S_0n(t_j))^2)
        })
      )
      piece_3 <- as.numeric(piece_3a %*% l_tilde(z_i,delta_i,t_i))
      return(piece_4*(piece_1-piece_2+piece_3))
    }))
    
  })()
  
  # Calculate Breslow estimator
  breslow <- Vectorize(memoise(function(x) {
    sum(sapply(t_ev, function(t_i) {
      as.integer(t_i<=x) / (sum(exp(lin)*as.integer(dat$y_star>=t_i)))
    }))
  }))
  
  # Calculate component estimator Q_n
  Q_n <- memoise(function(z_i,delta_i,t_i) {
    piece_1 <- (delta_i*as.integer(t_i<=t)) / S_0n(t_i)
    piece_2 <- exp(sum(z_i*theta_hat))
    piece_3 <- (1/N) * sum(sapply(t_ev, function(t_j) {
      as.integer(t_j<=min(t,t_i)) / (S_0n(t_j))^2
    }))
    return(piece_1-piece_2*piece_3)
  })
  
  # Note: omega_mean equals zero
  omega_mean <- (1/N) * sum(sapply(c(1:N), function(i) {
    omega_n(Z_[,i],D_[i],T_[i])
  }))
  var_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (omega_n(Z_[,i],D_[i],T_[i]))^2
  })) - ((1/N)*omega_mean^2)
  
  # NEW: Influence function of Breslow estiamtor
  infl_fn_2 <- function(z_i,delta_i,t_i) {
    pc_1 <- Q_n(z_i,delta_i,t_i)
    pc_2 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
      Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_hat))*breslow(T_[j])) *
        Q_n(Z_[,j],D_[j],T_[j])
    }))
    # pc_2 <- Reduce("+", lapply(c(1:N), function(j) { # !!!!!
    #   Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_hat))*breslow(T_[j])) * # !!!!! Removed "(1/N) * "
    #     Q_n(Z_[,j],D_[j],T_[j]) # !!!!!
    # })) # !!!!!
    pc_3 <- l_tilde(z_i,delta_i,t_i)
    return(pc_1-sum(pc_2*pc_3))
  }
  
  # # TEMP: Influence function components of Breslow estiamtor
  # infl_fn_2a_Bres <- function(z_i,delta_i,t_i) {
  #   pc_2 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
  #     Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_hat))*breslow(T_[j])) *
  #       Q_n(Z_[,j],D_[j],T_[j])
  #   }))
  #   # pc_2 <- Reduce("+", lapply(c(1:N), function(j) { # !!!!! Removed "(1/N) * "
  #   #   Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_hat))*breslow(T_[j])) * # !!!!!
  #   #     Q_n(Z_[,j],D_[j],T_[j]) # !!!!!
  #   # })) # !!!!!
  #   pc_3 <- l_tilde(z_i,delta_i,t_i)
  #   return(sum(pc_2*pc_3))
  # }
  # infl_fn_2b_Bres <- function(z_i,delta_i,t_i) {
  #   pc_1 <- Q_n(z_i,delta_i,t_i)
  #   return(pc_1)
  # }
  
  # NEW: Influence function of cumulative hazard estiamtor
  infl_fn_3 <- (function() {
    pc_1 <- exp(sum(theta_hat*z_0))
    pc_3 <- z_0 * breslow(t)
    return(function(z_i,delta_i,t_i) {
      pc_2 <- Q_n(z_i,delta_i,t_i)
      pc_4 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
        Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_hat))*breslow(T_[j])) *
          Q_n(Z_[,j],D_[j],T_[j])
      }))
      pc_5 <- l_tilde(z_i,delta_i,t_i)
      return(pc_1*(pc_2+sum((pc_3-pc_4)*pc_5)))
    })
  })()
  
  # NEW: variance estimate using influence function
  var_est2 <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_3(Z_[,i],D_[i],T_[i]))^2
  }))
  
  # Calculate variance of Breslow estimator (Lin formula)
  piece_1 <- (1/N) * Reduce("+", lapply(t_ev, function(t_i) {
    (as.integer(t_i<=t)*m_n(t_i)) / S_0n(t_i)
  }))
  piece_2 <- (1/N) * sum(sapply(t_ev, function(t_i) {
    as.integer(t_i<=t) / (S_0n(t_i))^2
  }))
  # sigma2n <- as.numeric(N*(t(piece_1)%*%I_tilde_inv%*%piece_1)+piece_2)/N
  sigma2n <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1)+piece_2)/N
  
  # # TEMP: variance component estimate
  # var_est2a <- (1/N^2) * sum(sapply(c(1:N), function(i) {
  #   (infl_fn_2a_Bres(Z_[,i],D_[i],T_[i]))^2
  # }))
  # var_est2b <- (1/N^2) * sum(sapply(c(1:N), function(i) {
  #   (infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2
  #   # (N*infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2 # !!!!!
  # }))
  # # sigma2na <- as.numeric(N*(t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  # sigma2na <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  # sigma2nb <- as.numeric(piece_2)/N
  
  # # NEW: Variance estimate of Breslow estimator using influence function
  # mean2 <- (1/N^2) * sum(sapply(c(1:N), function(i) {
  #   infl_fn_2(Z_[,i],D_[i],T_[i])
  # }))
  # sigma2n2 <- (1/N^2) * sum(sapply(c(1:N), function(i) {
  #   (infl_fn_2(Z_[,i],D_[i],T_[i]))^2
  # }))
  
  return(list(
    I_tilde_inv = I_tilde_inv,
    var_est = var_est,
    var_est2 = var_est2,
    sigma2n = sigma2n
    # sigma2n2 = sigma2n2,
    # mean2 = mean2,
    # var_est2a = var_est2a, # !!!!!
    # var_est2b = var_est2b, # !!!!!
    # sigma2na = sigma2na, # !!!!!
    # sigma2nb = sigma2nb # !!!!!
  ))
  
}



###############################.
##### Process sim results #####
###############################.

if (F) {
  
  sim %>% SimEngine::summarize(
    mean = list(
      list(name="se_bshz_MC", x="se_bshz_MC"),
      list(name="se_cmhz_Cox", x="se_cmhz_Cox"),
      list(name="se_cmhz_MC", x="se_cmhz_MC")
    ),
    sd = list(
      list(name="se_bshz_empr", x="est_bshz"),
      list(name="se_cmhz_empr", x="est_cmhz")
    ),
    coverage = list(
      list(name="cov_lp_MC", truth="true_lp", estimate="est_lp", se="se_lp_MC"),
      list(name="cov_exp_MC", truth="true_exp", estimate="est_exp", se="se_exp_MC"),
      list(name="cov_bshz_MC", truth="true_bshz", estimate="est_bshz", se="se_bshz_MC"),
      list(name="cov_cmhz_Cox", truth="true_cmhz", estimate="est_cmhz", se="se_cmhz_Cox"),
      list(name="cov_cmhz_MC", truth="true_cmhz", estimate="est_cmhz", se="se_cmhz_MC"),
      list(name="cov_w1_Cox", truth=sim$constants$alpha_1, estimate="est_w1_Cox", se="se_w1_Cox"),
      list(name="cov_w1_MC", truth=sim$constants$alpha_1, estimate="est_w1_Cox", se="se_w1_MC"),
      list(name="cov_w2_Cox", truth=sim$constants$alpha_2, estimate="est_w2_Cox", se="se_w2_Cox"),
      list(name="cov_w2_MC", truth=sim$constants$alpha_2, estimate="est_w2_Cox", se="se_w2_MC"),
      list(name="cov_a_Cox", truth=sim$levels$alpha_3, estimate="est_a_Cox", se="se_a_Cox"),
      list(name="cov_a_MC", truth=sim$levels$alpha_3, estimate="est_a_Cox", se="se_a_MC")
    )
  )
  
}



###################.
##### Archive #####
###################.

if (F) {
  
  # Estimate information matrix using score; should be consistent
  I_tilde2 <- matrix(0, nrow=d, ncol=d)
  for (i in c(1:N)) {
    score <- l_star(Z_[,i],D_[i],T_[i])
    I_tilde2 <- I_tilde + score %*% t(score)
  }
  I_tilde_inv2 <- solve(I_tilde2)
  return(list(
    I_tilde_inv = I_tilde_inv,
    I_tilde_inv2 = I_tilde_inv2,
    var_est = var_est
  ))
  
  # Calculate the baseline cumulative hazard manually
  cuml_haz_n <- function(t, w, a) {
    return(base_hz_n(t)*exp(cfs[["w1"]]*w[1]+cfs[["w2"]]*w[2]+cfs[["a"]]*a))
  }
  
  # Calculate the cumulative hazard via predict()
  newdata <- data.frame(y_star=14, delta_star=1, a=0, w1=0, w2=0)
  predict(model,newdata=newdata, type="expected", se.fit=F)
  x=survfit(model, newdata=newdata, stype=2)
  
  # !!!!! Testing new int function
  h <- function(x) { sqrt(x) }
  int_appx(h, G=base_hz_n, 0, 100)
  int_step(h, G=data.frame(x=bh$time, y=bh$hazard), 0, 100)
  
  # !!!!! Use true baseline hazard function
  H_0_true <- function(t) {
    # L$sc_params$lmbd * t^L$sc_params$v
    L$sc_params$lmbd*exp(-1.7) * t^L$sc_params$v # !!!!!
  }
  H_0_true(c(10,50,100,200))
  base_hz_n(c(10,50,100,200))
  
  t_ev <- sort(dat$y_star[dat$delta_star==1])
  breslow2 <- Vectorize(function(x) {
    sum(sapply(t_ev, function(t_i) {
      as.integer(t_i<=x) / (sum(exp(lin)*as.integer(dat$y_star>=t_i)))
    }))
    # sum(sapply(t_ev, function(x) {
    #   den <- sum(as.integer(dat$y_star>=s)*exp(lin)) # !!!!! Pre-calculate
    #   as.integer(x<=s) / den
    # }))
  })
  H_0_true(c(10,50,100,200,400,800))
  breslow(c(10,50,100,200,400,800))
  breslow2(c(10,50,100,200,400,800))
  base_hz_n(c(10,50,100,200,400,800))
  
  
}
