############################.
##### Define functions #####
############################.

# Estimates the variance of the cumulative hazard at time t for a subject with
#     covariate vector z_0
cox_var <- function(dat, dat_orig, theta_n, z_0, t, tau=NA, return_extras=F,
                    calc_bshz=F, calc_surv=F) {
  
  # Alias random variables
  WT <- dat$weights
  N <- round(sum(WT))
  n <- length(WT)
  Z_ <- t(matrix(c(dat$w$w1, dat$w$w2, dat$a), ncol=3))
  T_ <- dat$y_star
  D_ <- dat$delta_star
  lin <- as.numeric(t(theta_n)%*%Z_)
  d <- dim(Z_)[1]
  
  S_2n <- memoise(function(x) {
    res <- matrix(NA, nrow=d, ncol=d)
    for (i in c(1:d)) {
      for (j in c(1:d)) {
        if (!is.na(res[j,i])) {
          res[i,j] <- res[j,i]
        } else {
          res[i,j] <- (1/N)*sum(WT*as.integer(T_>=x)*Z_[i,]*Z_[j,]*exp(lin))
        }
      }
    }
    return(res)
  })
  
  S_0n <- memoise(function(x) {
    (1/N) * sum(WT*as.integer(T_>=x)*exp(lin))
  })
  
  S_1n <- memoise(function(x) {
    (1/N)*as.numeric(Z_ %*% (WT*as.integer(T_>=x)*exp(lin)))
  })
  
  m_n <- function(x) {
    S_1n(x) / S_0n(x)
  }
  
  h <- function(x) {
    (S_2n(x)/S_0n(x)) - m_n(x) %*% t(m_n(x))
  }
  
  # Create set of event times
  i_ev <- which(D_==1)
  t_ev <- T_[which(D_==1)]
  
  # Create estimated information matrix (for an individual)
  I_tilde <- Reduce("+", lapply(i_ev, function(i) {
    WT[i] * h(T_[i])
  }))
  I_tilde <- (1/N)*I_tilde
  I_tilde_inv <- solve(I_tilde)
  
  # Create score function
  l_star <- function(z_i,delta_i,t_i) {
    delta_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*exp(sum(z_i*theta_n))*as.integer(T_[j]<=t_i)*(z_i-m_n(T_[j]))) /
        S_0n(T_[j])
    }))
    # vec_list <- lapply(t_ev, function(t_j) {
    #   (exp(sum(z_i*theta_n))*as.integer(t_j<=t_i)*(z_i-m_n(t_j)))/S_0n(t_j)
    # })
    # delta_i*(z_i-m_n(t_i)) - (1/N)*Reduce("+",vec_list)
  }
  l_tilde <- memoise(function(z_i,delta_i,t_i) {
    I_tilde_inv %*% l_star(z_i,delta_i,t_i)
  })
  
  # Create omega function
  omega_n <- (function(){
    
    piece_3b <- (1/N) * sum(sapply(i_ev, function(j) {
      (WT[j]*as.integer(T_[j]<=t)) / S_0n(T_[j])
    }))
    piece_3c <- (1/N) * Reduce("+", lapply(i_ev, function(j) {
      (WT[j]*as.integer(T_[j]<=t)*m_n(T_[j])) / S_0n(T_[j])
    }))
    
      return(memoise(function(z_i,delta_i,t_i,z) {
      piece_1 <- (delta_i*as.integer(t_i<=t)) / S_0n(t_i)
      piece_2 <- (1/N) * exp(sum(z_i*theta_n)) * sum(
        sapply(i_ev, function(j) {
          (WT[j]*as.integer(T_[j]<=min(t,t_i))) / ((S_0n(T_[j]))^2)
        })
      )
      piece_3a <- t(z*piece_3b-piece_3c)
      piece_3 <- as.numeric(piece_3a %*% l_tilde(z_i,delta_i,t_i))
      piece_4 <- exp(sum(z*theta_n))
      return(piece_4*(piece_1-piece_2+piece_3))
    }))
    
  })()
  
  # Calculate Breslow estimator
  Lambda_n <- Vectorize(memoise(function(x) {
    (1/N) * sum(sapply(i_ev, function(i) {
      (WT[i] * as.integer(T_[i]<=x)) / S_0n(T_[i])
    }))
  }))
  
  # Influence function of marginalized survival
  S_n <- memoise(function(w1,w2,a) {
    exp(-exp(sum(c(w1,w2,a)*theta_n))*Lambda_n(t))
  })
  infl_fn_marg <- memoise(function(w_i,a_i,delta_i,t_i,wt_i,a) {
    
    piece_1 <- S_n(w_i[1],w_i[2],a)
    piece_2 <- (1/N) * sum(sapply(c(1:N), function(j) {
      w_j <- c(dat_orig$w$w1[j], dat_orig$w$w2[j])
      S_n(w_j[1],w_j[2],a) * (
        wt_i * ifelse(wt_i==0, 0, omega_n(c(w_i,a_i),delta_i,t_i,c(w_j,a))) - 1
      )
    }))
    
    return(piece_1+piece_2)
    
  })
  
  # Variance estimate of marginalized survival
  var_marg_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_marg(
      w_i = c(dat_orig$w$w1[i], dat_orig$w$w2[i]),
      a_i = replace_na(dat_orig$a[i],0),
      delta_i = dat_orig$delta_star[i],
      t_i = dat_orig$y_star[i],
      wt_i = dat_orig$weight[i],
      a = 0.5 # !!!!! A=0.5 hard-coded for now
    ))^2
    # (infl_fn_marg(Z_[1:2,i],Z_[3,i],D_[i],T_[i],WT[i],0.5))^2 # !!!!! A=0.5 hard-coded for now
  }))
  
  res <- list(
    I_tilde_inv = I_tilde_inv,
    var_marg_est = var_marg_est
  )
  
  if (return_extras) {
    
    res$S_0n <- S_0n
    res$S_1n <- S_1n
    res$S_2n <- S_2n
    res$I_tilde <- I_tilde
    res$l_tilde <- l_tilde
    res$omega_n <- omega_n
    res$Lambda_n <- Lambda_n
    
  }
  
  if (calc_bshz) {
    
    # Calculate component estimator Q_n
    # !!!!! Adapt to two-phase sampling
    Q_n <- memoise(function(z_i,delta_i,t_i) {
      piece_1 <- (delta_i*as.integer(t_i<=t)) / S_0n(t_i)
      piece_2 <- exp(sum(z_i*theta_n))
      piece_3 <- (1/N) * sum(sapply(t_ev, function(t_j) {
        as.integer(t_j<=min(t,t_i)) / (S_0n(t_j))^2
      }))
      return(piece_1-piece_2*piece_3)
    })
    
    # Influence function of Breslow estiamtor
    # !!!!! Adapt to two-phase sampling
    infl_fn_2 <- function(z_i,delta_i,t_i) {
      pc_1 <- Q_n(z_i,delta_i,t_i)
      pc_2 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
        Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_n))*Lambda_n(T_[j])) *
          Q_n(Z_[,j],D_[j],T_[j])
      }))
      pc_3 <- l_tilde(z_i,delta_i,t_i)
      return(pc_1-sum(pc_2*pc_3))
    }
    
    # Variance estimate of Breslow estimator using influence function
    # !!!!! Adapt to two-phase sampling
    var_bshz_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
      (infl_fn_2(Z_[,i],D_[i],T_[i]))^2
    }))
    
    res$var_bshz_est <- var_bshz_est
    
  }
  
  if (calc_surv) {
    
    # Variance of cumulative hazard estimator
    # !!!!! Adapt to two-phase sampling
    var_cmhz_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
      (omega_n(Z_[,i],D_[i],T_[i],z_0))^2
    }))
    
    # Influence function of cumulative hazard estiamtor
    # !!!!! Adapt to two-phase sampling
    infl_fn_3 <- (function() {
      pc_1 <- exp(sum(theta_n*z_0))
      pc_3 <- z_0 * Lambda_n(t)
      return(function(z_i,delta_i,t_i) {
        pc_2 <- Q_n(z_i,delta_i,t_i)
        pc_4 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
          Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_n))*Lambda_n(T_[j])) *
            Q_n(Z_[,j],D_[j],T_[j])
        }))
        pc_5 <- l_tilde(z_i,delta_i,t_i)
        return(pc_1*(pc_2+sum((pc_3-pc_4)*pc_5)))
      })
    })()
    
    # Variance estimate (cumulative hazard) using influence function
    # !!!!! Adapt to two-phase sampling
    var_cmhz_est2 <- (1/N^2) * sum(sapply(c(1:N), function(i) {
      (infl_fn_3(Z_[,i],D_[i],T_[i]))^2
    }))
    
    # Variance estimate (survival) using influence function
    # !!!!! Adapt to two-phase sampling
    surv_est <- exp(-exp(sum(z_0*theta_n))*Lambda_n(t))
    var_surv_est <- (1/N^2) * sum(sapply(c(1:N), function(i) {
      (surv_est*infl_fn_3(Z_[,i],D_[i],T_[i]))^2
    }))
    
    res$var_cmhz_est <- var_cmhz_est # !!!!! var_cmhz_est or var_cmhz_est2 ?????
    res$var_surv_est <- var_surv_est
    
  }
  
  return(res)
  
}



###############################.
##### Process sim results #####
###############################.

if (F) {
  
  sim %>% SimEngine::summarize(
    mean = list(
      # list(name="se_bshz_MC", x="se_bshz_MC"),
      # list(name="se_cmhz_MC", x="se_cmhz_MC"),
      # list(name="se_surv_MC", x="se_surv_MC"),
      list(name="se_marg_MC", x="se_marg_MC")
    ),
    sd = list(
      # list(name="se_bshz_empr", x="est_bshz"),
      # list(name="se_cmhz_empr", x="est_cmhz"),
      # list(name="se_surv_empr", x="est_surv"),
      list(name="se_marg_empr", x="est_marg")
    ),
    coverage = list(
      # list(name="cov_bshz_MC", truth="true_bshz", estimate="est_bshz", se="se_bshz_MC"),
      # list(name="cov_cmhz_MC", truth="true_cmhz", estimate="est_cmhz", se="se_cmhz_MC"),
      # list(name="cov_surv_MC", truth="true_surv", estimate="est_surv", se="se_surv_MC"),
      list(name="cov_marg_MC", truth="true_marg", estimate="est_marg", se="se_marg_MC")
    )
  )
  
}



###################.
##### Archive #####
###################.

if (F) {
  
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
  
  # TEMP: Influence function components of Breslow estiamtor
  infl_fn_2a_Bres <- function(z_i,delta_i,t_i) {
    pc_2 <- (1/N) * Reduce("+", lapply(c(1:N), function(j) {
      Z_[,j] * (D_[j]-exp(sum(Z_[,j]*theta_n))*Lambda_n(T_[j])) *
        Q_n(Z_[,j],D_[j],T_[j])
    }))
    pc_3 <- l_tilde(z_i,delta_i,t_i)
    return(sum(pc_2*pc_3))
  }
  infl_fn_2b_Bres <- function(z_i,delta_i,t_i) {
    pc_1 <- Q_n(z_i,delta_i,t_i)
    return(pc_1)
  }
  
  # Calculate variance of Breslow estimator (Lin formula)
  piece_1 <- (1/N) * Reduce("+", lapply(t_ev, function(t_i) {
    (as.integer(t_i<=t)*m_n(t_i)) / S_0n(t_i)
  }))
  piece_2 <- (1/N) * sum(sapply(t_ev, function(t_i) {
    as.integer(t_i<=t) / (S_0n(t_i))^2
  }))
  var_bshz_est2 <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1)+piece_2)/N
  
  # TEMP: variance component estimate
  var_cmhz_est2a <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_2a_Bres(Z_[,i],D_[i],T_[i]))^2
  }))
  var_cmhz_est2b <- (1/N^2) * sum(sapply(c(1:N), function(i) {
    (infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2
    # (N*infl_fn_2b_Bres(Z_[,i],D_[i],T_[i]))^2 # !!!!!
  }))
  # var_bshz_est_a <- as.numeric(N*(t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  var_bshz_est_a <- as.numeric((t(piece_1)%*%I_tilde_inv%*%piece_1))/N
  var_bshz_est_b <- as.numeric(piece_2)/N
  
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
    var_cmhz_est = var_cmhz_est
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
  Lambda_n2 <- Vectorize(function(x) {
    sum(sapply(t_ev, function(t_i) {
      as.integer(t_i<=x) / (sum(exp(lin)*as.integer(dat$y_star>=t_i)))
    }))
  })
  H_0_true(c(10,50,100,200,400,800))
  Lambda_n(c(10,50,100,200,400,800))
  Lambda_n2(c(10,50,100,200,400,800))
  base_hz_n(c(10,50,100,200,400,800))
  
  
}
