
###############################.
##### Process sim results #####
###############################.

if (F) {
  
  sim %>% SimEngine::summarize(
    mean = list(
      list(name="se_w1_MC", x="se_w1_MC"),
      list(name="se_w2_MC", x="se_w2_MC"),
      # list(name="se_bshz_MC", x="se_bshz_MC"),
      # list(name="se_cmhz_MC", x="se_cmhz_MC"),
      # list(name="se_surv_MC", x="se_surv_MC"),
      list(name="se_marg_MC", x="se_marg_MC")
    ),
    sd = list(
      list(name="se_w1_empr", x="est_w1"),
      list(name="se_w2_empr", x="est_w2"),
      # list(name="se_bshz_empr", x="est_bshz"),
      # list(name="se_cmhz_empr", x="est_cmhz"),
      # list(name="se_surv_empr", x="est_surv"),
      list(name="se_marg_empr", x="est_marg")
    ),
    coverage = list(
      list(name="cov_w1_MC", truth="true_w1", estimate="est_w1", se="se_w1_MC"),
      list(name="cov_w2_MC", truth="true_w2", estimate="est_w2", se="se_w2_MC"),
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
