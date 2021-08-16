
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
  p4 <- (1-hz3)*(1-hz2)*(1-hz1) # !!!!! Modification of Diaz and VDL
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
      wts_j <- dat$wts[indices]
      new_y_val <- sum(wts_j) / (n*s)
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
