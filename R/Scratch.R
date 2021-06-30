
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
    dat <- cbind(dat, wts=wts(dat))
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    s <- sum(dat$wts) / n
    Vectorize(function(x) { (1/(n*s))*sum(dat$wts*as.integer(dat$a<=x)) })
  }
  construct_Phi_n3 <- function (dat, type="ecdf") {
    # Adaptation of ecdf() source code
    dat <- cbind(dat, wts=wts(dat))
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
  
  n <- 500
  shape <- c(0.3,0.3)
  sample <- rbeta(n, shape1=shape[1], shape2=shape[2])
  grid <- seq(0,1,0.01)
  
  wts <- c(
    rep(1,length(sample))
    # rep(0.5,length(sample)/2),
    # rep(1.5,length(sample)/2)
  )
  wts <- wts / sum(wts)
  kde_gauss2 <- density(
    x = sort(sample),
    kernel = "gaussian",
    from = 0,
    to = 1,
    weights = wts
  )
  kd_gauss2 <- function(x) {
    index <- which.min(abs(kde_gauss2$x - x))
    area <- mean(kde_gauss2$y)
    return(kde_gauss2$y[index]/area)
  }
  kd_gauss <- kdensity(
    x = sample,
    start = "gumbel",
    kernel = "gaussian",
    support = c(0,1)
  )
  kd_beta <- kdensity(
    x = sample,
    start = "gumbel",
    kernel = "beta",
    support = c(0,1)
  )
  
  # kd_beta2_ests <- kde.boundary(
  #   x = sample,
  #   boundary.kernel = "beta",
  #   w = wts
  #   # xmin = 0,
  #   # xmax = 1
  #   # eval.points = as.matrix(grid)
  # )
  # kd_beta2 <- function(x) {
  #   len <- length(kd_beta2_ests$eval.points)
  #   k_x <- kd_beta2_ests$eval.points[2:(len-1)]
  #   k_dens <- kd_beta2_ests$estimate[2:(len-1)]
  #   index <- which.min(abs(k_x - x))
  #   return(k_dens[index])
  # }

  # # Histogram of sample
  # ggplot(data.frame(x=sample), aes(x=x)) + geom_histogram()
  
  # Graph of true Beta distribution against estimates
  dens_true <- sapply(grid, function(x) {
    dbeta(x, shape1=shape[1], shape2=shape[2])
  })
  dens_est_gauss <- sapply(grid, kd_gauss)
  dens_est_gauss2 <- sapply(grid, kd_gauss2)
  dens_est_beta <- sapply(grid, kd_beta)
  # dens_est_beta2 <- sapply(grid, kd_beta2)
  dat_plot <- data.frame(
    x = rep(grid, 4),
    density = c(dens_true, dens_est_gauss, dens_est_beta, dens_est_gauss2),
    which = rep(c("True", "Est: Gaussian", "Est: Beta", "Est: Gaussian 2"),
                each=length(grid))
  )
  
  # # Check to see if densities are properly normalized (should appx equal 1)
  # mean(dens_est_beta[1:100])
  # # mean(dens_est_beta2[1:100])
  # mean(dens_est_3[1:100])
  # mean(dens_est_gauss[1:100])
  
  # Plot densities
  ggplot(dat_plot, aes(x=x, y=density, color=factor(which))) +
    geom_line() +
    labs(color="Which")
  
}
