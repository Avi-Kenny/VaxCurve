#' Logit function
#' 
#' @param x Numeric input
#' @return Numeric output
logit <- function(x) {
  log(x/(1-x))
}



#' Expit function
#' 
#' @param x Numeric input
#' @return Numeric output
expit <- function(x) {
  1 / (1+exp(-x))
}



#' Derivative of expit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_expit <- function(x) {
  exp(x) / ((1+exp(x))^2)
}



#' Derivative of logit function
#' 
#' @param x Numeric input
#' @return Numeric output
deriv_logit <- function(x) {
  1 / (x-x^2)
}



#' Construct IP weights
#' 
#' @param dat Dataset returned by generate_data()
#' @return A sum-to-one vector of weights
#' @notes
#'   - Pi is accessed globally
#'   - Weights are normalized to sum to one because some functions require this.
#'     However, this is NOT the same as the stabilized weights, since those
#'     account for the zero-weighted observations (i.e. the NAs in the `a`
#'     variable)
wts <- function(dat) {
  
  if (attr(dat, "sampling")=="iid") {
    weights <- rep(1, nrow(dat))
  } else if (attr(dat, "sampling")=="two-phase") {
    weights <- 1 / Pi(dat$y, dat$w1, dat$w2)
  }
  weights <- weights / sum(weights)
  
  return(weights)
  
}



#' Probability of sampling
#' 
#' @param y Vector `y` of dataset returned by generate_data()
#' @param w1 Vector `w1` of dataset returned by generate_data()
#' @param w2 Vector `w2` of dataset returned by generate_data()
#' @return A vector of probabilities of sampling
Pi <- function(y, w1, w2) {
  
  # Note: does not currently depend on w1
  return(expit(2*y-w2))
  
}



#' Construct regression function mu_n
#' 
#' @param dat Dataset returned by generate_data()
#' @param type Currently only "logistic"
#' @param moment If moment=k, the regression E[Y^k|A,W] is estimated
#' @return Regression function
# !!!!! Potentially memoise/grid-appx later
# !!!!! Add smoothing spline type and others (e.g. random forest)
construct_mu_n <- function(dat, type, moment=1) {
  
  if (moment!=1) {
    dat %<>% mutate(y=y^moment)
  }
  
  if (type=="logistic") {
    model <- glm(y~w1+w2+a, data=dat, family="binomial", weights=wts(dat))
    coeffs <- as.integer(summary(model)$coefficients[,1])
    
    return(function(a, w1, w2){
      expit( coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2 + coeffs[4]*a )
    })
  }
  
  if (type=="smoothing spline") {
    # !!!!!
    # !!!!! vectorize and memoise
  }
  
}



#' Construct derivative estimator theta'_n
#' 
#' @param dat Dataset returned by generate_data()
#' @param theta_n Estimator of theta_n
#' @param grid Grid over which to approximate function
#' # !!!!! This is a placeholder estimator; replace eventually
construct_deriv_theta_n <- function(dat, theta_n, grid) {
  
  # Estimate entire function on grid
  theta_ns <- sapply(grid, theta_n)
  if (theta_ns[1]==theta_ns[length(grid)]) {
    stop("theta_n is flat")
  }
  
  grid_width <- grid[2] - grid[1]
  points_x <- c(grid[1])
  points_y <- c(theta_ns[1])
  for (i in 2:length(grid)) {
    if (theta_ns[i]-theta_ns[i-1]!=0) {
      points_x <- c(points_x, grid[i]-(grid_width/2))
      points_y <- c(points_y, mean(c(theta_ns[i],theta_ns[i-1])))
    }
  }
  points_x <- c(points_x, grid[length(grid)])
  points_y <- c(points_y, theta_ns[length(grid)])
  points_sl <- c()
  for (i in 2:length(points_x)) {
    slope <- (points_y[i]-points_y[i-1]) /
      (points_x[i]-points_x[i-1])
    points_sl <- c(points_sl, slope)
  }
  
  deriv_theta_n <- Vectorize(function(x) {
    if (x==0) {
      index <- 1
    } else {
      index <- which(x<=points_x)[1]-1
    }
    points_sl[index]
  })
  
  return(deriv_theta_n)
  
}



#' Construct conditional variance estimator function sigma^2_n(Y|A,W)
#' 
#' @param mu_n Dataset returned by generate_data()
#' @param mu2_n Currently only "logistic"
#' @return Conditional variane estimator function
# !!!!! Potentially grid-appx later
construct_sigma2_n <- function(mu_n, mu2_n) {
  
  return(memoise(Vectorize(function(a, w1, w2){
    mu2_n(a,w1,w2) - (mu_n(a,w1,w2))^2
  })))
  
}



#' Construct estimator of marginal density f_A
#' 
#' @param dat Dataset returned by generate_data()
#' @param type Currently unused
#' @return Density estimator function
# !!!!! Potentially grid-appx later
# !!!!! Add other "types"
construct_f_a_n <- function(dat, type) {
  
  # Run weighted KDE
  dat %<>% filter(!is.na(a))
  kde <- density(
    x = dat$a,
    kernel = "gaussian",
    weights = wts(dat),
    from = 0,
    to = 1
  )
  f_a_n <- function(x) {
    index <- which.min(abs(kde$x - x))
    area <- mean(kde$y)
    return(kde$y[index]/area)
  }
  
  # f_a_n <- kdensity(
  #   x = dat$a,
  #   start = "gumbel",
  #   kernel = "beta", # gaussian
  #   support = c(0,1)
  # )
  
  return(memoise(Vectorize(f_a_n)))
  
}



#' Construct estimator of conditional density f_{A|W}
#' 
#' @param dat Dataset returned by generate_data()
#' @param type Currently only "simple"
#' @return Density estimator function
# !!!!! Potentially grid-appx later
# !!!!! Add other "types"
# !!!!! Modify using cde()
construct_f_aIw_n <- function(dat, type) {
  
  if (type=="simple") {
    
    dat_0 <- dplyr::filter(dat,w2==0 & !is.na(a))
    dat_1 <- dplyr::filter(dat,w2==1 & !is.na(a))
    
    # print(paste("sim_uid:", L$sim_uid))
    # print(paste("dat contains", nrow(dat), "rows.")) # !!!!!
    # print(paste("dat contains", nrow(filter(dat, !is.na(a))), "complete rows.")) # !!!!!
    # print(paste("dat_0 contains", nrow(dat_0), "rows.")) # !!!!!
    # print(paste("dat_1 contains", nrow(dat_1), "rows.")) # !!!!!

    # kd_0 <- kdensity(
    #   x = dat_0$a,
    #   start = "gumbel",
    #   kernel = "gaussian"
    # )
    kde_0 <- density(
      x = dat_0$a,
      kernel = "gaussian",
      weights = wts(dat_0),
      from = 0,
      to = 1
    )
    kd_0 <- function(x) {
      index <- which.min(abs(kde_0$x - x))
      area <- mean(kde_0$y)
      return(kde_0$y[index]/area)
    }
    
    # kd_1 <- kdensity(
    #   x = dat_1$a,
    #   start = "gumbel",
    #   kernel = "gaussian"
    # )
    kde_1 <- density(
      x = dat_1$a,
      kernel = "gaussian",
      weights = wts(dat_1),
      from = 0,
      to = 1
    )
    kd_1 <- function(x) {
      index <- which.min(abs(kde_1$x - x))
      area <- mean(kde_1$y)
      return(kde_1$y[index]/area)
    }
    
    f_aIw_n <- function(a,w1,w2) {
      if (w2==0) { return(kd_0(a)) }
      if (w2==1) { return(kd_1(a)) }
    }
    
    return(memoise(Vectorize(f_aIw_n)))
    
  }
  
}



#' Construct density ratio estimator g_n
#' 
#' @param f_aIw_n Conditional density estimator returned by construct_f_aIw_n
#' @param f_a_n Density estimator returned by construct_f_a_n
#' @return Density ratio estimator function
# !!!!! Potentially grid-appx later
# !!!!! Add other "types"
# !!!!! Modify f_a_given_w using cde()
construct_g_n <- function(f_aIw_n, f_a_n) {
  
  return(memoise(Vectorize(function(a,w1,w2) {
    f_aIw_n(a,w1,w2) / f_a_n(a)
  })))
  
}



#' Construct Gamma_n primitive estimator
#' 
#' @param dat Dataset returned by generate_data()
#' @param mu_n A regression function returned by construct_mu_n()
#' @param g_n A density ratio estimator function returned by construct_g_n()
#' @return Gamma_n estimator
#' @notes This is the one-step estimator from Westling & Carone 2020
# !!!!! Potentially grid-appx later
construct_Gamma_n <- function(dat, mu_n, g_n) {
  
  if (attr(dat, "sampling")=="iid") {
    
    n <- nrow(dat)
    i_long <- rep(c(1:n), each=n)
    j_long <- rep(c(1:n), times=n)
    a_long <- dat$a[i_long]
    w1_long <- dat$w1[j_long]
    w2_long <- dat$w2[j_long]
    
    subpiece_1a <- (dat$y - mu_n(dat$a,dat$w1,dat$w2)) /
      g_n(dat$a,dat$w1,dat$w2)
    subpiece_2a <- mu_n(a_long,w1_long,w2_long)
    
    return(
      memoise(Vectorize(function(x) {
        
        subpiece_1b <- as.integer(dat$a<=x)
        piece_1 <- mean(subpiece_1a*subpiece_1b)
        
        subpiece_2b <- as.integer(a_long<=x)
        piece_2 <- mean(subpiece_2a*subpiece_2b)
        
        return(piece_1+piece_2)
        
      }))
    )
    
  } else if (attr(dat, "sampling")=="two-phase") {
    
    n_orig <- nrow(dat)
    dat %<>% filter(!is.na(a))
    n <- nrow(dat)
    s <- (1/n_orig) * sum(1/Pi(dat$y, dat$w1, dat$w2))
    
    i_long1 <- rep(c(1:n), each=n)
    j_long1 <- rep(c(1:n), times=n)
    i_long2 <- rep(i_long1, each=n)
    j_long2 <- rep(j_long1, each=n)
    k_long2 <- rep(j_long1, times=n)
    
    a_i_l1 <- dat$a[i_long1]
    y_i_l1 <- dat$y[i_long1]
    y_j_l1 <- dat$y[j_long1]
    w1_i_l1 <- dat$w1[i_long1]
    w1_j_l1 <- dat$w1[j_long1]
    w2_i_l1 <- dat$w2[i_long1]
    w2_j_l1 <- dat$w2[j_long1]
    
    a_i_l2 <- dat$a[i_long2]
    y_i_l2 <- dat$y[i_long2]
    y_j_l2 <- dat$y[j_long2]
    w1_i_l2 <- dat$w1[i_long2]
    w1_j_l2 <- dat$w1[j_long2]
    w1_k_l2 <- dat$w1[k_long2]
    w2_i_l2 <- dat$w2[i_long2]
    w2_j_l2 <- dat$w2[j_long2]
    w2_k_l2 <- dat$w2[k_long2]
    
    # Subpiece 1 (one-level sum)
    subpiece_1a <- (dat$y - mu_n(dat$a,dat$w1,dat$w2)) /
      (s*Pi(dat$y, dat$w1, dat$w2) * g_n(dat$a,dat$w1,dat$w2))
    
    # Subpiece two (two-level sum)
    subpiece_2a <- (
      mu_n(a_i_l1,w1_j_l1,w2_j_l1) * (2 + 1/(s*Pi(y_j_l1, w1_j_l1, w2_j_l1)))
    ) / (s*Pi(y_i_l1, w1_i_l1, w2_i_l1))
    
    # Subpiece three (three-level sum)
    subpiece_3a <- ( 2 * mu_n(a_i_l2,w1_k_l2,w2_k_l2) ) / (
      s*Pi(y_i_l2, w1_i_l2, w2_i_l2) * Pi(y_j_l2, w1_j_l2, w2_j_l2)
    )

    return(
      memoise(Vectorize(function(x) {
        
        subpiece_1b <- as.integer(dat$a<=x)
        piece_1 <- mean(subpiece_1a*subpiece_1b)
        subpiece_2b <- as.integer(a_i_l1<=x)
        piece_2 <- mean(subpiece_2a*subpiece_2b)
        subpiece_3b <- as.integer(a_i_l2<=x)
        piece_3 <- mean(subpiece_3a*subpiece_3b)
        
        return(piece_1+piece_2+piece_3)
        
      }))
    )
    
  }
  
}



#' Construct Phi_n and Phi_n^{-1}
#' 
#' @param dat Dataset returned by generate_data()
#' @param type One of c("ecdf", "inverse").
#' @return CDF or inverse CDF estimator function
#' @notes
#'   - This uses stabilized IP weights
#'   - This accesses wts() globally
construct_Phi_n <- function (dat, type="ecdf") {
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



#' Hypothesis test based on logistic regression
#' 
#' @param dat Data returned by generate_data_dr()
#' @param params Unused
#' @return Binary; is null rejected (1) or not (0)
test_wald <- function(dat, alt_type="incr", params) {
  
  model <- glm(y~w1+w2+a, data=dat, family="binomial")
  one_sided_p <- pnorm(summary(model)$coefficients["a",3], lower.tail=F)
  reject <- as.integer(one_sided_p<0.05)
  
  return (reject)
  
}
