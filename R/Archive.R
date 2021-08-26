
# Wald comparator test
if (F) {
  
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
  
}

# Old conditional variance function
if (F) {
  
  #' Construct conditional variance estimator function sigma^2_n(Y|A,W)
  #' 
  #' @param mu_n Dataset returned by generate_data()
  #' @param mu2_n Currently only "logistic"
  #' @return Conditional variane estimator function
  construct_sigma2_n <- function(mu_n, mu2_n) {
    
    return(memoise(Vectorize(function(a, w1, w2){
      mu2_n(a,w1,w2) - (mu_n(a,w1,w2))^2
    })))
    
  }
  
}

# Old visualizations (Bias, MSE, Coverage in one graph)
if (F) {
  
  # Export: 6" x 4"
  distr_A_ <- "Beta(0.9,1.1+0.4*w2)" # Unif(0,1) Beta(0.9,1.1+0.4*w2)
  estimand_ <- "Midpoint" # Midpoint Endpoint
  ggplot(
    filter(p_data, distr_A==distr_A_ & estimand==estimand_),
    aes(x=n, y=value, color=Estimator)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_point() +
    geom_line() +
    # geom_bar(stat="identity", position=position_dodge(),
    #          width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
                  scales=list(y=list(
                    Bias = scale_y_continuous(labels = percent_format()),
                    Coverage = scale_y_continuous(labels = percent_format()),
                    MSE = scale_y_continuous()
                  ))) +
    theme(legend.position="bottom") +
    # scale_color_manual(values=m_colors) +
    labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
         y=NULL, x=NULL, color="Estimator")
  
  
  # !!!!! Temp (export 7.5 x 4.5)
  distr_A_ <- "Beta(0.9,1.1+0.4*w2)" # Unif(0,1) Beta(0.9,1.1+0.4*w2)
  estimand_ <- "Endpoint" # Midpoint Endpoint
  ggplot(
    filter(p_data, distr_A==distr_A_ & estimand==estimand_),
    aes(x=sampling, y=value, fill=Estimator)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_bar(stat="identity", position=position_dodge(),
             width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
                  scales=list(y=list(
                    Bias = scale_y_continuous(labels = percent_format()),
                    Coverage = scale_y_continuous(labels = percent_format()),
                    MSE = scale_y_continuous()
                  ))) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors) +
    labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
         y=NULL, x=NULL, color="Estimator")
  
}

# Testing of regression estimators
if (F) {
  
  ##########################################.
  ##### TESTING: Regression estimators #####
  ##########################################.
  
  if (FALSE) {
    
    # Set levels here
    n <- 5000
    reg_true <- "Logistic" # Logistic GAM Complex
    sampling <- "iid" # iid two-phase
    
    # Generate data
    L <- list(alpha_3=0.7)
    C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3)
    dat <- generate_data(
      n = n,
      alpha_3 = 0.7,
      distr_A = "Unif(0,1)", # Unif(0,1) Beta(0.9,1.1+0.4*w2)
      reg_true = reg_true,
      sampling = sampling
    )
    
    # True regression function
    mu_0 <- function(a,w1,w2) {
      if (reg_true=="Logistic") {
        expit(-1.5 + 0.3*w1 + 0.7*w2 + 0.7*a)
      } else if (reg_true=="GAM") {
        expit(-1.5 + 0.3*w1 + 0.7*w2 + 0.7*sqrt(a))
      } else if (reg_true=="Complex") {
        expit(-1.5 + 0.3*sin(2*pi*w1) + 0.7*w2 + 0.7*sqrt(a) + -0.3*w1*w2)
      }
    }
    
    # Construct regression functions
    mu_n_logistic <- construct_mu_n(dat=dat, type="Logistic")
    mu_n_gam <- construct_mu_n(dat=dat, type="GAM")
    mu_n_rf <- construct_mu_n(dat=dat, type="Random forest")
    
    # Generate plot data
    grid <- seq(0,1,0.01)
    mu_models <- c("Truth", "Logistic", "Random forest", "GAM")
    n_models <- length(mu_models)
    len <- length(grid)
    plot_data <- data.frame(
      a = rep(grid, 4*n_models),
      mu = c(
        sapply(grid, function(a) { mu_0(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_0(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_0(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_0(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_logistic(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_rf(a, w1=0.8, w2=1) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.2, w2=0) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.8, w2=0) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.2, w2=1) }),
        sapply(grid, function(a) { mu_n_gam(a, w1=0.8, w2=1) })
      ),
      which = rep(mu_models, each=len*4),
      covariates = rep(c(
        rep("W1=0.2, W2=0",len),
        rep("W1=0.8, W2=0",len),
        rep("W1=0.2, W2=1",len),
        rep("W1=0.8, W2=1",len)
      ), n_models)
    )
    ggplot(plot_data, aes(x=a, y=mu, color=factor(which))) +
      geom_line() +
      facet_wrap(~covariates, ncol=2) +
      labs(color="Estimator", title="Estimation of regression: E[Y|W,A]")
    
  }

}

# Regression estimator
if (F) {
  
  #' Construct regression function mu_n
  #' 
  #' @param dat Dataset returned by generate_data(); accepts either full data or
  #'     truncated data
  #' @param type One of c("Logistic", "GAM", "Random forest")
  #' @param moment If moment=k, the regression E[Y^k|A,W] is estimated
  #' @return Regression function
  construct_mu_n <- function(dat, type, moment=1) {
    
    dat %<>% filter(!is.na(a))
    
    if (moment!=1) {
      dat %<>% mutate(y=y^moment)
    }
    
    if (type=="Logistic") {
      
      model <- glm(
        y~w1+w2+a,
        data = dat,
        family = "binomial",
        weights = wts(dat, scale="mean 1")
      )
      coeffs <- as.numeric(summary(model)$coefficients[,1])
      
      return(Vectorize(function(a, w1, w2){
        expit( coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2 + coeffs[4]*a )
      }))
      
    }
    
    if (type=="GAM") {
      
      k <- 4
      knots <- seq(0,1,1/(k+1))[2:(k+1)]
      ns_basis <- ns(dat$a, knots=knots, Boundary.knots=c(0,1))
      formula <- "y ~ w1 + w2"
      for (i in 1:(k+1)) {
        dat[paste0("b_",i)] <- ns_basis[,i]
        formula <- paste0(formula," + b_",i)
      }
      
      dat$wts <- wts(dat, scale="mean 1")
      model <- glm(
        formula,
        data = dat,
        family = "binomial",
        weights = wts
      )
      coeffs <- as.numeric(summary(model)$coefficients[,1])
      
      # Construct natural spline basis in advance over grid
      grid <- seq(0,1,0.01)
      nsb_grid <- ns(grid, knots=knots, Boundary.knots=c(0,1))
      
      return(Vectorize(function(a, w1, w2){
        index <- which.min(abs(a-grid))
        nsb_a <- as.numeric(nsb_grid[index,])
        lin_pred <- coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2
        for (i in 1:(k+1)) {
          lin_pred <- lin_pred + coeffs[i+3]*nsb_a[i]
        }
        return(expit(lin_pred))
      }))
      
    }
    
    if (type=="Random forest") {
      
      model <- ranger(
        y~w1+w2+a,
        data = dat,
        num.trees = 500,
        case.weights = wts(dat, scale="mean 1")
      )
      
      return(memoise(Vectorize(function(a, w1, w2){
        predict(model, data.frame(a=a, w1=w1, w2=w2))$predictions
      })))
    }
    
  }
  
}



# Chernoff realizations
if (F) {
  
  library(twostageTE)
  data(chernoff_realizations)
  
}

# Old deriv_theta_n estimator
if (F) {
  
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
  
}

# generate_data old regression code
if (F) {
  
  #' @param reg_true True functional form of the regression; one of the following:
  #'     - "Logistic": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*A)
  #'     - "GAM": E[Y|W,A]=expit(a0+a1*W1+a2*W2+a3*sqrt(A))
  #'     - "Complex": E[Y|W,A]=expit(a0+a1*sin(2*pi*W1)+a2*W2+a3*sqrt(A)+a4*W1*W2)
  
  # Compute true regression function (i.e. P(Y=1|W,A))
  if (reg_true=="Logistic") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*a)
  } else if (reg_true=="GAM") {
    probs <- expit(C$alpha_0 + C$alpha_1*w1 + C$alpha_2*w2 + alpha_3*sqrt(a))
  } else if (reg_true=="Complex") {
    probs <- expit(C$alpha_0 + C$alpha_1*sin(2*pi*w1) + C$alpha_2*w2 +
                     alpha_3*sqrt(a) + C$alpha_4*w1*w2)
  }
  
  # Sample outcome
  y <- rbinom(n, size=1, prob=probs)
  
}

# Old GAM code (using mgcv; slow predictions)
if (F) {
  
  model <- gam(
    y~w1+w2+s(a, bs="cr"),
    # y~w1+w2+s(a, bs="tp"),
    # y~w1+w2+s(a, k=5, bs="cr", fx=TRUE),
    data = dat,
    family = "binomial",
    weights = wts(dat, scale="mean 1")
  )

  return(memoise(Vectorize(function(a, w1, w2){
    predict.gam(model, list(a=a, w1=w1, w2=w2), type="response")
  })))
  
}

# Old construct_f_aIw_n function contents
if (F) {
  
  #' @param type One of c("KDE (Beta)", "Beta"). The former fits a KDE with Beta
  #'     kernels (Chen 1999). The latter fits Beta density via MLE.
  #' @notes
  #'   - We trim the first and last evaluation points (i.e. zero and one) because
  #'     the kde.boundary() function estimates the density as zero there; the
  #'     neighboring points are used as the estimates instead
  
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
      weights = wts(dat_0, scale="sum 1"),
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
      weights = wts(dat_1, scale="sum 1"),
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

# Old construct_f_a_n function contents
if (F) {
  
  # Run weighted KDE
  dat %<>% filter(!is.na(a))
  
  # KDE (Beta kernels)
  if (type=="KDE (Beta)") {
    
    kde <- kde.boundary(
      x = dat$a,
      boundary.kernel = "beta",
      w = wts(dat, scale="sum 1"),
      xmin = 0,
      xmax = 1
    )
    f_a_n <- function(x) {
      len <- length(kde$eval.points)
      k_x <- kde$eval.points[2:(len-1)] # Trimming off first and last point
      k_dens <- kde$estimate[2:(len-1)] # Trimming off first and last point
      index <- which.min(abs(k_x - x))
      return(k_dens[index])
    }
    
  }
  
  if (type=="Beta") {
    
    # # !!!!!
    # dat <- list()
    # dat$a <- rbeta(100, shape1=0.3, shape2=0.5)
    # # !!!!!
    
    # !!!!! Testing
    {
      # Beta(0.9,1.1+0.4*w2)
      # Generate true marginal distribution of A
      n <- 10000
      beta_samp_1 <- rbeta(n, shape1=0.9, shape2=1.1)
      beta_samp_2 <- rbeta(n, shape1=0.9, shape2=1.5)
      beta_samp <- c(beta_samp_1, beta_samp_2)
      ggplot(data.frame(x=beta_samp_1), aes(x=x)) + geom_histogram(bins=50)
      ggplot(data.frame(x=beta_samp_2), aes(x=x)) + geom_histogram(bins=50)
      ggplot(data.frame(x=beta_samp), aes(x=x)) + geom_histogram(bins=50)
      
    }
    
    dat %<>% filter(!is.na(a))
    n <- length(dat$a)
    
    # !!!!! Comparison
    Rfast::beta.mle(dat$a)
    
    # weights <- wts(dat, scale="sum 1") # !!!!!
    
    # Set the objective function (weighted likelihood)
    #   par[1] is alpha and par[2] is beta
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        loglik <- dbeta(dat$a[i], shape1=par[1], shape2=par[2], log=TRUE)
        wt <- 1 # !!!!! Testing
        # wt <- weights[i]
        return(-1*loglik*wt)
      }))
      
      return(sum_loglik)
      
    }
    optim(par=c(alpha=1, beta=1), fn=wlik)

  }
  
}

# Old GAM specification function
if (F) {
  
  if (mono_form=="identity") {
    mono_f <- function(x) {x}
  } else if (mono_form=="square") {
    mono_f <- function(x) {x^2}
  } else if (mono_form=="sqrt") {
    mono_f <- function(x) {sqrt(x)}
  } else if (mono_form=="step_0.2") {
    mono_f <- function(x) {as.integer(x>0.2)}
  } else if (mono_form=="step_0.8") {
    mono_f <- function(x) {as.integer(x>0.8)}
  } else {
    stop("mono_form incorrectly specified")
  }
  
}

# Old kernel density function
if (F) {
  
  # kde <- density(
  #   x = dat$a,
  #   kernel = "gaussian",
  #   weights = wts(dat),
  #   from = 0,
  #   to = 1
  # )
  # f_a_n <- function(x) {
  #   index <- which.min(abs(kde$x - x))
  #   area <- mean(kde$y)
  #   return(kde$y[index]/area)
  # }
  
  # f_a_n <- kdensity(
  #   x = dat$a,
  #   start = "gumbel",
  #   kernel = "beta", # gaussian
  #   support = c(0,1)
  # )
  
}

# Debugging code in construct_fns()
if (F) {
  mu_n_iid <- construct_mu_n(dat_iid, type="logistic")
  f_aIw_n_iid <- construct_f_aIw_n(dat_iid, type="simple")
  f_a_n_iid <- construct_f_a_n(dat_iid, type=NULL)
  g_n_iid <- construct_g_n(f_aIw_n_iid, f_a_n_iid)
  mu_n_tp <- construct_mu_n(dat_tp, type="logistic")
  f_aIw_n_tp <- construct_f_aIw_n(dat_tp, type="simple")
  f_a_n_tp <- construct_f_a_n(dat_tp, type=NULL)
  g_n_tp <- construct_g_n(f_aIw_n_tp, f_a_n_tp)
  Gamma_n_iid <- construct_Gamma_n(dat_iid, mu_n_iid, g_n_iid)
  Gamma_n_tp <- construct_Gamma_n(dat_tp, mu_n_tp, g_n_tp)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=function(a) {
    Gamma_n_iid(a)
  }) + ylim(0,0.5)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=function(a) {
    Gamma_n_tp(a)
  }) + ylim(0,0.5)
}



# Constructor functions for estimator componetns (hypothesis testing)
{
  
  # Set up component functions
  {
    # Set up G constructor
    construct_G <- function(params, dat) {
      if (params$G == "identity") {
        return(
          function(x) { x }
        )
      }
      if (params$G == "marginal") {
        ecdf_G <- ecdf(dat$a)
        return(
          function(x) { ecdf_G(x) }
        )
      }
    }
    
    # lambda function (returns a constant)
    lambda <- function(k, G, dat) {
      return( mean((G(dat$a, dat))^k) )
    }
    
    # Construct influence function
    # !!!!! Vectorize? Memoise?
    # !!!!! Also option to use one-step or plug-in ?????
    construct_infl_fn <- function(sub_x, x) {
      
      # !!!!!
      
      return(999)
      
    }
    
  }
    
}



# Test empirical CDF and empirical quantile functions
{
  dat <- list(a=c(0.6,0.8,0.9))
  Phi_n <- ecdf(dat$a)
  ggplot(data.frame(x=grid), aes(x=x)) +
    stat_function(fun = Phi_n) + ylim(0,1)
  Phi_n_inv <- function(x) {
    qemp(p=x, obs=dat$a, discrete=T)
  }
  ggplot(data.frame(x=grid), aes(x=x)) +
    stat_function(fun = Phi_n_inv) + ylim(0,1)
}



# Old functions
{
  #' Integral of expit function
  #' 
  #' @param x Numeric input
  #' @return Numeric output
  int_expit <- function(x) {log(1+exp(x))}
}



#' Hypothesis testing approach 2: slope (!!!!! OLD CODE !!!!!)
#'
#' @param dat Data returned by generate_data_dr()
#' @param params A list; `est` is the estimator used for Theta_hat; one of
#'     c("glm", "sm spline")
#' @return Binary; is null rejected (1) or not (0)
test_app2_dr <- function(dat, params) {

  Theta_hat_constr <- function(dat, subtype) {

    if (subtype=="glm") {

      # Run model and extract coefficients
      model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
      coeff <- summary(model)$coefficients
      alpha_0_hat <- coeff["(Intercept)",1]
      alpha_1_hat <- coeff["bmi",1]
      alpha_2_hat <- coeff["sex",1]
      alpha_3_hat <- coeff["antib",1]

      # Create Theta_hat function
      Theta_hat <- Vectorize(function(x) {

        t_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            log( (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] +
                  alpha_2_hat*r[["sex"]] + alpha_3_hat*x)) /
                (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]]))
            )
          }
        )

        return ((alpha_3_hat^-1)*mean(t_i))

      })

      return (Theta_hat)

    }

    if (subtype=="ss") {

      # Run model and extract coefficients
      model <- gam(
        infected ~ bmi + sex + s(antib, fx=FALSE, bs="cr", m=2, pc=0),
        data = dat,
        family = "binomial"
      )

      coeff <- model$coefficients
      alpha_0_hat <- coeff[["(Intercept)"]]
      alpha_1_hat <- coeff[["bmi"]]
      alpha_2_hat <- coeff[["sex"]]
      spline_vals <- as.numeric(predict(
        model,
        newdata = list(bmi=rep(25,101), sex=rep(0,101), antib=seq(0,1,0.01)),
        type = "terms"
      )[,3])

      # Construct theta_hat from GAM formula
      theta_hat <- Vectorize(function(x) {
        E_hat_i <- apply(
          X = dat,
          MARGIN = 1,
          FUN = function(r) {
            expit(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]] +
                  spline_vals[1:(round(100*round(x,2)+1,0))])
          }
        )
        return (mean(E_hat_i))
      })

      # Construct Theta_hat by integrating theta_hat
      # Approximating integral using a Riemann sum with ten intervals
      Theta_hat <- Vectorize(
        function (a) { a * mean(theta_hat(seq(a/10,a,a/10))) }
      )

      return (Theta_hat)

    }

  }

  # Construct Theta_hat function
  Theta_hat <- Theta_hat_constr(dat, params$subtype)

  # Calculate value of test statistic
  x <- dat$antib
  mu_2n <- mean(x^2)
  mu_3n <- mean(x^3)
  beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))

  # Define the statistic to bootstrap
  bootstat <- function(dat,indices) {
    d <- dat[indices,]
    Theta_hat <- Theta_hat_constr(d, params$subtype)
    x <- dat$antib
    mu_2n <- mean(x^2)
    mu_3n <- mean(x^3)
    return (mean((mu_2n*x^2 - mu_3n*x)*(Theta_hat(x)-x)))
  }

  # Run bootstrap
  boot_obj <- boot(data=dat, statistic=bootstat, R=100) # 1000

  # Calculate critical value
  # crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
  # !!!!! This is a one-sided test; either make this two-sided or do the Wald test above as a one-sided test
  crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))

  return(as.numeric(crit_val>0))

}

