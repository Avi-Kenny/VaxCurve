
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

