#' Testing approach 2: regression slope
#' 
#' @param dat Data returned by generate_data()
#' @param alt_type Type of alternative hypothesis; either "incr" or "decr"
#' @param params A list, containing the following:
#'   - `G` Domain transformation; one of c("identity","marginal")
#'   - `bootreps` Number of bootstrap replicates to run
#' @return Binary; is null rejected (1) or not (0)
#' @notes
#'   - Note

test2 <- function(dat, alt_type="incr", params) {
  
  # !!!!! Options for the following:
  #   domain transfer by marginal of A (yes/no)
  #   iid vs two-phase sampling
  
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
    
    # Fit a regression and return an estimator function mu_n(a,w)
    # !!!!! Potentially memoise more complex regressions later
    construct_mu_n <- function(dat, type) {
      
      if (type=="logistic") {
        model <- glm(y~w1+w2+a, data=dat, family="binomial")
        coeffs <- as.numeric(summary(model)$coefficients[,1])
        
        return(function(a, w1, w2){
          expit( coeffs[1] + coeffs[2]*w1 + coeffs[3]*w2 + coeffs[4]*a )
        })
      }
      
      if (type=="smoothing spline") {
        # !!!!!
      }
      
    }
    
    # Estimate density ratio and return an estimator function g_n(a,w)
    # !!!!! type currently unused
    construct_g_n <- function(dat, type) {
      
      f_a <- kdensity(
        x = dat$a,
        start = "gumbel",
        kernel = "gaussian"
      )
      
      k0 <- kdensity(
        x = dplyr::filter(dat,w2==0)$a,
        start = "gumbel",
        kernel = "gaussian"
      )
      
      k1 <- kdensity(
        x = dplyr::filter(dat,w2==1)$a,
        start = "gumbel",
        kernel = "gaussian"
      )
      
      # !!!!! Modify this estimator
      # f_a_given_w <- cde(
      #   x = cbind(dat$w1,dat$w2),
      #   y = dat$a
      # )
      f_a_given_w <- function(a,w1,w2) {
        if (w2==0) { return(k0(a)) }
        if (w2==1) { return(k1(a)) }
      }
      
      return(
        memoise(Vectorize(function(a,w1,w2) {
          f_a_given_w(a,w1,w2) / f_a(a)
        }))
      )
      
    }
    
    # Construct Gamma_n estimator
    # This is the one-step estimator from Westling & Carone 2020
    construct_Gamma_n <- function(dat, mu_n, g_n) {
      
      subpiece_1a <- (dat$y - mu_n(dat$a,dat$w1,dat$w2)) /
        g_n(dat$a,dat$w1,dat$w2)
      
      n <- nrow(dat)
      i_long <- rep(c(1:n), each=n)
      j_long <- rep(c(1:n), times=n)
      a_long <- dat$a[i_long]
      w1_long <- dat$w1[j_long]
      w2_long <- dat$w2[j_long]
      subpiece_2a <- mu_n(a_long,w1_long,w2_long)
      
      # !!!!! Memoise?
      return(
        Vectorize(function(x) {
          
          subpiece_1b <- as.numeric(dat$a<=x)
          piece_1 <- mean(subpiece_1a*subpiece_1b)
          
          subpiece_2b <- as.numeric(a_long<=x)
          piece_2 <- mean(subpiece_2a*subpiece_2b)
          
          return(piece_1+piece_2)
          
        })
      )
      
    }
    
    # Construct influence function
    # !!!!! Vectorize? Memoise?
    # !!!!! Also option to use one-step or plug-in ?????
    construct_infl_fn <- function(sub_x, x) {
      
      # !!!!!
      
      return(999)
      
    }
    
  }
  
  # Run bootstrap
  {
    
    # Pre-calculate values
    {
      piece_3 <- mean(
        (
          lambda(2,G,dat)*(G(dat$a,dat))^2 -
            lambda(3,G,dat)*G(dat$a,dat)
        ) * Gamma_n(dat$a, dat)
      )
      
      pre <- list(
        piece_3 = piece_3
      )
      
    }
    
    # Define the statistic to bootstrap
    bootstat <- function(dat,indices) {
      
      dat_boot <- dat[indices,]
      
      # Run regression
      
      
      x <- array(1:24, c(4,3,2))
      dim1 <- c(1.1,2.2,3.3,4.4)
      dim2 <- c(11,22,33)
      dim3 <- c(0,1)
      x[which(dim1==2.2),which(dim2==33),which(dim3==0)]
      
      # Pre-calculate requisite functions
      {
        # !!!!!
      }
      
      # Expectation
      # !!!!! Pre-calculate this; should be a single function that looks up values from a lookup table
      piece_1 <- mean(
        apply(
          X = expand.grid(x=dat$a, x_prime=dat_boot$a),
          MARGIN = 1,
          FUN = function(r) {
            x <- r[["x"]]
            x_prime <- r[["x_prime"]]
            return(
              infl_fn(sub_x=x, x=x_prime) * (
                lambda(2,G,dat_boot)*(G(x,dat_boot))^2 -
                  lambda(3,G,dat_boot)*G(x,dat_boot)
              )
            )
          }
        )
      )
      
      # Psi(P_n^#, Gamma_n)
      piece_2 <- mean(
        (
          lambda(2,G,dat_boot)*(G(dat_boot$a,dat_boot))^2 -
            lambda(3,G,dat_boot)*G(dat_boot$a,dat_boot)
        ) * Gamma_n(dat_boot$a, dat)
      )
      
      # Psi(P_n, Gamma_n)
      # Pre-calculated and accessed via parent environment
      piece_3 <- pre$piece_3
      
      return (piece_1+piece_2+piece_3)
      
    }
    
    # Run bootstrap
    boot_obj <- boot(data=dat, statistic=bootstat, R=params$boot_reps)
    
    # Calculate critical value (for a one-sided test)
    crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
    # crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
    
  }
  
  # !!!!! OLD CODE; RECYCLE
  
  #' #' Hypothesis testing approach 2: slope
  #' #' 
  #' #' @param dat Data returned by generate_data_dr()
  #' #' @param params A list; `est` is the estimator used for Theta_hat; one of
  #' #'     c("glm", "sm spline")
  #' #' @return Binary; is null rejected (1) or not (0)
  #' test_app2_dr <- function(dat, params) {
  #'   
  #'   Theta_hat_constr <- function(dat, subtype) {
  #'     
  #'     if (subtype=="glm") {
  #'       
  #'       # Run model and extract coefficients
  #'       model <- glm(infected~bmi+sex+antib, data=dat, family="binomial")
  #'       coeff <- summary(model)$coefficients
  #'       alpha_0_hat <- coeff["(Intercept)",1]
  #'       alpha_1_hat <- coeff["bmi",1]
  #'       alpha_2_hat <- coeff["sex",1]
  #'       alpha_3_hat <- coeff["antib",1]
  #'       
  #'       # Create Theta_hat function
  #'       Theta_hat <- Vectorize(function(x) {
  #'         
  #'         t_i <- apply(
  #'           X = dat,
  #'           MARGIN = 1,
  #'           FUN = function(r) {
  #'             log( (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] +
  #'                   alpha_2_hat*r[["sex"]] + alpha_3_hat*x)) /
  #'                 (1+exp(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]]))
  #'             )
  #'           }
  #'         )
  #'         
  #'         return ((alpha_3_hat^-1)*mean(t_i))
  #'         
  #'       })
  #'       
  #'       return (Theta_hat)
  #'       
  #'     }
  #'     
  #'     if (subtype=="ss") {
  #'       
  #'       # Run model and extract coefficients
  #'       model <- gam(
  #'         infected ~ bmi + sex + s(antib, fx=FALSE, bs="cr", m=2, pc=0),
  #'         data = dat,
  #'         family = "binomial"
  #'       )
  #'       
  #'       coeff <- model$coefficients
  #'       alpha_0_hat <- coeff[["(Intercept)"]]
  #'       alpha_1_hat <- coeff[["bmi"]]
  #'       alpha_2_hat <- coeff[["sex"]]
  #'       spline_vals <- as.numeric(predict(
  #'         model,
  #'         newdata = list(bmi=rep(25,101), sex=rep(0,101), antib=seq(0,1,0.01)),
  #'         type = "terms"
  #'       )[,3])
  #'       
  #'       # Construct theta_hat from GAM formula
  #'       theta_hat <- Vectorize(function(x) {
  #'         E_hat_i <- apply(
  #'           X = dat,
  #'           MARGIN = 1,
  #'           FUN = function(r) {
  #'             expit(alpha_0_hat + alpha_1_hat*r[["bmi"]] + alpha_2_hat*r[["sex"]] +
  #'                   spline_vals[1:(round(100*round(x,2)+1,0))])
  #'           }
  #'         )
  #'         return (mean(E_hat_i))
  #'       })
  #'       
  #'       # Construct Theta_hat by integrating theta_hat
  #'       # Approximating integral using a Riemann sum with ten intervals
  #'       Theta_hat <- Vectorize(
  #'         function (a) { a * mean(theta_hat(seq(a/10,a,a/10))) }
  #'       )
  #'       
  #'       return (Theta_hat)
  #'       
  #'     }
  #'     
  #'   }
  #'   
  #'   # Construct Theta_hat function
  #'   Theta_hat <- Theta_hat_constr(dat, params$subtype)
  #' 
  #'   # Calculate value of test statistic
  #'   x <- dat$antib
  #'   mu_2n <- mean(x^2)
  #'   mu_3n <- mean(x^3)
  #'   beta_n <- mean((mu_2n*x^2 - mu_3n*x)*Theta_hat(x))
  #'   
  #'   # Define the statistic to bootstrap
  #'   bootstat <- function(dat,indices) {
  #'     d <- dat[indices,]
  #'     Theta_hat <- Theta_hat_constr(d, params$subtype)
  #'     x <- dat$antib
  #'     mu_2n <- mean(x^2)
  #'     mu_3n <- mean(x^3)
  #'     return (mean((mu_2n*x^2 - mu_3n*x)*(Theta_hat(x)-x)))
  #'   }
  #'   
  #'   # Run bootstrap
  #'   boot_obj <- boot(data=dat, statistic=bootstat, R=100) # 1000
  #'   
  #'   # Calculate critical value
  #'   # crit_val <- as.numeric(quantile(boot_obj$t, 0.05))
  #'   # !!!!! This is a one-sided test; either make this two-sided or do the Wald test above as a one-sided test
  #'   crit_val <- qnorm(0.05, mean=mean(boot_obj$t), sd=sd(boot_obj$t))
  #'   
  #'   return(as.numeric(crit_val>0))
  #'   
  #' }
  
  return(as.numeric(crit_val>0))
  
}
