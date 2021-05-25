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
