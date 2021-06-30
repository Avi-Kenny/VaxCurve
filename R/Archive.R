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
