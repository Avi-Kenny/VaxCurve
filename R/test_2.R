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

test_2 <- function(dat, alt_type="incr", params) {
  
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
        
        ecdf_G <- construct_Phi_n(dat)
        
        return(
          function(x) { ecdf_G(x) }
        )
      }
    }
    
    # lambda function (returns a constant)
    lambda <- function(k, G, dat) {
      
      if (attr(dat, "sampling")=="iid") {
        
        return( mean((G(dat$a, dat))^k) )
        
      } else if (attr(dat, "sampling")=="two-phase") {
        
        n_orig <- nrow(dat)
        dat %<>% filter(!is.na(a))
        sum((G(dat$a, dat))^k) / n_orig
        # weights <- 1 / Pi(dat$y, dat$w1, dat$w2)
        
      }
      
      
      
      
      
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
  
  return(as.numeric(crit_val>0))
  
}
