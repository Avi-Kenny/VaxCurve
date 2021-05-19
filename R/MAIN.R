# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  setting = "pointwise est", # testing
  level_set_which = "level_set_1",
  run_or_update = "run",
  num_sim = 1000,
  pkgs = c("simba", "dplyr", "boot", "car", "mgcv", "kdensity",
           "memoise", "twostageTE"),
  pkgs_nocluster = c("ggplot2", "viridis"),
  parallel = "none",
  stop_at_error = FALSE
)

# Set cluster config
cluster_config <- list(
  js = "slurm",
  dir = "/home/akenny/VaxCurve"
)



#################.
##### SETUP #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Marco Carone/Project - monotone testing/z.monotest/R"))
  load_pkgs_local <- TRUE
} else {
  # Cluster
  setwd("z.monotest/R")
  load_pkgs_local <- FALSE
}

# Load packages (if running locally)
if (load_pkgs_local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    do.call("library", list(pkg))
  }
}

# Load simba + functions
library(simba)
source("generate_data.R")
source("test1.R")
source("test2.R")
if (cfg$setting %in% c("density", "doseresp")) {
  source(paste0("fns_", cfg$setting,".R"))
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("run") %in% c("first", "")) {
  
  if (cfg$setting=="density") {

    # Compare all methods
    level_set_1 <- list(
      n = c(10,20,30,40,50),
      test = list(
        "Z test" = list(type="test_sw_z", params=list(crit_val="simulated")),
        "L test" = list(type="test_sw_l", params=list(crit_val="simulated")),
        "App 2, p_star=U(0,1)" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="U(0,1)")
        ),
        "App 2, p_star=P_0" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="P_0")
        )
        # "slope (SS-adapted)" = list(type = "slope", params = list(subtype="SS-adapted")),
        # "slope (full bootstrap; quantile)" = list(type = "slope", params = list(subtype="full bootstrap", ci_type="quantile")),
        # "slope (mixed bootstrap; quantile)" = list(type = "slope", params = list(subtype="mixed bootstrap", ci_type="quantile")),
        # "cuml_incr (delta=1/2)" = list(type = "cuml_incr", params = list(delta=(1/2), wts=1)),
        # "cuml_incr (delta=1/3, equal weights)" = list(type = "cuml_incr",params = list(delta=(1/3), wts=c(0.5,0.5)))
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Different weights for cuml_incr
    level_set_2 <- list(
      n = c(10,30,50),
      test = list(
        "0.0" = list(type="test1", params=list(delta=1/3,wts=c(1.0,0.0))),
        "0.1" = list(type="test1", params=list(delta=1/3,wts=c(0.9,0.1))),
        "0.2" = list(type="test1", params=list(delta=1/3,wts=c(0.8,0.2))),
        "0.3" = list(type="test1", params=list(delta=1/3,wts=c(0.7,0.3))),
        "0.4" = list(type="test1", params=list(delta=1/3,wts=c(0.6,0.4))),
        "0.5" = list(type="test1", params=list(delta=1/3,wts=c(0.5,0.5))),
        "0.6" = list(type="test1", params=list(delta=1/3,wts=c(0.4,0.6))),
        "0.7" = list(type="test1", params=list(delta=1/3,wts=c(0.3,0.7))),
        "0.8" = list(type="test1", params=list(delta=1/3,wts=c(0.2,0.8))),
        "0.9" = list(type="test1", params=list(delta=1/3,wts=c(0.1,0.9))),
        "1.0" = list(type="test1", params=list(delta=1/3,wts=c(0.0,1.0)))
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Different delta values for cuml_incr
    level_set_3 <- list(
      n = c(10,30,50),
      test = list(
        "1/2" = list(type="test1", params=list(delta=1/2,wts="equal")),
        "1/3" = list(type="test1", params=list(delta=1/3,wts="equal")),
        "1/4" = list(type="test1", params=list(delta=1/4,wts="equal")),
        "1/5" = list(type="test1", params=list(delta=1/5,wts="equal")),
        "1/6" = list(type="test1", params=list(delta=1/6,wts="equal")),
        "1/7" = list(type="test1", params=list(delta=1/7,wts="equal")),
        "1/8" = list(type="test1", params=list(delta=1/8,wts="equal")),
        "1/9" = list(type="test1", params=list(delta=1/9,wts="equal")),
        "1/10" = list(type="test1", params=list(delta=1/10,wts="equal"))
      ),
      true_density = c("f(x)=1", "f(x)=2x", "f(x)=ke^x")
    )
    
    # Variants of approach 2 test
    level_set_4 <- list(
      n = c(5,10,15,20),
      test = list(
        "App 2, SS-adapted" = list(
          type = "test2",
          params = list(subtype="SS-adapted")
        ),
        "App 2, p_star=U(0,1)" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="U(0,1)")
        ),
        "App 2, p_star=P_0" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="P_0")
        )
      ),
      true_density = c("Uniform", "Linear", "Exponential", "Spline") # "Half-Unif"
    )
    
    # Illustration of "reversing" problem
    level_set_5 <- list(
      n = c(10,20,30),
      test = list(
        "App 2, p_star=U(0,1)" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="U(0,1)")
        ),
        "App 2, p_star=P_0" = list(
          type = "test2",
          params = list(subtype="asymptotic", p_star="P_0")
        )
      ),
      true_density = c("Exponential", "Exponential (decr)",
                       "Step", "Step (decr)")
    )
    
  }
  
  if (cfg$setting=="hazard") {
    #
  }
  
  if (cfg$setting=="regression") {
    
    # Compare four variants
    level_set_1a <- list(
      n = 30,
      alpha_3 = c(0,0.15,0.3),
      sigma = 0.1,
      mono_form = c("identity", "square", "step_0.2"),
      a_distr = list(
        "U" = list(shape1=0.3, shape2=0.3),
        "decr" = list(shape1=0.7, shape2=1.3),
        # "incr" = list(shape1=1.3, shape2=0.7),
        "unif" = list(shape1=1, shape2=1),
        "spike" = list(shape1=40, shape2=40)
      ),
      test = list(
        "App_2: var 1" = list(
          type = "test2",
          params = list(G="identity", P_star="uniform",
                          var="boot", bootreps=100)),
        "App_2: var 2" = list(
          type = "test2",
          params = list(G="marginal", P_star="uniform",
                        var="boot", bootreps=100)),
        "App_2: var 3" = list(
          type = "test2",
          params = list(G="identity", P_star="marginal",
                        var="boot", bootreps=100)),
        "App_2: var 4" = list(
          type = "test2",
          params = list(G="marginal", P_star="marginal",
                        var="boot", bootreps=100))
      )
    )
    
    # Compare four variants
    level_set_1b <- list(
      n = 30,
      alpha_3 = c(0,0.15,0.3),
      sigma = 0.1,
      mono_form = c("step_0.2", "step_0.5", "step_0.8"),
      a_distr = list(
        "decr" = list(shape1=0.7, shape2=1.3),
        "incr" = list(shape1=1.3, shape2=0.7),
        "unif" = list(shape1=1, shape2=1)
      ),
      test = list(
        "App_2: var 1" = list(
          type = "test2",
          params = list(G="identity", P_star="uniform",
                        var="boot", bootreps=100)),
        "App_2: var 2" = list(
          type = "test2",
          params = list(G="marginal", P_star="uniform",
                        var="boot", bootreps=100)),
        "App_2: var 3" = list(
          type = "test2",
          params = list(G="identity", P_star="marginal",
                        var="boot", bootreps=100)),
        "App_2: var 4" = list(
          type = "test2",
          params = list(G="marginal", P_star="marginal",
                        var="boot", bootreps=100))
      )
    )
    
    # Diagnose inflated alpha problem
    level_set_2 <- list(
      n = c(10,30,100,300,1000),
      alpha_3 = 0,
      sigma = 0.1,
      mono_form = "step_0.2",
      a_distr = list(
        "incr" = list(shape1=1.3, shape2=0.7)
      ),
      test = list(
        "App_2: var 1" = list(
          type = "test2",
          params = list(G="identity", P_star="uniform",
                        var="boot", bootreps=100))
      )
    )
    
  }
  
  if (cfg$setting=="doseresp") {
    
    # Compare all methods
    level_set_1 <- list(
      n = c(20,40),
      alpha_3 = c(0,0.7),
      mono_form = c("identity", "square", "step_0.2"),
      test = list(
        "Wald" = list(type="test_wald", params=NULL),
        "App_2: glm" = list(type="test2",
                            params=list(subtype="glm", boot_reps=100)), # 1000
        "App_2: ss" = list(type="test2",
                           params=list(subtype="ss", boot_reps=100)) # 1000
      ),
      sampling = list(
        "iid" = list(type="iid")
      )
    )
    
  }
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Use these commands to run on Slurm:
# sbatch --export=run='first',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-24000 --export=run='main',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.monotest' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

if (cfg$run_or_update=="run") {
  
  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        stop_at_error = cfg$stop_at_error,
        seed = 3,
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      
      # Add functions to simulation object
      sim %<>% add_creator(generate_data)
      sim %<>% add_method(test1)
      sim %<>% add_method(test2)
      
      # Code specific to "density" setting
      if (cfg$setting=="density") {
        
        # Add functions
        sim %<>% add_method(test_sw_z)
        sim %<>% add_method(test_sw_l)
        
        # Generate beta_n distribution and add as a simulation constant
        sim %<>% add_constants("beta_n_distr"=beta_n_distr(sim$levels$n))
        
        # Simulation script
        sim %<>% set_script(function() {
          dat <- generate_data(L$n, L$true_density)
          decr_densities <- c("Exponential (decr)", "Step (decr)")
          alt_type <- ifelse(L$true_density %in% decr_densities, "decr", "incr")
          reject <- do.call(L$test$type, list(dat, alt_type, L$test$params))
          return (list("reject"=reject))
        })
        
      }
      
      # Code specific to "regression" setting
      if (cfg$setting=="regression") {
        
        # Add functions
        # sim %<>% add_method(test_sw_z)
        
        # Simulation script
        sim %<>% set_script(function() {
          dat <- generate_data(L$n, L$alpha_3, L$sigma, L$mono_form, L$a_distr)
          reject <- do.call(L$test$type, list(dat, "incr", L$test$params))
          # return (list("reject"=reject))
          return (list(reject=reject$reject, z=reject$z,
                       beta_n=reject$beta_n, var_est=reject$var_est))
        })
        
      }
      
      # Code specific to "doseresp" setting
      if (cfg$setting=="doseresp") {
        
        # Add functions
        sim %<>% add_method(expit)
        sim %<>% add_method(intexpit)
        sim %<>% add_method(test_wald)
        
        # Add constants
        data(chernoff_realizations)
        sim %<>% add_constants(
          chern = chernoff_realizations
        )
        
        # Simulation script
        sim %<>% set_script(function() {
          dat <- generate_data_dr(L$n, L$alpha_3, L$mono_form, L$sampling)
          alt_type <- "incr" # Can change later
          reject <- do.call(L$test$type, list(dat, alt_type, L$test$params))
          return (list("reject"=reject))
        })
        
      }
      
    },
    
    main = {
      sim %<>% run()
    },
    
    last = {
      sim %>% summary() %>% print()
    },
    
    cluster_config = cluster_config
    
  )
  
}

if (cfg$run_or_update=="update") {
  
  update_on_cluster(
    
    first = {
      sim <- readRDS('/home/akenny/z.monotest/sim.simba')
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = {
      sim %<>% update()
    },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}



##########################################################.
##### VIZ (REGRESSION): Sim #1 (compare all methods) #####
##########################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary(mean=list(name="power",x="reject"))
  summ %<>% mutate(
    test = case_when(
      test=="App_2: var 1" ~ "V1: G0=id, Pstar=unif",
      test=="App_2: var 2" ~ "V2: G0=marg, Pstar=unif",
      test=="App_2: var 3" ~ "V3: G0=id, Pstar=marg",
      test=="App_2: var 4" ~ "V4: G0=marg, Pstar=marg"
    ),
    a_distr = case_when(
      # a_distr=="U" ~ "U-shaped",
      a_distr=="decr" ~ "Decreasing",
      a_distr=="unif" ~ "Uniform",
      # a_distr=="spike" ~ "Spiked"
      a_distr=="incr" ~ "Increasing"
    )
  )
  summ %<>% mutate(
    a_distr = factor(a_distr, levels=c("Uniform","Increasing","Decreasing"))
    # a_distr = factor(a_distr, levels=c("U-shaped","Decreasing","Uniform",
    #                                    "Spiked"))
  )
  
  # Visualize results (alpha_3==0)
  # Export 10" x 5"
  ggplot(summ, aes(x=alpha_3, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_grid(rows=vars(mono_form), cols=vars(a_distr)) +
    scale_color_manual(values=c("#D55E00", "#E69F00", "turquoise", "#009E73")) +
    labs(
      x = "Effect size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone regression ",
                     "(n=30; 1000 sims per level)")
    )
  
  # Plot of densities used
  # Export at 8"x3"
  grid <- seq(0,1,0.01)
  density_labels <- c("U-shaped","Decreasing","Uniform","Spiked")
  ggplot(
    data.frame(
      x = rep(grid,4),
      y = c(dbeta(grid, shape1=0.3, shape2=0.3),
            dbeta(grid, shape1=0.7, shape2=1.3),
            dbeta(grid, shape1=1.0, shape2=1.0),
            dbeta(grid, shape1=40, shape2=40)),
      density = factor(rep(density_labels, each=101), levels=density_labels)
    ),
    aes(x=x, y=y)
  ) +
    geom_line() +
    facet_wrap(~density, ncol=4) +
    labs(x="X", y="Density")
  
}



#######################################################.
##### VIZ (DENSITY): Sim #1 (compare all methods) #####
#######################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")

  summ <- sim %>% summary(mean=list(name="power",x="reject"))
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone density ",
                     "(1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}



###################################################################.
##### VIZ (DENSITY): Sim #2 (different weights for cuml_incr) #####
###################################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary() %>% rename("power"=`mean_reject`)

  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "w_2 (weight)", y = "Power", color = "n", group = "n",
      title = paste0("Power of cuml_incr test for constant vs. monotone ",
                     "density (10,000 sims per level)")
    )
  
}



########################################################################.
##### VIZ (DENSITY): Sim #3 (different delta values for cuml_incr) #####
########################################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim_3 %>% summary() %>% rename("power"=`mean_reject`)
  summ$test <- factor(summ$test, levels=c("1/2","1/3","1/4","1/5","1/6","1/7",
                                          "1/8","1/9","1/10"))
  
  # Visualize results
  ggplot(summ, aes(x=test, y=power, group=(n), color=factor(n))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=3, scales="free_y") +
    labs(
      x = "delta", y = "Power", color = "n", group = "n",
      title = paste0("Power of cuml_incr test for constant vs. monotone ",
                     "density (10,000 sims per level)")
    )
  
}



###############################################################.
##### VIZ (DENSITY): Sim #4 (Variants of approach 2 test) #####
###############################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary(mean=list(name="power",x="reject"))
  summ$true_density %<>% car::Recode("
    'Uniform' = '(a) Uniform';
    'Linear' = '(b) Linear';
    'Exponential' = '(c) Exponential';
    'Spline' = '(d) Spline';
  ")

  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=4, scales="free_y") +
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone density ",
                     "(1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}



#######################################################################.
##### VIZ (DENSITY): Sim #5 (Illustration of "reversing" problem) #####
#######################################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary(mean=list(name="power",x="reject"))
  summ$true_density %<>% car::Recode("
    'Exponential' = '(a) Exponential';
    'Exponential (decr)' = '(b) Exponential (decr)';
    'Step' = '(c) Step';
    'Step (decr)' = '(d) Step (decr)';
  ")
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_wrap(~true_density, ncol=4) +
    # facet_wrap(~true_density, ncol=4, scales="free_y") +
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone density ",
                     "(1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}



########################################################.
##### VIZ (DOSERESP): Sim #1 (compare all methods) #####
########################################################.

if (FALSE) {
  
  # sim <- readRDS("sim.simba")
  
  summ <- sim %>% summary(mean=list(name="power", x="reject"))
  
  # Visualize results
  ggplot(summ, aes(x=n, y=power, color=factor(test))) +
    geom_line() + geom_point() +
    facet_grid(rows=vars(alpha_3), cols=vars(mono_form)) + # scales="free_y"
    labs(
      x = "sample size", y = "Power", color = "Test type",
      title = paste0("Power of tests for constant vs. monotone causal ",
                     "dose-response curve (1,000 sims per level)")
    ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
                                    "green3", "darkorchid2", "orange"))
  
}



###########################################.
##### OTHER (DENSITY): Plot densities #####
###########################################.

if (FALSE) {
  
  # Plot 1: Main densities used
  
  # Generate data
  grid <- seq(0,1,0.01)
  d1 <- sapply(grid, function(x) {1})
  d2 <- sapply(grid, function(x) {2*x})
  d3 <- sapply(grid, function(x) {
    c <- 5
    c*(exp(c)-1-c)^-1 * (exp(c*x)-1)
  })
  d4 <- sapply(grid, function(x) {
    (x<=0.5) * 0.1 +
    (0.5<x)*(x<=0.6) * (20*x-9.9) +
    (x>0.6) * 2.1
  })
  curve_labels <- c("(a) Uniform", "(b) Linear", "(c) Exponential",
                    "(d) Spline")
  
  # Plot functions
  # Export: PDF 8"x3"
  ggplot(
    data.frame(
      x = rep(grid,4),
      y = c(d1,d2,d3,d4),
      fn = factor(rep(curve_labels, each=101), levels=curve_labels)
    ),
    aes(x=x, y=y)
  ) +
    geom_line() +
    facet_wrap(~fn, ncol=4) +
    labs(x="X", y="Density")
  
  # Plot 2: Densities corresponding to sim level set 5
  
  # Generate data
  grid <- seq(0,1,0.01)
  d1 <- sapply(grid, function(x) {
    c <- 5
    c*(exp(c)-1-c)^-1 * (exp(c*x)-1)
  })
  d2 <- sapply(grid, function(x) {
    c <- 5
    c*(exp(c)-1-c)^-1 * (exp(c*(1-x))-1)
  })
  d3 <- sapply(grid, function(x) {
    (x<=0.5) * 0.1 +
    (x>0.5) * 1.9
  })
  d4 <- sapply(grid, function(x) {
    (x<=0.5) * 1.9 +
      (x>0.5) * 0.1
  })
  curve_labels <- c("(a) Exponential", "(b) Exponential (decr)",
                    "(c) Step", "(d) Step (decr)")
  
  # Plot functions
  # Export: PDF 8"x3"
  ggplot(
    data.frame(
      x = rep(grid,4),
      y = c(d1,d2,d3,d4),
      fn = factor(rep(curve_labels, each=101), levels=curve_labels)
    ),
    aes(x=x, y=y)
  ) +
    geom_line() +
    facet_wrap(~fn, ncol=4) +
    labs(x="X", y="Density")
  
  # Plot 2: CDFs corresponding to sim level set 5
  
  # Generate data
  grid <- seq(0,1,0.01)
  d3 <- sapply(grid, function(x) {
    (x<=0.5) * 0.1*x +
    (x>0.5) * (1.9*x-0.9)
  })
  d4 <- sapply(grid, function(x) {
    (x<=0.5) * 1.9*x +
    (x>0.5) * (0.1*x+0.9)
  })
  curve_labels <- c("(c) Step", "(d) Step (decr)")
  n <- 50
  sample_3 <- runif(n=n, min=0.5, max=1) - 0.5*rbinom(n=n, size=1, prob=0.05)
  sample_4 <- runif(n=n, min=0, max=0.5) + 0.5*rbinom(n=n, size=1, prob=0.05)
  ecdf_3 <- ecdf(sample_3)(sample_3)
  ecdf_4 <- ecdf(sample_4)(sample_4)
  
  # Plot functions
  # Export: PDF 8"x3"
  ggplot(
    data.frame(
      x = rep(grid,2),
      y = c(d3,d4),
      fn = factor(rep(curve_labels, each=101), levels=curve_labels)
    ),
    aes(x=x, y=y)
  ) +
    geom_line() +
    facet_wrap(~fn, ncol=4) +
    labs(x="X", y="CDF") +
    geom_point(
      data = data.frame(
        x = c(sample_3, sample_4),
        y = c(ecdf_3, ecdf_4),
        # y = rep(0,2*n),
        fn = factor(rep(curve_labels, each=n), levels=curve_labels)
      ),
      alpha=0.2
    )
  
}



###########################################.
##### OTHER (DENSITY): Plot densities #####
###########################################.

if (FALSE) {
  
  # Graph for proof that beta_0=0 iff H_0 holds
  library(latex2exp)
  
  # Graph of big Theta
  # Export at 400x300
  Theta <- sapply(seq(0,1,0.01), function(x) {x^3})
  L <- sapply(seq(0,1,0.01), function(x) {0.6*x})
  Q <- sapply(seq(0,1,0.01), function(x) {
    b1 <- 0.4
    b0 <- 0.6 - b1*sqrt(0.6)
    return(b0*x + b1*x^2)
  })
  df <- data.frame(
    x = rep(seq(0,1,0.01), 3),
    y = c(Theta, L, Q),
    fn = rep(c("Theta_0", "L", "Q"), each=101)
  )
  ggplot(df, aes(x=x, y=y, color=fn)) + geom_line() +
    geom_point(data=data.frame(x=0.6^0.5,y=0.6^1.5), color="orange") +
    geom_point(data=data.frame(x=0,y=0), color="blue") +
    labs(y="f(x)", color="Function") +
    scale_color_discrete(labels = c("L", "Q", unname(TeX("$\\Theta_0"))))
  
  # # Graph of little theta
  # # Export at 400x300
  # theta <- sapply(seq(0,1,0.01), function(x) {3*x^2})
  # l <- rep(0.6, 101)
  # q <- sapply(seq(0,1,0.01), function(x) {0.3 + (0.3/(0.2^0.5))*x})
  # df2 <- data.frame(
  #   x = rep(seq(0,1,0.01), 3),
  #   y = c(theta, l, q),
  #   fn = rep(c("theta", "L'", "Q'"), each=101)
  # )
  # ggplot(df2, aes(x=x, y=y, color=fn)) + geom_line() +
  #   geom_point(data=data.frame(x=0.2^0.5,y=0.6), color="purple")
  
}



########################################################.
##### TESTING (DOSERESP): Generate testing dataset #####
########################################################.

if (FALSE) {
  
  dat <- generate_data_dr(
    n = 10,
    alpha_3 = 0.5,
    mono_form = "identity",
    sampling = list(type="iid")
  )
  
}
