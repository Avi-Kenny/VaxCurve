# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  which_sim = "estimation", # estimation testing
  level_set_which = "level_set_estimation_1",
  run_or_update = "run",
  num_sim = 1,
  pkgs = c("dplyr", "boot", "car", "mgcv", "kdensity", "memoise",
           "twostageTE", "EnvStats", "fdrtool"),
  pkgs_nocluster = c("ggplot2", "viridis", "sqldf", "facetscales", "scales",
                     "data.table", "latex2exp"), # devtools::install_github("zeehio/facetscales")
  parallel = "none", # none outer
  stop_at_error = FALSE
)

# Set cluster config
cluster_config <- list(
  js = "slurm",
  dir = "/home/akenny/z.VaxCurve"
)



#################.
##### SETUP #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Marco Carone/Project - VaxCurve/VaxCurve/R"))
  load_pkgs_local <- TRUE
} else {
  # Cluster
  setwd("z.VaxCurve/R")
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
source("one_simulation.R")
source("generate_data.R")
source("est_curve.R")
source("test_2.R")
source("fns_doseresp.R")



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("run") %in% c("first", "")) {
  
  # Estimation: compare all methods
  level_set_estimation_1 <- list(
    n = 200, # c(50,100,200)
    alpha_3 = 0.7,
    distr_A = c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)"),
    mono_form = c("identity", "step_0.2"), # "square"
    sampling = c("iid", "two-phase"),
    estimator = list(
      "G-comp (logistic)" = list(
        est = "G-comp (logistic)",
        params = list(boot_reps=100)),
      "Grenander (logit CIs)" = list(
        est = "Generalized Grenander",
        params = list(ci_type="logit")),
      "Grenander (split CIs; m=5)" = list(
        est = "Generalized Grenander",
        params = list(ci_type="sample split", m=5)),
      "Grenander (regular CIs)" = list(
        est = "Generalized Grenander",
        params = list(ci_type="regular"))
    )
  )
  
  # Testing: compare all methods
  level_set_testing_1 <- list(
    n = 200, # c(50,100,200)
    alpha_3 = 0.7,
    distr_A = c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)"),
    mono_form = c("identity", "step_0.2"), # "square"
    sampling = c("iid", "two-phase"),
    test = list(
      "One" = list(
        a = "a",
        b = b),
      "Two" = list(
        a = "a",
        b = b)
    )
  )
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Use these commands to run on Slurm:
# sbatch --export=run='first',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-16000 --export=run='main',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=run='last',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

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
      methods <- c(
        "est_curve", "test_2", "expit", "deriv_expit", "logit", "deriv_logit",
        "construct_mu_n", "construct_deriv_theta_n", "construct_sigma2_n",
        "construct_f_a_n", "construct_f_aIw_n", "construct_g_n",
        "construct_Gamma_n", "construct_Phi_n", "Pi", "wts"
      )
      for (method in methods) {
        sim %<>% add_method(method, eval(as.name(method)))
      }
      
      # Add constants
      data(chernoff_realizations)
      sim %<>% add_constants(chern=chernoff_realizations)
      
      # Simulation script
      sim %<>% set_script(one_simulation)
      
    },
    
    main = {
      sim %<>% run()
    },
    
    last = {
      
      sim %>% summarize() %>% print()

    },
    
    cluster_config = cluster_config
    
  )
  
}

if (cfg$run_or_update=="update") {
  
  update_sim_on_cluster(
    
    first = {
      sim <- readRDS('/home/akenny/z.VaxCurve/sim.simba')
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = {
      sim %<>% update_sim()
    },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}



#######################.
##### VIZ: Sim #1 #####
#######################.

if (FALSE) {
  
  # # Read in simulation object
  # sim <- readRDS("../simba.out/sim123_20210506.simba")
  
  # Summarize resuls
  summ <- summarize(
    sim_obj = sim,
    mean = list(
      list(x="est_mp"),
      list(x="est_ep")
    ),
    bias = list(
      list(name="bias_mp", estimate="est_mp", truth="theta_mp"),
      list(name="bias_ep", estimate="est_ep", truth="theta_ep")
    ),
    mse = list(
      list(name="mse_mp", estimate="est_mp", truth="theta_mp"),
      list(name="mse_ep", estimate="est_ep", truth="theta_ep")
    ),
    coverage = list(
      list(name="cov_mp", truth="theta_mp", estimate="est_mp",
           lower="ci_lo_mp", upper="ci_hi_mp"),
      list(name="cov_ep", truth="theta_ep", estimate="est_ep",
           lower="ci_lo_ep", upper="ci_hi_ep")
    )
  )
  
  summ %<>% rename("Estimator"=estimator)
  
  p_data <- sqldf("
    SELECT Estimator, mono_form, distr_A, 'Midpoint' AS estimand,
      bias_mp AS Bias, n, cov_mp AS Coverage, mse_mp AS MSE FROM summ
    UNION SELECT Estimator, mono_form, distr_A, 'Endpoint',
      bias_ep, n, cov_ep, mse_ep FROM summ
  ")
  p_data <- sqldf("
    SELECT Estimator, mono_form, distr_A, estimand, 'Bias' AS stat,
    Bias AS value, n FROM p_data
    UNION SELECT Estimator, mono_form, distr_A, estimand, 'Coverage',
    Coverage, n FROM p_data
    UNION SELECT Estimator, mono_form, distr_A, estimand, 'MSE',
    MSE, n FROM p_data
  ")
  
  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    `G-comp (logistic)` = cb_colors[2],
    `Generalized Grenander` = cb_colors[3]
  )

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
    facet_grid_sc(cols=vars(mono_form), rows=vars(stat), scales=list(y=list(
      Bias = scale_y_continuous(labels = percent_format()),
      Coverage = scale_y_continuous(labels = percent_format()),
      MSE = scale_y_continuous()
    ))) +
    theme(legend.position="bottom") +
    scale_color_manual(values=m_colors) +
    labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
         y=NULL, x=NULL, color="Estimator")
  
}



#####################################.
##### VIZ: Method demonstration #####
#####################################.

if (FALSE) {
  
  # Set parameters
  n <- 200
  alpha_0 <- -1.5
  alpha_1 <- 0.3
  alpha_2 <- 0.7
  alpha_3 <- 0.7
  mono_form <- "square"
  # mono_form <- "step_0.2"
  # mono_form <- "identity"
  mono_f <- function(x) {x^2}
  # mono_f <- function(x) {as.integer(x>0.2)}
  # mono_f <- function(x) {x}
  grid <- seq(0,1,0.01)
  
  # Generate dataset
  dat <- generate_data(
    n = n,
    alpha_3 = alpha_3,
    distr_A = "Beta(0.9,1.1+0.4*w2)",
    mono_form = mono_form,
    sampling = "iid"
  )
  
  # Approximate true values of theta_0
  m <- 10000
  w1 <- rnorm(m)
  w2 <- rbinom(m, size=1, prob=0.5)
  true_vals <- sapply(grid, function(x) {
    mean(expit(
      alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*mono_f(x)
    ))
  })
  
  # Estimate curve: G-comp (logistic)
  ests <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    # estimator = "G-comp (logistic)",
    estimator = "Generalized Grenander",
    params = list(ci_type="logit"),
    # params = NULL,
    points = grid
  )))
  ests$which <- "Grenander (logit CIs)"
  # ests$which <- "G-comp (logistic)"
  
  # Estimate curve: Generalized Grenander
  ests2 <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    estimator = "Generalized Grenander",
    params = list(ci_type="regular"),
    points = grid
  )))
  ests2$which <- "Grenander (regular CIs)"
  ests <- rbind(ests, ests2)
  
  # Attach true values
  ests <- rbind(ests, data.frame(
    point = grid,
    est = true_vals,
    ci_lo = true_vals,
    ci_hi = true_vals,
    # se = rep(0, length(grid)),
    which = rep("theta_0", length(grid))
  ))
  
  # Plot results
  # Export: 8" x 3"
  ggplot(
    ests,
    aes(x=point, y=est, color=which, fill=which)) +
    geom_line() +
    # geom_point() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_hi), # fill=as.factor(group_var)
      alpha = 0.2,
      linetype = "dotted"
    ) +
    facet_wrap(~which, ncol=3) +
    labs(
      x = "x",
      y = unname(latex2exp::TeX("$\\hat{\\theta}_n(x)$")),
      color = "Which",
      fill = "Which"
    )
  
  # end_true <- (filter(ests, point==1 & which=="theta_0"))$est
  # end_gcomp <- (filter(ests, point==1 & which=="G-comp (logistic)"))$est
  # print((end_true-end_gcomp)^2)
  
}



#######################################.
##### VIZ: Chernoff vs. N(0,0.52) #####
#######################################.

if (FALSE) {
  
  data(chernoff_realizations)
  dat <- chernoff_realizations
  dat$which <- "Chernoff"
  dat2 <- data.frame(
    xcoor = dat$xcoor,
    DF = pnorm(dat$xcoor, sd=0.52),
    density = dnorm(dat$xcoor, sd=0.52),
    which = "N(0,0.52^2)"
  )
  
  # Export: 6" x 4"
  ggplot(rbind(dat,dat2), aes(x=xcoor, y=density, color=factor(which))) +
    geom_line() +
    labs(x="x", y="Density", color="Distribution")
  ggplot(rbind(dat,dat2), aes(x=xcoor, y=DF, color=factor(which))) +
    geom_line() +
    labs(x="x", y="CDF", color="Distribution")
  
}
