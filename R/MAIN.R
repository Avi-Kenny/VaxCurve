# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
cfg <- list(
  which_sim = "estimation", # estimation testing
  level_set_which = "level_set_1",
  run_or_update = "run",
  num_sim = 100,
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
               "ch/Marco Carone/Project - VaxCurve/VaxCurve/R"))
  load_pkgs_local <- TRUE
} else {
  # Cluster
  setwd("VaxCurve/R")
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
source("test2.R")
source("fns_doseresp.R")



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("run") %in% c("first", "")) {
  
  # Compare all methods
  level_set_1 <- list(
    n = c(30,60),
    distr_A = c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)"),
    alpha_3 = 0.7,
    mono_form = c("identity", "square", "step_0.2"),
    estimator = c("G-comp (logistic)"),
    sampling = list(
      "iid" = list(type="iid")
    )
  )
  
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
      sim %<>% add_method(est_curve)
      sim %<>% add_method(test2)
      sim %<>% add_method(expit)
      sim %<>% add_method(deriv_expit)
      
      # Add constants
      data(chernoff_realizations)
      sim %<>% add_constants(
        chern = chernoff_realizations
      )
      
      # Simulation script
      sim %<>% set_script(one_simulation)
      
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
      sim <- readRDS('/home/akenny/VaxCurve/sim.simba')
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = {
      sim %<>% update()
    },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}



###############################################################.
##### VIZ (DOSERESP): Sim #1 (compare all methods; power) #####
###############################################################.

if (FALSE) {
  
  # # Read in simulation object
  # sim <- readRDS("../simba.out/sim123_20210506.simba")
  
  # Summarize resuls
  summ <- summary(
    sim_obj = sim,
    mean = list(
      list(name="est_mp", x="est_mp"),
      list(name="est_ep", x="est_ep")
    ),
    bias_pct = list(
      list(name="bias_mp", estimate="est_mp", truth="theta_mp"),
      list(name="bias_ep", estimate="est_ep", truth="theta_ep")
    ),
    mse = list(
      list(name="mse_mp", estimate="est_mp", truth="theta_mp"),
      list(name="mse_ep", estimate="est_ep", truth="theta_ep")
    ),
    coverage = list(
      list(name="cov_mp", truth="theta_mp", estimate="est_mp", se="se_mp"),
      list(name="cov_ep", truth="theta_ep", estimate="est_ep", se="se_ep")
    )
  )
  
  # !!!!! Adapt
  
  # summ %<>% mutate(
  #   method = factor(method, levels=s_methods),
  #   delay_model = factor(delay_model, levels=s_d_models),
  #   power_ate = 1 - beta_ate,
  #   power_lte = 1 - beta_lte
  # )
  # summ %<>% rename("Method"=method, "n_extra"=n_extra_time_points)
  
  p_data <- sqldf("
    SELECT Method, n_extra, delay_model, 'TATE' AS which, bias_ate AS bias,
      theta, cov_ate AS Coverage, power_ate AS Power, mse_ate AS MSE FROM summ
    UNION SELECT Method, n_extra, delay_model, 'LTE', bias_lte,
      theta, cov_lte, power_lte, mse_lte FROM summ
  ")
  p_data <- sqldf("
    SELECT Method, n_extra, delay_model, which, 'Bias' AS stat, Bias AS value,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'Coverage', Coverage,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'Power', Power,
      theta FROM p_data
    UNION SELECT Method, n_extra, delay_model, which, 'MSE', MSE,
      theta FROM p_data
  ")
  p_data %<>% mutate(which = factor(which, levels=c("TATE", "LTE")))
  p_data %<>% mutate(n_extra = as.character(n_extra))
  p_data %<>% rename("Extra time points"=n_extra)
  
  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    IT = cb_colors[2],
    ETI = cb_colors[4],
    `NCS (4df)` = cb_colors[3],
    `MEC (P=0.1)` = cb_colors[6],
    # `CUBIC-4df` = cb_colors[6],
    # `NCS-2df` = cb_colors[7],
    # `NCS-3df` = cb_colors[8],
    # `NCS-4df` = cb_colors[5],
    `RETI (3 steps)` = cb_colors[6],
    `RETI (4 steps)` = cb_colors[7],
    `ETI (RTE; height)` = cb_colors[6],
    `ETI (RTE; height+shape)` = cb_colors[7],
    `0` = cb_colors[4],
    `1` = cb_colors[7],
    `2` = cb_colors[6]
  )
  # viridis(5)
  
  # Export: 8: x 4"
  ggplot(
    p_data,
    aes(x=which, y=value, fill=Method)
  ) +
    geom_hline(
      aes(yintercept=y),
      data=data.frame(y=0.95, stat="Coverage"),
      linetype="longdash", color="grey"
    ) +
    geom_bar(stat="identity", position=position_dodge(),
             width=0.8, color="white", size=0.35) +
    facet_grid_sc(cols=vars(delay_model), rows=vars(stat), scales=list(y=list(
      Bias = scale_y_continuous(labels = percent_format()),
      Coverage = scale_y_continuous(labels = percent_format()),
      MSE = scale_y_continuous()
    ))) +
    theme(legend.position="bottom") +
    scale_fill_manual(values=m_colors) +
    labs(y=NULL, x=NULL, fill="Analysis model")
  
}



###############################################################.
##### VIZ (DOSERESP): Sim #1 (compare all methods; power) #####
###############################################################.

if (FALSE) {
  
  # # sim <- readRDS("sim.simba")
  # 
  # summ <- sim %>% summary(mean=list(name="power", x="reject"))
  # 
  # # Visualize results
  # ggplot(summ, aes(x=n, y=power, color=factor(test))) +
  #   geom_line() + geom_point() +
  #   facet_grid(rows=vars(alpha_3), cols=vars(mono_form)) + # scales="free_y"
  #   labs(
  #     x = "sample size", y = "Power", color = "Test type",
  #     title = paste0("Power of tests for constant vs. monotone causal ",
  #                    "dose-response curve (1,000 sims per level)")
  #   ) + scale_color_manual(values=c("turquoise", "salmon", "dodgerblue2",
  #                                   "green3", "darkorchid2", "orange"))
  
}
