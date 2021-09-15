# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
# GitHub packages:
#   - tedwestling/ctsCausal
#   - tedwestling/CFsurvival
#   - tedwestling/survSuperLearner
#   - zeehio/facetscales
cfg <- list(
  which_sim = "estimation", # estimation testing
  level_set_which = "level_set_estimation_1", # level_set_estimation_1 level_set_testing_1
  run_or_update = "run",
  num_sim = 300,
  pkgs = c("dplyr", "boot", "car", "mgcv", "memoise", "EnvStats", "fdrtool",
           "splines", "survival", "SuperLearner", "survSuperLearner",
           "randomForestSRC", "CFsurvival", "Rsolnp"),
  pkgs_nocluster = c("ggplot2", "viridis", "sqldf", "facetscales", "scales",
                     "data.table", "latex2exp", "tidyr"),
  parallel = "none",
  stop_at_error = FALSE,
  appx = list(t_e=1, w1=0.01, w1b=0.1, a=0.01)
  # appx = list(t_e=1, w1=0.001, w1b=0.01, a=0.001)
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

# Load SimEngine + functions
{
  library(SimEngine)
  source("one_simulation.R", local=TRUE)
  source("generate_data.R", local=TRUE)
  source("est_curve.R", local=TRUE)
  source("test_2.R", local=TRUE)
  source("fns_doseresp.R", local=TRUE)
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("sim_run") %in% c("first", "")) {
  
  # Level bank
  if (F) {
    level_bank <- list(
      alpha_3 = c(0, 0.4, 0.8),
      distr_A = c("Unif(0,1)", "Beta(1.5+w1,1.5+w2)"),
      edge = c("none", "expit", "complex"),
      surv_true = c("Cox PH", "complex"),
      sampling = c("iid", "two-phase"),
      deriv_type = c("linear", "gcomp", "spline"),
      estimator = list(
        "G-comp" = list(
          est = "G-comp", cf_folds = 1,
          params = list(S_n_type="Cox PH", boot_reps=100)
        ),
        "Grenander" = list(
          est = "Grenander",
          params = list(
            ci_type = c("regular", "logit", "sample split", "none"),
            cf_folds = c(1,10),
            S_n_type = c("Cox PH", "Random Forest"),
            g_n_type = c("parametric", "binning"),
            edge_corr = c("none", "point", "spread", "max"),
            deriv_type = c("line", "spline", "linear", "gcomp")
          )
        )
      ),
      test = list(
        "Slope (one-step Gamma_n)" = list(
          type = "test_2",
          params = list(
            var = "asymptotic",
            cf_folds = c(1,10),
            S_n_type = c("Cox PH", "Random Forest"),
            g_n_type = c("parametric", "binning"),
            est_known_nuis = c(TRUE,FALSE)
          )
        )
      )
    )
  }
  
  # Estimation: compare all methods
  # Not currently using (ci_type="sample split", m=5) or (ci_type="regular")
  level_set_estimation_1 <- list(
    n = 5000,
    alpha_3 = 0.8,
    distr_A = "Beta(1.5+w1,1.5+w2)",
    # distr_A = c("Unif(0,1)", "Beta(1.5+w1,1.5+w2)"),
    edge = "none",
    # edge = c("none", "expit"),
    surv_true = "complex",
    # surv_true = c("Cox PH", "complex"),
    sampling = "two-phase",
    estimator = list(
      "Grenander (spline)" = list(
        est = "Grenander",
        params = list(S_n_type="Random Forest", g_n_type="binning",
                      ci_type="regular", cf_folds=1, edge_corr="max",
                      deriv_type="spline")
      ),
      "Grenander (line)" = list(
        est = "Grenander",
        params = list(S_n_type="Random Forest", g_n_type="binning",
                      ci_type="regular", cf_folds=1, edge_corr="max",
                      deriv_type="line")
      )
    )
  )
  
  # Testing: compare all methods
  # Not currently using (var="boot")
  level_set_testing_1 <- list(
    n = c(1000,2000),
    alpha_3 = c(0,0.3,0.6),
    distr_A = c("Unif(0,1)","Beta(1.5+w1,1.5+w2)"),
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase",
    test = list(
      "Slope (one-step Gamma_n)" = list(
        type = "test_2",
        params = list(var="asymptotic", S_n_type="Cox PH", cf_folds=1,
                      g_n_type="parametric", est_known_nuis=FALSE)
      )
    )
  )
  
  # !!!!! Estimation/testing but with different nuisance estimators
  # S_n_type = c("Cox PH", "...")
  # g_n_type = c("parametric", "binning")
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Use these commands to run on Slurm:
# sbatch --export=sim_run='first',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-1200 --export=sim_run='main',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=sim_run='last',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

if (cfg$run_or_update=="run") {

  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        stop_at_error = cfg$stop_at_error,
        seed = 5,
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      
      # Add functions to simulation object
      sim %<>% add_creator(generate_data)
      methods <- c(
        "logit", "expit", "deriv_logit", "deriv_expit", "create_htab", "Pi",
        "wts", "construct_S_n", "construct_gcomp_n", "construct_deriv_theta_n",
        "construct_tau_n", "construct_gamma_n", "construct_f_aIw_n",
        "construct_f_a_n", "construct_g_n", "construct_omega_n",
        "construct_eta_n", "construct_Gamma_n", "construct_Phi_n",
        "construct_rho_n", "construct_xi_n", "construct_infl_fn_1",
        "construct_infl_fn_Gamma", "construct_infl_fn_2", "beta_n_var_hat",
        "create_val_list", "construct_Gamma_cf_k", "construct_Gamma_cf",
        "construct_pi_n", "theta_os_n", "sigma2_os_n",
        
        "est_curve", "generate_data",
        "lambda", "one_simulation", "test_2"
      )
      for (method in methods) {
        sim %<>% add_method(method, eval(as.name(method)))
      }
      
      # Add constants
      # lambda and v are the Weibull parameters for the survival distribution
      # lambda2 and v2 are the Weibull parameters for the censoring distribution
      # For the approximations, w1b is used for S_n and w1 is used elsewhere
      sim %<>% add_constants(
        lambda = 10^-4,
        v = 1.5,
        lambda2 = 0.5 * 10^-4,
        v2 = 1.5,
        points = seq(0,1,0.02),
        alpha_1 = 0.3,
        alpha_2 = 0.7,
        t_e = 200,
        appx = cfg$appx
      )
      
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
      sim <- readRDS('/home/akenny/z.VaxCurve/sim.rds')
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = {
      sim %<>% update_sim()
    },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
}



###################################.
##### VIZ: Estimation (sim 1) #####
###################################.

if (FALSE) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_est_20210819.rds")
  
  # Summarize results
  # !!!!! Count the number of NAs in coverage
  summ_bias <- list()
  summ_mse <- list()
  summ_cov <- list()
  for (i in c(1:51)) { # for (i in c(1:11)) {
    m <- format(round(i/50-0.02,2), nsmall=1) # m <- format(round(i/10-0.1,1), nsmall=1)
    summ_bias[[i]] <- list(
      name = paste0("bias_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
    summ_mse[[i]] <- list(
      name = paste0("mse_",m),
      estimate = paste0("est_",m),
      truth = paste0("theta_",m)
    )
    summ_cov[[i]] <- list(
      name = paste0("cov_",m),
      truth = paste0("theta_",m),
      lower = paste0("ci_lo_",m),
      upper = paste0("ci_hi_",m),
      na.rm = TRUE
    )
  }
  summ <- summarize(sim, bias_pct=summ_bias, mse=summ_mse, coverage=summ_cov)
  
  summ %<>% rename("Estimator"=estimator)
  
  # summ %<>% filter(distr_A=="Unif(0,1)") # !!!!!
  # summ %<>% filter(distr_A=="Beta(1.5+w1,1.5+w2)") # !!!!!
  
  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,alpha_3,distr_A,edge,surv_true,sampling,Estimator),
    # cols = -c(level_id,n,alpha_3,distr_A,edge,surv_true,sampling,Estimator),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point = as.numeric(point))
  
  cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999")
  m_colors <- c(
    # `G-comp` = cb_colors[1],
    `Grenander (logit CIs)` = cb_colors[2]
    # `Grenander (split CIs, m=5)` = cb_colors[3]
    # `Grenander (regular CIs)` = cb_colors[4],
  )
  
  # Bias plot
  # Export: 8" x 5"
  # p_data %<>% filter(sampling=="two-phase") # !!!!!
  ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    # geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(surv_true)) +
    # scale_y_continuous(labels=percent, limits=c(-0.4,0.4)) +
    scale_y_continuous(labels=percent, limits=c(-0.1,0.1)) +
    # scale_y_continuous(labels=percent) +
    # scale_color_manual(values=m_colors) +
    labs(title="Bias (%)", x="A", y=NULL, color="Estimator")
  
  # Coverage plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="cov"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_hline(aes(yintercept=0.95), linetype="longdash", color="grey") +
    # geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(surv_true)) +
    scale_y_continuous(labels=percent, limits=c(0.7,1)) +
    # scale_y_continuous(labels=percent) +
    # scale_color_manual(values=m_colors) +
    labs(title="Coverage (%)", x="A", y=NULL, color="Estimator")
  
  # MSE plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="mse"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    # geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(surv_true)) +
    # scale_color_manual(values=m_colors) +
    ylim(0,0.003) +
    labs(title="MSE", x="A", y=NULL, color="Estimator")
  
}



################################.
##### VIZ: Testing (sim 1) #####
################################.

if (FALSE) {
  
  # # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_testing_20210820.rds")
  
  # !!!!! Modify everything below
  #   color should be n-value
  #   x-axis should be alpha_3
  
  # Summarize resuls
  summ <- summarize(sim)
  
  summ %<>% rename(
    "Power" = mean_reject
  )
  
  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    `500` = cb_colors[2],
    `1000` = cb_colors[3]
    # `Slope (boot)` = cb_colors[2],
    # `Slope (mixed boot)` = cb_colors[3],
    # `Wald` = cb_colors[4]
  )
  
  # Export: 7" x 4.5"
  # distr_A_ <- "Unif(0,1)" # "Beta(1.5+w1,1.5+w2)"
  ggplot(
    summ,
    aes(x=alpha_3, y=Power, color=factor(n))
  ) +
    geom_point() +
    geom_line() +
    facet_grid(cols=dplyr::vars(test), rows=dplyr::vars(distr_A)) +
    # facet_grid(cols=dplyr::vars(surv_true), rows=dplyr::vars(distr_A)) +
    scale_y_continuous(labels=percent) +
    theme(legend.position="bottom") +
    scale_color_manual(values=m_colors) +
    labs(title="Testing", color="Sample size")
    # labs(title = paste0("alpha_3: ",alpha_3_,"; Sampling: ",sampling_),
    #      color = "Test")
  
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
