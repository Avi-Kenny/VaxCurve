# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
# GitHub packages: tedwestling/ctsCausal, tedwestling/CFsurvival, 
#                  tedwestling/survSuperLearner, zeehio/facetscales
cfg <- list(
  main_task = "run", # run update analysis_moderna.R analysis_705.R
  which_sim = "testing", # estimation edge testing
  level_set_which = "level_set_testing_1", # level_set_estimation_1 level_set_testing_1
  # keep = c(1:3,7:9,16:18,22:24),
  num_sim = 1000,
  pkgs = c("dplyr", "boot", "car", "mgcv", "memoise", "EnvStats", "fdrtool",
           "splines", "survival", "SuperLearner", "survSuperLearner",
           "randomForestSRC", "CFsurvival", "Rsolnp", "truncnorm"),
  pkgs_nocluster = c("ggplot2", "viridis", "sqldf", "facetscales", "scales",
                     "data.table", "latex2exp", "tidyr"),
  parallel = "none",
  stop_at_error = FALSE,
  appx = list(t_e=1, w1=0.1, w1b=0.1, a=0.01) # !!!!! a=0.001
)

# Set cluster config
cluster_config <- list(
  # js = "ge",
  # dir = "/home/users/avikenny/Desktop/z.VaxCurve"
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
  if (cfg$main_task %in% c("run", "update")) {
    load_pkgs_local <- FALSE
  } else {
    load_pkgs_local <- TRUE
  }
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
      n = 5000,
      alpha_3 = c(0,-2,-4),
      sc_params = list("sc_params"=list(lmbd=9e-7, v=1.5, lmbd2=1e-5, v2=1.5)),
      distr_A = c("Unif(0,1)", "N(0.5,0.01)", "N(0.5,0.04)"),
      edge = c("none", "expit", "complex"),
      surv_true = c("Cox PH", "complex"),
      sampling = c("iid", "two-phase (6%)", "two-phase (72%)"),
      estimator = list(
        "G-comp" = list(
          est = "gcomp", cf_folds = 1,
          params = list(S_n_type="Cox PH", boot_reps=100)
        ),
        "Grenander" = list(
          est = "Grenander",
          params = list(
            ci_type = c("regular", "logit", "sample split", "none"),
            cf_folds = c(1,10),
            S_n_type = c("Cox PH", "Random Forest"),
            g_n_type = c("parametric", "binning"),
            edge_corr = c("none", "point", "spread", "min"),
            deriv_type = c("line", "spline", "m-spline", "linear", "gcomp"),
            ecdf_type = c("step", "linear", "true")
          )
        )
      ),
      test = list(
        "Slope (one-step Gamma_os_n)" = list(
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
  
  # Estimation: ideal params
  level_set_temp <- list(
    n = 2000,
    alpha_3 = -2,
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_A = "Unif(0,1)",
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase (72%)",
    estimator = list(
      "Qbins (3)" = list(
        est = "Qbins",
        params = list(S_n_type="Cox PH", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel", n_bins=3)
      ),
      "Qbins (6)" = list(
        est = "Qbins",
        params = list(S_n_type="Cox PH", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel", n_bins=6)
      ),
      "Qbins (9)" = list(
        est = "Qbins",
        params = list(S_n_type="Cox PH", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel", n_bins=9)
      ),
      "Grenander" = list(
        est = "Grenander",
        params = list(S_n_type="Cox PH", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel", n_bins=0)
      )
    )
  )
  
  # Estimation: ideal params
  level_set_estimation_1 <- list(
    n = 5000,
    alpha_3 = -2,
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_A = c("Unif(0,1)", "N(0.5,0.01)", "N(0.5,0.04)"),
    edge = "none",
    # edge = c("none", "expit"),
    surv_true = c("Cox PH", "complex"),
    sampling = "two-phase (72%)",
    estimator = list(
      # "Grenander (no edge corr)" = list(
      #   est = "Grenander",
      #   params = list(S_n_type="Random Forest", g_n_type="binning",
      #                 ci_type="regular", cf_folds=1, edge_corr="none",
      #                 ecdf_type="linear (mid)", deriv_type="m-spline",
      #                 gamma_type="kernel")
      # ),
      # "Grenander (edge corr)" = list(
      #   est = "Grenander",
      #   params = list(S_n_type="Random Forest", g_n_type="binning",
      #                 ci_type="regular", cf_folds=1, edge_corr="min",
      #                 ecdf_type="linear (mid)", deriv_type="m-spline",
      #                 gamma_type="kernel")
      # )
      "Grenander (estimated nuisances)" = list(
        est = "Grenander",
        params = list(S_n_type="Random Forest", g_n_type="binning",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel")
      ),
      "Grenander (true nuisances)" = list(
        est = "Grenander",
        params = list(S_n_type="true", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="true", deriv_type="m-spline",
                      gamma_type="kernel")
      )
    )
  )
  
  # Estimation: trial params
  level_set_estimation_2 <- list(
    n = 14000,
    alpha_3 = -4,
    sc_params = list("sc_params"=list(lmbd=3e-5, v=1.5, lmbd2=3e-5, v2=1.5)),
    distr_A = "N(0.5,0.01)",
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase (6%)",
    estimator = list(
      "Grenander (true nuisances)" = list(
        est = "Grenander",
        params = list(S_n_type="true", g_n_type="true",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="true", deriv_type="m-spline",
                      gamma_type="kernel")
      ),
      "Grenander (est nuisances)" = list(
        est = "Grenander",
        params = list(S_n_type="Cox PH", g_n_type="binning",
                      ci_type="regular", cf_folds=1, edge_corr="none",
                      ecdf_type="linear (mid)", deriv_type="m-spline",
                      gamma_type="kernel")
      )
    )
  )
  
  # Testing: compare all methods
  level_set_testing_1 <- list(
    n = 2000,
    # n = c(1000,2000),
    alpha_3 = 0,
    # alpha_3 = c(0,-0.25,-0.5),
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_A = "Unif(0,1)",
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase (72%)",
    test = list(
      "Slope (true nuisances)" = list(
        type = "test_2",
        params = list(var="asymptotic", ecdf_type="true",
                      g_n_type="true", S_n_type="true", cf_folds=1,
                      est_known_nuis=FALSE)
      )
      # "Slope (est nuisances)" = list(
      #   type = "test_2",
      #   params = list(var="asymptotic", ecdf_type="linear (mid)",
      #                 g_n_type="binning", S_n_type="Cox PH", cf_folds=1,
      #                 est_known_nuis=FALSE)
      # )
    )
  )
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



######################################################.
##### MAIN: Setup and run simulation (or script) #####
######################################################.

# Use these commands to run on Slurm:
# sbatch --export=sim_run='first',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-600 --export=sim_run='main',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=sim_run='last',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./z.VaxCurve/sim_output.out --constraint=gizmok run_r.sh
# sbatch --export=cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./z.VaxCurve/slurm-%A_%a.out --constraint=gizmok run_r.sh

# Commands for job sumbission on SGE:
# qsub -v sim_run='first',cluster='bayes',type='R',project='z.VaxCurve' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992344 -t 1-600 -v sim_run='main',cluster='bayes',type='R',project='z.VaxCurve' -cwd -e ./io/ -o ./io/ run_r.sh
# qsub -hold_jid 1992345 -v sim_run='last',cluster='bayes',type='R',project='z.VaxCurve' -cwd -e ./io/ -o ./z.VaxCurve/ run_r.sh

if (cfg$main_task=="run") {
  
  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        stop_at_error = cfg$stop_at_error,
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      if (!is.null(cfg$keep)) { sim %<>% set_levels(.keep=cfg$keep) }

      # Add functions to simulation object
      sim %<>% add_creator(generate_data)
      methods <- c(
        "logit", "expit", "deriv_logit", "deriv_expit", "construct_superfunc",
        "Pi", "wts", "construct_S_n", "construct_gcomp_n",
        "construct_deriv_theta_n", "construct_tau_n", "construct_gamma_n",
        "construct_f_aIw_n", "construct_f_a_n", "construct_g_n",
        "construct_omega_n", "construct_eta_n", "construct_Gamma_os_n",
        "construct_Phi_n", "construct_rho_n", "construct_xi_n",
        "construct_infl_fn_1", "construct_infl_fn_Gamma", "construct_infl_fn_2",
        "beta_n_var_hat", "create_val_list", "construct_Gamma_cf_k",
        "construct_Gamma_cf", "construct_pi_n", "theta_os_n", "sigma2_os_n",
        "ss",
        
        "est_curve", "generate_data",
        "lambda", "one_simulation", "test_2"
      )
      for (method in methods) {
        sim %<>% add_method(method, eval(as.name(method)))
      }
      
      # Add constants
      # For the approximations, w1b is used for S_n and w1 is used elsewhere
      sim %<>% add_constants(
        points = seq(0,1,0.02),
        alpha_1 = 0.5,
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
  
} else if (cfg$main_task=="update") {
  
  update_sim_on_cluster(
    
    first = {
      sim <- readRDS(paste0(cluster_config$dir,"/sim.rds"))
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = {
      sim %<>% update_sim()
    },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
} else {
  
  source(cfg$main_task, local=TRUE)
  
}



###################################.
##### VIZ: Estimation (sim 1) #####
###################################.

if (FALSE) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_est_20211004.rds")
  
  # Summarize results
  summ_bias <- list()
  summ_mse <- list()
  summ_cov <- list()
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=1)
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
  
  # summ %<>% filter(distr_A=="N(0.5,0.01)")
  # summ %<>% filter(distr_A=="N(0.5,0.04)")
  # summ %<>% filter(distr_A=="N(0.4+0.2w1+0.1w2,0.01)")
  
  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,alpha_3,sc_params,distr_A,edge,
              surv_true,sampling,Estimator),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point = as.numeric(point))
  
  # cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
  #                "#0072B2", "#D55E00", "#CC79A7", "#999999")
  # m_colors <- c(cb_colors[2], cb_colors[5], cb_colors[6]
    # `Grenander (Est S_n/g_n)` = cb_colors[2],
  # )
  
  df_distr_A <- data.frame(
    x = rep(seq(0,1,0.01),3),
    ymax = c(dnorm(seq(0,1,0.01), mean=0.5, sd=0.1),
             dnorm(seq(0,1,0.01), mean=0.5, sd=0.2),
             rep(1,101)),
    distr_A = rep(c("N(0.5,0.01)", "N(0.5,0.04)", "Unif(0,1)"), each=101),
    value = 0
  )
  df_distr_A1 <- mutate(df_distr_A, ymin=-0.2, ymax=(ymax/15-0.2))
  df_distr_A2 <- mutate(df_distr_A, ymin=0.7, ymax=(ymax/20+0.7))
  df_vlines <- data.frame(
    x = c(qnorm(0.1,0.5,0.1),qnorm(0.1,0.5,0.2),0.1,
          qnorm(0.9,0.5,0.1),qnorm(0.9,0.5,0.2),0.9),
    distr_A = rep(c("N(0.5,0.01)", "N(0.5,0.04)", "Unif(0,1)"),2)
  )
  
  # distr_A_ <- "N(0.5,0.04)" # "Unif(0,1)"
  # df_distr_A1 %<>% filter(distr_A==distr_A_)
  # df_distr_A2 %<>% filter(distr_A==distr_A_)
  # df_vlines %<>% filter(distr_A==distr_A_)
  
  # Bias plots
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value)
    # aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
                data=df_distr_A1, fill="grey", color=NA, alpha=0.4) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_line() +
    facet_grid(rows=dplyr::vars(surv_true), cols=dplyr::vars(distr_A)) +
    scale_y_continuous(labels=percent, limits=c(-0.2,0.2)) +
    # scale_y_continuous(labels=percent) +
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(title="Bias (%)", x="A", y=NULL, color="Estimator")
  
  # Coverage plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="cov"),
    aes(x=point, y=value)
    # aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
                data=df_distr_A2, fill="grey", color=NA, alpha=0.4) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_hline(aes(yintercept=0.95), linetype="longdash", color="grey") +
    geom_line() +
    facet_grid(rows=dplyr::vars(surv_true), cols=dplyr::vars(distr_A)) +
    scale_y_continuous(labels=percent, limits=c(0.7,1)) +
    # scale_y_continuous(labels=percent) +
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(title="Coverage (%)", x="A", y=NULL, color="Estimator")
  
  # MSE plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="mse"),
    # aes(x=point, y=value)
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_hline(aes(yintercept=0.95), linetype="longdash", color="grey") +
    geom_line() +
    facet_grid(rows=dplyr::vars(surv_true), cols=dplyr::vars(distr_A)) +
    ylim(0,0.003) +
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
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
  # distr_A_ <- "Unif(0,1)"
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



##################################################.
##### MISC: Process edge estimation mini-sim #####
##################################################.

if (FALSE) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_edge_20210921.rds")
  
  # Summarize results
  summ <- sim %>% summarize(
    bias_pct = list(
      name = "bias_pct",
      estimate = "theta_est",
      truth = "theta_true"
    ),
    coverage = list(
      name = "coverage",
      truth = "theta_true",
      lower = "ci_lo",
      upper = "ci_hi"
      # na.rm = TRUE
    )
  )
  summ
  
}



###################################.
##### VIZ: Sample paths (CIs) #####
###################################.

if (FALSE) {
  
  # Read in simulation object
  # sim <- readRDS("../SimEngine.out/sim_est_20210921.rds") # edge_corr="none"
  sim <- readRDS("../SimEngine.out/sim_est_20210922.rds") # edge_corr="max"
  
  # Filter data
  distr_A_ <- "N(0.5,0.04)"
  surv_true_ <- "complex"
  d <- sim$results %>% filter(
    n==5000 & distr_A==distr_A_ & surv_true==surv_true_
  )
  d %<>% filter(estimator=="Grenander (Est S_n/g_n)")
  if (length(unique(d$level_id))!=1) { stop("adjust filtering") }
  
  # Set up vector containers
  theta_true <- c()
  theta_est <- c()
  ci_lo <- c()
  ci_hi <- c()
  which <- c()
  
  # Extract simulation data into vectors
  n_paths <- 6
  row_offset <- 10
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=1)
    theta_true <- c(theta_true, d[1,paste0("theta_",m)])
  }
  for (i in 1:n_paths) {
    for (j in c(1:51)) {
      m <- format(round(j/50-0.02,2), nsmall=1)
      theta_est <- c(theta_est, d[i+row_offset,paste0("est_",m)])
      ci_lo <- c(ci_lo, d[i+row_offset,paste0("ci_lo_",m)])
      ci_hi <- c(ci_hi, d[i+row_offset,paste0("ci_hi_",m)])
      which <- c(which, i)
    }
  }
  
  plot_data <- data.frame(
    x = rep(sim$constants$points, n_paths),
    y = theta_est,
    ci_lo = ci_lo,
    ci_hi = ci_hi,
    which = which
  )
  ggplot(
    plot_data,
    aes(x=x, y=y, fill=factor(which), color=factor(which))
  ) +
    geom_vline(
      xintercept = c(qnorm(0.1,0.5,0.2),qnorm(0.9,0.5,0.2)),
      color = "orange",
      linetype = "dashed"
    ) +
    geom_line(
      data = data.frame(x=sim$constants$points, y=theta_true),
      color = "#222222"
    ) +
    geom_line() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_hi),
      alpha = 0.2,
      linetype = "dotted"
    ) +
    facet_wrap(~which, ncol=3) +
    # xlim(c(0.2,0.8)) +
    ylim(c(0,0.01)) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
}



###################################.
##### VIZ: Sample paths (all) #####
###################################.

if (FALSE) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_est_20211004.rds")
  
  # Filter data
  # distr_A_ <- "N(0.5,0.04)"
  # surv_true_ <- "complex"
  # d <- sim$results %>% filter(
  #   n==5000 & distr_A==distr_A_ & surv_true==surv_true_
  # )
  # d %<>% filter(estimator=="Grenander (Est S_n/g_n)")
  d <- sim$results
  if (length(unique(d$level_id))!=1) { stop("adjust filtering") }
  
  # Set up vector containers
  theta_true <- c()
  theta_est <- c()
  which <- c()
  
  # Extract simulation data into vectors
  n_paths <- 500 # nrow(sim$results)
  row_offset <- 0
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=1)
    theta_true <- c(theta_true, d[1,paste0("theta_",m)])
  }
  for (i in 1:n_paths) {
    for (j in c(1:51)) {
      m <- format(round(j/50-0.02,2), nsmall=1)
      theta_est <- c(theta_est, d[i+row_offset,paste0("est_",m)])
      which <- c(which, i)
    }
  }
  
  plot_data <- data.frame(
    x = rep(sim$constants$points, n_paths),
    y = theta_est,
    which = which
  )
  ggplot(
    plot_data,
    aes(x=x, y=y, group=which)
  ) +
    geom_line(alpha=0.05) + # 0.4
    geom_line(
      data = data.frame(x=sim$constants$points, y=theta_true),
      aes(x=x, y=y),
      color = "white",
      inherit.aes = FALSE
    ) +
    geom_vline(
      xintercept = c(qnorm(0.1,0.5,0.2),qnorm(0.9,0.5,0.2)),
      color = "orange",
      linetype = "dashed"
    ) +
    # scale_x_continuous(breaks=seq(0.5,0.8,0.01), limits=c(0.5,0.8)) +
    ylim(c(0,1)) + # 0.4
    theme(legend.position="none")
  
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
