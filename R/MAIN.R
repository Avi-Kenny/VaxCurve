# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# To run multiple sims/analyses concurrently, ONLY change Slurm commands

# Set global config
# GitHub packages: tedwestling/ctsCausal, tedwestling/CFsurvival,
#                  tedwestling/survSuperLearner, zeehio/facetscales
cfg <- list(
  main_task = "run", # run update analysis.R
  which_sim = "estimation", # "estimation" "edge" "testing" "Cox" "debugging"
  level_set_which = "level_set_estimation_1", # level_set_estimation_1 level_set_testing_1 level_set_Cox_1
  # keep = c(1:3,7:9,16:18,22:24),
  num_sim = 500,
  pkgs = c("dplyr", "boot", "car", "mgcv", "memoise", "EnvStats", "fdrtool",
           "splines", "survival", "SuperLearner", "survSuperLearner",
           "randomForestSRC", "CFsurvival", "Rsolnp", "truncnorm", "tidyr",
           "ranger", "xgboost", "survey", "pbapply", "compiler", "simest"),
  pkgs_nocluster = c("ggplot2", "viridis", "sqldf", "facetscales", "scales",
                     "data.table", "latex2exp"),
  parallel = "none",
  stop_at_error = F,
  appx = list(t_0=1, x_tol=25, s=0.01) # !!!!! t_0=1, a=0.001
)

# Set cluster config
cluster_config <- list(
  # js = "ge",
  # dir = paste0("/home/users/avikenny/Desktop/", Sys.getenv("project"))
  js = "slurm",
  dir = paste0("/home/akenny/", Sys.getenv("project"))
)



#################.
##### SETUP #####
#################.

# Set local vs. cluster variables
if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
  # Local
  setwd(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Resear",
               "ch/Marco Carone/Project - VaxCurve/VaxCurve/R"))
  load_pkgs_local <- T
} else {
  # Cluster
  setwd(paste0(cluster_config$dir, "/R"))
  if (cfg$main_task %in% c("run", "update")) {
    load_pkgs_local <- F
  } else {
    load_pkgs_local <- T
  }
}

# Load packages (if running locally)
if (load_pkgs_local) {
  for (pkg in c(cfg$pkgs,cfg$pkgs_nocluster)) {
    suppressMessages({ do.call("library", list(pkg)) })
  }
}

# Load SimEngine + functions
{
  library(SimEngine)
  source("one_simulation.R", local=T)
  source("generate_data.R", local=T)
  source("est_curve.R", local=T)
  source("test_2.R", local=T)
  source("fns_doseresp.R", local=T)
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("sim_run") %in% c("first", "")) {
  
  # Estimation: ideal params
  level_set_estimation_1 <- list(
    n = 1000, # 1000
    alpha_3 = -2,
    dir = c("decr"), # "incr"
    # sc_params = list("no cens"=list(lmbd=1e-3, v=1.5, lmbd2=5e-7, v2=1.5)),
    # sc_params = list("exp"=list(lmbd=1e-3, v=1.5, lmbd2=5e-4, v2=1.5)),
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_S = c("Unif(0,1)", "N(0.5,0.04)", "N(0.3+0.4x2,0.04)"),
    edge = c("none"), # "expit 0.4"
    # surv_true = c("Complex"),
    surv_true = c("Cox PH", "Complex"), # "Cox PH" "Complex" "exp"
    sampling = c("two-phase (50%)"), # "iid"
    wts_type = c("estimated"), # c("true", "estimated")
    estimator = list(
      "Grenander (Cox)" = list(
        est = "Grenander",
        params = list(q_n_type="zero", Q_n_type="Cox PH", # Q_n_type="Random Forest", Q_n_type="true"
                      convex_type="GCM", ecdf_type="linear (mid)", # convex_type="LS"
                      edge_corr="none", # edge_corr="min"
                      deriv_type="m-spline", g_n_type="parametric") # g_n_type="binning", g_n_type="true"
      )
      # "Cox PH" = list(est="Cox gcomp")
      # "Qbins (true)" = list(
      #   est = "Qbins",
      #   params = list(n_bins=8, Q_n_type="Cox PH")
      # )
    )
  )
  
  # Estimation: trial params
  level_set_estimation_2 <- list(
    n = 15000,
    alpha_3 = -4,
    dir = c("incr", "decr"),
    sc_params = list("sc_params"=list(lmbd=3e-5, v=1.5, lmbd2=3e-5, v2=1.5)),
    distr_S = c("N(0.5,0.01)", "N(0.5,0.04)", "Unif(0,1)"),
    edge = "none",
    surv_true = "Cox PH",
    sampling = "two-phase (6%)",
    wts_type = "estimated",
    estimator = list(
      "Qbins (5)" = list(est="Qbins", params=list(n_bins=8, Q_n_type="Cox PH")),
      "Grenander" = list(
        est = "Grenander",
        params = list(q_n_type="new", Q_n_type="Cox PH")
      )
    )
  )
  
  # Testing: compare all methods
  level_set_testing_1 <- list(
    n = 1000,
    # n = c(100,200,400,800), # 1000
    # n = c(1000,2000),
    # alpha_3 = 0,
    alpha_3 = c(0,-0.25,-0.5),
    dir = "decr",
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_S = c("Unif(0,1)", "N(0.5,0.04)"), # "N(0.5,0.01)"
    edge = c("none", "expit 0.5"),
    surv_true = "Cox PH",
    sampling = c("iid", "two-phase (72%)"),
    # sampling = c("iid", "x1", "x2", "two-phase (50%)", "two-phase (50% random)"),
    # sampling = c("x1", "x2", "two-phase (50%)"),
    # sampling = c("iid", "cycle", "two-phase (72%)", "two-phase (70% random)"),
    wts_type = "true",
    test = list(
      "Slope (two-tailed, MC)" = list(
        type = "test_2",
        alt_type = "two-tailed", # decr
        params = list(g_n_type="binning", Q_n_type="Cox PH",
                      omega_n_type="estimated"),
        test_stat_only = F
      )
      # "Slope (two-tailed, boot)" = list(
      #   type = "test_2",
      #   alt_type = "two-tailed", # decr
      #   test_stat_only = F,
      #   params = list(g_n_type="binning", Q_n_type="Cox PH",
      #                 omega_n_type="estimated", var="boot", boot_reps=100)
      # )
    )
  )
  
  # Estimation: ideal params
  level_set_Cox_1 <- list(
    n = 500,
    alpha_3 = -2,
    dir = "decr",
    # wts_type = "estimated",
    sc_params = list("sc_params"=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5)),
    distr_S = c("Unif(0,1)"),
    # distr_S = c("Unif(0,1)", "N(0.5,0.01)", "N(0.5,0.04)"),
    edge = "none",
    sampling = "two-phase (50%)",
    wts_type = c("true", "estimated")
    # sampling = c("iid", "two-phase (72%)", "two-phase (50%)", "two-phase (25%)")
  )
  
  level_set <- get(cfg$level_set_which)
  
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
  
  # Set global constants
  C <- list(
    # points = round(seq(0,1,0.1),2), # !!!!!
    points = round(seq(0,1,0.02),2),
    alpha_1 = 0.5,
    alpha_2 = 0.7,
    t_0 = 200,
    appx = cfg$appx
  )
  
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
      
      # Simulation script
      sim %<>% set_script(one_simulation)
      
    },
    
    main = { sim %<>% run() },
    
    last = { sim %>% summarize() %>% print() },
    
    cluster_config = cluster_config
    
  )
  
} else if (cfg$main_task=="update") {
  
  update_sim_on_cluster(
    
    first = {
      sim <- readRDS(paste0(cluster_config$dir,"/sim.rds"))
      sim <- do.call(set_levels, c(list(sim), level_set))
    },
    
    main = { sim %<>% update_sim() },
    
    last = {},
    
    cluster_config = cluster_config
    
  )
  
} else {
  
  source(cfg$main_task, local=T)
  
}



###################################.
##### VIZ: Estimation (sim 1) #####
###################################.

if (F) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_est_20211004.rds")
  
  # Summarize results
  summ_bias <- list()
  if (F) {
    summ_biasG <- list()
    summ_biasP <- list()
  } # DEBUG: Gamma/Phi
  summ_mse <- list()
  summ_var <- list()
  summ_cov <- list()
  for (i in c(1:51)) {
  # for (i in c(1:5)) {
    m <- format(round(i/50-0.02,2), nsmall=2)
    # m <- format(round(i/4-0.25,2), nsmall=2)
    summ_bias[[i]] <- list(
      name = paste0("bias_",m),
      estimate = paste0("r_Mn_",m),
      truth = paste0("r_M0_",m)
    )
    if (F) {
      summ_biasG[[i]] <- list(
        name = paste0("biasG_",m),
        estimate = paste0("estG_",m),
        truth = paste0("Gamma_",m)
      )
      summ_biasP[[i]] <- list(
        name = paste0("biasP_",m),
        estimate = paste0("estP_",m),
        truth = paste0("Phi_",m)
      )
    } # DEBUG: Gamma/Phi
    summ_var[[i]] <- list(
      name = paste0("var_",m),
      x = paste0("r_Mn_",m)
    )
    summ_mse[[i]] <- list(
      name = paste0("mse_",m),
      estimate = paste0("r_Mn_",m),
      truth = paste0("r_M0_",m)
    )
    summ_cov[[i]] <- list(
      name = paste0("cov_",m),
      truth = paste0("r_M0_",m),
      lower = paste0("ci_lo_",m),
      upper = paste0("ci_hi_",m),
      na.rm = T
    )
  }
  summ <- summarize(sim, bias=summ_bias, mse=summ_mse, var=summ_var, coverage=summ_cov)
  if (F) {
    summ <- summarize(sim, bias=c(summ_bias,summ_biasG,summ_biasP), mse=summ_mse, coverage=summ_cov)
  } # DEBUG: Gamma/Phi
  
  summ %<>% rename("Estimator"=estimator)
  
  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,alpha_3,sc_params,distr_S,edge,
              surv_true,sampling,Estimator,dir,wts_type),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point = as.numeric(point))
  
  # cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
  #                "#0072B2", "#D55E00", "#CC79A7", "#999999")
  # m_colors <- c(cb_colors[2], cb_colors[5], cb_colors[6]
    # `Grenander (Est Q_n/g_n)` = cb_colors[2],
  # )
  
  # PLot Y-axis limits
  plot_lims <- list(b=c(-0.25,0.25), c=c(0,1), m=c(0,0.02), v=c(0,0.01)) # c=c(0.7,1)
  
  # Set faceting vectors
  distr_Ss <- c("Unif(0,1)", "N(0.5,0.04)", "N(0.3+0.4x2,0.04)")
  surv_trues <- c("Cox PH", "Complex")
  
  # Orange 10/90 quantile lines
  x2 <- rbinom(10^5, size=1, prob=0.5)
  q3 <- as.numeric(quantile(
    rtruncnorm(10^5, a=0, b=1, mean=0.3+0.4*x2, sd=0.2), c(0.1,0.9)
  ))
  df_vlines <- data.frame(
    # x = c(qnorm(0.1,0.5,0.1),qnorm(0.1,0.5,0.2),qunif(0.1,0,1),
    #       qnorm(0.9,0.5,0.1),qnorm(0.9,0.5,0.2),qunif(0.9,0,1)),
    # distr_S = rep(c("N(0.5,0.01)", "N(0.5,0.04)", "Unif(0,1)"),2)
    x = c(qunif(0.1,0,1), qtruncnorm(0.1, a=0, b=1, mean=0.5, sd=0.2), q3[1],
          qunif(0.9,0,1), qtruncnorm(0.9, a=0, b=1, mean=0.5, sd=0.2), q3[2]),
    distr_S = rep(distr_Ss,2)
  )
  
  # Grey background densities
  df_distr_S <- data.frame(
    x = rep(seq(0,1,0.01),3),
    # ymax = c(dnorm(seq(0,1,0.01), mean=0.5, sd=0.1),
    #          dnorm(seq(0,1,0.01), mean=0.5, sd=0.2),
    #          rep(1,101)),
    # distr_S = rep(c("N(0.5,0.01)", "N(0.5,0.04)", "Unif(0,1)"), each=101),
    ymax = c(rep(1,101),
             dtruncnorm(seq(0,1,0.01), a=0, b=1, mean=0.5, sd=0.2),
             0.5*dtruncnorm(seq(0,1,0.01), a=0, b=1, mean=0.3, sd=0.2) +
             0.5*dtruncnorm(seq(0,1,0.01), a=0, b=1, mean=0.7, sd=0.2)),
    distr_S = rep(distr_Ss, each=101),
    value = 0
  )
  if (sim$levels$edge!="none") {
    mass <- as.numeric(strsplit(sim$levels$edge," ",fixed=T)[[1]][2])
    height <- 10 * mass
    df_distr_S %<>% mutate(
      ymax = ifelse(x<0.1, 10*mass, (1-mass)*df_distr_S$ymax)
    )
    df_vlines <- df_vlines[4:6,]
  }
  df_distr_b <- mutate(df_distr_S, ymin=plot_lims$b[1],
                       ymax=((ymax*diff(plot_lims$b))/6+plot_lims$b[1]))
  df_distr_c <- mutate(df_distr_S, ymin=plot_lims$c[1],
                       ymax=((ymax*diff(plot_lims$c))/6+plot_lims$c[1]))
  df_distr_m <- mutate(df_distr_S, ymin=plot_lims$m[1],
                       ymax=((ymax*diff(plot_lims$m))/6+plot_lims$m[1]))
  df_distr_v <- mutate(df_distr_S, ymin=plot_lims$v[1],
                       ymax=((ymax*diff(plot_lims$v))/6+plot_lims$v[1]))
  
  # Bias plot
  # Export: 10" x 6"
  # Note: change "bias" to "biasG" for Gamma and "biasP" for Phi
  ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value, color=factor(Estimator), group=factor(Estimator))
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
                data=df_distr_b, fill="grey", color=NA, alpha=0.4) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_line() +
    facet_grid(rows = dplyr::vars(factor(surv_true, levels=surv_trues)),
               cols = dplyr::vars(factor(distr_S, levels=distr_Ss))) +
    scale_y_continuous(limits=plot_lims$b) + # labels=percent
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(y="Bias", x="S", color="Estimator")
  
  # Coverage plot
  # Export: 10" x 6"
  ggplot(
    filter(p_data, stat=="cov"),
    aes(x=point, y=value, color=factor(Estimator), group=factor(Estimator))
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
                data=df_distr_c, fill="grey", color=NA, alpha=0.4) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_hline(aes(yintercept=0.95), linetype="longdash", color="grey") +
    geom_line() +
    facet_grid(rows = dplyr::vars(factor(surv_true, levels=surv_trues)),
               cols = dplyr::vars(factor(distr_S, levels=distr_Ss))) +
    scale_y_continuous(labels=percent, limits=plot_lims$c) +
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(y="Coverage (%)", x="S", color="Estimator")
  
  # Variance plot
  # Export: 10" x 6"
  ggplot(
    filter(p_data, stat=="var"),
    aes(x=point, y=value, color=factor(Estimator), group=factor(Estimator))
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
                data=df_distr_v, fill="grey", color=NA, alpha=0.4) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_line() +
    facet_grid(rows = dplyr::vars(factor(surv_true, levels=surv_trues)),
               cols = dplyr::vars(factor(distr_S, levels=distr_Ss))) +
    scale_y_continuous(limits=plot_lims$v) +
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(y="Variance", x="S", color="Estimator")
  
  # # MSE plot
  # # Export: 10" x 6"
  # ggplot(
  #   filter(p_data, stat=="mse"),
  #   aes(x=point, y=value, color=factor(Estimator), group=factor(Estimator))
  # ) +
  #   geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax, color=NA, group=NA),
  #               data=df_distr_m, fill="grey", color=NA, alpha=0.4) +
  #   geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
  #              linetype="dashed") +
  #   geom_line() +
  #   facet_grid(rows = dplyr::vars(factor(surv_true, levels=surv_trues)),
  #              cols = dplyr::vars(factor(distr_S, levels=distr_Ss))) +
  #   scale_y_continuous(limits=plot_lims$m) +
  #   # scale_color_manual(values=m_colors) +
  #   theme(legend.position="bottom") +
  #   labs(y="MSE", x="S", color="Estimator")
  
}



################################.
##### VIZ: Testing (sim 1) #####
################################.

if (F) {
  
  # # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_testing_20210820.rds")
  
  # !!!!! Modify everything below
  #   color should be n-value
  #   x-axis should be alpha_3
  
  # Summarize resuls
  summ <- summarize(sim)
  
  summ %<>% rename("Power"=mean_reject)
  
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
  # distr_S_ <- "Unif(0,1)"
  ggplot(
    summ,
    aes(x=alpha_3, y=Power, color=factor(n))
  ) +
    geom_point() +
    geom_line() +
    facet_grid(cols=dplyr::vars(test), rows=dplyr::vars(distr_S)) +
    # facet_grid(cols=dplyr::vars(surv_true), rows=dplyr::vars(distr_S)) +
    scale_y_continuous(labels=percent) +
    theme(legend.position="bottom") +
    scale_color_manual(values=m_colors) +
    labs(title="Testing", color="Sample size")
    # labs(title = paste0("alpha_3: ",alpha_3_,"; Sampling: ",sampling_),
    #      color = "Test")
  
  # Diagnostics
  if (F) {
    
    # Filter (if necessary)
    r <- filter(sim$results, level_id==1)
    r1 <- filter(sim$results, level_id==1)
    r2 <- filter(sim$results, level_id==2)
    r3 <- filter(sim$results, level_id==3)
    
    # Print rejection rate
    mean(r$reject)
    
    # Histograms
    ggplot(data.frame(x=r$beta_n), aes(x=x)) +
      geom_histogram(bins=50) +
      labs(title="beta_n")
    ggplot(data.frame(x=r$Gamma_n_5), aes(x=x)) +
      geom_histogram(bins=50) +
      labs(title="Gamma_n_5") +
      xlim(0.24,0.38)
    # ggplot(data.frame(x=r$p_vals), aes(x=x)) +
    #   geom_histogram(bins=50) +
    #   labs(title="p_vals")
    ggplot(data.frame(x=r2$sd_n), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=sd(r2$beta_n), color="forestgreen", linetype="dashed") +
      labs(title="sd(beta_n)")
    ggplot(data.frame(x=sqrt(r2$Gamma_var_n)), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=sd(r2$Gamma_n_5), color="forestgreen", linetype="dashed") +
      labs(title="sd(Gamma_n(0.5))")
    ggplot(data.frame(x=r$Phi_n_5), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=0.5, color="forestgreen", linetype="dashed") +
      labs(title="Phi_n(0.5)") # + xlim(0.4,0.6)
    
    # Comparing component function distributions
    names <- c("two-phase (72%)", "two-phase (70% random)")
    cmp <- "sd_n" # Gamma_n_5 Gamma_var_n sd_n beta_n
    ggplot(data.frame(x=c(r1[[cmp]],r2[[cmp]]),
                      grp=c(rep(names[1], nrow(r1)),rep(names[2], nrow(r2)))),
           aes(x=x)) +
      facet_wrap(~grp, ncol=2) +
      geom_histogram(bins=25, alpha=0.7) +
      geom_vline(aes(xintercept=x), data.frame(
        x=c(quantile(r1[[cmp]])[2:4],quantile(r2[[cmp]])[2:4]),
        grp=rep(names,each=3)),
        color="grey", linetype="dashed") +
      labs(title=cmp)
    
    # Comparisons of SDs
    print(paste0("Actual SD(beta_n): ", sd(r1$beta_n)))
    print(paste0("Estimated SD(beta_n): ", mean(r1$sd_n)))
    print(paste0("Actual SD(Gamma_n): ", sd(r1$Gamma_n_5)))
    print(paste0("Estimated SD(Gamma_n): ", mean(sqrt(r1$Gamma_var_n))))
    
    # Checking infl_fn_1
    r1 <- filter(sim$results, level_id==1)
    sd(r1$Psi_1_n)
    mean(r1$Psi_1_sd)
    ggplot(data.frame(x=r1$Psi_1_sd), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=sd(r1$Psi_1_n), color="forestgreen", linetype="dashed") +
      labs(title="sd(partial_est)")
    
    # Checking infl_fn_4
    r1 <- filter(sim$results, level_id==1)
    sd(r1$Psi_4_n)
    mean(r1$Psi_4_sd)
    ggplot(data.frame(x=r1$Psi_4_sd), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=sd(r1$Psi_4_n), color="forestgreen", linetype="dashed") +
      labs(title="sd(partial_est)")
    
    # Checking infl_fn_5
    # True variance is (1/3)/n; SD is sqrt(1/(3*5000))
    sd(r1$Psi_5_n)
    mean(r1$Psi_5_sd)
    ggplot(data.frame(x=r1$Psi_5_sd), aes(x=x)) +
      geom_histogram(bins=50, alpha=0.7) +
      geom_vline(xintercept=sd(r1$Psi_5_n), color="forestgreen", linetype="dashed") +
      labs(title="sd(partial_est)")
    
  }
  
}



##################################################.
##### MISC: Process edge estimation mini-sim #####
##################################################.

if (F) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_edge_20210921.rds")
  
  # Summarize results
  summ <- sim %>% summarize(
    bias_pct = list(
      name = "bias_pct",
      estimate = "r_Mn",
      truth = "r_M0"
    ),
    coverage = list(
      name = "coverage",
      truth = "r_M0",
      lower = "ci_lo",
      upper = "ci_hi"
      # na.rm = T
    )
  )
  summ
  
}



#######################################.
##### MISC: Process Cox gcomp sim #####
#######################################.

if (F) {
  
  sim %>% SimEngine::summarize(
    mean = list(
      list(name="mean_runtime", x="runtime"),
      list(name="se_x1", x="se_x1"),
      list(name="se_x2", x="se_x2"),
      list(name="se_s", x="se_s"),
      list(name="se_est_bshz", x="se_est_bshz"),
      list(name="se_est_surv", x="se_est_surv"),
      list(name="se_est_marg", x="se_est_marg")
    ),
    sd = list(
      list(name="se_x1_empr", x="est_x1"),
      list(name="se_x2_empr", x="est_x2"),
      list(name="se_s_empr", x="est_s"),
      list(name="se_bshz_empr", x="est_bshz"),
      list(name="se_surv_empr", x="est_surv"),
      list(name="se_marg_empr", x="est_marg")
    ),
    coverage = list(
      list(name="cov_x1", truth="true_x1", estimate="est_x1", se="se_x1"),
      list(name="cov_x2", truth="true_x2", estimate="est_x2", se="se_x2"),
      list(name="cov_s", truth="true_s", estimate="est_s", se="se_s"),
      list(name="cov_bshz", truth="true_bshz", estimate="est_bshz", se="se_est_bshz"),
      list(name="cov_surv", truth="true_surv", estimate="est_surv", se="se_est_surv"),
      list(name="cov_marg", truth="true_marg", estimate="est_marg", se="se_est_marg")
    )
  )
  
}



###################################.
##### VIZ: Sample paths (CIs) #####
###################################.

if (F) {
  
  # Read in simulation object
  # sim <- readRDS("../SimEngine.out/sim_est_20210921.rds") # edge_corr="none"
  sim <- readRDS("../SimEngine.out/sim_est_20210922.rds") # edge_corr="max"
  
  # Filter data
  d <- sim$results
  d %<>% filter(level_id==1)
  
  # Set up vector containers
  r_M0 <- c()
  r_Mn <- c()
  ci_lo <- c()
  ci_hi <- c()
  which <- c()
  
  # Extract simulation data into vectors
  n_paths <- 6
  row_offset <- 10
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=1)
    r_M0 <- c(r_M0, d[1,paste0("r_M0_",m)])
  }
  for (i in 1:n_paths) {
    for (j in c(1:51)) {
      m <- format(round(j/50-0.02,2), nsmall=1)
      r_Mn <- c(r_Mn, d[i+row_offset,paste0("r_Mn_",m)])
      ci_lo <- c(ci_lo, d[i+row_offset,paste0("ci_lo_",m)])
      ci_hi <- c(ci_hi, d[i+row_offset,paste0("ci_hi_",m)])
      which <- c(which, i)
    }
  }
  
  plot_data <- data.frame(
    x = rep(sim$constants$points, n_paths),
    y = r_Mn,
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
      data = data.frame(x=sim$constants$points, y=r_M0),
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
    # ylim(c(0,0.01)) +
    theme(
      legend.position = "none",
      strip.background = element_blank(),
      strip.text.x = element_blank()
    )
  
}



###################################.
##### VIZ: Sample paths (all) #####
###################################.

if (F) {
  
  # Read in simulation object
  sim <- readRDS("../SimEngine.out/sim_est_20211004.rds")
  
  # Filter data
  d <- sim$results
  d %<>% filter(level_id==1)
  
  # Set up vector containers
  r_M0 <- c()
  # Theta_true <- c()
  r_Mn <- c()
  # Theta_est <- c()
  which <- c()
  
  # Extract simulation data into vectors
  n_paths <- nrow(d) # 500
  row_offset <- 0
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=2)
    r_M0 <- c(r_M0, d[1,paste0("r_M0_",m)])
    # Theta_true <- c(Theta_true, d[1,paste0("Theta_",m)])
  }
  for (i in 1:n_paths) {
    for (j in c(1:51)) {
      m <- format(round(j/50-0.02,2), nsmall=2)
      r_Mn <- c(r_Mn, d[i+row_offset,paste0("r_Mn_",m)])
      # Theta_est <- c(Theta_est, d[i+row_offset,paste0("esT_",m)])
      which <- c(which, i)
    }
  }
  
  # points <- get("points", envir=sim$vars$env)
  points <- C$points
  ggplot(
    data.frame(x=rep(points, n_paths), y=r_Mn, which=which),
    aes(x=x, y=y, group=which)
  ) +
    geom_line(alpha=0.05) +
    geom_line(
      data = data.frame(x=points, y=r_M0),
      aes(x=x, y=y),
      color = "white",
      inherit.aes = F
    ) +
    geom_vline(
      xintercept = c(qnorm(0.1,0.5,0.1),qnorm(0.9,0.5,0.1)),
      # xintercept = c(qunif(0.1,0,1),qunif(0.9,0,1)),
      # xintercept = c(qunif(0.1,0.3,0.7),qunif(0.9,0.3,0.7)),
      color = "orange",
      linetype = "dashed"
    ) +
    ylim(c(0,0.8)) + # 0.02
    theme(legend.position="none")
  
}



#######################################.
##### VIZ: Chernoff vs. N(0,0.52) #####
#######################################.

if (F) {
  
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
