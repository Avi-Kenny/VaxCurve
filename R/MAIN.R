# Title: "Nonparametric inference for controlled vaccine efficacy curves"
# Author: Avi Kenny, Marco Carone

##################.
##### CONFIG #####
##################.

# Set global config
# devtools::install_github("tedwestling/ctsCausal")
# devtools::install_github("zeehio/facetscales")
cfg <- list(
  which_sim = "estimation", # estimation testing
  level_set_which = "level_set_estimation_1", # level_set_estimation_1 level_set_testing_1
  run_or_update = "run",
  num_sim = 1000,
  pkgs = c("dplyr", "boot", "car", "mgcv", "memoise", "EnvStats",
           "fdrtool", "splines", "survival"), # "ranger", "ctsCausal", "SuperLearner", "earth", "Rsolnp", "sets"
  pkgs_nocluster = c("ggplot2", "viridis", "sqldf", "facetscales", "scales",
                     "data.table", "latex2exp", "tidyr"),
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
{
  library(simba)
  source("one_simulation.R", local=TRUE)
  source("generate_data.R", local=TRUE)
  source("est_curve.R", local=TRUE)
  source("test_2.R", local=TRUE)
  source("fns_doseresp.R", local=TRUE)
}



#############################################.
##### TESTING: Dataset for testing code #####
#############################################.

if (FALSE) {
  
  L <- list(n=1000, alpha_3=0.75, distr_A="Unif(0,1)",
            reg_true="Logistic", sampling="iid",
            # test=list(type="test_2",params=list(boot_reps=2))
            estimator=list(
              # est="Generalized Grenander",params=list(ci_type="logit")))
              est="G-comp",params=list(mu_n_type="Random forest", boot_reps=5)))
  C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3,
            points=seq(0,1,0.1))
  
  # Generate dataset
  dat <- generate_data(L$n, L$alpha_3, L$distr_A, L$reg_true, L$sampling)
  
}



##########################################################.
##### MAIN: Set level sets for different simulations #####
##########################################################.

if (Sys.getenv("simba_run") %in% c("first", "")) {
  
  # !!!!! Temp/testing !!!!!
  level_set_temp <- list(
    n = 500,
    alpha_3 = 0.75,
    distr_A = "Unif(0,1)",
    reg_true = "Logistic",
    sampling = "iid", # two-phase iid
    estimator = list(
      "G-comp" = list(
        est = "G-comp",
        params = list(mu_n_type="Logistic", boot_reps=10) # 100
      )
      # "Grenander (logit CIs)" = list(
      #   est = "Generalized Grenander",
      #   params = list(mu_n_type="Logistic", g_n_type="parametric", ci_type="logit") # mu_n_type="GAM"
      # )
      # "Grenander (split CIs, m=5)" = list(
      #   est = "Generalized Grenander",
      #   params = list(mu_n_type="GAM", g_n_type="parametric",
      #                 ci_type="sample split", m=5)
      # )
    )
  )
  
  # Estimation: compare all methods
  # Not currently using (ci_type="sample split", m=5) or (ci_type="regular")
  level_set_estimation_1 <- list(
    n = 5000,
    alpha_3 = 0.75,
    distr_A = c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)"), # "Beta(0.8+0.9*w1,0.8+0.4*w2)"
    reg_true = c("Logistic", "GAM"),
    sampling = c("iid", "two-phase"), # iid two-phase
    estimator = list(
      "G-comp" = list(
        est = "G-comp",
        params = list(mu_n_type="GAM", boot_reps=100)
      ),
      "Grenander (logit CIs)" = list(
        est = "Generalized Grenander",
        params = list(mu_n_type="GAM", g_n_type="parametric", ci_type="logit")
      )
      # "Grenander (split CIs, m=5)" = list(
      #   est = "Generalized Grenander",
      #   params = list(mu_n_type="GAM", g_n_type="parametric",
      #                 ci_type="sample split", m=5)
      # )
    )
  )
  
  # Testing: compare all methods
  # Not currently using (var="boot")
  level_set_testing_1 <- list(
    n = 200,
    # n = 1000,
    alpha_3 = 0.75,
    # alpha_3 = c(0,0.25,0.5), # c(0,0.25,0.5,0.75)
    distr_A = "Unif(0,1)",
    # distr_A = c("Unif(0,1)", "Beta(0.9,1.1+0.4*w2)"), # "Beta(0.8+0.9*w1,0.8+0.4*w2)"
    reg_true = "GAM",
    # reg_true = c("Logistic", "GAM"),
    sampling = "two-phase", # iid
    test = list(
      "Slope (mixed boot)" = list(
        type = "test_2",
        params = list(
          var = "mixed boot",
          mu_n_type = "GAM",
          g_n_type = "parametric",
          boot_reps = 3
          # boot_reps = 100
        )
      )
      # "Wald" = list(type="test_wald", params=list())
      # "Westling 2020" = list(type="test_causalnull", params=list()) # !!!!! Not yet ready
    )
  )
  
  # !!!!! Estimation/testing but with different nuisance estimators
  # mu_n_type = c("Logistic", "GAM", "Random forest")
  # g_n_type = c("parametric", "binning")
  
  level_set <- eval(as.name(cfg$level_set_which))
  
}



##########################################.
##### MAIN: Setup and run simulation #####
##########################################.

# Use these commands to run on Slurm:
# sbatch --export=simba_run='first',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:11 --array=1-1200 --export=simba_run='main',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh
# sbatch --depend=afterok:12 --export=simba_run='last',cluster='bionic',type='R',project='z.VaxCurve' -e ./io/slurm-%A_%a.out -o ./io/slurm-%A_%a.out --constraint=gizmok run_r.sh

if (cfg$run_or_update=="run") {

  run_on_cluster(
    
    first = {
      
      # Simulation setup
      sim <- new_sim()
      sim %<>% set_config(
        num_sim = cfg$num_sim,
        parallel = cfg$parallel,
        stop_at_error = cfg$stop_at_error,
        seed = 4,
        packages = cfg$pkgs
      )
      sim <- do.call(set_levels, c(list(sim), level_set))
      
      # Add functions to simulation object
      sim %<>% add_creator(generate_data)
      methods <- c(
        "est_curve", "test_2", "test_wald", "expit", "deriv_expit", "logit",
        "deriv_logit", "construct_S_n", "construct_deriv_theta_n",
        "construct_sigma2_n", "construct_f_a_n", "construct_f_aIw_n",
        "construct_g_n", "construct_Gamma_n", "construct_Phi_n", "Pi", "wts",
        "construct_gcomp", "ss"
      )
      for (method in methods) {
        sim %<>% add_method(method, eval(as.name(method)))
      }
      
      # Add constants
      # lambda and v are the Weibull parameters for the survival distribution
      # lambda2 and v2 are the Weibull parameters for the censoring distribution
      {
        sim %<>% add_constants(
          lambda = 10^(-4),
          v = 1.5,
          lambda2 = 0.5 * 10^(-4),
          v2 = 1.5,
          points = seq(0,1,0.1),
          # alpha_0 = -1.5,
          alpha_1 = 0.3,
          alpha_2 = 0.7,
          # alpha_4 = -0.3,
          t_e = 200
        )
      }
      
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



###################################.
##### VIZ: Estimation (sim 1) #####
###################################.

if (FALSE) {
  
  # Read in simulation object
  sim <- readRDS("../simba.out/sim_est_20210722_mod.simba")
  
  # Summarize results
  # !!!!! Count the number of NAs in coverage
  summ_bias <- list()
  summ_mse <- list()
  summ_cov <- list()
  for (i in c(1:11)) {
    m <- format(round(i/10-0.1,1), nsmall=1)
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
  
  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,alpha_3,distr_A,reg_true,sampling,Estimator),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  
  cb_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                 "#0072B2", "#D55E00", "#CC79A7", "#999999")
  m_colors <- c(
    `G-comp` = cb_colors[1],
    `Grenander (logit CIs)` = cb_colors[2],
    `Grenander (split CIs, m=5)` = cb_colors[3]
    # `Grenander (regular CIs)` = cb_colors[4],
  )
  
  # Bias plot
  # Export: 8" x 5"
  p_data %<>% filter(sampling=="two-phase") # !!!!!
  ggplot(
    filter(p_data, stat=="bias"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(reg_true)) +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values=m_colors) +
    labs(title="Bias (%)", x="A", y=NULL, color="Estimator")
  
  # Coverage plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="cov"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_hline(aes(yintercept=0.95), linetype="longdash", color="grey") +
    geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(reg_true)) +
    scale_y_continuous(labels=percent) +
    scale_color_manual(values=m_colors) +
    labs(title="Coverage (%)", x="A", y=NULL, color="Estimator")

  # MSE plot
  # Export: 8" x 5"
  ggplot(
    filter(p_data, stat=="mse"),
    aes(x=point, y=value, color=Estimator, group=Estimator)
  ) +
    geom_point() +
    geom_line() +
    facet_grid(rows=dplyr::vars(distr_A), cols=dplyr::vars(reg_true)) +
    scale_color_manual(values=m_colors) +
    labs(title="MSE", x="A", y=NULL, color="Estimator") +
    ylim(0,0.0115)
  
  # # Export: 6" x 4"
  # distr_A_ <- "Beta(0.9,1.1+0.4*w2)" # Unif(0,1) Beta(0.9,1.1+0.4*w2)
  # estimand_ <- "Midpoint" # Midpoint Endpoint
  # ggplot(
  #   filter(p_data, distr_A==distr_A_ & estimand==estimand_),
  #   aes(x=n, y=value, color=Estimator)
  # ) +
  #   geom_hline(
  #     aes(yintercept=y),
  #     data=data.frame(y=0.95, stat="Coverage"),
  #     linetype="longdash", color="grey"
  #   ) +
  #   geom_point() +
  #   geom_line() +
  #   # geom_bar(stat="identity", position=position_dodge(),
  #   #          width=0.8, color="white", size=0.35) +
  #   facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
  #                 scales=list(y=list(
  #                   Bias = scale_y_continuous(labels = percent_format()),
  #                   Coverage = scale_y_continuous(labels = percent_format()),
  #                   MSE = scale_y_continuous()
  #                 ))) +
  #   theme(legend.position="bottom") +
  #   # scale_color_manual(values=m_colors) +
  #   labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
  #        y=NULL, x=NULL, color="Estimator")
  
  
  # # !!!!! Temp (export 7.5 x 4.5)
  # distr_A_ <- "Beta(0.9,1.1+0.4*w2)" # Unif(0,1) Beta(0.9,1.1+0.4*w2)
  # estimand_ <- "Endpoint" # Midpoint Endpoint
  # ggplot(
  #   filter(p_data, distr_A==distr_A_ & estimand==estimand_),
  #   aes(x=sampling, y=value, fill=Estimator)
  # ) +
  #   geom_hline(
  #     aes(yintercept=y),
  #     data=data.frame(y=0.95, stat="Coverage"),
  #     linetype="longdash", color="grey"
  #   ) +
  #   geom_bar(stat="identity", position=position_dodge(),
  #            width=0.8, color="white", size=0.35) +
  #   facet_grid_sc(cols=dplyr::vars(reg_true), rows=dplyr::vars(stat),
  #                 scales=list(y=list(
  #                   Bias = scale_y_continuous(labels = percent_format()),
  #                   Coverage = scale_y_continuous(labels = percent_format()),
  #                   MSE = scale_y_continuous()
  #                 ))) +
  #   theme(legend.position="bottom") +
  #   scale_fill_manual(values=m_colors) +
  #   labs(title=paste0("Estimand: ",estimand_,"; MargDist(A): ",distr_A_),
  #        y=NULL, x=NULL, color="Estimator")
  
}



#########################################.
##### VIZ: Estimation (method demo) #####
#########################################.

if (FALSE) {
  
  # Set parameters
  n <- 1000
  alpha_0 <- -1.5
  alpha_1 <- 0.3
  alpha_2 <- 0.7
  alpha_3 <- 0.7
  alpha_4 <- -0.3
  reg_true <- "Logistic" # Logistic GAM Complex
  grid <- seq(0,1,0.1)
  # grid <- seq(0,1,0.01)
  
  # Generate dataset
  dat <- generate_data(
    n = n,
    alpha_3 = alpha_3,
    distr_A = "Beta(0.9,1.1+0.4*w2)",
    reg_true = reg_true,
    sampling = "two-phase"
  )
  
  # Approximate true values of theta_0
  m <- 10000
  w1 <- rnorm(m)
  w2 <- rbinom(m, size=1, prob=0.5)
  true_vals <- sapply(grid, function(x) {
    if (reg_true=="Logistic") {
      mean(expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*x))
    } else if (reg_true=="GAM") {
      mean(expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*sqrt(x)))
    } else if (reg_true=="Complex") {
      mean(expit(alpha_0 + alpha_1*sin(2*pi*w1) + alpha_2*w2 +
                   alpha_3*sqrt(x) + alpha_4*w1*w2))
    }
  })
  
  # Estimate curve: G-comp
  ests <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    estimator = "G-comp",
    list(boot_reps=100),
    points = grid
  )))
  ests$which <- "G-comp"
  
  # Estimate curve: Grenander (logit CIs)
  ests2 <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    estimator = "Generalized Grenander",
    params = list(ci_type="logit"),
    points = grid
  )))
  ests2$which <- "Grenander (logit CIs)"
  
  # Estimate curve: Grenander (regular CIs)
  ests3 <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    estimator = "Generalized Grenander",
    params = list(ci_type="regular"),
    points = grid
  )))
  ests3$which <- "Grenander (regular CIs)"
  
  # Estimate curve: Grenander (split CIs; m=5)
  ests4 <- as.data.frame(rbindlist(est_curve(
    dat = dat,
    estimator = "Generalized Grenander",
    params = list(ci_type="sample split", m=5),
    points = grid
  )))
  ests4$which <- "Grenander (split CIs; m=5)"
  
  ests <- rbind(ests, ests2, ests3, ests4)
  
  # Attach true values
  ests <- rbind(ests, data.frame(
    point = grid,
    est = true_vals,
    ci_lo = true_vals,
    ci_hi = true_vals,
    which = rep("theta_0", length(grid))
  ))
  
  # # Plot results
  # # Export: 8" x 3"
  # ggplot(
  #   ests,
  #   aes(x=point, y=est, color=which, fill=which)) +
  #   geom_line() +
  #   # geom_point() +
  #   geom_ribbon(
  #     aes(ymin=ci_lo, ymax=ci_hi), # fill=as.factor(group_var)
  #     alpha = 0.2,
  #     linetype = "dotted"
  #   ) +
  #   facet_wrap(~which, ncol=3) +
  #   labs(
  #     x = "x",
  #     y = unname(latex2exp::TeX("$\\hat{\\theta}_n(x)$")),
  #     color = "Which",
  #     fill = "Which"
  #   )
  
  # Plot results
  # Export: 8" x 3"
  d2 <- filter(ests, which=="theta_0")
  d2_nrows <- nrow(d2)
  d2 <- rbind(d2,d2,d2,d2)
  d2$which <- rep(unique(filter(ests, which!="theta_0")$which), each=d2_nrows)
  ggplot(
    filter(ests, which!="theta_0"),
    aes(x=point, y=est, color=which, fill=which)) +
    geom_line() +
    # geom_point() +
    geom_ribbon(
      aes(ymin=ci_lo, ymax=ci_hi), # fill=as.factor(group_var)
      alpha = 0.2,
      linetype = "dotted"
    ) +
    facet_wrap(~which, ncol=4) +
    theme(legend.position="bottom") +
    labs(
      x = "x",
      y = unname(latex2exp::TeX("$\\hat{\\theta}_n(x)$")),
      color = "Which",
      fill = "Which"
    ) +
    geom_line(
      data = d2,
      color = "black",
      linetype = "dashed"
    )
  
  # end_true <- (filter(ests, point==1 & which=="theta_0"))$est
  # end_gcomp <- (filter(ests, point==1 & which=="G-comp"))$est
  # print((end_true-end_gcomp)^2)
  
}



################################.
##### VIZ: Testing (sim 1) #####
################################.

if (FALSE) {
  
  # # Read in simulation object
  sim <- readRDS("../simba.out/sim_testing_20210705.simba")
  
  # !!!!! Modify everything below
  
  # Summarize resuls
  summ <- summarize(sim)
  
  summ %<>% rename(
    "Power" = mean_reject,
    "Test" = test
  )
  
  cb_colors <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
                 "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  m_colors <- c(
    `Slope (boot)` = cb_colors[2],
    `Slope (mixed boot)` = cb_colors[3],
    `Wald` = cb_colors[4]
  )
  
  # Export: 7" x 4.5"
  # distr_A_ <- "Beta(0.9,1.1+0.4*w2)" # Unif(0,1) Beta(0.9,1.1+0.4*w2)
  alpha_3_ <- 0.7 # 0 0.7
  sampling_ <- "two-phase" # iid two-phase
  ggplot(
    filter(summ, alpha_3==alpha_3_ & sampling==sampling_),
    aes(x=n, y=Power, color=Test)
  ) +
    geom_point() +
    geom_line() +
    facet_grid(cols=dplyr::vars(reg_true), rows=dplyr::vars(distr_A)) +
    scale_y_continuous(labels=percent) +
    theme(legend.position="bottom") +
    scale_color_manual(values=m_colors) +
    labs(title = paste0("alpha_3: ",alpha_3_,"; Sampling: ",sampling_),
         color = "Test")
  
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



#############################################################.
##### TESTING: Conditional survival function estimators #####
#############################################################.

if (FALSE) {
  
  # Generate data
  C <- list(points=c(0.5,0.8))
  dat <- generate_data(n=5000, alpha_3=0.7, distr_A="Unif(0,1)",
                       surv_true="CoxPH", sampling="two-phase")
  
  
  S_n <- construct_S_n(dat, type="Cox")
  
  S_0 <- Vectorize(function(t, w1, w2, a) {
    lambda <- 10^(-4)
    v <- 1.5
    alpha <- c(0.3,0.7,0.5)
    alpha_3 <- 0.7
    lin <- alpha[1]*w1 + alpha[2]*w2 + alpha_3*a
    return(exp(-1*lambda*(t^v)*exp(lin)))
  })
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 8),
    survival = c(
      S_n(t=times, w1=0, w2=1, a=0),
      S_n(t=times, w1=1, w2=1, a=0),
      S_n(t=times, w1=0, w2=1, a=1),
      S_n(t=times, w1=1, w2=1, a=1),
      S_0(t=times, w1=0, w2=1, a=0),
      S_0(t=times, w1=1, w2=1, a=0),
      S_0(t=times, w1=0, w2=1, a=1),
      S_0(t=times, w1=1, w2=1, a=1)
    ),
    which = rep(c("Cox","True S_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0,a=0","w1=1,a=0","w1=0,a=1","w1=1,a=1"),2),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: S_0[t|W,A]",
         color="Estimator")
    
}



##############################################################.
##### TESTING: Conditional censoring function estimators #####
##############################################################.

if (FALSE) {
  
  # Generate data
  C <- list(points=c(0.5,0.8))
  dat <- generate_data(n=5000, alpha_3=0.7, distr_A="Unif(0,1)",
                       surv_true="CoxPH", sampling="two-phase")
  
  Sc_n <- construct_S_n(dat, type="Cox", csf=TRUE)
  
  Sc_0 <- Vectorize(function(t, w1, w2, a) {
    lambda2 <- 0.5 * 10^(-4)
    v2 <- 1.5
    alpha <- c(0.3,0.7,0.5)
    alpha_3 <- 0.7
    lin <- alpha[1]*w1 + alpha[2]*w2 + alpha_3*a
    return(exp(-1*lambda2*(t^v2)*exp(lin)))
  })
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 8),
    survival = c(
      Sc_n(t=times, w1=0, w2=1, a=0),
      Sc_n(t=times, w1=1, w2=1, a=0),
      Sc_n(t=times, w1=0, w2=1, a=1),
      Sc_n(t=times, w1=1, w2=1, a=1),
      Sc_0(t=times, w1=0, w2=1, a=0),
      Sc_0(t=times, w1=1, w2=1, a=0),
      Sc_0(t=times, w1=0, w2=1, a=1),
      Sc_0(t=times, w1=1, w2=1, a=1)
    ),
    which = rep(c("Cox","True S^C_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0,a=0","w1=1,a=0","w1=0,a=1","w1=1,a=1"),2),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: S_0[t|W,A]",
         color="Estimator")
  
}



######################################.
##### TESTING: theta_n estimator #####
######################################.

if (FALSE) {
  
  # Generate data
  C <- list(points=c(0.5,0.8))
  dat <- generate_data(n=5000, alpha_3=0.7, distr_A="Unif(0,1)",
                       surv_true="CoxPH", sampling="two-phase")
  
  Sc_n <- construct_S_n(dat, type="Cox", csf=TRUE)
  
  Sc_0 <- Vectorize(function(t, w1, w2, a) {
    lambda2 <- 0.5 * 10^(-4)
    v2 <- 1.5
    alpha <- c(0.3,0.7,0.5)
    alpha_3 <- 0.7
    lin <- alpha[1]*w1 + alpha[2]*w2 + alpha_3*a
    return(exp(-1*lambda2*(t^v2)*exp(lin)))
  })
  
  # Plot true curve against estimated curve
  times <- c(1:200)
  df <- data.frame(
    time = rep(times, 8),
    survival = c(
      Sc_n(t=times, w1=0, w2=1, a=0),
      Sc_n(t=times, w1=1, w2=1, a=0),
      Sc_n(t=times, w1=0, w2=1, a=1),
      Sc_n(t=times, w1=1, w2=1, a=1),
      Sc_0(t=times, w1=0, w2=1, a=0),
      Sc_0(t=times, w1=1, w2=1, a=0),
      Sc_0(t=times, w1=0, w2=1, a=1),
      Sc_0(t=times, w1=1, w2=1, a=1)
    ),
    which = rep(c("Cox","True S^C_0"), each=4*length(times)),
    covs = rep(rep(c("w1=0,a=0","w1=1,a=0","w1=0,a=1","w1=1,a=1"),2),
               each=length(times))
  )
  ggplot(df, aes(x=time, y=survival, color=which)) +
    geom_line() +
    facet_wrap(~covs, ncol=2) +
    labs(title="Estimation of conditional survival: S_0[t|W,A]",
         color="Estimator")
  
}



###################################################.
##### TESTING: Conditional density estimators #####
###################################################.

if (FALSE) {
  
  # Set levels here
  n <- 5000
  distr_A <- "Beta(0.8+0.9*w1,0.8+0.4*w2)"  # Unif(0,1) Beta(0.9,1.1+0.4*w2) Beta(0.8+0.9*w1,0.8+0.4*w2)
  sampling <- "two-phase"                     # iid two-phase
  
  # Generate data
  C <- list(alpha_0=-1.5, alpha_1=0.3, alpha_2=0.7, alpha_4=-0.3)
  dat <- generate_data(
    n = n,
    alpha_3 = 0.7,
    distr_A = distr_A,
    reg_true = "Logistic",
    sampling = sampling
  )
  
  # True conditional density function
  f_aIw_0 <- function(a,w1,w2) {
    if (distr_A=="Unif(0,1)") {
      return(1)
    } else if (distr_A=="Beta(0.9,1.1+0.4*w2)") {
      shape1 <- 0.9
      shape2 <- 1.1 + 0.4*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    } else if (distr_A=="Beta(0.8+0.9*w1,0.8+0.4*w2)") {
      shape1 <- 0.8 + 0.9*w1
      shape2 <- 0.8 + 0.4*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    }
  }
  
  # Parametric estimate
  {
    dat_trunc <- filter(dat, !is.na(a))
    n_trunc <- nrow(dat_trunc)
    weights <- wts(dat_trunc, scale="none")
    wlik <- function(par) {
      
      sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
        shape1 <- par[1] + par[2]*dat_trunc$w1[i]
        shape2 <- par[3] + par[4]*dat_trunc$w2[i]
        loglik <- dbeta(dat_trunc$a[i], shape1=shape1, shape2=shape2, log=TRUE)
        return(loglik*weights[i])
      }))
      
      return(-1*sum_loglik)
      
    }
    opt <- optim(par=c(a1=0.5, a2=0.1, a3=0.5, a4=0.1), fn=wlik)
    f_aIw_n_para <- Vectorize(function(a, w1, w2){
      shape1 <- opt$par[1] + opt$par[2]*w1
      shape2 <- opt$par[3] + opt$par[4]*w2
      return(dbeta(a, shape1=shape1, shape2=shape2))
    })
  }
  
  # Nonparametric estimate (adapted from Diaz & VDL)
  {
    
    # k is the number of bins
    construct_f_aIw_n_nonpar <- function(dat, k) {
      
      alphas <- seq(0, 1, length.out=k+1)
      dat_trunc <- filter(dat, !is.na(a))
      n_trunc <- nrow(dat_trunc)
      weights <- wts(dat_trunc, scale="none")
      
      dens <- Vectorize(function(a, w1, w2, par) {
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        hz <- sapply(c(1:(ifelse(bin==k,k-1,bin))), function(j) {
          expit(par[j] + par[k]*w1 + par[k+1]*w2)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        dens <- k*p1*p2
        return(dens)
      }, vectorize.args=c("a","w1","w2"))
      
      wlik <- function(par) {
        
        # par[1] through par[k-1] are the hazard components for the bins 1 to k-1
        # par[k] and par[k+1] correspond to W1 and W2
        sum_loglik <- sum(sapply(c(1:n_trunc), function(i) {
          lik <- dens(a=dat_trunc$a[i], w1=dat_trunc$w1[i], w2=dat_trunc$w2[i],
                      par)
          return(weights[i]*log(lik))
        }))
        
        return(-1*sum_loglik)
        
      }
      
      opt <- optim(par=rep(0,k+1), fn=wlik, method="CG")
      if (opt$convergence!=0) {
        warning("Nonpar conditional density: optim() did not converge")
      }
      
      return(Vectorize(memoise(function(a, w1, w2){
        
        bin <- ifelse(a==1, k, which.min(a>=alphas)-1)
        par <- opt$par
        hz <- sapply(c(1:(k-1)), function(j) {
          expit(par[j] + par[k]*w1 + par[k+1]*w2)
        })
        p1 <- ifelse(bin==k, 1, hz[bin])
        p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
        
        return(k*p1*p2)
        
      })))
      
    }
    
    f_aIw_n_nonpar <- construct_f_aIw_n_nonpar(dat, k=10)
    
  }
  
  # Generate plot data
  grid <- seq(0.01,0.99,0.01)
  f_aIw_models <- c("Truth", "Parametric", "Nonparametric")
  n_models <- length(f_aIw_models)
  len <- length(grid)
  plot_data <- data.frame(
    a = rep(grid, 4*n_models),
    density = c(
      sapply(grid, function(a) { f_aIw_0(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_0(a, w1=0.8, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_para(a, w1=0.8, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.2, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.8, w2=0) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.2, w2=1) }),
      sapply(grid, function(a) { f_aIw_n_nonpar(a, w1=0.8, w2=1) })
    ),
    which = rep(f_aIw_models, each=len*4),
    covariates = rep(c(
      rep("W1=0.2, W2=0",len),
      rep("W1=0.8, W2=0",len),
      rep("W1=0.2, W2=1",len),
      rep("W1=0.8, W2=1",len)
    ), n_models)
  )
  ggplot(plot_data, aes(x=a, y=density, color=factor(which))) +
    geom_line() +
    facet_wrap(~covariates, ncol=4) +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of conditional density: f(A|W)") +
    ylim(c(0,NA))
  
}
