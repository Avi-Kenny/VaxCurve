
# !!!!! Need to update this file based on analysis_janssen.R

#################.
##### Setup #####
#################.

# Set configuration
cfg2 <- list(
  which_day = 29, # 29 57
  which_marker = "ID50", # "spike" "ID50"
  run_import_data = F,
  run_recomp_1 = T,
  run_recomp_2 = T,
  run_recomp_3 = T,
  run_recomp_4 = T,
  run_recomp_5 = T,
  run_dqa = F,
  run_graphs = T
)
cfg2$folder = paste0("moderna (SL, ",cfg2$which_marker,"-",cfg2$which_day,")")

# Create directory if it doesn't exist
if (!dir.exists(cfg2$folder)) { dir.create(cfg2$folder) }

# Helper function: save or load functions
save_or_load <- function(fns, folder, which) {
  if (which=="load") {
    for (fn in fns) { assign(fn, readRDS(paste0(folder,"/",fn,".rds"))) }
  } else if (which=="save") {
    for (fn in fns) { saveRDS(eval(as.name(fn)), paste0(folder,"/",fn,".rds")) }
  }
}



###########################.
##### Data processing #####
###########################.

if (cfg2$run_import_data) {
  
  # Read in raw data; DO NOT STORE THIS LOCALLY !!!!!
  df_raw <- read.csv(paste0("Z:/covpn/p3001/analysis/correlates/Part_A_Blinded",
                            "_Phase_Data/adata/P3001ModernaCOVEimmunemarkerdat",
                            "a_correlates_processed_v1.0_Oct28_2021.csv"))
  # df_raw <- read.csv(paste0("Z:/covpn/p3001/download_data/Moderna COVE mRNA 12",
  #                           "73P301_immune_20210915/moderna_real_data_processe",
  #                           "d_with_riskscore_Weiping.csv"))
  
  # Save control group data
  dat_ctrl <- filter(df_raw, Trt==0)
  
  # Filter to include only treatment group
  df_raw %<>% filter(Trt==1)
  
  # Create data structure for analysis
  dat_orig <- list(
    "id" = df_raw[["Ptid"]],
    "y_star" = df_raw[[paste0("EventTimePrimaryD",cfg2$which_day)]],
    "delta_star" = df_raw[[paste0("EventIndPrimaryD",cfg2$which_day)]],
    "w" = data.frame(
      "w1" = df_raw[["MinorityInd"]],
      "w2" = df_raw[["standardized_risk_score"]],
      "w3" = df_raw[["HighRiskInd"]]
    ),
    "weights" = df_raw[[paste0("wt.D",cfg2$which_day)]],
    "delta" = as.integer(df_raw[[paste0("ph2.D",cfg2$which_day)]])
  )
  
  if (cfg2$which_marker=="ID50") {
    dat_orig[["a"]] <- df_raw[[paste0("Day",cfg2$which_day,"pseudoneutid50")]]
  } else if (cfg2$which_marker=="spike") {
    dat_orig[["a"]] <- df_raw[[paste0("Day",cfg2$which_day,"bindSpike")]]
  }
  
  # Stabilize weights (rescale to sum to sample size)
  dat_orig$weights <- ifelse(dat_orig$delta==1, dat_orig$weights, 0)
  s <- sum(dat_orig$weights) / length(dat_orig$delta)
  dat_orig$weights <- dat_orig$weights / s
  
  saveRDS(dat_orig, file=paste0(cfg2$folder,"/dat_orig_moderna-",
                                cfg2$which_marker,"-",cfg2$which_day,".rds"))
  saveRDS(dat_ctrl, file=paste0(cfg2$folder,"/dat_ctrl_moderna.rds"))
  
} else {
  
  dat_orig <- readRDS(paste0(cfg2$folder,"/dat_orig_moderna-",
                             cfg2$which_marker,"-",cfg2$which_day,".rds"))
  dat_ctrl <- readRDS(paste0(cfg2$folder,"/dat_ctrl_moderna.rds"))
  
}



###############################.
##### Data quality checks #####
###############################.

if (cfg2$run_dqa) {
  
  # print(paste("# with missing covariate (MinorityInd):",
  #             sum(is.na(df_raw$MinorityInd))))
  # print(paste("# with missing covariate (standardized_risk_score):",
  #             sum(is.na(df_raw$standardized_risk_score))))
  
  # print(paste("# patients without missing weights:", nrow(df_raw)))
  
  # Missing data checks
  sum(is.na(dat_orig$id))
  sum(is.na(dat_orig$y_star))
  sum(is.na(dat_orig$delta_star))
  sum(is.na(dat_orig$a))
  sum(is.na(dat_orig$w$w1))
  sum(is.na(dat_orig$w$w2))
  sum(is.na(dat_orig$w$w3))
  sum(is.na(dat_orig$weights))
  sum(is.na(dat_orig$delta))
  
  # Missing data checks
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  sum(is.na(dat$id))
  sum(is.na(dat$y_star))
  sum(is.na(dat$delta_star))
  sum(is.na(dat$a))
  sum(is.na(dat$w$w1))
  sum(is.na(dat$w$w2))
  sum(is.na(dat$w$w3))
  sum(is.na(dat$weights))
  sum(is.na(dat$delta))
  
  # Cross-check `delta` with `a`
  # !!!!! Right now, this is showing 47 records for which delta==0 and !is.na(a)
  xtabs(~dat_orig$delta+is.na(dat_orig$a))
  xtabs(~df_raw$SubcohortInd+df_raw$TwophasesampIndD57)
  dat_orig$id[which(dat_orig$delta==0 & !is.na(dat_orig$a))]
  
  # Compare minority status variables
  xtabs(~df_raw$URMforsubcohortsampling, addNA=T)
  xtabs(~df_raw$MinorityInd, addNA=T)
  
  # Additional quality checks
  c(min(df$a),max(df$a))
  hist <- function(x){ ggplot(data.frame(x=x),aes(x=x))+geom_histogram() }
  sum(is.na(df$y_star))
  sum(is.na(df$delta_star))
  sum(is.na(df$a))
  sum(is.na(df$w1))
  sum(is.na(df$w2))
  df_raw$Day29pseudoneutid50 %>% hist()
  df_raw$Day57pseudoneutid50 %>% hist()
  
}



#########################.
##### Data analysis #####
#########################.

{
  # !!!!! Eventually, replace this with a call to est_curve()
  
  # Set up end time of interest
  C <- list(appx=cfg$appx)
  C$t_e <- ifelse(cfg2$which_day==29, 126, ifelse(cfg2$which_day==57, 100, NA))
  
  # Rescale A to lie in [0,1]
  a_shift <- -1 * floor(10*min(dat_orig$a, na.rm=T))/10
  dat_orig$a <- dat_orig$a + a_shift
  a_scale <- ceiling(10*max(dat_orig$a, na.rm=T))/10
  dat_orig$a <- dat_orig$a / a_scale
  
  # Round A, W2
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$w$w2 <- round(dat_orig$w$w2, 0) # !!!!!
  
  # Create truncated data object
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  
  # Create val_list
  vlist <- create_val_list(dat, C$appx)
  
  # Set estimation tuning parameters
  params <- list(S_n_type="Super Learner", g_n_type="binning",
                 ecdf_type="linear (mid)", deriv_type="m-spline",
                 gamma_type="kernel", ci_type="trunc")
  
  # Calculate control group survival
  var_time <- paste0("EventTimePrimaryD",cfg2$which_day)
  var_delta <- paste0("EventIndPrimaryD",cfg2$which_day)
  srv <- survfit(
    formula(paste0("Surv(",var_time,",",var_delta,")~1")),
    data = dat_ctrl
  )
  rate_ctrl <- 1 - srv$surv[which.min(abs(srv$time-C$t_e))]
  
}

# Construct/save component functions (SECTION 1)
fns <- c("Phi_n","Phi_n_inv","S_n","Sc_n","f_aIw_n","f_a_n","g_n","omega_n")
if (cfg2$run_recomp_1) {
  Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
  Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
  S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
  Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
  f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15) # !!!!! Also do k=0 for cross-validated selection of k
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  g_n <- construct_g_n(f_aIw_n, f_a_n)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  save_or_load(fns, cfg2$folder, "save")
} else {
  save_or_load(fns, cfg2$folder, "load")
}

# Construct/save component functions (SECTION 2)
fns <- c("Gamma_os_n")
if (cfg2$run_recomp_2) {
  Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n) # type="plug-in"
  save_or_load(fns, cfg2$folder, "save")
} else {
  save_or_load(fns, cfg2$folder, "load")
}

# Construct/save component functions (SECTION 3)
fns <- c("Psi_n","gcm","dGCM","theta_n_Gr","theta_n")
if (cfg2$run_recomp_3) {
  Psi_n <- Vectorize(function(x) {
    Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
  })
  gcm <- gcmlcm(x=seq(0,1,C$appx$a), y=Psi_n(seq(0,1,C$appx$a)), type="lcm")
  dGCM <- Vectorize(function(x) {
    index <- which(round(x,5)<=gcm$x.knots)[1]-1
    if (index==0) { index <- 1 }
    return(gcm$slope.knots[index])
  })
  theta_n_Gr <- Vectorize(function(x) { min(dGCM(Phi_n(x)),1) })
  theta_n <- theta_n_Gr
  save_or_load(fns, cfg2$folder, "save")
} else {
  save_or_load(fns, cfg2$folder, "load")
}

# Construct/save component functions (SECTION 4)
fns <- c("f_aIw_delta1_n","f_a_delta1_n","gamma_n","deriv_theta_n","tau_n")
if (cfg2$run_recomp_4) {
  f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type,
                                      k=15, delta1=TRUE) # !!!!! Also do k=0 for cross-validated selection of k
  f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                  f_aIw_delta1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, vlist$A_grid,
                               type=params$gamma_type, omega_n, f_aIw_n,
                               f_a_n, f_a_delta1_n)
  deriv_theta_n <- construct_deriv_theta_n(theta_n, type=params$deriv_type,
                                           dir="decr")
  tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
  save_or_load(fns, cfg2$folder, "save")
} else {
  save_or_load(fns, cfg2$folder, "load")
}

# Construct/save component functions (SECTION 5)
fns <- c("gcomp")
if (cfg2$run_recomp_5) {
  S_n2 <- construct_S_n(dat, vlist$S_n, type="Cox PH")
  gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n2)
  save_or_load(fns, cfg2$folder, "save")
} else {
  save_or_load(fns, cfg2$folder, "load")
}

if (cfg2$run_graphs) {
  
  # CVE graph
  {
    
    # Generate point estimates
    grid <- seq(0,1,0.01)
    cve <- Vectorize(function(x) { 1 - x/rate_ctrl })
    ests_gcomp <- gcomp(grid)
    ests_gren <- theta_n(grid)
    
    # Generate CIs
    tau_ns <- tau_n(grid)
    qnt <- 1.00
    n_orig <- length(dat_orig$delta)
    if (params$ci_type=="regular") {
      ci_lo <- cve(ests_gren - qnt * (tau_ns/(n_orig^(1/3))))
      ci_hi <- cve(ests_gren + qnt * (tau_ns/(n_orig^(1/3))))
    } else if (params$ci_type=="logit") {
      ci_lo <- expit(
        logit(cve(ests_gren)) -
          qnt * (tau_ns*deriv_logit(cve(ests_gren)))/(rate_ctrl*n_orig^(1/3))
      )
      ci_hi <- expit(
        logit(cve(ests_gren)) +
          qnt * (tau_ns*deriv_logit(cve(ests_gren)))/(rate_ctrl*n_orig^(1/3))
      )
    } else if (params$ci_type=="trunc") {
      ci_lo <- cve(ests_gren - qnt * (tau_ns/(n_orig^(1/3)))) %>%
        pmax(0) %>% pmin(1)
      ci_hi <- cve(ests_gren + qnt * (tau_ns/(n_orig^(1/3)))) %>%
        pmax(0) %>% pmin(1)
    }
    ests_gcomp <- cve(ests_gcomp)
    ests_gren <- cve(ests_gren)
    
    # Marginal distribution of A
    df_marg <- data.frame(
      x = grid*a_scale-a_shift,
      ymin = 0,
      ymax = f_a_n(grid)*0.12
    )
    
    # Truncate at 5/95 quantiles
    which <- grid>=Phi_n_inv(0.05) & grid<=Phi_n_inv(0.95)
    ests_gren <- ifelse(which,ests_gren,NA)
    ci_lo <- ifelse(which,ci_lo,NA)
    ci_hi <- ifelse(which,ci_hi,NA)
    
    # Labels
    if (cfg2$which_marker=="ID50") {
      x_lab <- "Pseudovirus-nAb cID50 (=s)"
    } else if (cfg2$which_marker=="spike") {
      x_lab <- "Anti Spike IgG (=s)"
    }
    if (cfg2$which_day==29) {
      y_lab <- "Controlled VE against COVID-19 by day 126"
    } else if (cfg2$which_day==57) {
      y_lab <- "Controlled VE against COVID-19 by day 100"
    }
    
    # Plot
    # Export: PNG 600 x 500
    print(ggplot(
      data.frame(
        x = rep(grid*a_scale-a_shift,2),
        y = c(ests_gcomp,ests_gren),
        which = rep(c("G-comp", "Grenander"), each=length(grid)),
        ci_lo = c(ests_gcomp,ci_lo),
        ci_hi = c(ests_gcomp,ci_hi)
      ),
      aes(x=x, y=y, color=which)
    ) +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.1,
        linetype = "dotted",
        fill = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_color_manual(values=c("purple", "darkblue")) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1.05),
                         breaks=seq(0,1,0.1)) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color="Estimator") +
      geom_line())
      
  }
  
  if (F) {
    
    # Gamma_os_n graph
    # Export: 6" x 5"
    grid <- seq(0,1,C$appx$a)
    df_marg <- data.frame(
      x = grid*a_scale-a_shift,
      ymin = 0,
      ymax = f_a_n(grid)*0.001
    )
    ggplot(
      data.frame(
        x = grid*a_scale-a_shift,
        y = Gamma_os_n(grid)
      ),
      aes(x=x, y=y)
    ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      theme(legend.position="bottom") +
      labs(x="Day 57 pseudoneut-ID50", y="Gamma_os_n",
           color="Estimator") +
      geom_line()
    
    # theta graph
    # Export: 6" x 5"
    grid <- seq(0,1,0.01)
    df_marg <- data.frame(
      x = grid*a_scale-a_shift,
      ymin = 0,
      ymax = f_a_n(grid)*0.01
    )
    ests_gcomp <- gcomp(grid)
    ests_gren <- theta_n(grid)
    ggplot(
      data.frame(
        x = rep(grid*a_scale-a_shift,2),
        y = c(ests_gcomp,ests_gren),
        which = rep(c("G-comp", "Grenander"), each=length(grid))
      ),
      aes(x=x, y=y, color=which)
    ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_color_manual(values=c("purple", "darkblue")) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,0.05),
                         breaks=seq(0,0.05,0.01)) +
      theme(legend.position="bottom") +
      labs(x="Day 57 pseudoneut-ID50", y="theta_n",
           color="Estimator") +
      geom_line()
    
  }
  
}

# Component function check: conditional survival function
if (F) {
  
  # !!!!! Look at what happens with/without supplying vlist$S_n
  # S_n_model <- construct_S_n(dat, vals=NA, type="Cox PH", return_model=T)
  # S_n <- construct_S_n(dat, vals=NA, type="Cox PH")
  # S_n <- construct_S_n(dat, vals=vlist$S_n, type="Cox PH")
  S_n(t=100,w=c(0,0,0),a=0.1)
  S_n(t=100,w=c(0,0,0),a=0.9)
  S_n(t=100,w=c(1,1,1),a=0.1)
  S_n(t=100,w=c(1,1,1),a=0.9)
  Sc_n(t=99,w=c(0,0,0),a=0.1)
  Sc_n(t=99,w=c(0,0,0),a=0.9)
  Sc_n(t=99,w=c(1,1,1),a=0.1)
  Sc_n(t=99,w=c(1,1,1),a=0.9)
  
}

# Component function check: conditional distribution
if (F) {

  grid <- seq(0,1,0.01)
  len <- length(grid)
  plot_data <- data.frame(
    a = rep(grid, 4),
    density = c(
      sapply(grid, function(a) { f_aIw_n(a, w=c(0.2,0,1)) }),
      sapply(grid, function(a) { f_aIw_n(a, w=c(0.8,0,1)) }),
      sapply(grid, function(a) { f_aIw_n(a, w=c(0.2,1,1)) }),
      sapply(grid, function(a) { f_aIw_n(a, w=c(0.8,1,1)) })
    ),
    covariates = c(
      rep("W1=0.2, W2=0, w3=1",len),
      rep("W1=0.8, W2=0, w3=1",len),
      rep("W1=0.2, W2=1, w3=1",len),
      rep("W1=0.8, W2=1, w3=1",len)
    )
  )
  ggplot(plot_data, aes(x=a, y=density)) +
    geom_line() +
    facet_wrap(~covariates, ncol=4) +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of conditional density: f(A|W)") +
    ylim(c(0,NA))

}

# Component function check: marginal distribution
if (F) {

  grid <- seq(0,1,0.01)
  len <- length(grid)
  plot_data <- data.frame(a=grid, density=f_a_n(grid))
  ggplot(plot_data, aes(x=a, y=density)) +
    geom_line() +
    theme(legend.position="bottom") +
    labs(color="Estimator", title="Estimation of marginal density: f(A)") +
    ylim(c(0,NA))

}
