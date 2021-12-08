# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

# !!!!! EventIndPrimaryD29 --> EventIndPrimaryIncludeNotMolecConfirmedD29
# !!!!! EventTimePrimaryD29 --> EventTimePrimaryIncludeNotMolecConfirmedD29

# !!!!! Standardize trt vs. tx, ctrl vs. ct, etc.

# 2*(10^0.7340943) # Should equal LOD
# 2*(10^0.847753) # Should equal LOD

# # !!!!!
# ggplot(data.frame(x=dat_orig_1$a_list[[1]]), aes(x=x)) + # page 20
#   geom_histogram(bins=40) +
#   xlim(0,3) +
#   labs(title="Day29bindSpike")
# ggplot(data.frame(x=dat_orig_1$a_list[[2]]), aes(x=x)) + # page 27
#   geom_histogram(bins=40) +
#   xlim(-1,3) +
#   labs(title="Day29bindRBD")

#################.
##### Setup #####
#################.

# Vector of markers
markers <- c("Day29bindSpike", "Day29bindRBD", "Day29ADCP")

# Set configuration
cfg2 <- list(
  # tid = 1,
  # marker_num = 1,
  tid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  marker_num = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  S_n_type = "Cox PH", # !!!!! "Cox PH" "Super Learner"
  run_import_data = F,
  run_recomp_1 = T,
  run_recomp_2 = T,
  run_recomp_3 = T,
  run_recomp_4 = T,
  run_recomp_5 = T,
  save_fns = F,
  run_dqa = F,
  run_graphs = T,
  edge_spread = F
)
cfg2$folder = paste0("Janssen functions/Janssen (tid ",cfg2$tid,")")
# cfg2$edge_spread <- ifelse(cfg2$marker_num==3, TRUE, FALSE) # !!!!!

# Create directory if it doesn't exist
if (!dir.exists(cfg2$folder)) { dir.create(cfg2$folder) }

# Helper function: save or load functions
save_or_load <- function(fns, folder, which) {
  if (which=="load") {
    for (fn in fns) {
      assign(fn, readRDS(paste0(folder,"/",fn,".rds")), envir=parent.frame())
    }
  } else if (which=="save") {
    for (fn in fns) {
      saveRDS(eval(as.name(fn)), paste0(folder,"/",fn,".rds"))
    }
  }
}



###########################.
##### Data processing #####
###########################.

if (cfg2$run_import_data) {
  
  df_raw_1 <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Blind",
                            "ed_Phase_Data/adata/janssen_pooled_real_data_proc",
                            "essed_with_riskscore.csv"))
  df_raw_adcp <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Bl",
                            "inded_Phase_Data/adata/janssen_pooled_realADCP_da",
                            "ta_processed_with_riskscore.csv"))
  
  # Subset data frames
  df_ph1_1 <- filter(df_raw_1, ph1.D29==T)
  df_ctrl_1 <- filter(df_ph1_1, Trt==0)
  df_trt_1 <- filter(df_ph1_1, Trt==1)
  df_ph1_adcp <- filter(df_raw_adcp, ph1.D29==T)
  df_ctrl_adcp <- filter(df_ph1_adcp, Trt==0)
  df_trt_adcp <- filter(df_ph1_adcp, Trt==1)
  
  # Create data structures for analysis
  dat_orig_1 <- list(
    "id" = df_trt_1[["Ptid"]],
    "y_star" = df_trt_1[["EventTimePrimaryIncludeNotMolecConfirmedD29"]],
    "delta_star" = df_trt_1[["EventIndPrimaryIncludeNotMolecConfirmedD29"]],
    "w" = data.frame(
      "w1" = as.integer(df_trt_1[["Region"]]==1),
      "w2" = as.integer(df_trt_1[["Region"]]==2),
      "w3" = df_trt_1[["risk_score"]]
    ),
    "weights" = df_trt_1[["wt.D29"]],
    "a_list" = list(df_trt_1[[markers[1]]], df_trt_1[[markers[2]]], NA),
    "delta" = as.integer(df_trt_1[["ph2.D29"]])
  )
  dat_orig_adcp <- list(
    "id" = df_trt_adcp[["Ptid"]],
    "y_star" = df_trt_adcp[["EventTimePrimaryIncludeNotMolecConfirmedD29"]],
    "delta_star" = df_trt_adcp[["EventIndPrimaryIncludeNotMolecConfirmedD29"]],
    "w" = data.frame(
      "w1" = as.integer(df_trt_adcp[["Region"]]==1),
      "w2" = as.integer(df_trt_adcp[["Region"]]==2),
      "w3" = df_trt_adcp[["risk_score"]]
    ),
    "weights" = df_trt_adcp[["wt.D29"]],
    "a_list" = list(NA, NA, df_trt_adcp[[markers[3]]]),
    "delta" = as.integer(df_trt_adcp[["ph2.D29"]])
  )
  
  # Stabilize weights (rescale to sum to sample size)
  dat_orig_1$weights <- ifelse(dat_orig_1$delta==1,
                               dat_orig_1$weights, 0)
  s <- sum(dat_orig_1$weights) / length(dat_orig_1$delta)
  dat_orig_1$weights <- dat_orig_1$weights / s
  dat_orig_adcp$weights <- ifelse(dat_orig_adcp$delta==1,
                                  dat_orig_adcp$weights, 0)
  s <- sum(dat_orig_adcp$weights) / length(dat_orig_adcp$delta)
  dat_orig_adcp$weights <- dat_orig_adcp$weights / s
  
  # Save datasets
  saveRDS(dat_orig_1, file="Janssen data/dat_orig_1_Janssen.rds")
  saveRDS(df_ph1_1, file="Janssen data/df_ph1_1_Janssen.rds")
  saveRDS(df_ctrl_1, file="Janssen data/df_ctrl_1_Janssen.rds")
  saveRDS(df_trt_1, file="Janssen data/df_trt_1_Janssen.rds")
  saveRDS(dat_orig_adcp, file="Janssen data/dat_orig_adcp_Janssen.rds")
  saveRDS(df_ph1_adcp, file="Janssen data/df_ph1_adcp_Janssen.rds")
  saveRDS(df_ctrl_adcp, file="Janssen data/df_ctrl_adcp_Janssen.rds")
  saveRDS(df_trt_adcp, file="Janssen data/df_trt_adcp_Janssen.rds")
  
} else {

  # Read datasets
  dat_orig_1 <- readRDS("Janssen data/dat_orig_1_Janssen.rds")
  df_ph1_1 <- readRDS("Janssen data/df_ph1_1_Janssen.rds")
  df_ctrl_1 <- readRDS("Janssen data/df_ctrl_1_Janssen.rds")
  df_trt_1 <- readRDS("Janssen data/df_trt_1_Janssen.rds")
  dat_orig_adcp <- readRDS("Janssen data/dat_orig_adcp_Janssen.rds")
  df_ph1_adcp <- readRDS("Janssen data/df_ph1_adcp_Janssen.rds")
  df_ctrl_adcp <- readRDS("Janssen data/df_ctrl_adcp_Janssen.rds")
  df_trt_adcp <- readRDS("Janssen data/df_trt_adcp_Janssen.rds")
  
}



###############################################.
##### Summary stats / data quality checks #####
###############################################.

if (cfg2$run_dqa) {
  
  # Alias vectors
  ind_tx <- df_trt_1$EventIndPrimaryIncludeNotMolecConfirmedD29
  ind_ct <- df_ctrl_1$EventIndPrimaryIncludeNotMolecConfirmedD29
  time_tx <- df_trt_1$EventTimePrimaryIncludeNotMolecConfirmedD29
  time_ct <- df_ctrl_1$EventTimePrimaryIncludeNotMolecConfirmedD29
  
  # Number of cases in each group
  num_case_tx <- sum(ind_tx)
  num_case_ct <- sum(ind_ct)
  num_case_tx_195 <- sum(ind_tx[time_tx<=195])
  num_case_ct_195 <- sum(ind_ct[time_ct<=195])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste("Number of cases in vaccine group:", num_case_tx))
  print(paste("Number of cases in control group:", num_case_ct))
  print(paste("Number of cases by day 195 in vaccine group:", num_case_tx_195))
  print(paste("Number of cases by day 195 in control group:", num_case_ct_195))
  print(paste("Number at-risk in vaccine group:", num_atrisk_tx))
  print(paste("Number at-risk in control group:", num_atrisk_ct))
  print(paste("Naive P(COVID by day 195) in vaccine group:",
              round(num_case_tx_195/num_atrisk_tx,3)))
  print(paste("Naive P(COVID by day 195) in control group:",
              round(num_case_ct_195/num_atrisk_ct,3)))
  print(paste("Naive vaccine efficacy:",
              round(1 - (num_case_tx_195/num_atrisk_tx) /
                      (num_case_ct_195/num_atrisk_ct),3)))
  
  # !!!!!
  get.marginalized.risk.no.marker <- function(dat) {
    fit.risk <- coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score +
        as.factor(Region),
      dat,
      model = T
    )
    dat[["EventTimePrimaryIncludeNotMolecConfirmedD29"]] <- 195
    risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
    mean(risks)
  }
  res.plac.cont <- get.marginalized.risk.no.marker(df_ctrl_1)
  res.vacc.cont <- get.marginalized.risk.no.marker(df_trt_1)
  overall.ve <- 1 - res.vacc.cont/res.plac.cont
  print(overall.ve)
  
  # srv_ct <- survfit(
  #   Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
  #        EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  #   data = df_ctrl_1
  # )
  # rate195_ct <- 1-srv_ct$surv[which(srv_ct$time==195)]
  # srv_tx <- survfit(
  #   Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
  #        EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  #   data = df_trt_1
  # )
  # rate195_tx <- 1-srv_tx$surv[which(srv_tx$time==195)]
  # print(rate195_ct)
  # print(rate195_tx)
  # print(1-(rate195_tx/rate195_ct))
  
  # Kaplan-Meier plot
  # !!!!!
  
  # Distribution of event times
  ggplot(
    data.frame(
      x = c(time_tx[which(ind_tx==1)], time_ct[which(ind_ct==1)]),
      which = c(rep("Tx",num_case_tx), rep("Ct",num_case_ct))
      ),
    aes(x=x, fill=which)
  ) +
    facet_wrap(~which) +
    geom_vline(xintercept=195, linetype="dashed", color="grey") +
    geom_histogram() +
    labs("Distribution of event times")

  # Distribution of censoring times
  ggplot(
    data.frame(
      x = c(time_tx[which(ind_tx==0)], time_ct[which(ind_ct==0)]),
      which = c(rep("Tx",num_atrisk_tx-num_case_tx),
                rep("Ct",num_atrisk_ct-num_case_ct))
    ),
    aes(x=x, fill=which)
  ) +
    facet_wrap(~which) +
    geom_vline(xintercept=195, linetype="dashed", color="grey") +
    geom_histogram() +
    labs("Distribution of event times")
  
  
  
  
  # # !!!!!
  # library(ggfortify)
  # survfit(
  #   Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
  #        EventIndPrimaryIncludeNotMolecConfirmedD29)~Trt,
  #   data = df_ph1_1
  # ) %>% autoplot()
  
  # Calculate control group survival
  srv_ct <- survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = df_ctrl
  )
  
  survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = filter(df_trt, ph2.D29==1),
    weights = wt.D29
  ) %>% autoplot()
  
}



##############################.
##### Data analysis prep #####
##############################.

{
  
  # Set the analysis dataframes
  if (cfg2$marker_num<=2) {
    dat_orig <- dat_orig_1
    df_ctrl <- df_ctrl_1
    df_trt <- df_trt_1
  } else {
    dat_orig <- dat_orig_adcp
    df_ctrl <- df_ctrl_adcp
    df_trt <- df_trt_adcp
  }
  
  # Set up end time of interest `t_e`
  C <- list(appx=cfg$appx, t_e=195)
  
  # !!!!! Reset appx for `t_e` and `a`
  C$appx$t_e <- 10
  C$appx$a <- 0.01 # !!!!!
  
  # Set `a` value from `a_list`
  dat_orig$a <- dat_orig$a_list[[cfg2$marker_num]]
  dat_orig$a_list <- NULL
  
  # Rescale A to lie in [0,1]
  a2 <- 1/C$appx$a
  a_shift <- -1 * floor(a2*min(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a + a_shift
  a_scale <- ceiling(a2*max(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a / a_scale

  # Round A, W3
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$w$w3 <- round(dat_orig$w$w2, 0) # !!!!! change to round(...,1)
  
  # Perform "spread" edge correction
  # !!!!! Should be symmetric around LLOD/2
  if (cfg2$edge_spread) {
    # width <- 0.02
    # noise_0 <- round(runif(length(dat_orig$a))*width, -log10(C$appx$a))
    # noise_1 <- round(1 - runif(length(dat_orig$a))*width, -log10(C$appx$a))
    # dat_orig$a <- ifelse(dat_orig$a==0, noise_0, dat_orig$a)
    # dat_orig$a <- ifelse(dat_orig$a==1, noise_1, dat_orig$a)
    # # dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  }

  # Set estimation tuning parameters
  params <- list(S_n_type=cfg2$S_n_type, g_n_type="binning",
                 ecdf_type="linear (mid)", deriv_type="m-spline",
                 gamma_type="kernel", ci_type="trunc", edge_corr="none",
                 cf_folds=1, n_bins=0)
  
}


#########################.
##### Data analysis #####
#########################.

# Create truncated data object
dat <- ss(dat_orig, which(dat_orig$delta==1))

# Obtain estimates
p_grid <- seq(0,1,0.01)
ests <- est_curve(
  dat_orig = dat_orig,
  estimator = "Grenander",
  params = params,
  points = p_grid,
  dir = "decr",
  return_extra = c("gcomp", "f_a_n", "Phi_n_inv")
)

# Calculate control group survival
srv_ct <- survfit(
  Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
       EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  data = df_ctrl
)
rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_e))]
ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-C$t_e))]
ci_hi_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-C$t_e))]
var_ct <- ((ci_hi_ct-ci_lo_ct)/3.92)^2

# Calculate treatment group survival
srv_tx <- survfit(
  Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
       EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  data = df_trt
)
rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-C$t_e))]
ci_lo_tx <- 1 - srv_tx$upper[which.min(abs(srv_tx$time-C$t_e))]
ci_hi_tx <- 1 - srv_tx$lower[which.min(abs(srv_tx$time-C$t_e))]
var_tx <- ((ci_hi_tx-ci_lo_tx)/3.92)^2

# Calculate overall vaccine efficacy
ve_overall <- 1 - (rate_tx/rate_ct)
ve_se <- sqrt(rate_ct^-2*var_tx + rate_tx^2*rate_ct^-4*var_ct)
ve_overall_lo <- ve_overall - 1.96*ve_se
ve_overall_hi <- ve_overall + 1.96*ve_se
print(paste0("Overall VE: ", round(100*ve_overall,1), "% (",
             round(100*ve_overall_lo,1), "% -- ", round(100*ve_overall_hi,1),
             "%)"))

# Calculate overall vaccine efficacy (within subcohort)
srv_tx_sub <- survfit(coxph(
  Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
       EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  data = filter(df_trt, ph2.D29==1),
  weights = wt.D29
))
rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_tx_sub$time-C$t_e))]
ve_subcohort <- 1 - (rate_tx_sub/rate_ct)
print(paste0("Overall VE (subcohort): ", round(100*ve_subcohort,1), "%"))



##################.
##### Graphs #####
##################.

if (cfg2$run_graphs) {
  
  # CVE graph
  {
    
    # Extract results
    theta_ests_gren <- ests$est
    ci_lo <- ests$ci_lo
    ci_hi <- ests$ci_hi
    theta_ests_gcomp <- ests$gcomp(p_grid)
    cve <- Vectorize(function(x) { 1 - x/rate_ct })
    
    # Generate CIs
    if (params$ci_type=="regular") {
      ci_lo <- cve(ci_hi)
      ci_hi <- cve(ci_lo)
    } else if (params$ci_type=="trunc") {
      ci_lo <- cve(ci_hi) %>% pmax(0) %>% pmin(1) # !!!!!
      ci_hi <- cve(ci_lo) %>% pmax(0) %>% pmin(1) # !!!!!
    }
    ests_gcomp <- cve(theta_ests_gcomp)
    ests_gren <- cve(theta_ests_gren)
    
    # Marginal distribution of A
    df_marg <- data.frame(
      x = p_grid*a_scale-a_shift,
      ymin = 0,
      ymax = (0.6/max(ests$f_a_n(p_grid))) * ests$f_a_n(p_grid)
    )
    
    # Truncate at 5/95 or 10/90 quantiles
    which <- p_grid>=ests$Phi_n_inv(0.05) & p_grid<=ests$Phi_n_inv(0.95)
    # if (cfg2$marker_num==3) {
    #   which <- p_grid<=Phi_n_inv(0.9)
    # } else {
    #   which <- p_grid>=Phi_n_inv(0.1) & p_grid<=Phi_n_inv(0.9)
    # }
    ests_gren <- ifelse(which,ests_gren,NA)
    ci_lo <- ifelse(which,ci_lo,NA)
    ci_hi <- ifelse(which,ci_hi,NA)
    
    # Labels
    x_labs <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)",
                "Phagocytic Score (=s)")
    x_lab <- x_labs[cfg2$marker_num]
    y_lab <- "Controlled VE against COVID by day 195"
    
    # Plots
    plot_1 <- ggplot(
      data.frame(
        x = rep(p_grid*a_scale-a_shift,2),
        y = c(ests_gcomp,ests_gren),
        which = rep(c("G-comp", "Grenander"), each=length(p_grid)),
        ci_lo = c(ests_gcomp,ci_lo),
        ci_hi = c(ests_gcomp,ci_hi)
      ),
      aes(x=x, y=y, color=which)
    ) +
      geom_hline(yintercept=0.465, alpha=0.3, size=0.4) + # !!!!! Later calculate this manually
      geom_hline(yintercept=c(0.411,0.515), alpha=0.3, # !!!!! Later calculate this manually
                 size=0.4, linetype="dotted") +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.1,
        linetype = "dotted",
        fill = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      # !!!!!
      scale_color_manual(values=c("purple", "darkblue")) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      scale_x_continuous(label=math_format(10^.x)) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color="Estimator") +
      geom_line()
    
    # !!!!! Consolodite this with the plot above
    plot_2 <- ggplot(
      data.frame(
        x = p_grid*a_scale-a_shift,
        y = ests_gren,
        ci_lo = ci_lo,
        ci_hi = ci_hi
      ),
      aes(x=x, y=y)
    ) +
      geom_hline(yintercept=0.465, alpha=0.3, size=0.4) + # !!!!! Later calculate this manually
      geom_hline(yintercept=c(0.411,0.515), alpha=0.3, # !!!!! Later calculate this manually
                 size=0.4, linetype="dotted") +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.1,
        linetype = "dotted",
        fill = "darkblue",
        color = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      scale_x_continuous(label=math_format(10^.x)) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color="Estimator") +
      geom_line(color="darkblue")
    
    # Save plots
    name_1 <- paste0("Janssen plots/plot1_",cfg2$tid,".pdf")
    name_2 <- paste0("Janssen plots/plot2_",cfg2$tid,".pdf")
    ggsave(filename=name_1, plot=plot_1, device="pdf", width=6, height=4)
    ggsave(filename=name_2, plot=plot_2, device="pdf", width=6, height=4)
    
  }

  # Marginal distribution histogram
  if (F) {

    # > markers
    # [1] "Day210ELCZ"
    # [2] "Day210ADCPgp140C97ZAfib"
    # [3] "Day210IgG340mdw_V1V2"
    # [4] "Day210IgG340mdw_gp120_gp140_vm"
    # [5] "Day210mdw_xassay"

    which <- 5
    x_labs <- c("IgG Vx-VT-C", "ADCP Vx-C97ZA",
                "IgG3 V1V2 Breadth", "IgG3 Env Breadth",
                "Multi-Epitopes Functions")
    mrk <- paste0(x_labs[which], " (", markers[which], ")")
    mrk <- "IgG3 gp41 (Day210IgG3gp4140delta)" # !!!!!

    ggplot(
      data.frame(x=df_trt[["Day210IgG3gp4140delta"]]),
      aes(x=x)
    ) +
      geom_histogram(bins=50) +
      labs(title=mrk)

  }

  # Gamma_os_n graph
  if (F) {

    # Export: 6" x 5"
    df_marg <- data.frame(
      x = p_grid*a_scale-a_shift,
      ymin = 0,
      ymax = ests$f_a_n(p_grid)*0.001
    )
    ggplot(
      data.frame(
        x = p_grid*a_scale-a_shift,
        y = Gamma_os_n(round(p_grid, -log10(C$appx$a)))
      ),
      aes(x=x, y=y)
    ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      theme(legend.position="bottom") +
      labs(x="Day 57 pseudoneut-ID50", y="Gamma_os_n",
           color="Estimator") +
      geom_line()

  }

  # theta graph
  # Export: 6" x 5"
  if (F) {
    
    df_marg <- data.frame(
      x = p_grid*a_scale-a_shift,
      ymin = 0,
      ymax = ests$f_a_n(p_grid)*0.01
    )
    ests_gcomp <- gcomp(p_grid)
    ests_gren <- theta_n(p_grid)
    ggplot(
      data.frame(
        x = rep(p_grid*a_scale-a_shift,2),
        y = c(ests_gcomp,ests_gren),
        which = rep(c("G-comp", "Grenander"), each=length(p_grid))
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
# !!!!! Update

# Component function check: conditional distribution
# !!!!! Update

# Component function check: marginal distribution
# !!!!! Update
