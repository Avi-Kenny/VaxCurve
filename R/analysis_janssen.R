# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

# !!!!! Standardize trt vs. tx, ctrl vs. ct, etc.

#################.
##### Setup #####
#################.

# Setup
markers <- c("Day29bindSpike", "Day29bindRBD", "Day29ADCP")
cfg2 <- list(
  run_import_data = F,
  run_dqa = F,
  run_graphs = T
  # edge_spread = F
)

# Catch TID and set config accordingly
# !!!!! Temporary
cfg2$tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# cfg2$tid <- 1
cfg2_df <- data.frame(
  tid = c(1:24),
  S_n_type = rep(c("Cox PH", "Super Learner"), 12),
  marker_num = rep(rep(c(1,2,3), each=2), 4),
  jit_L = rep(0.25, 24),
  jit_R = rep(c(0.6,0.75,0.9,1), each=6)
)
cfg2$S_n_type <- cfg2_df[cfg2$tid, "S_n_type"]
cfg2$marker_num <- cfg2_df[cfg2$tid, "marker_num"]
cfg2$jit_L <- cfg2_df[cfg2$tid, "jit_L"]
cfg2$jit_R <- cfg2_df[cfg2$tid, "jit_R"]

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
  
  # Distribution of event times (Ph2=0 vs. Ph2=1)
  ggplot(
    data.frame(
      x = time_tx[which(ind_tx==1)],
      ph2.D29 = df_trt_1$ph2.D29[which(ind_tx==1)]
    ),
    aes(x=x, fill=ph2.D29)
  ) +
    facet_wrap(~ph2.D29) +
    geom_vline(xintercept=c(138,195), linetype="dashed", color="#333333") +
    geom_histogram() +
    labs(title="Distribution of event times, by Ph2 indicator", x="Time")
  
  # Distribution of event times (Tx vs. Ct)
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
  
  # Distribution of censoring times (Tx vs. Ct)
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
  
  # Treatment group survival (entire cohort)
  survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = filter(df_trt, ph2.D29==1),
    weights = wt.D29
  ) %>% autoplot()
  
  # Treatment group survival (subcohort)
  
  
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



###################.
##### Scratch #####
###################.

if (F) {
  # Alias markers
  a <- list(
    dat_orig_1$a_list[[1]],
    dat_orig_1$a_list[[2]],
    dat_orig_adcp$a_list[[3]]
  )
  a <- lapply(a, function(x) { x[!is.na(x)] }) # Remove NA values
  a <- lapply(a, function(x) { 10^x }) # Re-express on natural scale
  
  # Check that min marker values equal one-half the positivity cutoff or LOD
  lod <- c(
    2*min(a[[1]], na.rm=T), # 10.84237 = PosCutoff/2
    2*min(a[[2]], na.rm=T), # 14.08585 = PosCutoff/2
    2*min(a[[3]], na.rm=T)  # 11.57 = LOD/2
  )
  
  # Check marker quantiles
  quantile(a[[1]], probs=seq(0,1,0.1))
  quantile(a[[2]], probs=seq(0,1,0.1))
  quantile(a[[3]], probs=seq(0,1,0.1))
  
  # Check percent mass at left edge
  sum(a[[1]]==min(a[[1]]))/length(a[[1]])
  sum(a[[2]]==min(a[[2]]))/length(a[[2]])
  sum(a[[3]]==min(a[[3]]))/length(a[[3]])
  
  # Check percent mass at right edge
  sum(a[[1]]==max(a[[1]]))/length(a[[1]])
  sum(a[[2]]==max(a[[2]]))/length(a[[2]])
  sum(a[[3]]==max(a[[3]]))/length(a[[3]])
  
  # Explore jittering values
  sort(unique(a[[1]]))[1:5]
  sort(unique(a[[2]]))[1:5]
  sort(unique(a[[3]]))[1:5]
  lens <- c(length(a[[1]]), length(a[[2]]), length(a[[3]]))
  a_jit <- list()
  # jit_1 <- runif(n=lens[1], min=0.25*lod[1], max=0.75*lod[1])
  # jit_2 <- runif(n=lens[2], min=0.25*lod[2], max=0.75*lod[2])
  # jit_3 <- runif(n=lens[3], min=0.25*lod[3], max=0.75*lod[3])
  jit_1 <- runif(n=lens[1], min=0.4*lod[1], max=0.6*lod[1])
  jit_2 <- runif(n=lens[2], min=0.4*lod[2], max=0.6*lod[2])
  jit_3 <- runif(n=lens[3], min=0.4*lod[3], max=0.6*lod[3])
  a_jit[[1]] <- ifelse(a[[1]]==min(a[[1]]), jit_1, a[[1]])
  a_jit[[2]] <- ifelse(a[[2]]==min(a[[2]]), jit_1, a[[2]])
  a_jit[[3]] <- ifelse(a[[3]]==min(a[[3]]), jit_1, a[[3]])
  
  # Histogram (unjittered vs. jittered)
  ggplot(
    data.frame(
      x = log10(c(a[[1]],a[[2]],a[[3]],a_jit[[1]],a_jit[[2]],a_jit[[3]])),
      # x = c(a[[1]],a[[2]],a[[3]],a_jit[[1]],a_jit[[2]],a_jit[[3]]),
      which = rep(c("original","jittered"), each=sum(lens)),
      marker = rep(c(rep(1,lens[1]), rep(2,lens[2]), rep(3,lens[3])),2)
    ),
    aes(x=x)
  ) +
    geom_histogram(bins=50) +
    facet_grid(rows=dplyr::vars(which), cols=dplyr::vars(marker))
  
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
  C <- list(appx=cfg$appx, t_e=138)
  
  # !!!!! Reset appx for `t_e` and `a`
  C$appx$t_e <- 10
  C$appx$a <- 0.01 # !!!!!
  
  # Set `a` value from `a_list`
  dat_orig$a <- dat_orig$a_list[[cfg2$marker_num]]
  dat_orig$a_list <- NULL
  
  # Archive unjittered marker (for histogram)
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  
  # Jitter A values (left endpoint)
  # !!!!! Experiment with different width values
  dat_orig$a <- 10^dat_orig$a
  lod <- 2*min(dat_orig$a, na.rm=T)
  jit <- runif(n=length(dat_orig$a), min=cfg2$jit_L*lod, max=cfg2$jit_R*lod)
  dat_orig$a <- ifelse(dat_orig$a==min(dat_orig$a, na.rm=T), jit, dat_orig$a)
  dat_orig$a <- log10(dat_orig$a)
  
  # Jitter A values (right endpoint)
  # !!!!! TO DO !!!!!
  
  # Rescale A to lie in [0,1]
  a2 <- 1/C$appx$a
  a_shift <- -1 * floor(a2*min(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a + a_shift
  a_scale <- ceiling(a2*max(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a / a_scale

  # Round A, W3
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$w$w3 <- round(dat_orig$w$w3, 0) # !!!!! change to round(...,1)
  
  # Set estimation tuning parameters
  params <- list(S_n_type=cfg2$S_n_type, g_n_type="binning",
                 ecdf_type="linear (mid)", deriv_type="m-spline",
                 gamma_type="kernel", ci_type="regular", edge_corr="none",
                 cf_folds=1, n_bins=0)
  
}


#########################.
##### Data analysis #####
#########################.

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

# Calculate control/vaccine group marginalized survival
get.marginalized.risk.no.marker <- function(dat, tfinal.tpeak) {
  fit.risk <- coxph(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score +
      as.factor(Region),
    dat,
    model = T
  )
  dat[["EventTimePrimaryIncludeNotMolecConfirmedD29"]] <- tfinal.tpeak
  risks <- 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
  mean(risks)
}
rate_ct <- get.marginalized.risk.no.marker(df_ctrl, C$t_e)
rate_tx <- get.marginalized.risk.no.marker(df_trt, C$t_e)
round(1-(rate_tx/rate_ct),3)

# !!!!! Move/reorganize everything in this section
if (F) {
  
  # Calculate control group survival (OLD)
  srv_ct <- survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = df_ctrl
  )
  rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_e))]
  # ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-C$t_e))]
  # ci_hi_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-C$t_e))]
  # var_ct <- ((ci_hi_ct-ci_lo_ct)/3.92)^2
  
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
  
}



##################.
##### Graphs #####
##################.

if (cfg2$run_graphs) {
  
  # CVE graph
  {
    
    # ests <- readRDS("ests.rds")
    # ests_gcomp <- readRDS("ests_gcomp.rds")
    
    # Extract results
    theta_ests_gren <- ests$est
    ci_lo_gren <- ests$ci_lo
    ci_hi_gren <- ests$ci_hi
    theta_ests_gcomp <- ests$gcomp(p_grid)
    cve <- Vectorize(function(x) { 1 - x/rate_ct })
    
    # Generate CIs
    ci_lo <- cve(ci_hi_gren) %>% pmax(0) %>% pmin(1)
    ci_hi <- cve(ci_lo_gren) %>% pmax(0) %>% pmin(1)
    ests_gcomp <- cve(theta_ests_gcomp)
    ests_gren <- cve(theta_ests_gren)
    
    # Truncate at 5/95 quantiles and at histogram edges
    # !!!!! replace NA values with first non-NA value for point mass at left
    which1 <- p_grid>=ests$Phi_n_inv(0.05) & p_grid<=ests$Phi_n_inv(0.95)
    a_grid <- p_grid*a_scale-a_shift
    which2 <- a_grid>=min(a_orig) & a_grid<=max(a_orig)
    which <- which1 & which2
    ests_gren <- ifelse(which,ests_gren,NA)
    ests_gcomp <- ifelse(which,ests_gcomp,NA)
    ci_lo <- ifelse(which,ci_lo,NA)
    ci_hi <- ifelse(which,ci_hi,NA)
    x1 <- a_grid[min(which(which))]
    x2 <- a_grid[max(which(which))]
    
    # Labels
    x_labs <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)",
                "Phagocytic Score (=s)")
    x_lab <- x_labs[cfg2$marker_num]
    y_lab <- "Controlled VE against COVID by day 195"
    
    # VE overall est/ci
    ve_overall <- 0.501 # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    ve_ci <- c(0.437,0.564) # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    
    # Plots: setup
    n_bins <- 20
    xlim <- c(min(a_orig)-(max(a_orig)-min(a_orig))/n_bins,
              max(a_orig)+(max(a_orig)-min(a_orig))/n_bins)
    
    # Plots: data
    plot_data_1 <- data.frame(
      x = c(a_grid,x1,x2),
      # x = c(a_grid,xlim),
      y = c(ests_gren,rep(ve_overall,2)),
      which = c(rep("Controlled VE", length(p_grid)),rep("Overall VE",2)),
      ci_lo = c(ci_lo,rep(ve_ci[1],2)),
      ci_hi = c(ci_hi,rep(ve_ci[2],2))
    )
    plot_data_2 <- data.frame(
      x = c(rep(a_grid,2),x1,x2),
      # x = c(rep(a_grid,2),xlim),
      y = c(ests_gcomp,ests_gren,rep(ve_overall,2)),
      which = c(rep(c("Cox model","Controlled VE"), each=length(p_grid)),
                rep("Overall VE",2)),
      ci_lo = c(ests_gcomp,ci_lo,rep(ve_ci[1],2)),
      ci_hi = c(ests_gcomp,ci_hi,rep(ve_ci[2],2))
    )
    
    plot_1 <- ggplot(plot_data_1, aes(x=x, y=y, color=which)) +
      geom_ribbon(aes(ymin=ci_lo,ymax=ci_hi,fill=which), alpha=0.1,
                  linetype="dotted") +
      geom_histogram(mapping=aes(x=x,y=(0.6*..count..)/max(..count..)),
                     data=data.frame(x=a_orig), bins=n_bins, fill="forestgreen",
                     alpha=0.3, inherit.aes=F) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      # scale_x_continuous(label=math_format(10^.x)) +
      scale_x_continuous(label=math_format(10^.x),limits=xlim) +
      scale_color_manual(values=c("darkblue","darkgrey")) +
      scale_fill_manual(values=c("darkblue","darkgrey")) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color=NULL, fill=NULL) +
      geom_line()
    plot_1 # !!!!!
    plot_2 <- plot_1 %+% plot_data_2
    suppressMessages({
      plot_2 <- plot_2 +
        scale_color_manual(values=c("darkblue","purple","darkgrey")) +
        scale_fill_manual(values=c("darkblue","purple","darkgrey"))
    })
    plot_2 # !!!!!
    
    # Save plots
    name_1 <- paste0("Janssen plots/plot_",cfg2$tid,".pdf")
    name_2 <- paste0("Janssen plots/plot_w_Cox_",cfg2$tid,".pdf")
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
      x = a_grid,
      ymin = 0,
      ymax = ests$f_a_n(p_grid)*0.001
    )
    ggplot(
      data.frame(
        x = a_grid,
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
      x = a_grid,
      ymin = 0,
      ymax = ests$f_a_n(p_grid)*0.01
    )
    ests_gcomp <- gcomp(p_grid)
    ests_gren <- theta_n(p_grid)
    ggplot(
      data.frame(
        x = rep(a_grid,2),
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
