# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

#################.
##### Setup #####
#################.

# Set up configuration variables
C <- list(
  t_e=66,
  appx = list(t_e=1, w1=0.1, w1b=0.1, w3=1, a=0.01) # !!!!! a=0.001, w3=0.1
)
cfg2 <- list(
  markers = c("Day29bindSpike", "Day29bindRBD", "Day29ADCP",
              "Day29pseudoneutid50"),
  datasets = c("1", "adcp", "psv"),
  run_import_data = F,
  run_dqa = F,
  run_graphs = T
  # edge_spread = F
)

# Catch TID and set config accordingly
cfg2$tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
# cfg2$tid <- 1
cfg2_df <- data.frame(
  # tid = c(1:32),
  # S_n_type = rep(c("Cox PH", "Super Learner"), 16),
  # marker_num = rep(rep(c(1:4), each=2), 4),
  # jit_L = rep(0.25, 32),
  # jit_R = rep(c(0.6,0.75,0.9,1), each=8)
  tid = c(1:48),
  S_n_type = rep(c("Cox PH", "Super Learner", "Random Forest"), 16),
  marker_num = rep(rep(c(1:4), each=2), 6),
  jit_L = rep(0.25, 48),
  jit_R = rep(c(0.6,0.75,0.9,1), each=12)
)
cfg2_df <- rbind(cfg2_df,cfg2_df) # !!!!!
cfg2_df <- arrange(cfg2_df, marker_num, jit_L, jit_R, S_n_type) # !!!!!
cfg2_df$tid <- c(1:nrow(cfg2_df)) #!!!!!
cfg2$S_n_type <- cfg2_df[cfg2$tid, "S_n_type"]
cfg2$marker_num <- cfg2_df[cfg2$tid, "marker_num"]
cfg2$jit_L <- cfg2_df[cfg2$tid, "jit_L"]
cfg2$jit_R <- cfg2_df[cfg2$tid, "jit_R"]
cfg2$d <- case_when(
  cfg2$marker_num<=2 ~ cfg2$datasets[1],
  cfg2$marker_num==3 ~ cfg2$datasets[2],
  cfg2$marker_num==4 ~ cfg2$datasets[3],
)

cfg2$folder = paste0("Janssen functions/Janssen (tid ",cfg2$tid,")")
# cfg2$edge_spread <- ifelse(cfg2$marker_num==3, TRUE, FALSE) # !!!!!

# Create directory if it doesn't exist
if (!dir.exists(cfg2$folder)) { dir.create(cfg2$folder) }



###########################.
##### Data processing #####
###########################.

if (cfg2$run_import_data) {
  
  df_raw_1 <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Blind",
                              "ed_Phase_Data/adata/janssen_pooled_real_data_pr",
                              "ocessed_with_riskscore.csv"))
  df_raw_adcp <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Bl",
                                 "inded_Phase_Data/adata/janssen_pooled_realAD",
                                 "CP_data_processed_with_riskscore.csv"))
  df_raw_psv <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Bli",
                                "nded_Phase_Data/adata/janssen_pooled_realPsV_",
                                "data_processed_with_riskscore.csv"))
  
  for (d in cfg2$datasets) {
    
    # Subset data frames
    assign(paste0("df_ph1_",d),
           filter(eval(as.name(paste0("df_raw_",d))), ph1.D29==T))
    assign(paste0("df_ct_",d),
           filter(eval(as.name(paste0("df_ph1_",d))), Trt==0))
    assign(paste0("df_tx_",d),
           filter(eval(as.name(paste0("df_ph1_",d))), Trt==1))
    
    # Create data structures for analysis
    df <- get(paste0("df_tx_",d))
    if (d=="1") {
      a_list <- list(df[[cfg2$markers[1]]], df[[cfg2$markers[2]]], NA, NA)
    } else if (d=="adcp") {
      a_list <- list(NA, NA, df[[cfg2$markers[3]]], NA)
    } else if (d=="psv") {
      a_list <- list(NA, NA, NA, df[[cfg2$markers[4]]])
    }
    dorig <- list(
      "id" = df[["Ptid"]],
      "y_star" = df[["EventTimePrimaryIncludeNotMolecConfirmedD29"]],
      "delta_star" = df[["EventIndPrimaryIncludeNotMolecConfirmedD29"]],
      "w" = data.frame(
        "w1" = as.integer(df[["Region"]]==1),
        "w2" = as.integer(df[["Region"]]==2),
        "w3" = df[["risk_score"]]
      ),
      "weights" = df[["wt.D29"]],
      "a_list" = a_list,
      "delta" = as.integer(df[["ph2.D29"]])
    )
    
    # Stabilize weights (rescale to sum to sample size)
    dorig$weights <- ifelse(dorig$delta==1, dorig$weights, 0)
    s <- sum(dorig$weights) / length(dorig$delta)
    dorig$weights <- dorig$weights / s
    
    # Write to object
    assign(paste0("dat_orig_",d), dorig)
    rm(dorig)
    
    # Save datasets
    saveRDS(get(paste0("dat_orig_",d)),
            file=paste0("Janssen data/dat_orig_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ph1_",d)),
            file=paste0("Janssen data/df_ph1_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ct_",d)),
            file=paste0("Janssen data/df_ct_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_tx_",d)),
            file=paste0("Janssen data/df_tx_",d,"_Janssen.rds"))
    
  }
  
} else {

  for (d in cfg2$datasets) {
    
    # Read datasets
    assign(paste0("dat_orig_",d),
           readRDS(paste0("Janssen data/dat_orig_",d,"_Janssen.rds")))
    assign(paste0("df_ph1_",d),
           readRDS(paste0("Janssen data/df_ph1_",d,"_Janssen.rds")))
    assign(paste0("df_ct_",d),
           readRDS(paste0("Janssen data/df_ct_",d,"_Janssen.rds")))
    assign(paste0("df_tx_",d),
           readRDS(paste0("Janssen data/df_tx_",d,"_Janssen.rds")))
    
  }
  
}



###############################################.
##### Summary stats / data quality checks #####
###############################################.

if (cfg2$run_dqa) {
  
  # Alias vectors
  ind_tx <- df_tx_1$EventIndPrimaryIncludeNotMolecConfirmedD29
  ind_ct <- df_ct_1$EventIndPrimaryIncludeNotMolecConfirmedD29
  time_tx <- df_tx_1$EventTimePrimaryIncludeNotMolecConfirmedD29
  time_ct <- df_ct_1$EventTimePrimaryIncludeNotMolecConfirmedD29
  
  # Number of cases in each group
  num_case_tx <- sum(ind_tx)
  num_case_ct <- sum(ind_ct)
  num_case_tx_66 <- sum(ind_tx[time_tx<=66])
  num_case_ct_66 <- sum(ind_ct[time_ct<=66])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste("Number of cases in vaccine group:", num_case_tx))
  print(paste("Number of cases in control group:", num_case_ct))
  print(paste("Number of cases by day 66 in vaccine group:", num_case_tx_66))
  print(paste("Number of cases by day 66 in control group:", num_case_ct_66))
  print(paste("Number at-risk in vaccine group:", num_atrisk_tx))
  print(paste("Number at-risk in control group:", num_atrisk_ct))
  print(paste("Naive P(COVID by day 66) in vaccine group:",
              round(num_case_tx_66/num_atrisk_tx,3)))
  print(paste("Naive P(COVID by day 66) in control group:",
              round(num_case_ct_66/num_atrisk_ct,3)))
  print(paste("Naive vaccine efficacy:",
              round(1 - (num_case_tx_66/num_atrisk_tx) /
                      (num_case_ct_66/num_atrisk_ct),3)))
  
  # !!!!!
  get.marginalized.risk.no.marker <- function(dat) {
    fit.risk <- coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score +
        as.factor(Region),
      dat,
      model = T
    )
    dat[["EventTimePrimaryIncludeNotMolecConfirmedD29"]] <- 66
    risks = 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
    mean(risks)
  }
  res.plac.cont <- get.marginalized.risk.no.marker(df_ct_1)
  res.vacc.cont <- get.marginalized.risk.no.marker(df_tx_1)
  overall.ve <- 1 - res.vacc.cont/res.plac.cont
  print(overall.ve)
  
  # srv_ct <- survfit(
  #   Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
  #        EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  #   data = df_ct_1
  # )
  # rate66_ct <- 1-srv_ct$surv[which(srv_ct$time==66)]
  # srv_tx <- survfit(
  #   Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
  #        EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
  #   data = df_tx_1
  # )
  # rate66_tx <- 1-srv_tx$surv[which(srv_tx$time==66)]
  # print(rate66_ct)
  # print(rate66_tx)
  # print(1-(rate66_tx/rate66_ct))
  
  # Kaplan-Meier plot
  # !!!!!
  
  # Distribution of event times (Ph2=0 vs. Ph2=1)
  ggplot(
    data.frame(
      x = time_tx[which(ind_tx==1)],
      ph2.D29 = df_tx_1$ph2.D29[which(ind_tx==1)]
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
    data = filter(df_tx, ph2.D29==1),
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
    data = df_ct
  )
  
  survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = filter(df_tx, ph2.D29==1),
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
  dat_orig <- get(paste0("dat_orig_",cfg2$d))
  df_ct <- get(paste0("df_ct_",cfg2$d))
  df_tx <- get(paste0("df_tx_",cfg2$d))
  
  # Set `a` value from `a_list`
  dat_orig$a <- dat_orig$a_list[[cfg2$marker_num]]
  dat_orig$a_list <- NULL
  
  # Archive unjittered marker (for histogram)
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  
  # Jitter A values (left endpoint)
  dat_orig$a <- 10^dat_orig$a
  lod <- 2*min(dat_orig$a, na.rm=T)
  jit <- runif(n=length(dat_orig$a), min=cfg2$jit_L*lod, max=cfg2$jit_R*lod)
  dat_orig$a <- ifelse(dat_orig$a==min(dat_orig$a, na.rm=T), jit, dat_orig$a)
  dat_orig$a <- log10(dat_orig$a)
  
  # Rescale A to lie in [0,1]
  a2 <- 1/C$appx$a
  a_shift <- -1 * floor(a2*min(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a + a_shift
  a_scale <- ceiling(a2*max(dat_orig$a, na.rm=T))/a2
  dat_orig$a <- dat_orig$a / a_scale

  # Round A, W3
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$w$w3 <- round(dat_orig$w$w3, -log10(C$appx$w3))
  
  # Set estimation tuning parameters
  params <- list(S_n_type=cfg2$S_n_type, g_n_type="binning",
                 ecdf_type="linear (mid)", deriv_type="m-spline",
                 gamma_type="kernel", ci_type="regular", edge_corr="none",
                 omega_n_type="estimated", cf_folds=1, n_bins=0)
  
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
  return_extra = c("gcomp", "Phi_n_inv")
)

# Calculate control/vaccine group marginalized survival
get.marginalized.risk.no.marker <- function(dat, tfinal.tpeak, wts=F) {
  if (!wts) {
    fit.risk <- coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score +
        as.factor(Region),
      dat,
      model = T
    )
  } else {
    fit.risk <- coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29) ~ risk_score +
        as.factor(Region),
      dat,
      model = T,
      weights = wt.D29
    )
  }
  dat[["EventTimePrimaryIncludeNotMolecConfirmedD29"]] <- tfinal.tpeak
  risks <- 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
  mean(risks)
}
rate_ct <- get.marginalized.risk.no.marker(df_ct, C$t_e)
rate_tx <- get.marginalized.risk.no.marker(df_tx, C$t_e)
rate_tx_0 <- get.marginalized.risk.no.marker( # !!!!! Need to incorporate weights
  filter(df_tx, Day29bindSpike==min(df_tx$Day29bindSpike, na.rm=T)),
  C$t_e, wts=T) # !!!!!
rate_tx_1 <- get.marginalized.risk.no.marker( # !!!!! Need to incorporate weights
  filter(df_tx, Day29bindSpike!=min(df_tx$Day29bindSpike, na.rm=T)),
  C$t_e, wts=T) # !!!!!
round(1-(rate_tx/rate_ct),3) # Overall VE
# round(1-(rate_tx_0/rate_ct),3)
round(1-(rate_tx_1/rate_ct),3)

# !!!!! Move/reorganize everything in this section
if (F) {
  
  # Calculate control group survival (OLD)
  srv_ct <- survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = df_ct
  )
  rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_e))]
  # ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-C$t_e))]
  # ci_hi_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-C$t_e))]
  # var_ct <- ((ci_hi_ct-ci_lo_ct)/3.92)^2
  
  # Calculate treatment group survival
  srv_tx <- survfit(
    Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
         EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
    data = df_tx
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
    data = filter(df_tx, ph2.D29==1),
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
                "Phagocytic Score (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)")
    x_lab <- x_labs[cfg2$marker_num]
    y_lab <- paste("Controlled VE against COVID by day", C$t_e)
    
    # VE overall est/ci
    ve_overall <- 0.650 # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    ve_ci <- c(0.451,0.772) # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    
    # Plots: setup
    n_bins <- 20
    xlim <- c(min(a_orig)-(max(a_orig)-min(a_orig))/n_bins,
              max(a_orig)+(max(a_orig)-min(a_orig))/n_bins)
    
    # Plots: data
    plot_data_1 <- data.frame(
      x = c(a_grid,x1,x2), # x = c(a_grid,xlim),
      y = c(ests_gren,rep(ve_overall,2)),
      which = c(rep("Controlled VE", length(p_grid)),rep("Overall VE",2)),
      ci_lo = c(ci_lo,rep(ve_ci[1],2)),
      ci_hi = c(ci_hi,rep(ve_ci[2],2))
    )
    plot_data_2 <- data.frame(
      x = c(rep(a_grid,2),x1,x2), # x = c(rep(a_grid,2),xlim),
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
      scale_x_continuous(label=math_format(10^.x),limits=xlim) +
      scale_color_manual(values=c("darkblue","darkgrey")) +
      scale_fill_manual(values=c("darkblue","darkgrey")) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color=NULL, fill=NULL) +
      geom_line()
    # plot_1 # !!!!!
    plot_2 <- plot_1 %+% plot_data_2
    suppressMessages({
      plot_2 <- plot_2 +
        scale_color_manual(values=c("darkblue","purple","darkgrey")) +
        scale_fill_manual(values=c("darkblue","purple","darkgrey"))
    })
    # plot_2 # !!!!!
    
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
    mrk <- paste0(x_labs[which], " (", cfg2$markers[which], ")")
    mrk <- "IgG3 gp41 (Day210IgG3gp4140delta)" # !!!!!

    ggplot(
      data.frame(x=df_tx[["Day210IgG3gp4140delta"]]),
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
