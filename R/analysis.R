# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

#################.
##### Setup #####
#################.

{
  # Set seed
  set.seed(1)
  
  # Set configuration variables
  # Note: NA param values are set below based on the task ID
  # !!!!! In the correlates repo, if t_e=0, it is inferred from the data
  cfg2 <- list(
    analysis = "Janssen", # Janssen Moderna
    run_dqa = F,
    run_debug = list(gren_var=F, objs=F),
    run_graphs = T,
    save_plot_data = T
  )
  
  # Set up analysis-specific configuration variables. Each row in the cfg2$map
  #     dataframe represents the set of indices to use for a particular
  #     analysis.
  if (cfg2$analysis=="Janssen") {
    
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29ADCP",
                     "Day29pseudoneutid50")
    cfg2$x_lab <- c("Anti Spike IgG (BAU/ml) (=s)",
                    "Anti RBD IgG (BAU/ml) (=s)",
                    "Phagocytic Score (=s)",
                    "Pseudovirus-nAb ID50 (IU50/ml) (=s)")
    cfg2$title <- c("Binding Antibody to Spike",
                    "Binding Antibody to RBD",
                    "Phagocytic Score",
                    "PsV Neutralization 50% Titer")
    cfg2$day <- c(29)
    cfg2$t_e <- c(66)
    cfg2$dataset <- c(
      "janssen_pooled_realbAb_data_processed_with_riskscore.csv",
      "janssen_pooled_realADCP_data_processed_with_riskscore.csv",
      "janssen_pooled_realPsV_data_processed_with_riskscore.csv"
    )
    cfg2$edge_corr <- c("min")
    cfg2$v <- list(
      id = c("Ptid"),
      time = c("EventTimePrimaryIncludeNotMolecConfirmedD29"),
      event = c("EventIndPrimaryIncludeNotMolecConfirmedD29"),
      wt = c("wt.D29start1"),
      ph1 = c("ph1.D29start1"),
      ph2 = c("ph2.D29start1"),
      covariates = c("~. + risk_score + as.factor(Region)")
    )
    cfg2$qnt_lo <- c(0)
    cfg2$qnt_hi <- c(0.9)
    cfg2$ve_overall <- c(0.650, 0.451, 0.772) # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    cfg2$folder_local <- "Janssen data/"
    cfg2$folder_cluster <- paste0("Z:/covpn/p3003/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/")
    cfg2$params = list(
      S_n_type=NA, g_n_type="binning", ecdf_type="linear (mid)",
      deriv_type="m-spline", gamma_type="kernel", ci_type="regular", # !!!!! gamma_type="Super Learner"
      edge_corr=NA, omega_n_type="estimated", cf_folds=1, n_bins=0,
      marg=NA, lod_shift="none" # 3/4
    )
    C <- list(appx=list(t_e=1,w_tol=25,y_star=1,a=0.01)) # !!!!! a=0.001
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4),
      x_lab = c(1,2,3,4),
      title = c(1,2,3,4),
      day = c(1,1,1,1),
      t_e = c(1,1,1,1),
      dataset = c(1,1,2,3),
      edge_corr = c(1,1,1,1),
      v_id = c(1,1,1,1),
      v_time = c(1,1,1,1),
      v_event = c(1,1,1,1),
      v_wt = c(1,1,1,1),
      v_ph1 = c(1,1,1,1),
      v_ph2 = c(1,1,1,1),
      v_covariates = c(1,1,1,1)
    )
    
    # Secondary map for variations within a graph; map_row corresponds to which
    #     row of cfg2$map to use
    cfg2$map2 <- data.frame(
      tid = c(1:4),
      map_row = c(1:4),
      S_n_type = rep("Super Learner",4),
      marg = rep("Gamma_star",4),
      trim = rep(F,4)
    )
    
  } else if (cfg2$analysis=="Moderna") {
    
    cfg2$marker <- c("Day29bindSpike", "Day57bindSpike", "Day29bindRBD",
                     "Day57bindRBD", "Day29pseudoneutid50",
                     "Day57pseudoneutid50", "Day29pseudoneutid80",
                     "Day57pseudoneutid80", "Day29liveneutmn50",
                     "Day57liveneutmn50")
    cfg2$x_lab <- c("Anti Spike IgG (BAU/ml) (=s)",
                    "Anti RBD IgG (BAU/ml) (=s)",
                    "Pseudovirus-nAb ID50 (IU50/ml) (=s)",
                    "Pseudovirus-nAb ID80 (IU80/ml) (=s)",
                    "Live Virus-mnAb ID50 (IU50/ml) (=s)")
    cfg2$title <- c("Binding Antibody to Spike",
                    "Binding Antibody to RBD",
                    "PsV Neutralization 50% Titer",
                    "PsV Neutralization 80% Titer",
                    "Live Virus Micro Neut 50% Titer")
    cfg2$day <- c(29,57)
    cfg2$t_e <- c(126,100)
    # cfg2$dataset <- c(paste0("P3001ModernaCOVEimmunemarkerdata_correlates_proc",
    #                          "essed_v1.0_Oct28_2021.csv"))
    cfg2$dataset <- c(paste0("P3001ModernaCOVEimmunemarkerdata_correlates_proc",
                             "essed_v1.1_lvmn_added_Jan14_2022.csv"))
    cfg2$edge_corr <- c("none", "min")
    cfg2$v <- list(
      id = c("Ptid"),
      time = c("EventTimePrimaryD29", "EventTimePrimaryD57"),
      event = c("EventIndPrimaryD29", "EventIndPrimaryD57"),
      wt = c("wt.D29", "wt.D57"),
      ph1 = c("ph1.D29", "ph1.D57"),
      ph2 = c("ph2.D29", "ph2.D57"),
      covariates = c("~. + MinorityInd + HighRiskInd + risk_score")
    )
    cfg2$qnt_lo <- c(0.05)
    cfg2$qnt_hi <- c(0.95)
    cfg2$ve_overall <- c(0.923, 0.899, 0.945) # !!!!! Later calculate the overall VE numbers manually or pull from Youyi report
    cfg2$folder_local <- "Moderna data/"
    cfg2$folder_cluster <- paste0("Z:/covpn/p3001/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/")
    cfg2$params = list(
      S_n_type=NA, g_n_type="binning", ecdf_type="linear (mid)",
      deriv_type="m-spline", gamma_type="kernel", ci_type="regular", # !!!!! gamma_type="Super Learner"
      edge_corr=NA, omega_n_type="estimated", cf_folds=1, n_bins=0,
      marg=NA, lod_shift="none"
    )
    C <- list(appx=list(t_e=1,w_tol=25,y_star=1,a=0.01)) # !!!!! a=0.001
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4,5,6,7,8,9,10),
      x_lab = c(1,1,2,2,3,3,4,4,5,5),
      title = c(1,1,2,2,3,3,4,4,5,5),
      day = c(1,2,1,2,1,2,1,2,1,2),
      t_e = c(1,2,1,2,1,2,1,2,1,2),
      dataset = c(1,1,1,1,1,1,1,1,1,1),
      edge_corr = c(1,1,1,1,2,1,2,1,2,1),
      v_id = c(1,1,1,1,1,1,1,1,1,1),
      v_time = c(1,2,1,2,1,2,1,2,1,2),
      v_event = c(1,2,1,2,1,2,1,2,1,2),
      v_wt = c(1,2,1,2,1,2,1,2,1,2),
      v_ph1 = c(1,2,1,2,1,2,1,2,1,2),
      v_ph2 = c(1,2,1,2,1,2,1,2,1,2),
      v_covariates = c(1,1,1,1,1,1,1,1,1,1)
    )
    
    # Secondary map for variations within a graph; map_row corresponds to which
    #     row of cfg2$map to use
    cfg2$map2 <- data.frame(
      tid = c(1:10),
      map_row = c(1:10),
      S_n_type = rep("Super Learner",10),
      marg = rep("Gamma_star", 10),
      trim = rep(F,10)
    )
    
  } else if (cfg2$analysis=="HVTN 705") {
    
    # !!!!! TO DO !!!!!
    
  }
  
  # Set config based on local vs. cluster
  if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
    cfg2$tid <- 1
    cfg2$dataset <- paste0(cfg2$folder_cluster,cfg2$dataset)
  } else {
    cfg2$tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
    cfg2$dataset <- paste0(cfg2$folder_local,cfg2$dataset)
  }
  
  # Set config based on cfg2$map and cfg2$map2
  i <- cfg2$map2[cfg2$tid,"map_row"]
  for (x in c("marker", "x_lab", "title", "day", "dataset", "t_e")) {
    cfg2[[x]] <- cfg2[[x]][cfg2$map[i,x]]
  }
  for (x in c("id", "time", "event", "wt", "ph1", "ph2", "covariates")) {
    cfg2$v[[x]] <- cfg2$v[[x]][cfg2$map[i,paste0("v_",x)]]
  }
  cfg2$params$edge_corr <- cfg2$edge_corr[cfg2$map[i,"edge_corr"]]
  cfg2$v$covariates <- formula(cfg2$v$covariates)
  cfg2$params$S_n_type <- cfg2$map2[cfg2$tid,"S_n_type"]
  cfg2$params$marg <- cfg2$map2[cfg2$tid,"marg"]
  cfg2$params$trim <- cfg2$map2[cfg2$tid,"trim"] # !!!!!
  C$t_e <- cfg2$t_e
  if ((i %in% c(5,7,9)) && cfg2$analysis=="Moderna") { cfg2$qnt_lo <- 0 } # !!!!! hack
  
}



###########################.
##### Data processing #####
###########################.

{
  
  df_raw <- read.csv(cfg2$dataset)
  
  # Subset data frames
  df_ph1 <- dplyr::filter(df_raw, !!rlang::sym(cfg2$v$ph1)==T)
  df_ct <- dplyr::filter(df_ph1, Trt==0)
  df_tx <- dplyr::filter(df_ph1, Trt==1)
  rm(df_raw)
  
  # Parse covariate data frame
  f <- (function(f) {
    f <- deparse(rlang::f_rhs(f))
    f <- gsub(" ","",f)
    f <- strsplit(f, "+", fixed=T)[[1]]
    f <- f[f!="."]
    factors <- c()
    vars <- c()
    for (i in 1:length(f)) {
      c1 <- (substr(f[i],1,10)=="as.factor(")
      c2 <- (substr(f[i],1,7)=="factor(")
      factors[i] <- as.integer(c1 || c2)
      if (c1) {
        vars[i] <- substr(f[i],11,nchar(f[i])-1)
      } else if (c2) {
        vars[i] <- substr(f[i],8,nchar(f[i])-1)
      } else {
        vars[i] <- f[i]
      }
    }
    return(list("vars"=vars, "factors"=factors))
  })(cfg2$v$covariates)
  df_w <- data.frame(x=c(1:length(df_tx[[f$vars[1]]])))
  col <- 1
  for (i in c(1:length(f$vars))) {
    if (f$factors[i]==0) {
      df_w[[paste0("w",col)]] <- df_tx[[f$vars[i]]]
      col <- col + 1
    } else {
      w_col <- as.factor(df_tx[[f$vars[i]]])
      levs <- unique(w_col)
      if (length(levs)==1) {
        stop(paste("Covariate", f$vars[i], "only has one unique level"))
      } else {
        for (j in c(1:(length(levs)-1))) {
          df_w[[paste0("w",col)]] <- as.integer(df_tx[[f$vars[i]]]==levs[j])
          col <- col + 1
        }
      }
    }
  }
  df_w$x <- NULL
  
  # Create data structures for analysis
  dat_orig <- list(
    "id" = df_tx[[cfg2$v$id]],
    "y_star" = df_tx[[cfg2$v$time]],
    "delta_star" = df_tx[[cfg2$v$event]],
    "w" = df_w,
    "weights" = df_tx[[cfg2$v$wt]],
    "a" = df_tx[[cfg2$marker]],
    "delta" = as.integer(df_tx[[cfg2$v$ph2]])
  )
  
  # Stabilize weights (rescale to sum to sample size)
  dat_orig$weights <- ifelse(dat_orig$delta==1, dat_orig$weights, 0)
  s <- sum(dat_orig$weights) / length(dat_orig$delta)
  dat_orig$weights <- dat_orig$weights / s
  
}



###############################################.
##### Summary stats / data quality checks #####
###############################################.

if (cfg2$run_dqa) {
  
  library(ggfortify)
  
  # Alias vectors
  ind_tx <- df_tx[[cfg2$v$event]]
  ind_ct <- df_ct[[cfg2$v$event]]
  time_tx <- df_tx[[cfg2$v$time]]
  time_ct <- df_ct[[cfg2$v$time]]
  
  # Number of cases in each group
  num_case_tx <- sum(ind_tx)
  num_case_ct <- sum(ind_ct)
  num_case_tx_t_e <- sum(ind_tx[time_tx<=C$t_e])
  num_case_ct_t_e <- sum(ind_ct[time_ct<=C$t_e])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste0("Number of cases in vaccine group: ", num_case_tx))
  print(paste0("Number of cases in control group: ", num_case_ct))
  print(paste0("Number of cases by day ", C$t_e, " in vaccine group: ",
              num_case_tx_t_e))
  print(paste0("Number of cases by day ", C$t_e, " in control group: ",
              num_case_ct_t_e))
  print(paste0("Number at-risk in vaccine group: ", num_atrisk_tx))
  print(paste0("Number at-risk in control group: ", num_atrisk_ct))
  print(paste0("Naive P(COVID by day ", C$t_e, ") in vaccine group: ",
              round(num_case_tx_t_e/num_atrisk_tx,3)))
  print(paste0("Naive P(COVID by day ", C$t_e, ") in control group: ",
              round(num_case_ct_t_e/num_atrisk_ct,3)))
  print(paste0("Naive vaccine efficacy: ",
              round(1 - (num_case_tx_t_e/num_atrisk_tx) /
                      (num_case_ct_t_e/num_atrisk_ct),3)))
  
  # Fraction of point mass at edge
  a <- dat_orig$a
  round(sum(a==min(a,na.rm=T),na.rm=T) / sum(!is.na(a)), 3)
  
  # Distribution of event times (Ph2=0 vs. Ph2=1)
  ggplot(
    data.frame(
      x = time_tx[which(ind_tx==1)],
      ph2 = df_tx[[cfg2$v$ph2]][which(ind_tx==1)]
    ),
    aes(x=x, fill=ph2)
  ) +
    facet_wrap(~ph2) +
    # geom_vline(xintercept=c(138,195), linetype="dashed", color="#333333") +
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
    formula(paste0("Surv(",cfg2$v$time,",",cfg2$v$event,")~1")),
    data = df_tx
  ) %>% autoplot()
  
  # Treatment group survival (subcohort)
  survfit(
    formula(paste0("Surv(",cfg2$v$time,",",cfg2$v$event,")~1")),
    data = df_tx,
    weights = wt.D29start1
  ) %>% autoplot()
  
}



#########################.
##### Data analysis #####
#########################.

{
  
  # if (cfg2$params$trim) {
  #   q02 <- quantile(dat_orig$a, na.rm=T, probs=0.02)[[1]]
  #   if (length(which(dat_orig$a<q02))>0) {
  #     indicies_to_keep <- c(1:length(dat_orig$a))[-(which(dat_orig$a<q02))]
  #     dat_orig <- ss(dat_orig, indicies_to_keep)
  #   }
  # }
  
  # # !!!!! temporary
  # mrks <- c("Day29pseudoneutid50", "Day29pseudoneutid50", "Day29liveneutmn50")
  # if (cfg2$marker %in% mrks && cfg2$analysis=="Moderna") {
  #   q98 <- quantile(dat_orig$a, na.rm=T, probs=0.98)[[1]]
  #   indicies_to_keep <- c(1:length(dat_orig$a))[-(which(dat_orig$a>q98))]
  #   dat_orig <- ss(dat_orig, indicies_to_keep)
  # }
  
  # Archive original marker and lod*1/2 value (for plot)
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  lod12 <- min(dat_orig$a,na.rm=T)
  
  # Obtain estimates
  p_grid <- seq(0,1,0.01)
  return_extra <- c("Phi_n_inv", "deriv_theta_n", "f_a_n", "gamma_n", "Psi_n",
                    "omega_n", "f_aIw_n", "S_n", "gcm", "dGCM")
  if (cfg2$params$marg=="Theta") {
    return_extra <- c(return_extra, "etastar_n")
  }
  a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=101)
  ests <- est_curve(
    dat_orig = dat_orig,
    estimator = "Grenander",
    params = cfg2$params,
    points = a_grid,
    dir = "decr",
    return_extra = return_extra
  )
  
  # Debugging: Grenander variance scale factor components
  if (cfg2$run_debug$gren_var) {
    print("Grenander variance scale factor components")
    print("deriv_theta_n")
    print(ests$deriv_theta_n(seq(0,1,0.05)))
    print("f_a_n")
    print(ests$f_a_n(seq(0,1,0.05)))
    print("gamma_n")
    print(ests$gamma_n(seq(0,1,0.05)))
    print("deriv_theta_n*f_a_n")
    print(ests$deriv_theta_n(seq(0,1,0.05))*ests$f_a_n(seq(0,1,0.05)))
    print("deriv_theta_n*gamma_n")
    print(ests$deriv_theta_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
    print("f_a_n*gamma_n")
    print(ests$f_a_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
    print("deriv_theta_n*f_a_n*gamma_n")
    print(ests$deriv_theta_n(seq(0,1,0.05))*ests$f_a_n(seq(0,1,0.05))*
            ests$gamma_n(seq(0,1,0.05)))
  }
  
  # Calculate control/vaccine group marginalized survival
  get.marginalized.risk.no.marker <- function(dat, tfinal.tpeak) {
    fit.risk <- coxph(
      update(as.formula(paste0("Surv(",cfg2$v$time,",",cfg2$v$event,")~1")),
             cfg2$v$covariates),
      dat,
      # weights = cfg2$v$wt,
      model = T
    )
    dat[[cfg2$v$time]] <- tfinal.tpeak
    risks <- 1 - exp(-predict(fit.risk, newdata=dat, type="expected"))
    mean(risks)
  }
  rate_ct <- get.marginalized.risk.no.marker(df_ct, C$t_e)
  rate_tx <- get.marginalized.risk.no.marker(df_tx, C$t_e)
  # print(paste0("Overall VE: ", round(1-rate_tx/rate_ct, 3)))
  # rate_tx_0 <- get.marginalized.risk.no.marker(
  #   dplyr::filter(df_tx, Day29bindSpike==min(df_tx$Day29bindSpike, na.rm=T)),
  #   C$t_e, wts=T) # !!!!!
  # rate_tx_1 <- get.marginalized.risk.no.marker(
  #   dplyr::filter(df_tx, Day29bindSpike!=min(df_tx$Day29bindSpike, na.rm=T)),
  #   C$t_e, wts=T) # !!!!!
  
}



########################################.
##### Examine intermediate objects #####
########################################.

if (cfg2$run_debug$objs) {
  
  grid <- seq(0,1,0.01)
  gcm <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots, ties="ordered")
  omega_n <- Vectorize(function(a) {
    ests$omega_n(w=c(0,0),a,y_star=100,delta_star=0)
  })
  etastar_n <- Vectorize(function(a) {
    ests$etastar_n(a,w=c(0,0))
  })
  S_n <- Vectorize(function(a) {
    ests$S_n(t=C$t_e, w=c(0,0), a)
  })
  
  plot_data <- data.frame(
    x = rep(grid,3),
    y = c(ests$Psi_n(grid), gcm(grid), ests$dGCM(grid)),
    which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*theta_n)"), each=101)
  )
  plot1 <- ggplot(plot_data, aes(x=x, y=y, color=which)) +
    geom_line() +
    theme(legend.position="bottom")
  
  ggsave(
    filename = paste0(cfg2$analysis," plots/debug_",cfg2$tid,".pdf"),
    plot=plot1, device="pdf", width=6, height=4
  )
  
  if (F) {
    
    # S_n: conditional survival function (as a function of A)
    # !!!!! Continue
    # as.data.frame(cbind(w1=rep(0.2,n), w2=rep(1,n)))
    plot_data2 <- data.frame(
      x = grid,
      y = S_n(grid)
      # which = rep(c("",""),each=101)
    )
    ggplot(plot_data2, aes(x=x, y=y)) + # color=which
      geom_line() +
      theme(legend.position="bottom")
    
    # # omega_n
    # plot_data3 <- data.frame(
    #   x = grid,
    #   y = omega_n(grid)
    #   # which = rep(c("",""),each=101)
    # )
    # ggplot(plot_data3, aes(x=x, y=y)) + # color=which
    #   geom_line() +
    #   theme(legend.position="bottom")
    
    # etastar_n
    plot_data4 <- data.frame(
      x = grid,
      y = etastar_n(grid)
      # which = rep(c("",""),each=101)
    )
    ggplot(plot_data4, aes(x=x, y=y)) + # color=which
      geom_line() +
      theme(legend.position="bottom") +
      labs(title="etastar_n")
    
  }
  
}



##################.
##### Graphs #####
##################.

if (cfg2$run_graphs) {
  
  # CVE graph
  {
    
    # Extract results
    theta_ests_gren <- ests$est
    ci_lo_gren <- ests$ci_lo
    ci_hi_gren <- ests$ci_hi
    cve <- Vectorize(function(x) { 1 - x/rate_ct })
    
    # Generate CIs
    ci_lo <- cve(ci_hi_gren) %>% pmax(0) %>% pmin(1)
    ci_hi <- cve(ci_lo_gren) %>% pmax(0) %>% pmin(1)
    ests_gren <- cve(theta_ests_gren)
    
    # Truncate at 10/90 quantiles and at histogram edges
    which1 <- p_grid>=ests$Phi_n_inv(cfg2$qnt_lo) &
      p_grid<=ests$Phi_n_inv(cfg2$qnt_hi)
    which2 <- a_grid<=max(a_orig) # ?????
    # which2 <- a_grid>=min(a_orig) & a_grid<=max(a_orig) # ?????
    which12 <- which1 & which2
    ests_gren <- ifelse(which12,ests_gren,NA)
    ci_lo <- ifelse(which12,ci_lo,NA)
    ci_hi <- ifelse(which12,ci_hi,NA)
    # x1 <- a_grid[min(which(which12))]
    x1 <- lod12
    x2 <- max(a_grid) # !!!!!
    
    # Labels
    y_lab <- paste("Controlled VE against COVID by day", C$t_e)
    
    # Plots: setup
    n_bins <- 20
    xlim <- c(min(a_orig)-(max(a_orig)-min(a_orig))/(n_bins*(2/3)),
              max(a_orig)+(max(a_orig)-min(a_orig))/(n_bins*(2/3)))
    
    # Plots: data
    plot_data <- data.frame(
      x = c(x1,a_grid,x1,x2),
      y = c(ests_gren[1],ests_gren,rep(cfg2$ve_overall[1],2)),
      curve = factor(
        c(rep("Controlled VE", length(p_grid)+1),rep("Overall VE",2)),
        levels=c("Overall VE","Controlled VE")
      ),
      ci_lo = c(ci_lo[1],ci_lo,rep(cfg2$ve_overall[2],2)),
      ci_hi = c(ci_hi[1],ci_hi,rep(cfg2$ve_overall[3],2))
    )
    if (cfg2$save_plot_data) {
      saveRDS(
        list(plot_data=plot_data, a_orig=a_orig, n_bins=n_bins, xlim=xlim,
             cfg2=cfg2, y_lab=y_lab),
        paste0(cfg2$analysis," plots/plot_data_",cfg2$tid,".rds")
      )
    }
    
    plot <- ggplot(plot_data, aes(x=x, y=y, color=curve)) +
      geom_ribbon(aes(ymin=ci_lo,ymax=ci_hi,fill=curve), alpha=0.1,
                  linetype="dotted") +
      geom_histogram(mapping=aes(x=x,y=(0.6*..count..)/max(..count..)),
                     data=data.frame(x=a_orig), bins=n_bins, fill="forestgreen",
                     alpha=0.3, inherit.aes=F) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      scale_x_continuous(label=math_format(10^.x),limits=xlim) +
      scale_color_manual(values=c("darkgrey","darkblue")) +
      scale_fill_manual(values=c("darkgrey","darkblue")) +
      theme(legend.position="bottom") +
      labs(title=paste0(cfg2$title,": Day ",cfg2$day),
           x=cfg2$x_lab, y=y_lab, color=NULL, fill=NULL) +
      geom_line()
    
    # Save plot
    ggsave(
      filename = paste0(cfg2$analysis," plots/plot_",cfg2$tid,".pdf"),
      plot=plot, device="pdf", width=6, height=4
    )
    
    # Alternate plots (for David)
    if (F) {
      
      for (i in c(1:10)) {
        
        plot_data_ <- readRDS(paste0("Moderna plots/plot_data_",i,".rds"))
        attach(plot_data_)
        
        ind_lo <- min(which(!is.na(plot_data$y))) - 3
        if (ind_lo<0) { ind_lo <- 1 }
        ind_hi <- max(head(which(!is.na(plot_data$y)),-2)) + 3
        new_xlims <- c(plot_data$x[ind_lo], plot_data$x[ind_hi])
        
        plot <- ggplot(plot_data, aes(x=x, y=y, color=curve)) +
          geom_ribbon(aes(ymin=ci_lo,ymax=ci_hi,fill=curve), alpha=0.1,
                      linetype="dotted") +
          geom_histogram(
            mapping=aes(x=x,y=0.526+0.45*((0.6*..count..)/max(..count..))),
            data=data.frame(x=a_orig), bins=n_bins, fill="forestgreen",
            alpha=0.3, inherit.aes=F) +
          coord_cartesian(xlim=new_xlims, ylim=c(0.55,1)) +
          scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1),
                             breaks=seq(-1,1,0.1), minor_breaks=NULL) +
          theme(panel.grid.major=element_line(colour="white", size=0.3),
                panel.grid.minor=element_line(colour="white", size=0.3)) +
          scale_x_continuous(label=math_format(10^.x),limits=xlim) +
          scale_color_manual(values=c("darkgrey","darkblue")) +
          scale_fill_manual(values=c("darkgrey","darkblue")) +
          theme(legend.position="bottom") +
          labs(title=paste0(cfg2$title,": Day ",cfg2$day),
               x=cfg2$x_lab, y=y_lab, color=NULL, fill=NULL) +
          geom_line()
        
        # Save plot
        ggsave(
          filename = paste0(cfg2$analysis," plots/plot_zoomed_",i,".pdf"),
          plot=plot, device="pdf", width=6, height=4
        )
        
        detach(plot_data_)
        
      }
      
    }
    
  }
  
}



###################.
##### Archive #####
###################.

if (F) {
  
  # Data processing
  {
    # Read in raw data
    df_raw_1 <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_Bli",
                                "nded_Phase_Data/adata/janssen_pooled_real_dat",
                                "a_processed_with_riskscore.csv"))
    df_raw_adcp <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_",
                                   "Blinded_Phase_Data/adata/janssen_pooled_re",
                                   "alADCP_data_processed_with_riskscore.csv"))
    df_raw_psv <- read.csv(paste0("Z:/covpn/p3003/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/janssen_pooled_real",
                                  "PsV_data_processed_with_riskscore.csv"))
    
    # Save datasets
    saveRDS(get(paste0("dat_orig_",d)),
            file=paste0("Janssen data/dat_orig_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ph1_",d)),
            file=paste0("Janssen data/df_ph1_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_ct_",d)),
            file=paste0("Janssen data/df_ct_",d,"_Janssen.rds"))
    saveRDS(get(paste0("df_tx_",d)),
            file=paste0("Janssen data/df_tx_",d,"_Janssen.rds"))
    
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
  
  # Use Kaplan-Meier to calculate survival and efficacy with CI
  {
    # Calculate control group survival (KM; with SE)
    srv_ct <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_ct
    )
    rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_e))]
    ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-C$t_e))]
    ci_hi_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-C$t_e))]
    var_ct <- ((ci_hi_ct-ci_lo_ct)/3.92)^2
    
    # Calculate treatment group survival (KM; with SE)
    srv_tx <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_tx
    )
    rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-C$t_e))]
    ci_lo_tx <- 1 - srv_tx$upper[which.min(abs(srv_tx$time-C$t_e))]
    ci_hi_tx <- 1 - srv_tx$lower[which.min(abs(srv_tx$time-C$t_e))]
    var_tx <- ((ci_hi_tx-ci_lo_tx)/3.92)^2
    
    # Calculate overall vaccine efficacy (KM; delta method; with SE+CI)
    ve_overall <- 1 - (rate_tx/rate_ct)
    ve_se <- sqrt(rate_ct^-2*var_tx + rate_tx^2*rate_ct^-4*var_ct)
    ve_overall_lo <- ve_overall - 1.96*ve_se
    ve_overall_hi <- ve_overall + 1.96*ve_se
    print(paste0("Overall VE: ",round(100*ve_overall,1),"% (",
                 round(100*ve_overall_lo,1),"% -- ",round(100*ve_overall_hi,1),
                 "%)"))
    
    # Calculate overall vaccine efficacy (KM; within subcohort)
    srv_tx_sub <- survfit(coxph(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = dplyr::filter(df_tx, ph2.D29start1==1),
      weights = wt.D29start1
    ))
    rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_tx_sub$time-C$t_e))]
    ve_subcohort <- 1 - (rate_tx_sub/rate_ct)
    print(paste0("Overall VE (subcohort): ", round(100*ve_subcohort,1), "%"))
    
    # Plot KM curve
    library(ggfortify)
    survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~Trt,
      data = df_ph1_1
    ) %>% autoplot()
    
  }
  
  # Check marker quantiles and edge mass
  {
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
  }
  
  # Extra pieces to overlay Cox model gcomp estimator in plot
  {
    theta_ests_gcomp <- ests$gcomp(p_grid)
    ests_gcomp <- cve(theta_ests_gcomp)
    ests_gcomp <- ifelse(which,ests_gcomp,NA)
    
    plot_data_2 <- data.frame(
      x = c(x1,rep(a_grid,2),x1,x2),
      y = c(ests_gren[1],ests_gren,ests_gcomp,rep(ve_overall,2)),
      which = c("Controlled VE",rep(c("Controlled VE","Cox model"),
                                    each=length(p_grid)),
                rep("Overall VE",2)),
      ci_lo = c(ci_lo[1],ci_lo,ests_gcomp,rep(ve_ci[1],2)),
      ci_hi = c(ci_hi[1],ci_hi,ests_gcomp,rep(ve_ci[2],2))
    )
    plot_2 <- plot_1 %+% plot_data_2
    suppressMessages({
      plot_2 <- plot_2 +
        scale_color_manual(values=c("darkblue","purple","darkgrey")) +
        scale_fill_manual(values=c("darkblue","purple","darkgrey"))
    })
    name_2 <- paste0(cfg2$analysis," plots/plot_w_Cox_",cfg2$tid,".pdf")
    ggsave(filename=name_2, plot=plot_2, device="pdf", width=6, height=4)
    
  }
  
}
