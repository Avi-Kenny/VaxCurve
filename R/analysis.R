# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

# !!!!! Rewrite with inner parallelization via parLapply

#################.
##### Setup #####
#################.

{
  # Choose analysis
  which_analysis <- "Janssen" # Janssen Moderna AMP
  
  # # TEMP: AMP analysis
  # ..tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  # if (..tid<=10) {
  #   which_analysis <- "Moderna"
  # } else {
  #   which_analysis <- "Janssen"
  #   Sys.setenv("SLURM_ARRAY_TASK_ID"=as.character(..tid-10))
  # }
  
  # # TEMP: Uncomment this code to run both analyses (1=10=Moderna, 11-14=Janssen)
  # ..tid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  # if (..tid<=10) {
  #   which_analysis <- "Moderna"
  # } else {
  #   which_analysis <- "Janssen"
  #   Sys.setenv("SLURM_ARRAY_TASK_ID"=as.character(..tid-10))
  # }
  
  # Set seed
  set.seed(1)
  
  # Set configuration variables
  # Note: NA param values are set below based on the task ID
  # !!!!! In the correlates repo, if t_e=0, it is inferred from the data
  cfg2 <- list(
    analysis = which_analysis,
    run_analysis = T,
    run_dqa = F,
    run_debug = list(gren_var=F, objs=F),
    run_hyptest = F
  )
  
  # Set up analysis-specific configuration variables. Each row in the cfg2$map
  #     dataframe represents the set of indices to use for a particular
  #     analysis.
  if (cfg2$analysis=="Janssen") {
    
    cfg2$plot_cve <- list(overall="Cox", est=c("Grenander", "Cox")) # "Qbins", "Cox gcomp"
    cfg2$plot_risk <- list(overall="Cox", est=c("Grenander", "Cox")) # "Qbins", "Cox gcomp"
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50",
                     "Day29ADCP")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29",
                    "Binding Antibody to RBD: Day 29",
                    "PsV Neutralization 50% Titer: Day 29",
                    "Phagocytic Score: Day 29")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)",
                    "Anti RBD IgG (BAU/ml) (=s)",
                    "Pseudovirus-nAb ID50 (IU50/ml) (=s)",
                    "Phagocytic Score (=s)")
    cfg2$endpoint <- "COVID"
    cfg2$t_e <- c(54)
    cfg2$dataset <- c(
      "janssen_pooled_real_data_processed_with_riskscore.csv",
      "janssen_pooled_realADCP_data_processed_with_riskscore.csv"
    )
    cfg2$txct <- T
    cfg2$cr2_trial <- c("janssen_pooled_real", "janssen_pooled_realADCP")
    cfg2$cr2_COR <- c("D29IncludeNotMolecConfirmedstart1")
    cfg2$cr2_marker <- c(1,2,3)
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
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0,0.95),
      "CVE, nonparametric" = c(0,0.95),
      "Risk, Qbins" = c(0,0.95),
      "CVE, Qbins" = c(0,0.95),
      "Risk, Cox model" = c(0,0.975),
      "CVE, Cox model" = c(0,0.975)
    )
    cfg2$zoom_x <- NA
    cfg2$zoom_y <- NA
    cfg2$folder_local <- "Janssen data/"
    cfg2$folder_cluster <- paste0("Z:/covpn/p3003/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/")
    cfg2$params = list(
      g_n_type="binning", ecdf_type="linear (mid)", deriv_type="linear",
      gamma_type="Super Learner", gamma_which="new", ci_type="regular",
      omega_n_type="estimated", cf_folds=1, n_bins=5, lod_shift="none"
    )
    C <- list(appx=list(t_e=1,w_tol=25,a=0.01)) # !!!!! a=0.001
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4),
      lab_x = c(1,2,3,4),
      lab_title = c(1,2,3,4),
      day = c(1,1,1,1),
      t_e = c(1,1,1,1),
      dataset = c(1,1,1,2),
      cr2_trial = c(1,1,1,2),
      cr2_COR = c(1,1,1,1),
      cr2_marker = c(1,2,3,1),
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
      marg = rep("Gamma_star2",4)
    )
    
  } else if (cfg2$analysis=="Moderna") {
    
    cfg2$plot_cve <- list(overall="Cox", est=c("Grenander","Cox")) # "Qbins",
    cfg2$plot_risk <- list(overall="Cox", est=c("Grenander","Cox")) # "Qbins",
    cfg2$marker <- c("Day29bindSpike", "Day57bindSpike", "Day29bindRBD",
                     "Day57bindRBD", "Day29pseudoneutid50",
                     "Day57pseudoneutid50", "Day29pseudoneutid80",
                     "Day57pseudoneutid80", "Day29liveneutmn50",
                     "Day57liveneutmn50")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29",
                        "Binding Antibody to Spike: Day 57",
                        "Binding Antibody to RBD: Day 29",
                        "Binding Antibody to RBD: Day 57",
                        "PsV Neutralization 50% Titer: Day 29",
                        "PsV Neutralization 50% Titer: Day 57",
                        "PsV Neutralization 80% Titer: Day 29",
                        "PsV Neutralization 80% Titer: Day 57",
                        "Live Virus Micro Neut 50% Titer: Day 29",
                        "Live Virus Micro Neut 50% Titer: Day 57")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)",
                    "Anti RBD IgG (BAU/ml) (=s)",
                    "Pseudovirus-nAb ID50 (IU50/ml) (=s)",
                    "Pseudovirus-nAb ID80 (IU80/ml) (=s)",
                    "Live Virus-mnAb ID50 (IU50/ml) (=s)")
    cfg2$endpoint <- "COVID"
    cfg2$t_e <- c(126,100)
    cfg2$dataset <- c(paste0("P3001ModernaCOVEimmunemarkerdata_correlates_proc",
                             "essed_v1.1_lvmn_added_Jan14_2022.csv"))
    cfg2$txct <- T
    cfg2$cr2_trial <- c("moderna_real")
    cfg2$cr2_COR <- c("D29", "D57")
    cfg2$cr2_marker <- c(1,2,3,4,5)
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
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.05,0.95),
      "CVE, nonparametric" = c(0.05,0.95),
      "Risk, Qbins" = c(0.05,0.95),
      "CVE, Qbins" = c(0.05,0.95),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y <- "zoomed"
    cfg2$folder_local <- "Moderna data/"
    cfg2$folder_cluster <- paste0("Z:/covpn/p3001/analysis/correlates/Part_A_B",
                                  "linded_Phase_Data/adata/")
    cfg2$params = list(
      g_n_type="binning", ecdf_type="linear (mid)", deriv_type="linear",
      gamma_type="Super Learner", gamma_which="new", ci_type="regular",
      omega_n_type="estimated", cf_folds=1, n_bins=5, lod_shift="none"
    )
    C <- list(appx=list(t_e=1,w_tol=25,a=0.01)) # !!!!! a=0.001
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4,5,6,7,8,9,10),
      lab_x = c(1,1,2,2,3,3,4,4,5,5),
      lab_title = c(1,2,3,4,5,6,7,8,9,10),
      day = c(1,2,1,2,1,2,1,2,1,2),
      t_e = c(1,2,1,2,1,2,1,2,1,2),
      dataset = c(1,1,1,1,1,1,1,1,1,1),
      cr2_trial = c(1,1,1,1,1,1,1,1,1,1),
      cr2_COR = c(1,2,1,2,1,2,1,2,1,2),
      cr2_marker = c(1,1,2,2,3,3,4,4,5,5),
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
      marg = rep("Gamma_star2", 10)
    )
    
  } else if (cfg2$analysis=="AMP") {
    
    cfg2$plot_cve <- list(overall=FALSE, est=FALSE)
    cfg2$plot_risk <- list(overall="KM", est=c("Grenander")) # "Qbins"
    cfg2$marker <- "bweight"
    cfg2$lab_title <- c("HVTN703/HPTN081","HVTN704/HPTN085","Pooled AMP trials")
    cfg2$lab_x <- "Body Weight (kg)"
    cfg2$endpoint <- "HIV-1 infection"
    cfg2$amp_protocol <- c("HVTN 703", "HVTN 704", "Pooled")
    cfg2$amp_tx <- c("C3", "T1", "T2", "T1+T2")
    cfg2$t_e <- c(595,595) # c(601,609)
    cfg2$dataset <- "amp_survival_all.csv"
    cfg2$txct <- F
    cfg2$cr2_trial <- F
    cfg2$cr2_COR <- F
    cfg2$cr2_marker <- F
    cfg2$edge_corr <- "none"
    cfg2$v <- list(
      id = "pub_id",
      time = "hiv1survday",
      event = "hiv1event",
      wt = NA,
      ph1 = NA,
      ph2 = NA,
      covariates = c("~. + age + standardized_risk_score",
                     "~. + as.factor(protocol) + age + standardized_risk_score")
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.1,0.9),
      "CVE, nonparametric" = c(0.1,0.9),
      "Risk, Qbins" = c(0.1,0.9),
      "CVE, Qbins" = c(0.1,0.9),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y <- "zoomed"
    cfg2$folder_local <- "AMP data/"
    cfg2$folder_cluster <- paste0("Z:/vaccine/p704/analysis/datashare/avi_kenn",
                                  "y/adata/")
    cfg2$params = list(
      g_n_type="binning", ecdf_type="linear (mid)", deriv_type="linear",
      gamma_type="Super Learner", gamma_which="new", ci_type="regular",
      omega_n_type="estimated", cf_folds=1, n_bins=5, lod_shift="none"
    )
    C <- list(appx=list(t_e=10,w_tol=25,a=0.01)) # !!!!! a=0.001
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      amp_protocol = c(1,1,1,1,2,2,2,2,3,3,3,3),
      amp_tx = c(1,2,3,4,1,2,3,4,1,2,3,4),
      marker = rep(1,12),
      lab_x = rep(1,12),
      lab_title = c(1,1,1,1,2,2,2,2,3,3,3,3),
      day = rep(1,12),
      t_e = c(1,1,1,1,2,2,2,2,1,1,1,1),
      dataset = rep(1,12),
      cr2_trial = rep(1,12),
      cr2_COR = rep(1,12),
      cr2_marker = rep(1,12),
      edge_corr = rep(1,12),
      v_id = rep(1,12),
      v_time = rep(1,12),
      v_event = rep(1,12),
      v_wt = rep(1,12),
      v_ph1 = rep(1,12),
      v_ph2 = rep(1,12),
      v_covariates = c(1,1,1,1,1,1,1,1,2,2,2,2)
    )
    
    # Secondary map for variations within a graph; map_row corresponds to which
    #     row of cfg2$map to use
    cfg2$map2 <- data.frame(
      tid = c(1:12),
      map_row = c(1:12),
      S_n_type = rep("Super Learner",12),
      marg = rep("Gamma_star2", 12)
    )
    
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
  for (x in c("marker", "lab_x", "lab_title", "day", "dataset", "cr2_trial",
              "cr2_COR", "cr2_marker", "t_e")) {
    cfg2[[x]] <- cfg2[[x]][cfg2$map[i,x]]
  }
  for (x in c("id", "time", "event", "wt", "ph1", "ph2", "covariates")) {
    cfg2$v[[x]] <- cfg2$v[[x]][cfg2$map[i,paste0("v_",x)]]
  }
  if (cfg2$analysis=="AMP") {
    cfg2[["amp_protocol"]] <- cfg2[["amp_protocol"]][cfg2$map[i,"amp_protocol"]]
    cfg2[["amp_tx"]] <- cfg2[["amp_tx"]][cfg2$map[i,"amp_tx"]]
  }
  cfg2$params$edge_corr <- cfg2$edge_corr[cfg2$map[i,"edge_corr"]]
  cfg2$v$covariates <- formula(cfg2$v$covariates)
  cfg2$params$S_n_type <- cfg2$map2[cfg2$tid,"S_n_type"]
  cfg2$params$marg <- cfg2$map2[cfg2$tid,"marg"]
  C$t_e <- cfg2$t_e
  if ((i %in% c(5,7,9)) && cfg2$analysis=="Moderna") {
    cfg2$qnt <- lapply(cfg2$qnt, function(x) { c(0,x[2]) }) # !!!!! temp hack
  }
  
}



###########################.
##### Data processing #####
###########################.

{
  # Read in primary data
  df_raw <- read.csv(cfg2$dataset)
  
  # Subset data frames
  if (!is.na(cfg2$v$ph1)) {
    df_ph1 <- dplyr::filter(df_raw, !!rlang::sym(cfg2$v$ph1)==T)
  } else {
    df_ph1 <- df_raw
  }
  if (cfg2$analysis=="AMP") {
    if (cfg2$amp_protocol!="Pooled") {
      df_ph1 %<>% dplyr::filter(protocol==cfg2$amp_protocol)
    }
    if (cfg2$amp_tx=="T1+T2") {
      df_ph1 %<>% dplyr::filter(tx_pool=="T1+T2")
    } else {
      df_ph1 %<>% dplyr::filter(tx==cfg2$amp_tx)
    }
    # print(paste("# infections:", sum(df_ph1$hiv1event))) # !!!!! QA
  }
  if (cfg2$txct) {
    df_ct <- dplyr::filter(df_ph1, Trt==0)
    df_tx <- dplyr::filter(df_ph1, Trt==1)
    df_analysis <- df_tx
  } else {
    df_analysis <- df_ph1
  }
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
  df_w <- data.frame(x=c(1:length(df_analysis[[f$vars[1]]])))
  col <- 1
  for (i in c(1:length(f$vars))) {
    if (f$factors[i]==0) {
      df_w[[paste0("w",col)]] <- df_analysis[[f$vars[i]]]
      col <- col + 1
    } else {
      w_col <- as.factor(df_analysis[[f$vars[i]]])
      levs <- unique(w_col)
      if (length(levs)==1) {
        stop(paste("Covariate", f$vars[i], "has only one unique level"))
      } else {
        for (j in c(1:(length(levs)-1))) {
          df_w[[paste0("w",col)]] <- as.integer(
            df_analysis[[f$vars[i]]]==levs[j]
          )
          col <- col + 1
        }
      }
    }
  }
  df_w$x <- NULL
  
  # Create data structures for analysis
  if (is.na(cfg2$v$ph2)) {
    df_delta <- rep(1, nrow(df_analysis))
  } else {
    df_delta <- as.integer(df_analysis[[cfg2$v$ph2]])
  }
  if (is.na(cfg2$v$wt)) {
    df_weights <- rep(1, nrow(df_analysis))
  } else {
    df_weights <- df_analysis[[cfg2$v$wt]]
  }
  dat_orig <- list(
    "id" = df_analysis[[cfg2$v$id]],
    "y_star" = df_analysis[[cfg2$v$time]],
    "delta_star" = df_analysis[[cfg2$v$event]],
    "w" = df_w,
    "weights" = df_weights,
    "a" = df_analysis[[cfg2$marker]],
    "delta" = df_delta
  )
  rm(df_w,df_weights,df_delta)
  
  # Create data structure to hold results
  plot_data_risk <- data.frame(
    x = double(),
    y = double(),
    curve = character(),
    ci_lo = double(),
    ci_hi = double(),
    tag = character()
  )
  plot_data_cve <- plot_data_risk
  
  # Stabilize weights (rescale to sum to sample size)
  dat_orig$weights <- ifelse(dat_orig$delta==1, dat_orig$weights, 0)
  s <- sum(dat_orig$weights) / length(dat_orig$delta)
  dat_orig$weights <- dat_orig$weights / s
  
}



############################################################.
##### Import functions from correlates_reporting2 repo #####
############################################################.

{
  
  # Calculate marginalized risk (adapted)
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
  
  # Get histogram (adapted)
  get.marker.histogram <- function(marker, wt=NA, trial) {
    # first call hist to get breaks, then call weighted.hist
    tmp.1 <- hist(marker, breaks=ifelse(trial=="moderna_real", 25, 15), plot=F)
    if (is.na(wt[1])) { wt <- rep(1, length(marker)) }
    tmp <- plotrix::weighted.hist(marker, wt, breaks=tmp.1$breaks, plot=F)
    attr(tmp,"class") <- "histogram"
    return(tmp)
  }
  
  # New axis labels function
  draw.x.axis.cor <- function(xlim, llox) {
    xx <- seq(ceiling(xlim[1]), floor(xlim[2]))
    x_axis <- list(ticks=c(), labels=list())
    if (is.na(llox)) {
      for (x in xx) {
        if (x>=3) {
          # label <- scales::math_format(10^.x)
          label <- bquote(10^.(x))
        } else {
          label <- 10^x
        }
        x_axis$ticks[length(x_axis$ticks)+1] <- x
        x_axis$labels[[length(x_axis$labels)+1]] <- label
      }
    } else {
      # for (x in xx) {
      #   if (x>log10(llox*1.8)) {
      #     if (log10(llox)==x) { label <- "lod" } else if (x>=3) { label <- bquote(10^.(x)) } else { label <- 10^x }
      #     x_axis$ticks[length(x_axis$ticks)+1] <- x
      #     x_axis$labels[[length(x_axis$labels)+1]] <- label
      #   }
      # }
      # axis(1, at=log10(llox), labels=config$llox_label)
    }
    if (length(xx)<=3) {
      for (i in 2:length(xx)) {
        x=xx[i-1]
        if (x>=3) {
          label <- bquote(3%*%10^.(x))
        } else {
          label <- 3*10^x
        }
        x_axis$ticks[length(x_axis$ticks)+1] <- x+log10(3)
        x_axis$labels[[length(x_axis$labels)+1]] <- label
      }
    }
    return(x_axis)
  }  
  
}



######################################.
##### Import Cox model estimates #####
######################################.

if (cfg2$run_analysis && any(unlist(c(cfg2$plot_cve,cfg2$plot_risk))=="Cox")) {
  
  path1 <- paste0(cfg2$folder_local, "output/", cfg2$cr2_trial, "/",
                  cfg2$cr2_COR, "/marginalized.risk.Rdata")
  path2 <- gsub("risk.Rdata", "risk.no.marker.Rdata", path1)
  load(path1) # risks.all.1
  load(path2) # overall.ve, res.plac.cont, res.vacc.cont, prev.plac, prev.vacc
  rm(risks.all.2,risks.all.3,path1,path2)
  risks <- risks.all.1[[cfg2$cr2_marker]]
  cox_cve_boot <- 1 - t( t(risks$boot)/res.plac.cont[2:(1+ncol(risks$boot))] )
  cox_cve_cis <- apply(cox_cve_boot, 1, function (x) {
    quantile(x, c(0.025,0.975))
  })
  cox_risk_boot <- risks$boot
  cox_risk_cis <- apply(cox_risk_boot, 1, function (x) {
    quantile(x, c(0.025,0.975))
  })
  
  if (cfg2$plot_cve$overall=="Cox") {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = c(999,999),
      y = rep(overall.ve[[1]], 2),
      curve = rep("Overall VE", 2),
      ci_lo = rep(overall.ve[[2]], 2),
      ci_hi = rep(overall.ve[[3]], 2),
      tag = c("Overall L", "Overall R")
    ))
  }
  
  if ("Cox" %in% cfg2$plot_cve$est) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(1-risks$prob/res.plac.cont["est"]),
      curve = rep("CVE, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_cve_cis[1,]),
      ci_hi = as.numeric(cox_cve_cis[2,]),
      tag = rep("Cox", length(as.numeric(risks$marker)))
    ))
  }
  
  if (cfg2$plot_risk$overall=="Cox") {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = rep(999,4),
      y = c(rep(prev.plac["est"], 2), rep(prev.vacc["est"], 2)),
      curve = c(rep("Placebo overall",2), rep("Vaccine overall",2)),
      ci_lo = c(rep(prev.plac["2.5%"], 2), rep(prev.vacc["2.5%"], 2)),
      ci_hi = c(rep(prev.plac["97.5%"], 2), rep(prev.vacc["97.5%"], 2)),
      tag = rep(c("Overall L", "Overall R"),2)
    ))
  }
  
  if ("Cox" %in% cfg2$plot_risk$est) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(risks$prob),
      curve = rep("Risk, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_risk_cis[1,]),
      ci_hi = as.numeric(cox_risk_cis[2,]),
      tag = rep("Cox", length(as.numeric(risks$marker)))
    ))
  }
  
}



#################################################.
##### Calculate Kaplan-Meier risk estimates #####
#################################################.

if (cfg2$run_analysis && cfg2$plot_risk$overall=="KM") {
  
  srv_ov <- survfit(Surv(dat_orig$y_star,dat_orig$delta_star)~1)
  risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
  ci_lo_ov <- 1 - srv_ov$upper[which.min(abs(srv_ov$time-C$t_e))]
  ci_hi_ov <- 1 - srv_ov$lower[which.min(abs(srv_ov$time-C$t_e))]
  
  plot_data_risk <- rbind(plot_data_risk, data.frame(
    x = c(999,999),
    y = rep(risk_ov, 2),
    curve = rep("Overall risk", 2),
    ci_lo = rep(ci_lo_ov, 2),
    ci_hi = rep(ci_hi_ov, 2),
    tag = c("Overall L", "Overall R")
  ))
  
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
  if (!is.na(cfg2$v$ph2)) {
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
  }
  
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



####################################################.
##### Helper function to process CVE estimates #####
####################################################.

{
  
  process_ests <- function(ests, a_grid, run_cve=F, lab_risk=NA, lab_cve=NA,
                           tag="0") {
    
    # Extract risk estimates
    ests_risk <- ests$est
    ci_lo_risk <- ests$ci_lo %>% pmax(0) %>% pmin(1)
    ci_hi_risk <- ests$ci_hi %>% pmax(0) %>% pmin(1)
    
    # Compute CVE estimates
    if (run_cve) {
      if (!exists("df_ct")) { stop("df_ct does not exist") }
      rate_ct <- get.marginalized.risk.no.marker(df_ct, C$t_e)
      cve <- Vectorize(function(x) { 1 - x/rate_ct })
      ests_cve <- cve(ests_risk)
      ci_lo_cve <- cve(ci_hi_risk) %>% pmax(0) %>% pmin(1)
      ci_hi_cve <- cve(ci_lo_risk) %>% pmax(0) %>% pmin(1)
    }
    
    plot_data_risk <- data.frame(
      x = a_grid,
      y = ests_risk,
      curve = rep(lab_risk, length(ests_risk)),
      ci_lo = ci_lo_risk,
      ci_hi = ci_hi_risk,
      tag = rep(tag, length(ests_risk))
    )
    if (run_cve) {
      plot_data_cve <- data.frame(
        x = a_grid,
        y = ests_cve,
        curve = rep(lab_cve, length(ests_cve)),
        ci_lo = ci_lo_cve,
        ci_hi = ci_hi_cve,
        tag = rep(tag, length(ests_cve))
      )
    } else {
      plot_data_cve <- NA
    }
    
    return(list(risk=plot_data_risk, cve=plot_data_cve))
    
  }
  
}



#####################################.
##### Data analysis (Grenander) #####
#####################################.

if (cfg2$run_analysis &&
    any(unlist(c(cfg2$plot_cve$est,cfg2$plot_risk$est))=="Grenander")) {
  
  # Obtain estimates
  return_extra <- c("Phi_n_inv_notrans")
  if (cfg2$run_debug$objs) {
    return_extra <- c(return_extra, "Phi_n_inv", "Psi_n", "omega_n", "f_aIw_n",
                      "S_n", "gcm", "dGCM", "etastar_n")
  }
  if (cfg2$run_debug$gren_var) {
    return_extra <- c(return_extra, "deriv_theta_n", "f_a_n", "gamma_n")
  }
  if (cfg2$params$marg=="Theta") { return_extra <- c(return_extra,"etastar_n") }
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=101)
  ests <- est_curve(
    dat_orig = dat_orig,
    estimator = "Grenander",
    params = cfg2$params,
    points = a_grid,
    dir = "decr",
    return_extra = return_extra
  )
  
  saveRDS(ests, paste0(cfg2$analysis," plots/ests_g_",cfg2$tid,".rds")) # !!!!!
  
  run_cve <- as.logical("Grenander" %in% cfg2$plot_cve$est)
  ests2 <- process_ests(ests, a_grid, run_cve=run_cve,
                        lab_risk="Risk, nonparametric",
                        lab_cve="CVE, nonparametric", tag="Gren")
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



#################################.
##### Data analysis (Qbins) #####
#################################.

if (cfg2$run_analysis &&
    any(unlist(c(cfg2$plot_cve$est,cfg2$plot_risk$est))=="Qbins")) {
  
  # Obtain estimates
  return_extra <- c("Phi_n_inv_notrans")
  if (cfg2$run_debug$objs) {
    return_extra <- c(return_extra, "Phi_n_inv", "Psi_n", "omega_n", "f_aIw_n",
                      "S_n", "gcm", "dGCM", "etastar_n")
  }
  if (cfg2$run_debug$gren_var) {
    return_extra <- c(return_extra, "deriv_theta_n", "f_a_n", "gamma_n")
  }
  if (cfg2$params$marg=="Theta") { return_extra <- c(return_extra,"etastar_n") }
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=101)
  ests <- est_curve(
    dat_orig = dat_orig,
    estimator = "Qbins",
    params = cfg2$params,
    points = a_grid,
    dir = "decr",
    return_extra = return_extra
  )
  
  saveRDS(ests, paste0(cfg2$analysis," plots/ests_q_",cfg2$tid,".rds")) # !!!!!
  
  run_cve <- as.logical("Qbins" %in% cfg2$plot_cve$est)
  ests2 <- process_ests(ests, a_grid, run_cve=run_cve,
                        lab_risk="Risk, Qbins", lab_cve="CVE, Qbins",
                        tag="Qbins")
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



#####################################.
##### Data analysis (Cox gcomp) #####
#####################################.

if (cfg2$run_analysis &&
    any(unlist(c(cfg2$plot_cve$est,cfg2$plot_risk$est))=="Cox gcomp")) {
  
  # Obtain estimates
  return_extra <- c("Phi_n_inv_notrans")
  a_orig <- dat_orig$a[!is.na(dat_orig$a)]
  # a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=101)
  a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=5) # !!!!!
  ests <- est_curve(
    dat_orig = dat_orig,
    estimator = "Cox gcomp",
    params = cfg2$params,
    points = a_grid,
    dir = "decr",
    return_extra = return_extra
  )
  
  saveRDS(ests, paste0(cfg2$analysis," plots/ests_c_",cfg2$tid,".rds")) # !!!!!
  
  run_cve <- as.logical("Cox gcomp" %in% cfg2$plot_cve$est)
  ests2 <- process_ests(ests, a_grid, run_cve=run_cve,
                        lab_risk="Risk, Cox (analytic)",
                        lab_cve="CVE, Cox (analytic)", tag="Cox gcomp")
  
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



###########################################.
##### Data analysis (hypothesis test) #####
###########################################.

if (cfg2$run_hyptest) {
  
  res <- test_2(
    dat_orig = dat_orig,
    alt_type = "two-tailed", # decr
    params = list(),
    return_extras = T # !!!!!
  )
  res_df <- data.frame(reject=res$reject, p_val=res$p_val, beta_n=res$beta_n,
                       sd_n=res$sd_n, var_n=res$var_n)
  write.table(
    res_df,
    file = paste0(cfg2$analysis," plots/hyptest_",cfg2$tid,".csv"),
    sep = ",",
    row.names = FALSE
  )
  saveRDS(
    res$Theta_os_n,
    paste0(cfg2$analysis," plots/Theta_os_n_",cfg2$tid,".rds")
  )
  
}



#############################.
##### Plotting function #####
#############################.

{
  
  #' Create a controlled vaccine efficacy (CVE) plot
  #' 
  #' @param plot_data A dataframe containing the following columns:
  #'   - `x`: X-value, in log10(marker) units
  #'   - `y`: Y-value (CVE)
  #'   - `curve`: One of c("Overall VE", "CVE, Cox model", "CVE, nonparametric",
  #'     "Placebo overall", "Vaccine overall", "Risk, Cox model",
  #'     "Risk, nonparametric")
  #'   - `ci_lo`: Lower confidence bound for CVE
  #'   - `ci_hi`: Upper confidence bound for CVE
  #' @param which One of c("CVE", "Risk")
  #' @param zoom_x Either a numeric vector of length 2 representing the plot X
  #'     limits, or the string "zoomed", in which case the plot will be zoomed
  #'     into the Grenander cutoff quantiles (e.g. 5%/95%), with 10% padding on
  #'     each side. Defaults to the width of the histogram plus 5% padding on
  #'     each side.
  #' @param zoom_y Either a numeric vector of length 2 representing the plot Y
  #' limits or the string "zoomed", in which case the plot height will be zoomed
  #'     such that the upper/lower CIs are in the frame with 10% padding above
  #'     and below. Defaults to 0%--105%.
  #' @param labs A named list of plot labels containing; names include
  #'     c("title", "x", "y")
  #' @param hst A histogram object returned by get.marker.histogram()
  #' @param rr_y_axis Boolean; if true, a secondary risk ratio Y-axis will be
  #'     displayed
  #' @param log10_x_axis Boolean; if true, the X-axis values are log-10
  #' @return A ggplot2 plot object
  #' @notes
  #'   - This function plots pointwise estimates and confidence intervals from
  #'     both Cox model and nonparametric approaches
  create_plot <- function(plot_data, which, zoom_x=NA, zoom_y=NA, labs, hst,
                          rr_y_axis=F, log10_x_axis=F) {
    
    # # !!!!! Add lod/2 ?????
    # lod12 <- min(dat_orig$a,na.rm=T)
    
    # Change curve labels to factors and set color scale
    # goldenrod3, darkseagreen3, darkslategray3, darkorange3, darkgoldenrod3, cyan3
    curves <- c("Placebo overall", "Vaccine overall", "Overall VE",
                "Risk, Cox model", "CVE, Cox model", "Risk, Qbins",
                "CVE, Qbins", "Risk, nonparametric", "CVE, nonparametric",
                "Risk, Cox (analytic)", "CVE, Cox (analytic)", "Control",
                "VRC01 10mg/kg", "VRC01 30mg/kg", "VRC01 Pooled")
    curve_colors <- c("darkgrey", "darkgrey", "darkgrey", "darkorchid3",
                      "darkorchid3", "firebrick3", "firebrick3",
                      "deepskyblue3", "deepskyblue3", "deepskyblue3",
                      "deepskyblue3", "deepskyblue3", "darkorchid3",
                      "firebrick3", "darkolivegreen3")
    names(curve_colors) <- curves
    indices <- which(curves %in% unique(plot_data$curve))
    curve_colors <- curve_colors[indices]
    plot_data$curve <- factor(plot_data$curve, levels=curves[indices])
    
    # Replace placeholder "Overall" X-values
    plot_data[plot_data$tag=="Overall L","x"] <- min(hst$breaks)
    plot_data[plot_data$tag=="Overall R","x"] <- max(hst$breaks)
    
    # Set default zoom levels
    if (is.na(zoom_x[1])) {
      z_x_L <- min(plot_data$x)
      z_x_R <- max(plot_data$x)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, tag %in% c("Gren", "Qbins") & !is.na(y))$x
      z_x_L <- min(zz, na.rm=T)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.1*(z_x_R-z_x_L),
                  z_x_R + 0.1*(z_x_R-z_x_L))
    }
    if (is.na(zoom_y[1])) {
      zoom_y <- c(0,1)
      zoom_y[2] <- zoom_y[2] + 0.05*(zoom_y[2]-zoom_y[1])
    } else if (zoom_y[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, x>=z_x_L & x<=z_x_R)
      z_y_L <- min(zz$ci_lo, na.rm=T)
      z_y_U <- max(zz$ci_hi, na.rm=T)
      zoom_y <- c(z_y_L - 0.1*(z_y_U-z_y_L),
                  z_y_U + 0.1*(z_y_U-z_y_L))
    } else if (zoom_y[1]=="zoomed (risk)") {
      zz <- dplyr::filter(plot_data, x>=z_x_L & x<=z_x_R)
      z_y_L <- 0
      z_y_U <- max(plot_data_risk$ci_hi, na.rm=T)
      zoom_y <- c(z_y_L - 0.1*(z_y_U-z_y_L),
                  z_y_U + 0.1*(z_y_U-z_y_L))
    }
    
    # Generate histogram data
    ymax <- 0.6 * (hst$counts/max(hst$counts)) * (zoom_y[2]/1.05-zoom_y[1]) +
      zoom_y[1]
    hist_data <- data.frame(
      xmin = hst$breaks[-length(hst$breaks)],
      xmax = hst$breaks[-1],
      ymin = rep(zoom_y[1], length(hst$counts)),
      ymax = ymax
    )
    
    # Create and return ggplot2 object
    # Note: using geom_rect for the histogram so that it can be shifted
    #     vertically
    if (rr_y_axis) {
      sec.axis <- sec_axis(~1-., breaks=seq(-1,1,0.1),
                           name="Risk ratio (vaccine/placebo)")
    } else {
      sec.axis <- waiver()
    }
    y_ticks <- ifelse(which=="CVE", 0.1, 0.01)
    plot <- ggplot(plot_data, aes(x=x, y=y, color=curve)) +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi, fill=curve),
        alpha = 0.05,
        linetype = "dotted"
      ) +
      geom_rect(
        aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
        hist_data,
        alpha = 0.3,
        fill = "forestgreen",
        inherit.aes = F
      ) +
      geom_line(size=0.7) +
      scale_y_continuous(
        labels = scales::label_percent(accuracy=1),
        breaks = seq(-1,1,y_ticks),
        minor_breaks = NULL,
        sec.axis = sec.axis
      ) +
      theme(
        panel.border = element_rect(color="#bbbbbb", fill=NA),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      coord_cartesian(xlim=zoom_x, ylim=zoom_y, expand=F) +
      scale_color_manual(
        values = curve_colors,
        breaks = curves[!(curves %in% c("Placebo overall", "Vaccine overall"))]
      ) +
      scale_fill_manual(
        values = curve_colors,
        breaks = curves[!(curves %in% c("Placebo overall", "Vaccine overall"))]
      ) +
      theme(legend.position="bottom") +
      labs(title=labs$title, x=labs$x, y=labs$y, color=NULL, fill=NULL)
    if (log10_x_axis) {
      xlim <- c(min(dat_orig$a, na.rm=T), max(dat_orig$a, na.rm=T))
      x_axis <- draw.x.axis.cor(xlim, NA)
      
      # !!!!! Temp Janssen
      if (cfg2$analysis=="Janssen" && cfg2$marker=="Day29pseudoneutid50") {
        x_axis$ticks[4] <- log10(2.7426)
        x_axis$labels[[4]] <- "LLOQ"
      }
      
      plot <- plot + scale_x_continuous(
        labels = do.call(expression,x_axis$labels),
        breaks = x_axis$ticks
      )
    }
    if (which=="Risk") {
      x_r <- max(hst$breaks)
      y_plac <- filter(plot_data_risk, curve=="Placebo overall")[1,"y"]
      y_vacc <- filter(plot_data_risk, curve=="Vaccine overall")[1,"y"]
      plot <- plot + annotate("text", x=x_r, y=y_plac, label="Placebo overall",
                              size=2.5, hjust=1, vjust=-0.5)
      plot <- plot + annotate("text", x=x_r, y=y_vacc, label="Vaccine overall",
                              size=2.5, hjust=1, vjust=-0.5)
    }
    
    return(plot)
    
  }
  
}



###################################.
##### Generate and save plots #####
###################################.

if (nrow(plot_data_risk)>0 || nrow(plot_data_cve)>0) {
  
  # Create cutoff values corresponding to cfg2$qnt
  cutoffs <- lapply(cfg2$qnt, function(qnt) {
    as.numeric(quantile(dat_orig$a, na.rm=T, probs=qnt))
  })
  
  # Trim estimates at specified quantiles
  trim_plot_data <- function(plot_data) {
    
    which_curves <- names(cfg2$qnt)
    cut_lo <- sapply(plot_data$curve, function(curve) {
      ifelse(curve %in% which_curves, cutoffs[[curve]][1], NA)
    }, USE.NAMES=F)
    cut_hi <- sapply(plot_data$curve, function(curve) {
      ifelse(curve %in% which_curves, cutoffs[[curve]][2], NA)
    }, USE.NAMES=F)
    rows_1 <- which(plot_data$curve %in% which_curves) # !!!!! Might not be needed anymore
    rows_2 <- which(plot_data$x <= cut_lo | plot_data$x >= cut_hi)
    rows <- intersect(rows_1, rows_2)
    plot_data[rows, c("y", "ci_lo", "ci_hi")] <- NA
    return(plot_data)
  }
  
  hst <- get.marker.histogram(
    marker = dat_orig$a[!is.na(dat_orig$a)],
    wt = dat_orig$weights[!is.na(dat_orig$a)],
    trial = cfg2$cr2_trial
  )
  
  if (nrow(plot_data_risk)>0) {
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", cfg2$t_e)
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_risk),
      which = "Risk",
      zoom_x = cfg2$zoom_x, # This might mess up AMP plots
      zoom_y = "zoomed (risk)", # This might mess up AMP plots
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      log10_x_axis = T
    )
    ggsave(
      filename = paste0(cfg2$analysis," plots/plot_risk_",cfg2$tid,".pdf"),
      plot=plot, device="pdf", width=6, height=4
    )
  }
  
  if (nrow(plot_data_cve)>0) {
    cfg2$lab_y <- paste0("Controlled VE against ", cfg2$endpoint,
                         " by day ", cfg2$t_e)
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_cve),
      which = "CVE",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y,
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      rr_y_axis = T,
      log10_x_axis = T
    )
    ggsave(
      filename = paste0(cfg2$analysis," plots/plot_cve_",cfg2$tid,".pdf"),
      plot=plot, device="pdf", width=6, height=4
    )
  }
  
}



######################################################.
##### Generate and save AMP plots (run manually) #####
######################################################.

if (F) {
  
  # Temp
  cfg2$t_e <- 595
  cfg2$lab_title <- c("HVTN703/HPTN081", "HVTN704/HPTN085", "Pooled AMP trials")
  cfg2$amp_tx2 <- rep(c("Control", "VRC01 10mg/kg", "VRC01 30mg/kg",
                    "VRC01 Pooled"), 3)
  
  # Generate histograms and KM objects
  # !!!!! Eventually replace this section and instead save the objects above,
  #       accounting for the fact that 13:15 don't fall into the framework
  for (i in c(1:15)) {
    
    # TEMP: Generate dat_orig
    {
      df_analysis <- df_raw
      if (i %in% c(1:4,13)) {
        df_analysis %<>% dplyr::filter(protocol=="HVTN 703")
      }
      if (i %in% c(5:8,14)) {
        df_analysis %<>% dplyr::filter(protocol=="HVTN 704")
      }
      if (i %in% c(1,5,9)) { df_analysis %<>% dplyr::filter(tx=="C3") }
      if (i %in% c(2,6,10)) { df_analysis %<>% dplyr::filter(tx=="T1") }
      if (i %in% c(3,7,11)) { df_analysis %<>% dplyr::filter(tx=="T2") }
      if (i %in% c(4,8,12)) { df_analysis %<>% dplyr::filter(tx_pool=="T1+T2") }
      dat_orig <- list(
        y_star = df_analysis[["hiv1survday"]],
        delta_star = df_analysis[["hiv1event"]],
        weights = rep(1, nrow(df_analysis)),
        a = df_analysis[["bweight"]]
      )
    }
    
    # Generate a_orig
    a_orig <- dat_orig$a[!is.na(dat_orig$a)]
    saveRDS(a_orig, paste0(cfg2$analysis, " plots/a_orig_", i, ".rds"))
    
    if (i %in% c(1,5,9,13,14,15)) {
      
      # Generate histogram
      hst <- get.marker.histogram(
        marker = dat_orig$a[!is.na(dat_orig$a)],
        wt = dat_orig$weights[!is.na(dat_orig$a)],
        trial = F
      )
      saveRDS(hst, paste0(cfg2$analysis, " plots/hist_", i, ".rds"))
      
      # Generate KM object
      srv_ov <- survfit(Surv(dat_orig$y_star,dat_orig$delta_star)~1)
      risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
      ci_lo_ov <- 1 - srv_ov$upper[which.min(abs(srv_ov$time-C$t_e))]
      ci_hi_ov <- 1 - srv_ov$lower[which.min(abs(srv_ov$time-C$t_e))]
      km <- data.frame(
        x = c(999,999),
        y = rep(risk_ov, 2),
        curve = rep("Overall risk", 2),
        ci_lo = rep(ci_lo_ov, 2),
        ci_hi = rep(ci_hi_ov, 2),
        tag = c("Overall L", "Overall R")
      )
      saveRDS(km, paste0(cfg2$analysis, " plots/km_", i, ".rds"))
      
    }
    
  }
  
  # Mapping of plot objects
  plot_map <- list(
    lab_title = c(1,2,3,1,2,3,1,2,3),
    gren = list(1, 5, 9, c(1:3), c(5:7), c(9:11), c(1,4), c(5,8), c(9,12)),
    hist = c(1,5,9,13,14,15,13,14,15),
    overall = c(1,5,9,13,14,15,13,14,15)
  )
  
  for (i in c(1:9)) {
    
    # Read in data objects
    plot_data_risk <- readRDS(paste0(cfg2$analysis, " plots/km_",
                         plot_map$overall[i], ".rds"))
    hst <- readRDS(paste0(cfg2$analysis, " plots/hist_",
                          plot_map$hist[i], ".rds"))
    for (j in plot_map$gren[[i]]) {
      a_orig <- readRDS(paste0(cfg2$analysis, " plots/a_orig_",
                             j, ".rds"))
      a_grid <- seq(from=min(a_orig), to=max(a_orig), length.out=101)
      ests <- readRDS(paste0(cfg2$analysis, " plots/ests_g_",
                             j, ".rds"))
      ests2 <- process_ests(ests, a_grid, run_cve=F, lab_risk=cfg2$amp_tx2[j],
                            tag="Gren")
      plot_data_risk <- rbind(plot_data_risk, ests2$risk)
    }
    
    # Generate and save plot
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", cfg2$t_e)
    lab_title <- cfg2$lab_title[plot_map$lab_title[i]]
    plot <- create_plot(
      plot_data = plot_data_risk,
      which = "Risk",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y,
      labs = list(title=lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst
    )
    ggsave(
      filename = paste0(cfg2$analysis," plots/plot_risk_",i,".pdf"),
      plot=plot, device="pdf", width=6, height=4
    )
    
  }
  
}



###############################.
##### Grenander debugging #####
###############################.

if (cfg2$run_debug$objs) {
  
  # Debugging: Examine intermediate objects
  
  grid <- round(seq(0,1,0.01),2)
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
  
  int_data <- data.frame(
    x = rep(grid,3),
    y = c(ests$Psi_n(grid), gcm(grid), ests$dGCM(grid)),
    which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*theta_n)"), each=101)
  )
  plot1 <- ggplot(int_data, aes(x=x, y=y, color=which)) +
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
    int_data2 <- data.frame(
      x = grid,
      y = S_n(grid)
      # which = rep(c("",""),each=101)
    )
    ggplot(int_data2, aes(x=x, y=y)) + # color=which
      geom_line() +
      theme(legend.position="bottom")
    
    # # omega_n
    # int_data3 <- data.frame(
    #   x = grid,
    #   y = omega_n(grid)
    #   # which = rep(c("",""),each=101)
    # )
    # ggplot(int_data3, aes(x=x, y=y)) + # color=which
    #   geom_line() +
    #   theme(legend.position="bottom")
    
    # etastar_n
    int_data4 <- data.frame(
      x = grid,
      y = etastar_n(grid)
      # which = rep(c("",""),each=101)
    )
    ggplot(int_data4, aes(x=x, y=y)) + # color=which
      geom_line() +
      theme(legend.position="bottom") +
      labs(title="etastar_n")
    
  }
  
}

if (cfg2$run_debug$gren_var) {
  
  # Debugging: Grenander variance scale factor components
  
  print("Grenander variance scale factor components")
  print("deriv_theta_n")
  print(ests$deriv_theta_n(seq(0,1,0.05)))
  print("f_a_n")
  print(ests$f_a_n(seq(0,1,0.05)))
  # print("gamma_n")
  # print(ests$gamma_n(seq(0,1,0.05)))
  print("deriv_theta_n*f_a_n")
  print(ests$deriv_theta_n(seq(0,1,0.05))*ests$f_a_n(seq(0,1,0.05)))
  # print("deriv_theta_n*gamma_n")
  # print(ests$deriv_theta_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  # print("f_a_n*gamma_n")
  # print(ests$f_a_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  # print("deriv_theta_n*f_a_n*gamma_n")
  # print(ests$deriv_theta_n(seq(0,1,0.05))*ests$f_a_n(seq(0,1,0.05))*
  #         ests$gamma_n(seq(0,1,0.05)))
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
      y = c(ests_cve[1],ests_cve,ests_gcomp,rep(ve_overall,2)),
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
