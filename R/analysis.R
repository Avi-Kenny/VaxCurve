# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

#################.
##### Setup #####
#################.

{
  # Print timestamp
  print(paste("START:", Sys.time()))
  
  # Choose analysis
  # "Janssen" "Moderna" "AMP" "AZD1222" "Janssen (partA)" "Profiscov"
  # "HVTN 705 (primary)" "HVTN 705 (all)" "RV144" "HVTN 705 (second)"
  # "HVTN 705 (compare RV144)"
  cfg2 <- list(analysis="Janssen (partA)", seed=1)
  
  # Set proper task ID variable
  if (cluster_config$js=="slurm") {
    .tid_var <- "SLURM_ARRAY_TASK_ID"
  } else if (cluster_config$js=="ge") {
    .tid_var <- "SGE_TASK_ID"
  } else if (cluster_config$js=="") {
    .tid_var <- NA
  } else {
    stop("Invalid cluster_config$js")
  }
  
  # Run multiple analyses at once
  {
    
    # # RV144 vs HVTN 705: uncomment to run
    # ..tid <- as.integer(Sys.getenv(.tid_var))
    # if (..tid %in% c(1:6)) {
    #   cfg2$analysis <- "RV144"
    #   .tid_lst = list(as.character(round(..tid)))
    # } else if (..tid %in% c(7:12)) {
    #   cfg2$analysis <- "HVTN 705 (compare RV144)"
    #   .tid_lst = list(as.character(round(..tid-6)))
    # }
    # names(.tid_lst) = .tid_var
    # do.call(Sys.setenv, .tid_lst)
    
    # 128 plots across 4 analyses: uncomment to run
    ..tid <- as.integer(Sys.getenv(.tid_var))
    if (..tid %in% c(1:4)) { # 4 markers
      cfg2$analysis <- "Janssen"
      .tid_lst = list(as.character(round(..tid)))
    } else if (..tid %in% c(5:14)) { # 10 markers
      cfg2$analysis <- "Moderna"
      .tid_lst = list(as.character(round(..tid-4)))
    } else if (..tid %in% c(15:53)) { # 39 markers
      cfg2$analysis <- "HVTN 705 (all)"
      .tid_lst = list(as.character(round(..tid-14)))
    } else if (..tid %in% c(54:65)) { # 12 markers
      cfg2$analysis <- "HVTN 705 (ICS)"
      .tid_lst = list(as.character(round(..tid-53)))
    } else if (..tid %in% c(66:123)) { # 58 markers
      cfg2$analysis <- "Janssen (partA)"
      .tid_lst = list(as.character(round(..tid-65)))
    } else if (..tid %in% c(124:128)) { # 5 markers
      cfg2$analysis <- "RV144"
      .tid_lst = list(as.character(round(..tid-123)))
    }
    names(.tid_lst) = .tid_var
    do.call(Sys.setenv, .tid_lst)
    
  }
  
  # Set seed
  set.seed(cfg2$seed)
  
  # Set analysis-specific flags
  # Note: some flags are set at the end of the "Setup" block because they are
  #       dependent on cfg2 variables
  flags <- list(
    run_hyptest = F,
    run_mediation = T,
    hvtn705_abstract_fig = F,
    table_of_vals = F,
    save_data_objs = F,
    save_plot_objs = F,
    save_diagnostics = F,
    paper_npcve = F,
    paper_cox = F,
    hvtn124_plot = F,
    partA_mnscrpt2 = F
  )
  
  # Set default cfg2 values (+ those common to multiple analyses)
  {
    # Estimator types: "Cox import", "Grenander", "Cox gcomp", "Qbins",
    #     "Cox (spline 4 df)", "Cox edge"
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp"))
    cfg2$edge_corr <- c(FALSE,TRUE)
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.05,0.95),
      "CVE, nonparametric" = c(0.05,0.95),
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975),
      "Risk, Cox (spline 4 df)" = c(0.025,0.975),
      "CVE, Cox (spline 4 df)" = c(0.025,0.975),
      "Risk, Cox (analytic)" = c(0.025,0.975),
      "CVE, Cox (analytic)" = c(0.025,0.975),
      "Risk, Cox (basic)" = c(0.025,0.975),
      "CVE, Cox (basic)" = c(0.025,0.975),
      "Risk, Cox (spline 4 df)" = c(0.025,0.975),
      "CVE, Cox (spline 4 df)" = c(0.025,0.975),
      "Risk, Cox (edge)" = c(0.025,0.975),
      "CVE, Cox (edge)" = c(0.025,0.975)
    )
    cfg2$plots <- c("Risk", "CVE")
    cfg2$endpoint <- c("COVID", "HIV-1 infection", "HIV")
    cfg2$cr2_marker <- c(1:1000)
    cfg2$params = list(
      g_n_type = "binning",
      deriv_type = "m-spline",
      ci_type = "transformed 2",
      Q_n_type = "survSL"
    )
    cfg2$dir <- c("decr", "incr")
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y_cve <- "zoomed"
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_risk_max <- NA
    cfg2$more_ticks <- c(1,2)
    cfg2$llox_label <- NULL
    cfg2$llox <- NULL
    cfg2$density_type <- "histogram"
    
  }
  
  # Set up analysis-specific configuration variables. Each row in the cfg2$map
  #     dataframe represents the set of indices to use for a particular
  #     analysis.
  if (cfg2$analysis=="Janssen") {
    
    # Override default config
    # cfg2$params$ci_type <- "truncated" # Used historically
    cfg2$zoom_x <- NA
    cfg2$zoom_y_cve <- NA
    
    # Analysis-specific config
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29ADCP")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to RBD: Day 29", "PsV Neutralization 50% Titer: Day 29", "Phagocytic Score: Day 29")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Phagocytic Score (=s)")
    cfg2$t_0 <- 54
    cfg2$dataset <- c("janssen_pooled_real_data_processed_with_riskscore.csv", "janssen_pooled_realADCP_data_processed_with_riskscore.csv")
    cfg2$folder_local <- "Janssen data/"
    # cfg2$folder_cluster <- "Z:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$folder_cluster <- "C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Marco Carone/Project - VaxCurve/VaxCurve/R/Janssen data/"
    cfg2$cr2_trial <- c("janssen_pooled_real", "janssen_pooled_realADCP")
    cfg2$cr2_COR <- "D29IncludeNotMolecConfirmedstart1"
    cfg2$v <- list(
      id = "Ptid",
      time = "EventTimePrimaryIncludeNotMolecConfirmedD29",
      event = "EventIndPrimaryIncludeNotMolecConfirmedD29",
      wt = "wt.D29start1",
      ph1 = "ph1.D29start1",
      ph2 = "ph2.D29start1",
      covariates = "~. + risk_score + as.factor(Region)"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(1, 4),
      marker = c(1,2,3,4),
      lab_x = c(1,2,3,4),
      lab_title = c(1,2,3,4),
      t_0 = rep(1, 4),
      dataset = c(1,1,1,2),
      cr2_trial = c(1,1,1,2),
      cr2_COR = rep(1,4),
      cr2_marker = c(1,2,3,1),
      edge_corr = c(2,2,2,2),
      v_id = rep(1, 4),
      v_time = rep(1, 4),
      v_event = rep(1, 4),
      v_wt = rep(1, 4),
      v_ph1 = rep(1, 4),
      v_ph2 = rep(1, 4),
      v_covariates = rep(1, 4),
      dir = rep(1, 4),
      zoom_x = rep(1, 4),
      zoom_y_cve = rep(1, 4),
      zoom_y_risk = rep(1, 4),
      more_ticks = rep(1, 4)
    )
    
  }
  
  if (cfg2$analysis=="Moderna") {
    
    # Override default config
    # cfg2$params$ci_type <- "truncated" # Used historically
    cfg2$zoom_y_risk <- list(c(-0.002,0.072))
    cfg2$llox_label <- "LOD" # NEW
    cfg2$llox <- c(0.3076,1.594,2.42,15.02,22.66) # NEW
    
    # Analysis-specific config
    cfg2$marker <- c("Day29bindSpike", "Day57bindSpike", "Day29bindRBD", "Day57bindRBD", "Day29pseudoneutid50", "Day57pseudoneutid50", "Day29pseudoneutid80", "Day57pseudoneutid80", "Day29liveneutmn50", "Day57liveneutmn50")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to Spike: Day 57", "Binding Antibody to RBD: Day 29", "Binding Antibody to RBD: Day 57", "PsV Neutralization 50% Titer: Day 29", "PsV Neutralization 50% Titer: Day 57", "PsV Neutralization 80% Titer: Day 29", "PsV Neutralization 80% Titer: Day 57", "Live Virus Micro Neut 50% Titer: Day 29", "Live Virus Micro Neut 50% Titer: Day 57")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Pseudovirus-nAb ID80 (IU80/ml) (=s)", "Live Virus-mnAb ID50 (IU50/ml) (=s)")
    cfg2$t_0 <- c(126,100) # Try changing to 0
    cfg2$dataset <- "P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv"
    cfg2$folder_local <- "Moderna data/"
    cfg2$folder_cluster <- "Z:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$cr2_trial <- "moderna_real"
    cfg2$cr2_COR <- c("D29", "D57")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD29", "EventTimePrimaryD57"),
      event = c("EventIndPrimaryD29", "EventIndPrimaryD57"),
      wt = c("wt.D29", "wt.D57"),
      ph1 = c("ph1.D29", "ph1.D57"),
      ph2 = c("ph2.D29", "ph2.D57"),
      covariates = "~. + MinorityInd + HighRiskInd + risk_score"
    )
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(1, 10),
      marker = c(1,2,3,4,5,6,7,8,9,10),
      lab_x = c(1,1,2,2,3,3,4,4,5,5),
      lab_title = c(1,2,3,4,5,6,7,8,9,10),
      t_0 = c(1,2,1,2,1,2,1,2,1,2),
      dataset = rep(1, 10),
      cr2_trial = rep(1, 10),
      cr2_COR = c(1,2,1,2,1,2,1,2,1,2),
      cr2_marker = c(1,1,2,2,3,3,4,4,5,5),
      edge_corr = c(1,1,1,1,2,1,2,1,2,1),
      v_id = rep(1, 10),
      v_time = c(1,2,1,2,1,2,1,2,1,2),
      v_event = c(1,2,1,2,1,2,1,2,1,2),
      v_wt = c(1,2,1,2,1,2,1,2,1,2),
      v_ph1 = c(1,2,1,2,1,2,1,2,1,2),
      v_ph2 = c(1,2,1,2,1,2,1,2,1,2),
      v_covariates = rep(1, 10),
      dir = rep(1, 10),
      zoom_x = rep(1, 10),
      zoom_y_cve = rep(1, 10),
      zoom_y_risk = rep(1, 10),
      more_ticks = rep(1, 10),
      llox_label = rep(1, 10),
      llox = c(1,1,2,2,3,3,4,4,5,5)
    )
    
    # Flag-specific operation
    if (flags$paper_npcve || flags$paper_cox) {
      
      cfg2$zoom_x <- c("zoomed", "zoomed llox")
      cfg2$map$zoom_x <- c(1,1,1,1,2,1,2,1,2,1)
      
    }
    
  }
  
  if (cfg2$analysis=="AMP") {
    
    # Override default config
    cfg2$estimators <- list(overall="KM", cr=c("Grenander"))
    cfg2$plots <- c("Risk")
    cfg2$qnt[["Risk, nonparametric"]] <- c(0.1,0.9)
    cfg2$qnt[["CVE, nonparametric"]] <- c(0.1,0.9)
    cfg2$params$deriv_type <- "linear"
    cfg2$params$ci_type <- "regular"

    # Analysis-specific config
    cfg2$marker <- "bweight"
    cfg2$lab_title <- c("HVTN703/HPTN081", "HVTN704/HPTN085", "Pooled AMP trials")
    cfg2$lab_x <- "Body Weight (kg)"
    cfg2$t_0 <- 595
    cfg2$dataset <- "amp_survival_all.csv"
    cfg2$folder_local <- "AMP data/"
    cfg2$folder_cluster <- "Z:/vaccine/p704/analysis/datashare/avi_kenny/adata/"
    cfg2$cr2_trial <- F
    cfg2$cr2_COR <- F
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
    cfg2$amp_protocol <- c("HVTN 703", "HVTN 704", "Pooled")
    cfg2$amp_tx <- c("C3", "T1", "T2", "T1+T2")
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(2, 12),
      amp_protocol = c(1,1,1,1,2,2,2,2,3,3,3,3),
      amp_tx = c(1,2,3,4,1,2,3,4,1,2,3,4),
      marker = rep(1, 12),
      lab_x = rep(1, 12),
      lab_title = c(1,1,1,1,2,2,2,2,3,3,3,3),
      t_0 = rep(1, 12),
      dataset = rep(1, 12),
      cr2_trial = rep(1, 12),
      cr2_COR = rep(1, 12),
      cr2_marker = rep(1, 12),
      edge_corr = rep(1, 12),
      v_id = rep(1, 12),
      v_time = rep(1, 12),
      v_event = rep(1, 12),
      v_wt = rep(1, 12),
      v_ph1 = rep(1, 12),
      v_ph2 = rep(1, 12),
      v_covariates = c(1,1,1,1,1,1,1,1,2,2,2,2),
      dir = rep(1, 12),
      zoom_x = rep(1, 12),
      zoom_y_cve = rep(1, 12),
      zoom_y_risk = rep(1, 12),
      more_ticks = rep(1, 12)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (all)") {
    
    # Override default config
    cfg2$params$deriv_type <- "line"

    # Analysis-specific config
    cfg2$marker <- c("Day210ELCZ", "Day210ELMo", "Day210ADCPgp140C97ZAfib", "Day210ADCPgp140Mos1fib", "Day210IgG3gp140C97ZAfibritin40delta", "Day210IgG3gp140Mos1fibritin40delta", "Day210IgG340mdw_gp120", "Day210IgG340mdw_gp140", "Day210IgG340mdw_V1V2", "Day210IgG3gp4140delta", "Day210IgG340mdw_multi", "Day210IgG340mdw_gp120_gp140_vm", "Day210IgG50mdw_V1V2", "Day210mdw_xassay", "Day210ADCCCAP8_pk", "Day210ADCCCH58_pk", "Day210ADCCWITO_pk", "Day210ADCCCAP8_pAUC", "Day210ADCCCH58_pAUC", "Day210ADCCWITO_pAUC", "Day210ICS4AnyEnvIFNg_OR_IL2", "Day210ICS8AnyEnvIFNg_OR_IL2", "Day210IgG3AE.A244.V1V2.Tags_293F40delta", "Day210IgG3C.1086C.V1.V2.Tags40delta", "Day210IgG3gp70.001428.2.42.V1V240delta", "Day210IgG3gp70.1012.11.TC21.3257.V1V240delta", "Day210IgG3gp70.1394C9G1.V1V240delta", "Day210IgG3gp70.BF1266.431a.V1V240delta", "Day210IgG3gp70.Ce1086.B2.V1V240delta", "Day210IgG3gp70.B.CaseA.V1.V240delta", "Day210IgGAE.A244.V1V2.Tags_293F50delta", "Day210IgGC.1086C.V1.V2.Tags50delta", "Day210IgGgp70_001428.2.42.V1V250delta", "Day210IgGgp70_1012.11.TC21.3257.V1V250delta", "Day210IgGgp70_1394C9G1.V1V250delta", "Day210IgGgp70_9004SS.A3.4.V1V250delta", "Day210IgGgp70_BF1266.431a.V1V250delta", "Day210IgGgp70_Ce1086.B2.V1V250delta", "Day210IgGgp70.B.CaseA.V1.V250delta")
    cfg2$lab_title <- c("IgG to VT-C (EU/ml): Month 7", "IgG to VT-M (EU/ml): Month 7", "Average phagocytosis score to gp140 C97ZA: Month 7", "Average phagocytosis score to gp140 Mos1: Month 7", "IgG3 Net MFI to gp140 C97ZA: Month 7", "IgG3 Net MFI to gp140 Mosaic: Month 7", "IgG3 gp120 breadth (Weighted avg log10 Net MFI): Month 7", "IgG3 gp140 breadth (Weighted avg log10 Net MFI): Month 7", "IgG3 V1V2 breadth (Weighted avg log10 Net MFI): Month 7", "IgG3 Net MFI to gp41: Month 7", "IgG3 multi-epitope breadth (Wt avg log10 Net MFI): Month 7", "IgG3 gp120 + gp140 breadth (Wt avg log10 Net MFI): Month 7", "IgG V1V2 breadth (Wt avg log10 Net MFI): Month 7", "Overall maximal diversity score: Month 7", "Peak baseline-subtracted pct loss luc activity to CAP8: Month 7", "Peak baseline-subtracted pct loss luc activity to CH58: Month 7", "Peak baseline-subtracted pct loss luc activity to WITO: Month 7", "AUC baseline-subtracted pct loss luc activity to CAP8: Month 7", "AUC baseline-subtracted pct loss luc activity to CH58: Month 7", "AUC baseline-subtracted pct loss luc activity to WITO: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2: Month 7", "IgG3 Net MFI to AE.A244 V1V2 Tags 293F: Month 7", "IgG3 Net MFI to C.1086C V1V2 Tags: Month 7", "IgG3 Net MFI to gp70-001428.2.42 V1V2: Month 7", "IgG3 Net MFI to gp70-1012.11.TC21.3257 V1V2: Month 7", "IgG3 Net MFI to gp70-1394C9G1 V1V2: Month 7", "IgG3 Net MFI to gp70-BF1266 431a V1V2: Month 7", "IgG3 Net MFI to gp70-Ce1086 B2 V1V2: Month 7", "IgG3 Net MFI to gp70-B.CaseAV1V2: Month 7", "IgG Net MFI to AE.A244 V1V2 Tags 293F: Month 7", "IgG Net MFI to C.1086C V1V2 Tags: Month 7", "IgG Net MFI to gp70-001428.2.42 V1V2: Month 7", "IgG Net MFI to gp70-1012.11.TC21.3257 V1V2: Month 7", "IgG Net MFI to gp70-1394C9G1 V1V2: Month 7", "IgG Net MFI to gp70-9004SS.A3.4 V1V2: Month 7", "IgG Net MFI to gp70-BF1266.431a V1V2: Month 7", "IgG Net MFI to gp70-Ce1086.B2 V1V2: Month 7", "IgG Net MFI to gp70.B.CaseA V1V2: Month 7")
    cfg2$lab_x <- c("IgG to VT-C (=s)", "IgG to VT-M (=s)", "ADCP gp140 C97ZA (=s)", "ADCP gp140 Mos1 (=s)", "IgG3 gp140 C97ZA (=s)", "IgG3 gp140 Mosaic (=s)", "IgG3 gp120 breadth (=s)", "IgG3 gp140 breadth (=s)", "IgG3 V1V2 breadth (=s)", "IgG3 gp41 (=s)", "IgG3 multi-epitope breadth (=s)", "IgG3 gp120+gp140 breadth (=s)", "IgG V1V2 breadth (=s)", "Overall max diversity score (=s)", "ADCC Peak CAP8 (=s)", "ADCC Peak CH58 (=s)", "ADCC Peak WITO (=s)", "ADCC  AUC CAP8 (=s)", "ADCC AUC CH58 (=s)", "ADCC AUC WITO (=s)", "CD4+ T-cells IFN-g/IL-2 (=s)", "CD8+ T-cells IFN-g/IL-2 (=s)", "IgG3 AE.A244 V1V2 Tags 293F (=s)", "IgG3 C.1086C V1V2 Tags (=s)", "IgG3 gp70-001428.2.42 V1V2 (=s)", "IgG3 gp70-1012.11.TC21.3257 V1V2 (=s)", "IgG3 gp70-1394C9G1 V1V2 (=s)", "IgG3 gp70-BF1266 431a V1V2 (=s)", "IgG3 gp70-Ce1086 B2 V1V2 (=s)", "IgG3 gp70-B.CaseAV1V2 (=s)", "IgG AE.A244 V1V2 Tags 293F (=s)", "IgG C.1086C V1V2 Tags (=s)", "IgG gp70-001428.2.42 V1V2 (=s)", "IgG gp70-1012.11.TC21.3257 V1V2 (=s)", "IgG gp70-1394C9G1 V1V2 (=s)", "IgG gp70-9004SS.A3.4 V1V2 (=s)", "IgG gp70-BF1266.431a V1V2 (=s)", "IgG gp70-Ce1086.B2 V1V2 (=s)", "IgG gp70.B.CaseA V1V2 (=s)")
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_excludeELISpotmarkers.csv"
    cfg2$folder_local <- "HVTN 705 (all) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(3, 39),
      marker = c(1:39),
      lab_x = c(1:39),
      lab_title = c(1:39),
      t_0 = rep(1, 39),
      dataset = rep(1, 39),
      cr2_trial = rep(1, 39),
      cr2_COR = rep(1, 39),
      cr2_marker = c(1:39),
      edge_corr = c(1,1,1,1,1,1,1,1,2,1,1,1,1,1,2,2,2,2,2,2,
                    2,2,1,1,2,2,2,2,2,2,1,1,2,2,2,2,2,1,2),
      v_id = rep(1, 39),
      v_time = rep(1, 39),
      v_event = rep(1, 39),
      v_wt = rep(1, 39),
      v_ph1 = rep(1, 39),
      v_ph2 = rep(1, 39),
      v_covariates = rep(1, 39),
      dir = rep(1, 39),
      zoom_x = rep(1, 39),
      zoom_y_cve = rep(1, 39),
      zoom_y_risk = rep(1, 39),
      more_ticks = rep(1, 39)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (second)") {
    
    # Override default config
    cfg2$params$deriv_type <- "line"
    cfg2$density_type <- "kde edge"
    
    # Analysis-specific config
    cfg2$marker <- c("Day210ELCZ", "Day210ADCPgp140C97ZAfib", "Day210IgG340mdw_V1V2", "Day210IgG340mdw_gp120_gp140_vm", "Day210ELISpotPTEEnv", "Day210mdw_xassay_overall", "Day210ADCPgp140Mos1fib", "Day210IgG50mdw_V1V2", "Day210ADCCCAP8_pAUC", "Day210ADCCCH58_pAUC", "Day210ADCCWITO_pAUC", "Day210ICS4AnyEnvIFNg_OR_IL2", "Day210ICS8AnyEnvIFNg_OR_IL2", "Day210IgG3AE.A244.V1V2.Tags_293F40delta", "Day210IgG3C.1086C.V1.V2.Tags40delta", "Day210IgG3gp70.001428.2.42.V1V240delta", "Day210IgG3gp70.1012.11.TC21.3257.V1V240delta", "Day210IgG3gp70.1394C9G1.V1V240delta", "Day210IgG3gp70.BF1266.431a.V1V240delta", "Day210IgG3gp70.Ce1086.B2.V1V240delta", "Day210IgG3gp70.B.CaseA.V1.V240delta", "Day210mdw_xassay_select_igg3v1v2")
    cfg2$lab_title <- c("IgG gp140 C97ZA: Month 7", "ADCP gp140 C97ZA: Month 7", "IgG3 V1V2 breadth: Month 7", "IgG3 gp120+gp140 breadth: Month 7", "ELISPot PTE Env: Month 7", "Multi-epitope functions: Month 7", "ADCP gp140 Mos1: Month 7", "IgG V1V2 breadth: Month 7", "ADCC AUC CAP8: Month 7", "ADCC AUC CH58: Month 7", "ADCC AUC WITO: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2: Month 7", "IgG3 AE.A244 V1V2 Tags 293F: Month 7", "IgG3 C.1086C V1V2 Tags: Month 7", "IgG3 gp70-001428.2.42 V1V2: Month 7", "IgG3 gp70-1012.11.TC21.3257 V1V2: Month 7", "IgG3 gp70-1394C9G1 V1V2: Month 7", "IgG3 gp70-BF1266 431a V1V2: Month 7", "IgG3 gp70-Ce1086 B2 V1V2: Month 7", "IgG3 gp70-B.CaseAV1V2: Month 7", "IgG3 V1V2 A244/1086/CaseA: Month 7")
    cfg2$lab_x <- c("IgG gp140 C97ZA (EU/ml) (=s)", "Average phagocytosis score to gp140 C97ZA (=s)", "IgG3 V1V2 breadth (Wt avg log10 Net MFI) (=s)", "IgG3 gp120 + gp140 breadth (Wt avg log10 Net MFI) (=s)", "ELISPot PTE Env (=s)", "Multi-epitope functions (=s)", "Average phagocytosis score to gp140 Mos1 (=s)", "IgG V1V2 breadth (Wt avg log10 Net MFI) (=s)", "AUC baseline-subtracted CAP8 (% loss of luc activity) (=s)", "AUC baseline-subtracted CH58 (% loss of luc activity) (=s)", "AUC baseline-subtracted WITO (% loss of luc activity) (=s)", "CD4+ T cell responses to any Env peptide pools (=s)", "CD8+ T cell responses to any Env peptide pools (=s)", "IgG3 Net MFI to AE.A244 V1V2 Tags 293F (=s)", "IgG3 Net MFI to C.1086C V1V2 Tags (=s)", "IgG3 Net MFI to gp70-001428.2.42 V1V2 (=s)", "IgG3 Net MFI to gp70-1012.11.TC21.3257 V1V2 (=s)", "IgG3 Net MFI to gp70-1394C9G1 V1V2 (=s)", "IgG3 Net MFI to gp70-BF1266 431a V1V2 (=s)", "IgG3 Net MFI to gp70-Ce1086 B2 V1V2 (=s)", "IgG3 Net MFI to gp70-B.CaseAV1V2 (=s)", "IgG3 V1V2 breadth (AE.A244/C.1086/B.CaseA) (=s)")
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_v9.csv"
    cfg2$folder_local <- "HVTN 705 (second) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(3, 22),
      marker = c(1:22),
      lab_x = c(1:22),
      lab_title = c(1:22),
      t_0 = rep(1, 22),
      dataset = rep(1, 22),
      cr2_trial = rep(1, 22),
      cr2_COR = rep(1, 22),
      cr2_marker = c(1:22),
      edge_corr = c(1,1,2,1,2,1,1,1,2,2,2,2,2,1,1,2,2,2,2,2,2,1),
      v_id = rep(1, 22),
      v_time = rep(1, 22),
      v_event = rep(1, 22),
      v_wt = rep(1, 22),
      v_ph1 = rep(1, 22),
      v_ph2 = rep(1, 22),
      v_covariates = rep(1, 22),
      dir = rep(1, 22),
      zoom_x = rep(1, 22),
      zoom_y_cve = rep(1, 22),
      zoom_y_risk = rep(1, 22),
      more_ticks = rep(1, 22),
      llox_label = rep(1, 22),
      llox = rep(1, 22)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (compare RV144)") {
    
    # Override default config
    cfg2$params$deriv_type <- "line"
    cfg2$density_type <- "kde edge"
    cfg2$zoom_y_cve <- list(c(-1,1.05))
    cfg2$zoom_y_risk_max <- 0.1
    
    # Analysis-specific config
    cfg2$marker <- c("Day210IgG3AE.A244.V1V2.Tags_293F40delta", "Day210IgG3C.1086C.V1.V2.Tags40delta", "Day210IgG3gp70.Ce1086.B2.V1V240delta", "Day210IgG3gp70.B.CaseA.V1.V240delta", "Day210IgAA1.con.env03.140.CF10delta", "Day210mdw_xassay_select_igg3v1v2")
    cfg2$lab_title <- c("IgG3 Net MFI to AE.A244 V1V2 Tags 293F: Week 30", "IgG3 Net MFI to C.1086C V1V2 Tags: Week 30", "IgG3 Net MFI to gp70-Ce1086 B2 V1V2: Week 30", "IgG3 Net MFI to gp70-B.CaseAV1V2: Week 30", "IgA A1.con.env03 140 CF: Week 30", "IgG3 V1V2 breadth (AE.A244, C.1086, B.CaseA): Week 30")
    cfg2$lab_x <- c("IgG3 AE.A244 V1V2 Tags 293F (=s)", "IgG3 C.1086C V1V2 Tags (=s)", "IgG3 gp70-Ce1086 B2 V1V2 (=s)", "IgG3 gp70-B.CaseAV1V2 (=s)", "IgA A1.con.env03 140 CF (=s)", "IgG3 V1V2 A244, 1086, CaseA (=s)")
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_v9.csv"
    cfg2$folder_local <- "HVTN 705 (compare RV144) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(3, 6),
      marker = c(1:6),
      lab_x = c(1:6),
      lab_title = c(1:6),
      t_0 = rep(1, 6),
      dataset = rep(1, 6),
      cr2_trial = rep(1, 6),
      cr2_COR = rep(1, 6),
      cr2_marker = c(1:6),
      edge_corr = c(1,1,2,2,2,1),
      v_id = rep(1, 6),
      v_time = rep(1, 6),
      v_event = rep(1, 6),
      v_wt = rep(1, 6),
      v_ph1 = rep(1, 6),
      v_ph2 = rep(1, 6),
      v_covariates = rep(1, 6),
      dir = c(1,1,1,1,2,1),
      zoom_x = rep(1, 6),
      zoom_y_cve = rep(1, 6),
      zoom_y_risk = rep(1, 6),
      more_ticks = rep(1, 6),
      llox_label = rep(1, 6),
      llox = rep(1, 6)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (ICS)") {
    
    # Override default config
    cfg2$params$deriv_type <- "line"

    # Analysis-specific config
    cfg2$marker <- c("Day210ICS4JMos1gp120IFNg_OR_IL2", "Day210ICS4JMos1gp41IFNg_OR_IL2", "Day210ICS4JMos2GagIFNg_OR_IL2", "Day210ICS4JMos2RNAseIntIFNg_OR_IL2", "Day210ICS4JMos2Sgp120IFNg_OR_IL2", "Day210ICS4JMos2Sgp41IFNg_OR_IL2", "Day210ICS8JMos1gp120IFNg_OR_IL2", "Day210ICS8JMos1gp41IFNg_OR_IL2", "Day210ICS8JMos2GagIFNg_OR_IL2", "Day210ICS8JMos2RNAseIntIFNg_OR_IL2", "Day210ICS8JMos2Sgp120IFNg_OR_IL2", "Day210ICS8JMos2Sgp41IFNg_OR_IL2")
    cfg2$lab_title <- c("Pct CD4+ T-cells expressing IFN-g/IL-2 JMos1 gp120: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos1 gp41: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2 Gag: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2 RNAseInt: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2s gp120: Month 7", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2s gp41: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos1 gp120: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos1 gp41: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2 Gag: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2 RNAseInt: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2s gp120: Month 7", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2s gp41: Month 7")
    cfg2$lab_x <- c("CD4+ T-cells IFN-g/IL-2 JMos1 gp120 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos1 gp41 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2 Gag (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2 RNAseInt (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2s gp120 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2s gp41 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos1 gp120 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos1 gp41 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2 Gag (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2 RNAseInt (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2s gp120 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2s gp41 (=s)")
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_excludeELISpotmarkersaddICSmarkers.csv"
    cfg2$folder_local <- "HVTN 705 (ICS) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(3, 12),
      marker = c(1:12),
      lab_x = c(1:12),
      lab_title = c(1:12),
      t_0 = rep(1, 12),
      dataset = rep(1, 12),
      cr2_trial = rep(1, 12),
      cr2_COR = rep(1, 12),
      cr2_marker = c(1:12),
      edge_corr = c(2,2,2,2,2,2,2,2,2,2,2,2),
      v_id = rep(1, 12),
      v_time = rep(1, 12),
      v_event = rep(1, 12),
      v_wt = rep(1, 12),
      v_ph1 = rep(1, 12),
      v_ph2 = rep(1, 12),
      v_covariates = rep(1, 12),
      dir = rep(1, 12),
      zoom_x = rep(1, 12),
      zoom_y_cve = rep(1, 12),
      zoom_y_risk = rep(1, 12),
      more_ticks = rep(1, 12)
    )
    
  }
  
  if (cfg2$analysis=="AZD1222") {
    
    # Override default config
    cfg2$qnt[["Risk, nonparametric"]] <- c(0.1,0.9)
    cfg2$qnt[["CVE, nonparametric"]] <- c(0.1,0.9)
    cfg2$params$deriv_type <- "line"
    cfg2$params$ci_type <- "regular"
    cfg2$zoom_x <- NA
    cfg2$zoom_y_cve <- NA
    
    # Analysis-specific config
    cfg2$marker <- c("Day29pseudoneutid50", "Day57pseudoneutid50", "Day29bindSpike", "Day57bindSpike")
    cfg2$lab_title <- c("PsV Neutralization 50% Titer: Day 29", "PsV Neutralization 50% Titer: Day 57", "Binding Antibody to Spike: Day 29", "Binding Antibody to Spike: Day 57")
    cfg2$lab_x <- c("Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Anti Spike IgG (BAU/ml) (=s)")
    cfg2$t_0 <- c(117,92)
    cfg2$dataset <- c("azd1222_data_processed_with_riskscore.csv", "azd1222_bAb_data_processed_with_riskscore.csv")
    cfg2$folder_local <- "AZD1222 data/"
    cfg2$folder_cluster <- "Z:/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$cr2_trial <- c("azd1222", "azd1222_bAb")
    cfg2$cr2_COR <- c("D29", "D57")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD29", "EventTimePrimaryD57"),
      event = c("EventIndPrimaryD29", "EventIndPrimaryD57"),
      wt = c("wt.D29", "wt.D57"),
      ph1 = c("ph1.D29", "ph1.D57"),
      ph2 = c("ph2.D29", "ph2.D57"),
      covariates = "~. + Age + risk_score"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(1, 4),
      marker = c(1,2,3,4),
      lab_x = c(1,1,2,2),
      lab_title = c(1,2,3,4),
      t_0 = c(1,2,1,2),
      dataset = c(1,1,2,2),
      cr2_trial = c(1,1,2,2),
      cr2_COR = c(1,2,1,2),
      cr2_marker = rep(1, 4),
      edge_corr = c(2,1,1,1),
      v_id = rep(1, 4),
      v_time = c(1,2,1,2),
      v_event = c(1,2,1,2),
      v_wt = c(1,2,1,2),
      v_ph1 = c(1,2,1,2),
      v_ph2 = c(1,2,1,2),
      v_covariates = rep(1, 4),
      dir = rep(1, 4),
      zoom_x = rep(1, 4),
      zoom_y_cve = rep(1, 4),
      zoom_y_risk = rep(1, 4),
      more_ticks = rep(1, 4)
    )
    
  }
  
  if (cfg2$analysis=="Janssen (partA)") {
    
    # Initial analysis: 1-58
    # Second manuscript: 1-3,13-15,25-27,37-39,45-47,49-51,59-64
    
    # Override default config
    cfg2$qnt[["Risk, nonparametric"]] <- c(0,0.9)
    cfg2$qnt[["CVE, nonparametric"]] <- c(0,0.9)
    cfg2$params$g_n_type <- "parametric (edge) 2"
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk_max <- 0.15
    cfg2$more_ticks <- c(1,2)
    cfg2$llox_label <- c("LOD", "LLOQ")
    cfg2$llox <- c(NA, 4.8975)
    
    # Analysis-specific config
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29ADCP", "Day29pseudoneutid50la", "Day29pseudoneutid50sa")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to RBD: Day 29", "PsV Neutralization 50% Titer: Day 29", "Phagocytic Score: Day 29", "PsV Neutralization 50% Titer (LA): Day 29", "PsV Neutralization 50% Titer (SA): Day 29")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Phagocytic Score (=s)", "Pseudovirus-nAb ID50 LA (IU50/ml) (=s)", "Pseudovirus-nAb ID50 SA (IU50/ml) (=s)")
    cfg2$t_0 <- 0
    cfg2$dataset <- c("janssen_pooled_partA_data_processed_with_riskscore.csv",
                      "janssen_pooled_partAsenior_data_processed_with_riskscore.csv",
                      "janssen_pooled_partAnonsenior_data_processed_with_riskscore.csv",
                      "janssen_na_partA_data_processed_with_riskscore.csv",
                      "janssen_na_partAsenior_data_processed_with_riskscore.csv",
                      "janssen_na_partAnonsenior_data_processed_with_riskscore.csv",
                      "janssen_la_partA_data_processed_with_riskscore.csv",
                      "janssen_la_partAsenior_data_processed_with_riskscore.csv",
                      "janssen_la_partAnonsenior_data_processed_with_riskscore.csv",
                      "janssen_sa_partA_data_processed_with_riskscore.csv",
                      "janssen_sa_partAnonsenior_data_processed_with_riskscore.csv")
    cfg2$folder_local <- "Janssen (partA) data/"
    cfg2$folder_cluster <- "Z:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$cr2_trial <- c("janssen_pooled_partA",
                        "janssen_pooled_partAsenior",
                        "janssen_pooled_partAnonsenior",
                        "janssen_na_partA",
                        "janssen_na_partAsenior",
                        "janssen_na_partAnonsenior",
                        "janssen_la_partA",
                        "janssen_la_partAsenior",
                        "janssen_la_partAnonsenior",
                        "janssen_sa_partA",
                        "janssen_sa_partAnonsenior")
    cfg2$cr2_COR <- c("D29IncludeNotMolecConfirmed",
                      "D29SevereIncludeNotMolecConfirmed",
                      "D29ModerateIncludeNotMolecConfirmed")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryIncludeNotMolecConfirmedD29",
               "SevereEventTimePrimaryIncludeNotMolecConfirmedD29",
               "ModerateEventTimePrimaryIncludeNotMolecConfirmedD29"),
      event = c("EventIndPrimaryIncludeNotMolecConfirmedD29",
                "SevereEventIndPrimaryIncludeNotMolecConfirmedD29",
                "ModerateEventIndPrimaryIncludeNotMolecConfirmedD29"),
      wt = "wt.D29",
      ph1 = "ph1.D29",
      ph2 = "ph2.D29",
      covariates = c("~. + risk_score + as.factor(Region)", "~. + risk_score")
    )
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(1, 64), # C (same)
      marker = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      lab_x = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      lab_title = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      t_0 = rep(1, 64), # C (same)
      dataset = c(rep(c(1:11), each=4), rep(1,4), rep(7,4), c(7:11), 7, c(1,1,1,7,7,7)), # D (region X age)
      cr2_trial = c(rep(c(1:11), each=4), rep(1,4), rep(7,4), c(7:11), 7, c(1,1,1,7,7,7)), # D (region X age)
      cr2_COR = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      cr2_marker = c(rep(c(1:4),6), rep(c(1,2,3,5),5), c(1:4), c(1,2,3,5), rep(4,6), c(1,2,3,1,2,3)),
      edge_corr = rep(2, 64), # C (same)
      v_id = rep(1, 64), # C (same)
      v_time = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      v_event = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      v_wt = rep(1, 64), # C (same)
      v_ph1 = rep(1, 64), # C (same)
      v_ph2 = rep(1, 64), # C (same)
      v_covariates = c(rep(1,11), rep(2,33), rep(1,4), rep(2,4), rep(2,6), c(1,1,1,2,2,2)), # B (pooled or not)
      dir = rep(1, 64), # C (same)
      zoom_x = rep(1, 64), # C (same)
      zoom_y_cve = rep(1, 64), # C (same)
      zoom_y_risk = rep(1, 64), # C (same)
      # more_ticks = replace(rep(1,64),c(3,47,61),2), # G (custom)
      more_ticks = rep(2, 64), # G (custom)
      llox_label = c(rep(c(1,1,2,2), 13), rep(2,6), c(1,1,2,1,1,2)), # F (llox)
      llox = c(rep(c(1,1,2,2), 13), rep(2,6), c(1,1,2,1,1,2)) # F (llox)
    )
    
  }
  
  if (cfg2$analysis=="Profiscov") {
    
    # Override default config
    # cfg2$params$g_n_type <- "parametric" # Historical
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk_max <- 0.21
    cfg2$density_type <- "kde edge"
    cfg2$llox_label <- c("LLOQ", "LOD")
    cfg2$llox <- c(49*0.009, 70*0.009, 72*0.009, 32*0.009, 35*0.0272, 224*0.0272, 53*0.0272, 91*0.0272, 46*0.00236, 27.56)
    
    # Analysis-specific config
    cfg2$marker <- c("Day43bindSpike", "Day43bindSpike_B.1.1.7", "Day43bindSpike_B.1.351", "Day43bindSpike_P.1", "Day43bindRBD", "Day43bindRBD_B.1.1.7", "Day43bindRBD_B.1.351", "Day43bindRBD_P.1", "Day43bindN",
                     "Day91bindSpike", "Day91bindSpike_B.1.1.7", "Day91bindSpike_B.1.351", "Day91bindSpike_P.1", "Day91bindRBD", "Day91bindRBD_B.1.1.7", "Day91bindRBD_B.1.351", "Day91bindRBD_P.1", "Day91bindN",
                     "Day43liveneutmn50")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 43", "Binding Antibody to Spike B.1.1.7: Day 43", "Binding Antibody to Spike B.1.351: Day 43", "Binding Antibody to Spike P.1: Day 43", "Binding Antibody to RBD: Day 43", "Binding Antibody to RBD B.1.1.7: Day 43", "Binding Antibody to RBD B.1.351: Day 43", "Binding Antibody to RBD P.1: Day 43", "Binding Antibody to Nucleocapsid: Day 43",
                        "Binding Antibody to Spike: Day 91", "Binding Antibody to Spike B.1.1.7: Day 91", "Binding Antibody to Spike B.1.351: Day 91", "Binding Antibody to Spike P.1: Day 91", "Binding Antibody to RBD: Day 91", "Binding Antibody to RBD B.1.1.7: Day 91", "Binding Antibody to RBD B.1.351: Day 91", "Binding Antibody to RBD P.1: Day 91", "Binding Antibody to Nucleocapsid: Day 91",
                        "Live Virus Micro Neut 50% Titer: Day 43")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti Spike B.1.1.7 IgG (BAU/ml) (=s)", "Anti Spike B.1.351 IgG (BAU/ml) (=s)", "Anti Spike P.1 IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Anti RBD B.1.1.7 IgG (BAU/ml) (=s)", "Anti RBD B.1.351 IgG (BAU/ml) (=s)", "Anti RBD P.1 IgG (BAU/ml) (=s)", "Anti N IgG (BAU/ml) (=s)", "Live Virus-mnAb ID50 (IU50/ml) (=s)")
    cfg2$t_0 <- c(114,66)
    cfg2$dataset <- c("profiscov_data_processed_with_riskscore.csv",
                      "profiscov_lvmn_data_processed_with_riskscore.csv")
    cfg2$folder_local <- "Profiscov data/"
    cfg2$folder_cluster <- "Y:/cavd/Objective 4/GH-VAP/ID127-Gast/correlates/adata/"
    cfg2$cr2_trial <- c("profiscov", "profiscov_lvmn") # dummy; not using Cox estimates
    cfg2$cr2_COR <- c("D43", "D91")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD43", "EventTimePrimaryD91"),
      event = c("EventIndPrimaryD43", "EventIndPrimaryD91"),
      wt = c("wt.D43", "wt.D91"),
      ph1 = c("ph1.D43", "ph1.D91"),
      ph2 = c("ph2.D43", "ph2.D91"),
      covariates = "~.+ HighRiskInd + Sex + Age + BMI"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(1, 19), # F (constant)
      marker = c(1:19), # A (same)
      lab_x = c(rep(c(1:9),2),10), # B (marker type)
      lab_title = c(1:19), # A (same)
      t_0 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      dataset = c(rep(1, 18), 2), # D (regular vs. LVNT)
      cr2_trial = c(rep(1, 18), 2), # D (regular vs. LVNT)
      cr2_COR = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      cr2_marker = c(rep(c(1:9),2),1), # E (marker order cr2)
      edge_corr = rep(1, 19), # F (constant)
      v_id = rep(1, 19), # F (constant)
      v_time = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_event = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_wt = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_ph1 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_ph2 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_covariates = rep(1, 19), # F (constant)
      dir = rep(1, 19), # F (constant)
      zoom_x = rep(1, 19), # F (constant)
      zoom_y_cve = rep(1, 19), # F (constant)
      zoom_y_risk = rep(1, 19), # F (constant)
      more_ticks = c(2,2,2,2,2,2,2,2,1,2,
                     2,2,2,2,2,2,2,1,2), # F (constant)
      llox_label = c(rep(1, 18), 2), # D (regular vs. LVNT)
      llox = c(rep(c(1:9),2),10) # B (marker type)
    )
    
  }
  
  if (cfg2$analysis=="RV144") {
    
    # Override default config
    cfg2$zoom_y_cve <- list(c(-1,1.05))
    cfg2$zoom_y_risk_max <- 0.05
    cfg2$density_type <- "kde edge"
    
    # Analysis-specific config
    cfg2$marker <- c("Day182AEA244V1V2Tags293F", "Day182C1086C_V1_V2Tags", "Day182gp70_C1086CV1V2293F", "Day182gp70_BCaseA_V1_V2", "Day182iga_A1conenv03140CF", "mdw.select.igg3.v1v2")
    cfg2$lab_title <- c("IgG3 Net MFI to AE.A244 V1V2 Tags 293F: Week 26", "IgG3 Net MFI to C.1086C V1V2 Tags: Week 26", "IgG3 Net MFI to gp70-Ce1086 B2 V1V2: Week 26", "IgG3 Net MFI to gp70-B.CaseAV1V2: Week 26", "IgA A1.con.env03 140 CF: Week 26", "IgG3 V1V2 breadth (AE.A244, C.1086, B.CaseA): Week 26")
    cfg2$lab_x <- c("IgG3 AE.A244 V1V2 Tags 293F (=s)", "IgG3 C.1086C V1V2 Tags (=s)", "IgG3 gp70-Ce1086 B2 V1V2 (=s)", "IgG3 gp70-B.CaseAV1V2 (=s)", "IgA A1.con.env03 140 CF (=s)", "IgG3 V1V2 A244, 1086, CaseA (=s)")
    cfg2$t_0 <- c(578) # c(0)
    cfg2$dataset <- c("rv144_ank.csv")
    cfg2$folder_local <- "RV144 data/"
    cfg2$folder_cluster <- "C:/Users/avike/OneDrive/Desktop/Avi/Biostats + Research/Research/Marco Carone/Project - VaxCurve/VaxCurve/R/RV144 data/"
    cfg2$cr2_trial <- ""
    cfg2$cr2_COR <- ""
    cfg2$v <- list(
      id = "pin",
      time = "time_to_infect",
      event = "infect",
      wt = "wt",
      ph1 = "ph1",
      ph2 = "ph2",
      covariates = "~. + dem_sex + as.factor(BRA_risk)"
    )

    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      endpoint = rep(3, 6),
      marker = c(1:6),
      lab_x = c(1:6),
      lab_title = c(1:6),
      t_0 = rep(1, 6),
      dataset = rep(1, 6),
      cr2_trial = rep(1, 6),
      cr2_COR = rep(1, 6),
      cr2_marker = rep(1, 6),
      edge_corr = c(1,1,1,1,2,1),
      v_id = rep(1, 6),
      v_time = rep(1, 6),
      v_event = rep(1, 6),
      v_wt = rep(1, 6),
      v_ph1 = rep(1, 6),
      v_ph2 = rep(1, 6),
      v_covariates = rep(1, 6),
      dir = c(1,1,1,1,2,1),
      zoom_x = rep(1, 6),
      zoom_y_cve = rep(1, 6),
      zoom_y_risk = rep(1, 6),
      more_ticks = rep(1, 6),
      llox_label = rep(1, 6),
      llox = rep(1, 6)
    )
    
  }
  
  # Set config based on local vs. cluster
  if (Sys.getenv("USERDOMAIN")=="AVI-KENNY-T460") {
    cfg2$tid <- 1
    cfg2$dataset <- paste0(cfg2$folder_cluster,cfg2$dataset)
  } else {
    cfg2$tid <- as.integer(Sys.getenv(.tid_var))
    cfg2$dataset <- paste0(cfg2$folder_local,cfg2$dataset)
  }
  
  # Set config based on cfg2$map
  for (x in c("endpoint", "marker", "lab_x", "lab_title", "day", "dataset",
              "cr2_trial", "cr2_COR", "cr2_marker", "t_0", "dir", "zoom_x",
              "zoom_y_cve", "zoom_y_risk", "more_ticks", "llox_label", "llox",
              "edge_corr")) {
    cfg2[[x]] <- cfg2[[x]][[cfg2$map[cfg2$tid,x]]]
  }
  for (x in c("id", "time", "event", "wt", "ph1", "ph2", "covariates")) {
    cfg2$v[[x]] <- cfg2$v[[x]][[cfg2$map[cfg2$tid,paste0("v_",x)]]]
  }
  cfg2$v$covariates <- formula(cfg2$v$covariates)
  
  # AMP-specific code
  if (cfg2$analysis=="AMP") {
    cfg2[["amp_protocol"]] <-
      cfg2[["amp_protocol"]][[cfg2$map[cfg2$tid,"amp_protocol"]]]
    cfg2[["amp_tx"]] <- cfg2[["amp_tx"]][[cfg2$map[cfg2$tid,"amp_tx"]]]
  }
  
  # Moderna-specific code
  if (cfg2$analysis=="Moderna" && (cfg2$tid %in% c(5,7,9))) {
    cfg2$qnt <- lapply(cfg2$qnt, function(x) { c(0,x[2]) })
  }
  
  # Janssen-specific code
  flags$janssen_id50_lloq <- cfg2$analysis=="Janssen" &&
    cfg2$marker=="Day29pseudoneutid50"
  flags$bsero <- cfg2$analysis %in% c("Janssen (partA)", "Profiscov")
  
  # HVTN705-specific code
  flags$hvtn705_supress <- cfg2$analysis=="HVTN 705 (all)" && cfg2$tid==37
  
}



###########################.
##### Data processing #####
###########################.

{
  # Read in primary data
  df_raw <- read.csv(cfg2$dataset)
  
  # Subset to ph1 cohort
  if (!is.na(cfg2$v$ph1)) {
    df_ph1_full <- dplyr::filter(df_raw, !!rlang::sym(cfg2$v$ph1)==T)
  } else {
    df_ph1_full <- df_raw
  }
  rm(df_raw)
  
  # Subset to filter out baseline seropositive individuals
  if (flags$bsero) { df_ph1_full %<>% dplyr::filter(Bserostatus==0) }
  
  # Subset to filter out individuals without risk scores
  if (!is.null(df_ph1_full$risk_score)) {
    df_ph1_full %<>% dplyr::filter(!is.na(risk_score))
  }
  
  # AMP-specific code
  if (cfg2$analysis=="AMP") {
    if (cfg2$amp_protocol!="Pooled") {
      df_ph1_full %<>% dplyr::filter(protocol==cfg2$amp_protocol)
    }
    if (cfg2$amp_tx=="T1+T2") {
      df_ph1_full %<>% dplyr::filter(tx_pool=="T1+T2")
    } else {
      df_ph1_full %<>% dplyr::filter(tx==cfg2$amp_tx)
    }
  }
  
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
  df_x <- data.frame(tmp=c(1:length(df_ph1_full[[f$vars[1]]])))
  col <- 1
  for (i in c(1:length(f$vars))) {
    if (f$factors[i]==0) {
      df_x[[paste0("x",col)]] <- df_ph1_full[[f$vars[i]]]
      col <- col + 1
    } else {
      x_col <- as.factor(df_ph1_full[[f$vars[i]]])
      levs <- unique(x_col)
      if (length(levs)==1) {
        stop(paste("Covariate", f$vars[i], "has only one unique level"))
      } else {
        for (j in c(1:(length(levs)-1))) {
          df_x[[paste0("x",col)]] <- as.integer(
            df_ph1_full[[f$vars[i]]]==levs[j]
          )
          col <- col + 1
        }
      }
    }
  }
  df_x$tmp <- NULL
  
  # Create phase-two indicator
  if (is.na(cfg2$v$ph2)) {
    df_z <- rep(1, nrow(df_ph1_full))
  } else {
    df_z <- as.integer(df_ph1_full[[cfg2$v$ph2]])
  }
  
  # Create weights variable
  if (is.na(cfg2$v$wt)) {
    df_weights <- rep(1, nrow(df_ph1_full))
  } else {
    df_weights <- df_ph1_full[[cfg2$v$wt]]
  }
  
  # Create data structures to hold results
  plot_data_risk <- data.frame(
    x = double(),
    y = double(),
    curve = character(),
    ci_lo = double(),
    ci_up = double(),
    overall = character()
  )
  plot_data_cve <- plot_data_risk
  
  df_ph1 <- cbind(
    time = df_ph1_full[[cfg2$v$time]],
    event = df_ph1_full[[cfg2$v$event]],
    vacc = df_ph1_full$Trt,
    marker = df_ph1_full[[cfg2$marker]],
    df_x,
    weights = df_weights,
    ph2 = df_z
  )
  
  # Create data object needed by `vaccine` package functions
  dat <- vaccine::load_data(
    time = "time",
    event = "event",
    vacc = "vacc",
    marker = "marker",
    covariates = names(df_x),
    weights = "weights",
    ph2 = "ph2",
    data = df_ph1
  )
  
  # Create cfg2$t_0 variable if it is set to zero
  if (cfg2$t_0==0) {
    if (cfg2$analysis=="Janssen (partA)") {
      SubcohortInd <- df_ph1_full[df_ph1_full$Trt==1,]$SubcohortInd
      if (length(SubcohortInd)!=length(dat$v$z)) {
        stop("Error; lengths differ between SubcohortInd and dat$v.")
      }
      indices_1 <- which(dat$v$z==1 & dat$v$delta==1)
      indices_2 <- which(dat$v$z==1 & SubcohortInd==1)
      time_1 <- max(dat$v$y[indices_1])
      time_2 <- sort(dat$v$y[indices_2], decreasing=T)[15] - 1
      s_num <- sum(dat$v$s==min(dat$v$s, na.rm=T), na.rm=T)
      s_den <- sum(!is.na(dat$v$s))
      cfg2$t_0 <- min(time_1,time_2)
    } else {
      cfg2$t_0 <- max(dat$v$y[dat$v$z==1 & dat$v$delta==1])
    }
  }
  
  # Generate grid of points
  s_grid <- seq(
    from = min(dat$v$s, na.rm=T),
    to = max(dat$v$s, na.rm=T),
    length.out = 101
  )
  
  # Edge mass
  edge_mass <- sum(dat$v$s==min(dat$v$s, na.rm=T), na.rm=T) /
    sum(!is.na(dat$v$s))
  print(paste("Edge mass:", round(edge_mass, 2)))
  
}



#########################.
##### Summary stats #####
#########################.

{
  vaccine::summary_stats(dat=dat)
}



###################################.
##### Overall ests of risk/VE #####
###################################.

if (cfg2$estimators$overall %in% c("Cox gcomp", "KM")) {
  
  # !!!!!
  if (F) {
    
    # Create new covariate
    x1_new <- cut(dat$v$x$x1, breaks=quantile(dat$v$x$x1))
    x1_new_v2 <- cut(dat$v$x$x1, breaks=quantile(dat$v$x$x1, probs=c(0,0.5,1)))
    
    # Unadjusted KM curve (ph1)
    srv_p <- survfit(Surv(dat$v$y,dat$v$delta)~1)
    round(1 - srv_p$surv[which.min(abs(srv_p$time-cfg2$t_0))], 3)
    
    # Unadjusted KM curve (ph2)
    y_ph2 <- dat$v$y[dat$v$z==1]
    delta_ph2 <- dat$v$delta[dat$v$z==1]
    weight_ph2 <- dat$v$weights[dat$v$z==1]
    srv_p2 <- survfit(Surv(y_ph2,delta_ph2)~1, weights=weight_ph2)
    round(1 - srv_p2$surv[which.min(abs(srv_p2$time-cfg2$t_0))],3)
    
    srv_p3 <- survfit(Surv(dat$v$y,dat$v$delta)~factor(x1_new)+dat$v$x$x2+dat$v$x$x3)
    1 - srv_p3$surv[which.min(abs(srv_p3$time-cfg2$t_0))]
    
    model <- coxph(Surv(dat$v$y,dat$v$delta)~factor(x1_new)+dat$v$x$x2+dat$v$x$x3)
    srv_cox <- survfit(model)
    1 - srv_cox$surv[which.min(abs(srv_cox$time-cfg2$t_0))]
    
    
    srv_df <- data.frame(time=srv_p$time, n.risk=srv_p$n.risk,
                         n.event=srv_p$n.event, n.censor=srv_p$n.censor,
                         surv=srv_p$surv, inv_surv=(1-srv_p$surv))
    print(srv_df)
    
    # Plots
    reg_1 <- dat$v$x$x2
    reg_2 <- dat$v$x$x3
    
    plot(survfit(Surv(dat$v$y,dat$v$delta)~1), ylim=c(0.9,1))
    plot(survfit(Surv(dat$v$y,dat$v$delta)~reg_1+reg_2), ylim=c(0.9,1))
    
    plot(survfit(Surv(dat$v$y,dat$v$delta)~x1_new_v2), ylim=c(0.9,1))
    plot(survfit(Surv(dat$v$y,dat$v$delta)~x1_new_v2+reg_1+reg_2), ylim=c(0.9,1))
    
    
    # est <- 1 - srv_p$surv[which.min(abs(srv_p$time-t_0))]
    # ci_lo <- 1 - srv_p$upper[which.min(abs(srv_p$time-t_0))]
    # ci_up <- 1 - srv_p$lower[which.min(abs(srv_p$time-t_0))]
    # se <- srv_p$std.err[which.min(abs(srv_p$time-t_0))]
    
  }
  
  if (cfg2$estimators$overall=="Cox gcomp") {
    method <- "Cox"
  } else if (cfg2$estimators$overall=="KM") {
    method <- "KM"
  } else {
    stop("cfg2$estimators$overall incorrectly specified.")
  }
  
  ests_ov <- vaccine::est_overall(dat=dat, t_0=cfg2$t_0, method=method)
  
  if ("CVE" %in% cfg2$plots) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = c(999,999),
      y = rep(ests_ov[ests_ov$stat=="ve","est"], 2),
      curve = rep("Overall VE", 2),
      ci_lo = rep(ests_ov[ests_ov$stat=="ve","ci_lower"], 2),
      ci_up = rep(ests_ov[ests_ov$stat=="ve","ci_upper"], 2),
      overall = c("Overall L", "Overall R")
    ))
  }
  
  if ("Risk" %in% cfg2$plots) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = rep(999,4),
      y = c(rep(ests_ov[ests_ov$group=="placebo","est"], 2),
            rep(ests_ov[ests_ov$group=="vaccine","est"], 2)),
      curve = c(rep("Placebo overall",2),
                rep("Vaccine overall",2)),
      ci_lo = c(rep(ests_ov[ests_ov$group=="placebo","ci_lower"], 2),
                rep(ests_ov[ests_ov$group=="vaccine","ci_lower"], 2)),
      ci_up = c(rep(ests_ov[ests_ov$group=="placebo","ci_upper"], 2),
                rep(ests_ov[ests_ov$group=="vaccine","ci_upper"], 2)),
      overall = rep(c("Overall L", "Overall R"),2)
    ))
  }
    
}



############################################################.
##### Import functions from correlates_reporting2 repo #####
############################################################.

{
  
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
  draw.x.axis.cor <- function(xlim, llox, more_ticks=1) {
    xlim_scaled <- xlim*more_ticks
    xx <- seq(ceiling(xlim_scaled[1]), floor(xlim_scaled[2])) / more_ticks
    x_axis <- list(ticks=c(), labels=list())
    for (x in xx) {
      if (x>=3) {
        x_axis$ticks[length(x_axis$ticks)+1] <- x
        x_axis$labels[[length(x_axis$labels)+1]] <- bquote(10^.(x))
      } else {
        x_axis$ticks[length(x_axis$ticks)+1] <- log10(signif(10^x,1))
        x_axis$labels[[length(x_axis$labels)+1]] <- signif(10^x,1)
      }
    }
    if (!is.na(llox)) {
      x_axis$ticks[length(x_axis$ticks)+1] <- log10(cfg2$llox)
      x_axis$labels[[length(x_axis$labels)+1]] <- cfg2$llox_label
      # which(abs(x_axis$ticks-cfg2$llox)<0.1) # !!!!! TO DO: suppress label if it overlaps with LOD
    }
    return(x_axis)
  }  
  
}



######################################.
##### Import Cox model estimates #####
######################################.

if (cfg2$estimators$overall=="Cox import" ||
    "Cox import" %in% cfg2$estimators$cr) {
  
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
  
  if (cfg2$estimators$overall=="Cox import" && "CVE" %in% cfg2$plots) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = c(999,999),
      y = rep(overall.ve[[1]], 2),
      curve = rep("Overall VE", 2),
      ci_lo = rep(overall.ve[[2]], 2),
      ci_up = rep(overall.ve[[3]], 2),
      overall = c("Overall L", "Overall R")
    ))
  }
  
  if (cfg2$estimators$cr=="Cox import" && "CVE" %in% cfg2$plots) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(1-risks$prob/res.plac.cont["est"]),
      curve = rep("CVE, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_cve_cis[1,]),
      ci_up = as.numeric(cox_cve_cis[2,]),
      overall = rep("", length(as.numeric(risks$marker)))
    ))
  }
  
  if (cfg2$estimators$overall=="Cox import" && "Risk" %in% cfg2$plots) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = rep(999,4),
      y = c(rep(prev.plac["est"], 2), rep(prev.vacc["est"], 2)),
      curve = c(rep("Placebo overall",2), rep("Vaccine overall",2)),
      ci_lo = c(rep(prev.plac["2.5%"], 2), rep(prev.vacc["2.5%"], 2)),
      ci_up = c(rep(prev.plac["97.5%"], 2), rep(prev.vacc["97.5%"], 2)),
      overall = rep(c("Overall L", "Overall R"),2)
    ))
  }
  
  if (cfg2$estimators$cr=="Cox import" && "Risk" %in% cfg2$plots) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(risks$prob),
      curve = rep("Risk, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_risk_cis[1,]),
      ci_up = as.numeric(cox_risk_cis[2,]),
      overall = rep("", length(as.numeric(risks$marker)))
    ))
  }
  
}



####################################################.
##### Helper function to process CVE estimates #####
####################################################.

{
  
  process_ests <- function(ests, s_grid, run_cve=F, lab_risk=NA, lab_cve=NA) {
    
    # Extract risk estimates and CIs
    ests_risk <- ests$cr$est %>% pmax(0) %>% pmin(1)
    ci_lo_risk <- ests$cr$ci_lower %>% pmax(0) %>% pmin(1)
    ci_up_risk <- ests$cr$ci_upper %>% pmax(0) %>% pmin(1)

    # Compute CVE estimates
    if (run_cve) {
      ests_cve <- ests$cve$est
      ci_lo_cve <- ests$cve$ci_lower %>% pmin(1)
      ci_up_cve <- ests$cve$ci_upper %>% pmin(1)
    }
    
    plot_data_risk <- data.frame(
      x = s_grid,
      y = ests_risk,
      curve = rep(lab_risk, length(ests_risk)),
      ci_lo = ci_lo_risk,
      ci_up = ci_up_risk,
      overall = rep("", length(ests_risk))
    )
    if (run_cve) {
      plot_data_cve <- data.frame(
        x = s_grid,
        y = ests_cve,
        curve = rep(lab_cve, length(ests_cve)),
        ci_lo = ci_lo_cve,
        ci_up = ci_up_cve,
        overall = rep("", length(ests_cve))
      )
    } else {
      plot_data_cve <- NA
    }
    
    return(list(risk=plot_data_risk, cve=plot_data_cve))
    
  }
  
}



#####################################.
##### Data analysis (Mediation) #####
#####################################.

if (flags$run_mediation) {
  
  set.seed(cfg2$seed)
  ests_np_med <- vaccine::est_med(
    dat = dat,
    type = "NP",
    t_0 = cfg2$t_0,
    scale = "VE",
    params_np = vaccine::params_med_np(
      grid_size = list(y=101, s=101, x=5),
      surv_type = cfg2$params$Q_n_type,
      density_type = cfg2$params$g_n_type
      # density_bins = 0 # !!!!!
    )
  )
  saveRDS(ests_np_med, paste0(cfg2$analysis," plots/ests_med_",cfg2$tid,".rds"))
  
  print("MEDIATION ANALYSIS RESULTS")
  print("--------------------------")
  print(ests_np_med)
  # stop("Halted after mediation")
  
  # Process results
  if (F) {
    
    # Choose analysis
    # "Janssen" "Moderna" "AMP" "AZD1222" "Janssen (partA)" "Profiscov"
    # "HVTN 705 (primary)" "HVTN 705 (all)" "RV144" "HVTN 705 (second)"
    # "HVTN 705 (compare RV144)"
    analysis <- "Janssen (partA)"
    if (analysis=="Moderna") {
      inds <- c(5,7,9)
    } else if (analysis=="Janssen (partA)") {
      inds <- c(1:64)
      # inds <- c(1:58)
    } else if (analysis=="HVTN 705 (all)") {
      inds <- c(9,15:22,25:30,33:37,39)
    }
    
    df_med <- data.frame(
      "i"=integer(),
      "nde_est"=double(),"nde_se"=double(),"nde_lo"=double(),"nde_up"=double(),
      "nie_est"=double(),"nie_se"=double(),"nie_lo"=double(),"nie_up"=double(),
      "pm_est"=double(),"pm_se"=double(),"pm_lo"=double(),"pm_up"=double()
    )

    for (i in inds) {
      ests <- readRDS(paste0(analysis," plots/ests_med_",i,".rds"))
      e_nde <- as.list(ests[ests$effect=="NDE",][2:5])
      e_nie <- as.list(ests[ests$effect=="NIE",][2:5])
      e_pm <- as.list(ests[ests$effect=="PM",][2:5])
      e_all <- c(e_nde, e_nie, e_pm)
      e_all <- lapply(e_all, function(x) { round(x,3) })
      df_med[nrow(df_med)+1,] <- c(list(i=i), e_all)

    }
    write.table(
      df_med,
      file = paste0("mediation_results - ",analysis,".csv"),
      sep = ",",
      row.names = F
    )
    
  }
  
}



#################################.
##### Data analysis (NPCVE) #####
#################################.

if ("Grenander" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    set.seed(cfg2$seed)
    ests <- vaccine::est_ce(
      dat = dat,
      type = "NP",
      t_0 = cfg2$t_0,
      cve = T,
      s_out = s_grid,
      ci_type = cfg2$params$ci_type,
      placebo_risk_method = "Cox", # !!!!! Reevaluate this after finishing alternate estimator
      return_extras = T,
      params_np = vaccine::params_ce_np(
        dir = cfg2$dir,
        edge_corr = cfg2$edge_corr,
        grid_size = list(y=101, s=101, x=5),
        surv_type = cfg2$params$Q_n_type,
        # surv_type = "survML-G",
        density_type = cfg2$params$g_n_type,
        # density_bins = 0, # !!!!!
        deriv_type = cfg2$params$deriv_type
      )
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_g_",cfg2$tid,".rds"))
    }
    
    if (flags$save_diagnostics) {
      ggsave(
        filename = paste0(cfg2$analysis," plots/diagnostics_",cfg2$tid,".pdf"),
        plot = diagnostics(ests), device="pdf", width=12, height=8
      )
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_g_",cfg2$tid,".rds"))
    
  }
  
  run_cve <- as.logical("CVE" %in% cfg2$plots)
  ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                        lab_risk="Risk, nonparametric",
                        lab_cve="CVE, nonparametric")
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



#############################################.
##### Data analysis (Cox (spline 3 df)) #####
#############################################.

if ("Cox (spline 3 df)" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    set.seed(cfg2$seed)
    ests <- vaccine::est_ce(
      dat = dat,
      type = "Cox",
      t_0 = cfg2$t_0,
      cve = T,
      s_out = s_grid,
      ci_type = "transformed",
      placebo_risk_method = "Cox",
      params_cox = vaccine::params_ce_cox(spline_df=3)
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_z3_",cfg2$tid,".rds"))
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_z3_",cfg2$tid,".rds"))
    
  }
  
  run_cve <- as.logical("CVE" %in% cfg2$plots)
  ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                        lab_risk="Risk, Cox (spline 3 df)",
                        lab_cve="CVE, Cox (spline 3 df)")
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



#############################################.
##### Data analysis (Cox (spline 4 df)) #####
#############################################.

if ("Cox (spline 4 df)" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    set.seed(cfg2$seed)
    if (cfg2$analysis!="RV144") {
      ests <- vaccine::est_ce(
        dat = dat,
        type = "Cox",
        t_0 = cfg2$t_0,
        cve = T,
        s_out = s_grid,
        ci_type = "transformed",
        placebo_risk_method = "Cox",
        params_cox = vaccine::params_ce_cox(spline_df=4)
      )
    } else {
      ests <- vaccine::est_ce(
        dat = dat,
        type = "Cox",
        t_0 = cfg2$t_0,
        cve = T,
        s_out = s_grid,
        ci_type = "transformed",
        placebo_risk_method = "Cox",
        params_cox = vaccine::params_ce_cox(spline_knots=c(0,5,6.5,8))
      )
    }
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_z4_",cfg2$tid,".rds"))
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_z4_",cfg2$tid,".rds"))
    
  }
  
  run_cve <- as.logical("CVE" %in% cfg2$plots)
  ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                        lab_risk="Risk, Cox (spline 4 df)",
                        lab_cve="CVE, Cox (spline 4 df)")
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



#####################################.
##### Data analysis (Cox gcomp) #####
#####################################.

if ("Cox gcomp" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    set.seed(cfg2$seed)
    ests <- vaccine::est_ce(
      dat = dat,
      type = "Cox",
      t_0 = cfg2$t_0,
      cve = T,
      s_out = s_grid,
      ci_type = "transformed",
      placebo_risk_method = "Cox"
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_c_",cfg2$tid,".rds"))
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_c_",cfg2$tid,".rds"))
    
  }
  
  run_cve <- as.logical("CVE" %in% cfg2$plots)
  ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                        lab_risk="Risk, Cox model", # "Risk, Cox (basic)"
                        lab_cve="CVE, Cox model") # "CVE, Cox (basic)"
  
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



####################################.
##### Data analysis (Cox edge) #####
####################################.

if ("Cox edge" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    set.seed(cfg2$seed)
    ests <- vaccine::est_ce(
      dat = dat,
      type = "Cox",
      t_0 = cfg2$t_0,
      cve = T,
      s_out = s_grid,
      ci_type = "transformed",
      placebo_risk_method = "Cox",
      params_cox = vaccine::params_ce_cox(edge_ind=T)
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_e_",cfg2$tid,".rds"))
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_e_",cfg2$tid,".rds"))
    
  }
  
  run_cve <- as.logical("CVE" %in% cfg2$plots)
  ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                        lab_risk="Risk, Cox (edge)",
                        lab_cve="CVE, Cox (edge)")
  
  plot_data_risk <- rbind(plot_data_risk, ests2$risk)
  if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
  
}



###########################################.
##### Data analysis (Hypothesis test) #####
###########################################.

if (flags$run_hyptest) {
  
  set.seed(cfg2$seed)
  test_results <- test_2(
    dat_orig = dat$v,
    alt_type = "decr",
    # alt_type = "two-tailed",
    params = list(
      type = c("simple (with constant)", "edge", "combined", "combined 2"),
      Q_n_type = "survSL"
    )
  )
  
  if (F) {
    saveRDS(
      test_results,
      paste0(cfg2$analysis," plots/test_results_",cfg2$tid,".rds")
    )
  }
  
  test_results$extras <- NULL
  write.table(
    do.call(rbind, test_results),
    file = paste0(cfg2$analysis," plots/hyptest_",cfg2$tid,".csv"),
    sep = ",",
    row.names = FALSE
  )
  
  # Process hyp test results
  if (F) {
    
    # folder <- "Janssen (partA) plots/Run 11 (58 graphs, 0.90 cutoff)/Hyptest"
    folder <- "Profiscov plots/Run 5 (added hyptest)/Hyptest"
    files <- dir(folder)
    n_pvals <- length(files)
    p_vals <- rep(NA, n_pvals)
    
    for (i in c(1:n_pvals)) {
      
      file <- paste0(folder, "/hyptest_", i, ".csv")
      df <- read.csv(file)
      # df %<>% filter(type=="combined 2")
      p_vals[i] <- df$p_val
      
    }
    
    # Save results
    write.table(data.frame(i=c(1:n_pvals), p_val=round(p_vals,3)),
                file = paste0(folder, "/", "p_vals.csv"),
                sep = ",",
                row.names = F)
    
  }
  
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
  #'   - `ci_up`: Upper confidence bound for CVE
  #' @param which One of c("CVE", "Risk")
  #' @param zoom_x Either a numeric vector of length 2 representing the plot X
  #'     limits, or the string "zoomed", in which case the plot will be zoomed
  #'     into the cutoff quantiles (e.g. 5%/95%), with 10% padding on each side.
  #'     Defaults to the width of the histogram plus 5% padding on each side.
  #' @param zoom_y Either a numeric vector of length 2 representing the plot Y
  #' limits or one of the strings c("zoomed", "zo0med (risk)"). With "zoomed",
  #'     the plot height will be zoomed such that the upper/lower CIs are in the
  #'     frame with 10% padding above and below. With "zoomed (risk)", behavior
  #'     is like "zoomed" but the lower bound is set to zero (plus padding).
  #'     Defaults to 0%--105%.
  #' @param labs A named list of plot labels containing; names include
  #'     c("title", "x", "y")
  #' @param hst A histogram object returned by get.marker.histogram()
  #' @param rr_y_axis Boolean; if true, a secondary risk ratio Y-axis will be
  #'     displayed
  #' @param log10_x_axis Boolean; if true, the X-axis values are log-10
  #' @param case_dots Boolean; if true, dots are plotted on the X-axis
  #'     corresponding to vaccine cases
  #' @param density_type One of c("histogram", "kde", "none")
  #' @return A ggplot2 plot object
  #' @notes
  #'   - This function plots pointwise estimates and confidence intervals from
  #'     both Cox model and nonparametric approaches
  # !!!!! cfg2 and dat currently accessed globally
  create_plot <- function(
    plot_data, which, zoom_x=NA, zoom_y=NA, zoom_y_max=NA, labs, hst,
    rr_y_axis=F, log10_x_axis=F, log10_y_axis=F, case_dots=F,
    density_type="histogram"
  ) {
    
    # Change curve labels to factors and set color scale
    curves <- c(
      "Placebo overall", "Vaccine overall", "Overall VE",
      "Risk, Cox model", "CVE, Cox model", "Risk, Qbins", "CVE, Qbins",
      "Risk, nonparametric", "CVE, nonparametric",
      "Risk, Cox (analytic)", "CVE, Cox (analytic)",
      "Control", "VRC01 10mg/kg", "VRC01 30mg/kg", "VRC01 Pooled",
      "Risk, Cox (basic)", "CVE, Cox (basic)",
      "Risk, Cox (spline 4 df)", "CVE, Cox (spline 4 df)",
      "Risk, Cox (spline 3 df)", "CVE, Cox (spline 3 df)",
      "Risk, Cox (edge)", "CVE, Cox (edge)"
    )
    curve_colors <- c(
      "darkgrey", "darkgrey", "darkgrey",
      "darkorchid3", "darkorchid3", "firebrick3", "firebrick3", # !!!!!
      # "firebrick3", "firebrick3", "firebrick3", "firebrick3", # !!!!!
      "deepskyblue3", "deepskyblue3",
      "deepskyblue3", "deepskyblue3",
      "deepskyblue3", "darkorchid3", "firebrick3", "darkolivegreen3",
      "darkorchid3", "darkorchid3",
      "firebrick3", "firebrick3",
      "firebrick3", "firebrick3",
      "darkgreen", "darkgreen"
    )
    names(curve_colors) <- curves
    indices <- which(curves %in% unique(plot_data$curve))
    curve_colors <- curve_colors[indices]
    plot_data$curve <- factor(plot_data$curve, levels=curves[indices])
    
    # Replace placeholder "Overall" X-values
    plot_data[plot_data$overall=="Overall L","x"] <- min(hst$breaks)
    plot_data[plot_data$overall=="Overall R","x"] <- max(hst$breaks)
    
    # Set default zoom levels
    if (is.na(zoom_x[1])) {
      z_x_L <- min(plot_data$x)
      z_x_R <- max(plot_data$x)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- min(zz, na.rm=T)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed llox") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- log10(cfg2$llox/2)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    }
    # browser() # !!!!!
    if (is.na(zoom_y[1])) {
      zoom_y <- c(0,1.05)
    } else if (zoom_y[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- min(zz$ci_lo, na.rm=T)
      z_y_U <- max(zz$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    } else if (zoom_y[1]=="zoomed (risk)") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- 0
      z_y_U <- max(plot_data$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    }
    if (!is.na(zoom_y_max)) { zoom_y[2] <- min(zoom_y_max, zoom_y[2]) }
    
    # Hack to get zoom_y value for HVTN 124 plots
    if (flags$hvtn124_plot) { zoom_x <<- zoom_x; zoom_y <<- zoom_y; }
    
    # Generate histogram/KDE data
    dens_height <- 0.6 * (zoom_y[2]/1.05-zoom_y[1])
    
    # !!!!! For now, accessing data globally; change
    min_s <- min(dat$v$s, na.rm=T)
    p_edge <- mean(dat$v$s==min_s, na.rm=T) # !!!!! Make this weighted
    
    if (p_edge<0.03 & density_type=="kde edge") { density_type <- "kde" }
    if (density_type=="histogram") {
      
      hist_data <- data.frame(
        xmin = hst$breaks[-length(hst$breaks)],
        xmax = hst$breaks[-1],
        ymin = rep(zoom_y[1], length(hst$counts)),
        ymax = dens_height * (hst$counts/max(hst$counts)) + zoom_y[1]
      )
      
    } else if (density_type=="kde") {
      
      # !!!!! For now, accessing data globally; change
      df_dens <- data.frame(
        s = dat$v$s[!is.na(dat$v$s)],
        weights = dat$v$weights[!is.na(dat$v$s)]
      )
      df_dens$weights <- df_dens$weights / sum(df_dens$weights)
      dens <- stats::density(
        x = df_dens$s,
        bw = "ucv",
        # adjust = 2, # !!!!!
        weights = df_dens$weights
      )
      kde_data <- data.frame(
        x = dens$x,
        ymin = zoom_y[1],
        ymax = dens_height * (dens$y/max(dens$y)) + zoom_y[1]
      )
      
    } else if (density_type=="kde edge") {
      
      # !!!!! For now, accessing data globally; change
      df_dens <- data.frame(
        s = dat$v$s[!is.na(dat$v$s) & dat$v$s!=min_s],
        weights = dat$v$weights[!is.na(dat$v$s) & dat$v$s!=min_s]
      )
      df_dens$weights <- df_dens$weights / sum(df_dens$weights)
      dens <- stats::density(
        x = df_dens$s,
        bw = "ucv",
        # adjust = 2, # !!!!!
        weights = df_dens$weights
      )
      dens$y <- dens$y * (1-p_edge)
      
      plot_width <- zoom_x[2]-zoom_x[1]
      # rect_x <- c(zoom_x[1]+0.01*plot_width, zoom_x[1]+0.06*plot_width)
      rect_x <- c(min_s-0.025*plot_width, min_s+0.025*plot_width)
      rect_y <- p_edge / (rect_x[2]-rect_x[1])
      inds_to_remove <- dens$x>rect_x[2]
      dens$x <- dens$x[inds_to_remove]
      dens$y <- dens$y[inds_to_remove]
      dens$x[length(dens$x)+1] <- rect_x[1]
      dens$y[length(dens$y)+1] <- rect_y
      dens$x[length(dens$x)+1] <- rect_x[2]
      dens$y[length(dens$y)+1] <- rect_y
      dens$x[length(dens$x)+1] <- rect_x[2] + plot_width/10^5
      dens$y[length(dens$y)+1] <- zoom_y[1]
      
      kde_data <- data.frame(
        x = dens$x,
        ymin = zoom_y[1],
        ymax = dens_height * (dens$y/max(dens$y)) + zoom_y[1]
      )
      

    }
    
    # Create and return ggplot2 object
    # Note: using geom_rect for the histogram so that it can be shifted up/down
    if (rr_y_axis) {
      syc_sec.axis <- sec_axis(~1-., breaks=seq(0,2,0.1),
                           name="Risk ratio (vaccine/placebo)")
    } else {
      syc_sec.axis <- waiver()
    }
    y_ticks <- ifelse(which=="CVE", 0.1, 0.01)
    syc_breaks <- seq(-1,1,y_ticks)
    if (log10_y_axis) {
      # zoom_y[2] <- max(plot_data$y, na.rm=T)
      # zoom_y[2] <- min(zoom_y[2], 0.985)
      zoom_y[2] <- 0.99
      plot_data %<>% mutate(
        y = ifelse(y==1, 0.999, y),
        ci_lo = ifelse(ci_lo==1, 0.999, ci_lo),
        ci_up = ifelse(ci_up==1, 0.999, ci_up)
      )
      syc_trans <- scales::trans_new(
        name = "log10_RR",
        transform = function(x) { -log10(1-x) },
        inverse = function(x) { 1 - 10^(-x) }
      )
      # browser() # !!!!! asdf
      hist_data$ymax <- 1 - 10^(-(hist_data$ymax*(-log10(1-zoom_y[2]))))
      # hist_data$ymax <- 1 - 10^(-((hist_data$ymax-hist_data$ymin)/2*(-log10(1-zoom_y[2]))))
      if (which=="CVE") {
        syc_breaks <- c(-1,0,0.5,0.75,0.9,0.95)
        if (rr_y_axis) {
          syc_sec.axis <- sec_axis(~1-., breaks=(1-syc_breaks),
                                   name="Risk ratio (vaccine/placebo)")
        }
      }
    } else {
      syc_trans <- "identity"
    }
    
    # Assemble plot
    plot <- ggplot(plot_data, aes(x=x, y=y, color=curve))
    if (density_type=="histogram") {
      
       plot <- plot + geom_rect(
          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
          hist_data,
          alpha = 0.3,
          fill = "orange", # "forestgreen"
          inherit.aes = F
        )
      
    } else if (density_type %in% c("kde", "kde edge")) {
      
      plot <- plot + geom_ribbon(
        aes(x=x, ymin=ymin, ymax=ymax),
        data = kde_data,
        inherit.aes = F,
        color = "white",
        fill = "orange", # "forestgreen"
        alpha = 0.3
      )
      
    }
    plot <- plot +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_up, fill=curve),
        alpha = 0.05,
        linetype = "dotted"
      ) +
      geom_line(linewidth=0.7) +
      scale_y_continuous(
        labels = scales::label_percent(accuracy=1),
        breaks = syc_breaks,
        minor_breaks = NULL,
        trans = syc_trans,
        sec.axis = syc_sec.axis
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
      if (is.null(cfg2$llox)) {
        x_axis <- draw.x.axis.cor(zoom_x, NA, cfg2$more_ticks)
      } else {
        x_axis <- draw.x.axis.cor(zoom_x, cfg2$llox, cfg2$more_ticks)
      }
      
      if (flags$janssen_id50_lloq) {
        x_axis$ticks[4] <- log10(2.7426)
        x_axis$labels[[4]] <- "LLOQ"
      }
      
      plot <- plot + scale_x_continuous(
        labels = do.call(expression,x_axis$labels),
        breaks = x_axis$ticks
      )
      
    }
    if (which=="Risk") {
      y_plac <- filter(plot_data, curve=="Placebo overall")[1,"y"]
      y_vacc <- filter(plot_data, curve=="Vaccine overall")[1,"y"]
      plot <- plot + annotate("text", label="Placebo overall", x=zoom_x[2],
                              y=y_plac, size=2.5, hjust=1.05, vjust=-0.5)
      plot <- plot + annotate("text", label="Vaccine overall", x=zoom_x[2],
                              y=y_vacc, size=2.5, hjust=1.05, vjust=-0.5)
    }
    
    if (case_dots) {
      s_cases <- dat$v$s[dat$v$delta==1]
      s_cases <- s_cases[!is.na(s_cases)]
      plot <- plot + geom_jitter(
        mapping = aes(x=x,y=y),
        data = data.frame(x=s_cases, y=0.05),
        width = 0.01,
        height = 0.01,
        inherit.aes = F,
        alpha = 0.9,
        color = "darkred",
        size = 1
      )
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
    as.numeric(quantile(dat$v$s, na.rm=T, probs=qnt))
  })
  
  # Trim estimates at specified quantiles
  trim_plot_data <- function(plot_data, cutoffs, cfg2) {
    
    which_curves <- names(cfg2$qnt)
    cut_lo <- sapply(plot_data$curve, function(curve) {
      ifelse(curve %in% which_curves, cutoffs[[curve]][1], NA)
    }, USE.NAMES=F)
    cut_hi <- sapply(plot_data$curve, function(curve) {
      ifelse(curve %in% which_curves, cutoffs[[curve]][2], NA)
    }, USE.NAMES=F)
    rows_1 <- which(plot_data$curve %in% which_curves)
    rows_2 <- which(plot_data$x < cut_lo | plot_data$x > cut_hi)
    rows <- intersect(rows_1, rows_2)
    plot_data[rows, c("y", "ci_lo", "ci_up")] <- NA
    return(plot_data)
  }
  
  hst <- get.marker.histogram(
    marker = dat$v$s[!is.na(dat$v$s)],
    wt = dat$v$weights[!is.na(dat$v$s)],
    trial = cfg2$cr2_trial
  )
  
  if (flags$hvtn705_supress) {
    if (cfg2$marker=="Day210IgGgp70_BF1266.431a.V1V250delta") {
      ind_lo_cve <- min(which(plot_data_cve$curve=="CVE, nonparametric"))
      ind_up_cve <- max(which(plot_data_cve$ci_lo<=-2))
      plot_data_cve[c(ind_lo_cve:ind_up_cve), c("y", "ci_lo", "ci_up")] <- NA
      ind_lo_risk <- min(which(plot_data_risk$curve=="Risk, nonparametric"))
      ind_up_risk <- max(which(plot_data_risk$ci_up>=0.5))
      plot_data_risk[c(ind_lo_risk:ind_up_risk), c("y", "ci_lo", "ci_up")] <- NA
    } else {
      stop("HVTN 705 config changed")
    }
  }
  
  if (flags$paper_npcve || flags$paper_cox) { cfg2$lab_title <- NULL }
  
  if (nrow(plot_data_risk)>0) {
    
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", cfg2$t_0)
    
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_risk, cutoffs, cfg2),
      # plot_data = trim_plot_data(filter(plot_data_risk, curve!="Risk, Cox model"), cutoffs, cfg2),
      which = "Risk",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y_risk,
      zoom_y_max = cfg2$zoom_y_risk_max,
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      log10_x_axis = T,
      density_type = cfg2$density_type
    )
    
    if (cfg2$analysis=="Janssen (partA)") {
      svr <- case_when(
        cfg2$cr2_COR=="D29SevereIncludeNotMolecConfirmed" ~ "_severe",
        cfg2$cr2_COR=="D29ModerateIncludeNotMolecConfirmed" ~ "_moderate",
        TRUE ~ ""
      )
      filename <- paste0("cr_tid",sprintf("%02d", cfg2$tid),"_",cfg2$cr2_trial,
                         svr,"_mrk",cfg2$cr2_marker,".pdf")
    } else {
      filename <- paste0("plot_risk_",cfg2$tid,".pdf")
    }
    
    ggsave(
      filename = paste0(cfg2$analysis," plots/",filename),
      plot=plot, device="pdf", width=6, height=4
    )
    
    if (flags$table_of_vals) {
      write.table(trim_plot_data(plot_data_risk, cutoffs, cfg2),
                  file=paste0(cfg2$analysis," plots/risk_",cfg2$tid,".csv"),
                  sep=",",
                  row.names=FALSE)
    }
    
  }
  
  if (nrow(plot_data_cve)>0) {
    
    cfg2$lab_y <- paste0("Controlled VE against ", cfg2$endpoint,
                         " by day ", cfg2$t_0)
    
    if (flags$hvtn705_abstract_fig) {
      cfg2$lab_title <- "IgG3 V1V2 breadth (Weighted avg log10 Net MFI): Month 7"
      draw.x.axis.cor <- function(xlim, llox) { # llox currently unused
        xx <- seq(ceiling(xlim[1]), floor(xlim[2]))
        x_axis <- list(ticks=c(), labels=list())
        if (is.na(llox)) {
          for (x in xx) {
            label <- 10^x
            x_axis$ticks[length(x_axis$ticks)+1] <- x
            x_axis$labels[[length(x_axis$labels)+1]] <- label
          }
        }
        if (length(xx)<=3) {
          for (i in 2:length(xx)) {
            x=xx[i-1]
            label <- 3*10^x
            x_axis$ticks[length(x_axis$ticks)+1] <- x+log10(3)
            x_axis$labels[[length(x_axis$labels)+1]] <- label
          }
        }
        return(x_axis)
      }
    }
    
    log10_y_axis <- as.logical(flags$partA_mnscrpt2)
    case_dots <- as.logical(flags$partA_mnscrpt2)
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_cve, cutoffs, cfg2),
      # plot_data = trim_plot_data(filter(plot_data_cve, curve!="CVE, Cox model"), cutoffs, cfg2),
      which = "CVE",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y_cve,
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      rr_y_axis = T,
      log10_x_axis = T,
      log10_y_axis = log10_y_axis,
      case_dots = case_dots,
      density_type = cfg2$density_type
    )
    
    if (cfg2$analysis=="Janssen (partA)") {
      svr <- case_when(
        cfg2$cr2_COR=="D29SevereIncludeNotMolecConfirmed" ~ "_severe",
        cfg2$cr2_COR=="D29ModerateIncludeNotMolecConfirmed" ~ "_moderate",
        TRUE ~ ""
      )
      filename <- paste0("cve_tid",sprintf("%02d", cfg2$tid),"_",cfg2$cr2_trial,
                         svr,"_mrk",cfg2$cr2_marker,".pdf")
    } else {
      filename <- paste0("plot_cve_",cfg2$tid,".pdf")
    }
    
    ggsave(
      filename = paste0(cfg2$analysis," plots/",filename),
      plot=plot, device="pdf", width=6, height=4
    )
    
    if (flags$save_plot_objs) {
      saveRDS(plot_data_cve, paste0(cfg2$analysis," plots/plot_data_cve_",
                                    cfg2$tid,".rds"))
      saveRDS(hst, paste0(cfg2$analysis," plots/hst_",cfg2$tid,".rds"))
      saveRDS(cutoffs, paste0(cfg2$analysis," plots/cutoffs_",cfg2$tid,".rds"))
      saveRDS(cfg2, paste0(cfg2$analysis," plots/cfg2_",cfg2$tid,".rds"))
      dat_v <- list(s=dat$v$s, weights=dat$v$weights)
      saveRDS(dat_v, paste0(cfg2$analysis," plots/dat_v_",cfg2$tid,".rds"))
      saveRDS(cfg2, paste0(cfg2$analysis," plots/cfg2_",cfg2$tid,".rds"))
    }
    
    if (flags$table_of_vals) {
      write.table(trim_plot_data(plot_data_cve, cutoffs, cfg2),
                  file=paste0(cfg2$analysis," plots/cve_",cfg2$tid,".csv"),
                  sep=",",
                  row.names=FALSE)
    }
    
  }
  
  # Print timestamp
  print(paste("END:", Sys.time()))
  
}



##########################.
##### HVTN 124 plots #####
##########################.

if (F) {
  
  # Read in data
  dat_124 <- read.csv(paste0("C:/Users/avike/OneDrive/Desktop/Avi/Biostats + R",
                             "esearch/Research/Peter Gilbert/Project - HVTN 12",
                             "4/bama_igg3_data_v124_for_avi_v2.csv"))
  
  # Map marker names
  cfg2$marker2 <- case_when(
    cfg2$marker=="Day210IgG3gp70.001428.2.42.V1V240delta" ~ "gp70-001428.2.42 V1V2",
    cfg2$marker=="Day210IgG3gp70.BF1266.431a.V1V240delta" ~ "gp70-BF1266_431a_V1V2",
    cfg2$marker=="Day210IgG3gp70.Ce1086.B2.V1V240delta" ~ "gp70-Ce1086_B2 V1V2",
    cfg2$marker=="Day210IgG3gp70.B.CaseA.V1.V240delta" ~ "gp70_B.CaseA_V1_V2",
    cfg2$marker=="Day210IgG3AE.A244.V1V2.Tags_293F40delta" ~ "gp70-CM244.ec1 V1V2",
    cfg2$marker=="Day210IgG340mdw_V1V2" ~ "5-antigen V1V2 panel",
    TRUE ~ "ERROR"
  )
  
  # Extract marker values
  dat_124 %<>% filter(antigen==cfg2$marker2 & trt=="T2")
  mrk_vals <- log10(dat_124$delta)

  # Generate histogram
  get.marker.histogram2 <- function(marker) {
    tmp.1 <- hist(marker, breaks=15, plot=F)
    wt <- rep(1, length(marker))
    tmp <- plotrix::weighted.hist(marker, wt, breaks=tmp.1$breaks, plot=F)
    attr(tmp,"class") <- "histogram"
    return(tmp)
  }
  hst2 <- get.marker.histogram2(mrk_vals)
  # zoom_y comes from hack in plotting function using flags$hvtn124_plot==T
  ymax2 <- 0.6 * (hst2$counts/max(hst2$counts)) * (zoom_y[2]/1.05-zoom_y[1]) +
    zoom_y[1]
  hist_data <- data.frame(
    xmin = hst2$breaks[-length(hst2$breaks)],
    xmax = hst2$breaks[-1],
    ymin = rep(zoom_y[1], length(hst2$counts)),
    ymax = ymax2
  )
  if (zoom_x[2]<max(mrk_vals)) {
    zoom_x[2] <- max(mrk_vals) + 0.1*(max(mrk_vals)-zoom_x[1])
  }
  
  plot + geom_rect(
    aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax),
    hist_data,
    color = "brown4",
    size = 0.4,
    fill = "brown4",
    alpha = 0.1,
    inherit.aes = F
  ) + coord_cartesian(xlim=zoom_x, ylim=zoom_y, expand=F)
  
}



#######################################.
##### Super Learner risk analysis #####
#######################################.

if (F) {
  
  # Set library
  SL.library <- c("SL.mean", "SL.gam", "SL.ranger", "SL.earth", "SL.nnet",
                  "SL.svm", "SL.glmnet")
  
  # # Create data objects
  # newX <- distinct(dat$p$x)
  
  set.seed(cfg2$seed)
  model_sl <- SuperLearner(
    Y = dat$p$delta,
    X = dat$p$x,
    # newX = newX,
    family = "binomial",
    SL.library = SL.library,
    verbose = F
  )
  as.numeric(model_sl$SL.predict)
  
  set.seed(cfg2$seed)
  model_slcv <- CV.SuperLearner(
    Y = dat$p$delta,
    X = dat$p$x,
    # newX = newX,
    family = "binomial",
    SL.library = SL.library,
    verbose = F
  )
  as.numeric(model_sl$SL.predict)
  
}



######################################################.
##### Generate and save AMP plots (run manually) #####
######################################################.

if (F) {
  
  # Temp
  cfg2$t_0 <- 595
  cfg2$lab_title <- c("HVTN703/HPTN081", "HVTN704/HPTN085", "Pooled AMP trials")
  cfg2$amp_tx2 <- rep(c("Control", "VRC01 10mg/kg", "VRC01 30mg/kg",
                        "VRC01 Pooled"), 3)
  
  # Generate histograms and KM objects
  # !!!!! Eventually replace this section and instead save the objects above,
  #       accounting for the fact that 13:15 don't fall into the framework
  for (i in c(1:15)) {
    
    # TEMP: Generate dat_amp
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
      dat_amp <- list(
        y = df_analysis[["hiv1survday"]],
        delta = df_analysis[["hiv1event"]],
        weights = rep(1, nrow(df_analysis)),
        s = df_analysis[["bweight"]]
      )
    }
    
    # Generate s_orig
    s_orig <- dat_amp$s[!is.na(dat_amp$s)]
    saveRDS(s_orig, paste0(cfg2$analysis, " plots/s_orig_", i, ".rds"))
    
    if (i %in% c(1,5,9,13,14,15)) {
      
      # Generate histogram
      hst <- get.marker.histogram(
        marker = dat_amp$s[!is.na(dat_amp$s)],
        wt = dat_amp$weights[!is.na(dat_amp$s)],
        trial = F
      )
      saveRDS(hst, paste0(cfg2$analysis, " plots/hist_", i, ".rds"))
      
      # Generate KM object
      srv_ov <- survfit(Surv(dat_amp$y,dat_amp$delta)~1)
      risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-cfg2$t_0))]
      ci_lo_ov <- 1 - srv_ov$upper[which.min(abs(srv_ov$time-cfg2$t_0))]
      ci_up_ov <- 1 - srv_ov$lower[which.min(abs(srv_ov$time-cfg2$t_0))]
      km <- data.frame(
        x = c(999,999),
        y = rep(risk_ov, 2),
        curve = rep("Overall risk", 2),
        ci_lo = rep(ci_lo_ov, 2),
        ci_up = rep(ci_up_ov, 2),
        overall = c("Overall L", "Overall R")
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
      s_orig <- readRDS(paste0(cfg2$analysis, " plots/s_orig_",
                               j, ".rds"))
      s_grid <- seq(from=min(s_orig), to=max(s_orig), length.out=101)
      ests <- readRDS(paste0(cfg2$analysis, " plots/ests_g_",
                             j, ".rds"))
      ests2 <- process_ests(ests, s_grid, run_cve=F, lab_risk=cfg2$amp_tx2[j])
      plot_data_risk <- rbind(plot_data_risk, ests2$risk)
    }
    
    # Generate and save plot
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", cfg2$t_0)
    lab_title <- cfg2$lab_title[plot_map$lab_title[i]]
    plot <- create_plot(
      plot_data = plot_data_risk,
      which = "Risk",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y_risk,
      labs = list(title=lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst
    )
    ggsave(
      filename = paste0(cfg2$analysis," plots/plot_risk_",i,".pdf"),
      plot=plot, device="pdf", width=6, height=4
    )
    
  }
  
}



########################.
##### Violin plots #####
########################.

if (F) {
  
  mrk <- "Day210mdw_xassay_overall"
  mrk_lab <- "Expanded Multi-epitope functions: Month 7"
  title_lab <- "Figure 1-5. HVTN 705 Expanded Multi-epitope functions at Month 7"
  
  df_tx <- df_ph1_full[which(df_ph1_full$Trt==1),]
  df_ct <- df_ph1_full[which(df_ph1_full$Trt==0),]
  df_tx2 <- filter(df_tx, Ph2ptids.D210==1)
  df_ct2 <- filter(df_ct, Ph2ptids.D210==1)
  
  df_plot <- data.frame(
    casect = factor(ifelse(c(df_tx2[,"Delta.D210"], df_ct2[,"Delta.D210"]),
                           "Case", "Control")),
    val = c(df_tx2[,mrk], df_ct2[,mrk]),
    arm = factor(c(rep("Vaccine",nrow(df_tx2)), rep("Placebo",nrow(df_ct2))),
                 levels=c("Vaccine", "Placebo"))
  )
  medians <- round(10^c(
    median(filter(df_plot, arm=="Vaccine" & casect=="Case")$val),
    median(filter(df_plot, arm=="Vaccine" & casect=="Control")$val),
    median(filter(df_plot, arm=="Placebo" & casect=="Case")$val),
    median(filter(df_plot, arm=="Placebo" & casect=="Control")$val)
  ),3)
  
  # Export image: 700 x 400
  ggplot(df_plot, aes(x=val, y=casect, group=casect, color=casect)) +
    geom_violin() +
    geom_boxplot(width=0.3) +
    theme(
      panel.border = element_rect(color="#000000", fill=NA),
      panel.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      strip.background = element_rect(fill="#ffffff"),
      plot.title = element_text(size=10, hjust=0.5)
    ) +
    scale_x_continuous(labels=scales::math_format(10^.x), limits=c(-2,2.5)) +
    facet_wrap(~arm, ncol=2) +
    geom_jitter(height=0.1, width=0) +
    geom_text(
      data = data.frame(
        val = rep(2.2,4),
        lab = paste("Median:", medians),
        casect = c("Case", "Control", "Case", "Control"),
        arm = factor(c(rep("Vaccine",2), rep("Placebo",2)),
                     levels=c("Vaccine", "Placebo"))
      ),
      aes(label=lab),
      color = "#000000",
      size = 3
    ) +
    coord_flip() +
    labs(x=mrk_lab, y="Cohort",
         title=title_lab)
  
}



##################################.
##### Generate RV144 dataset #####
##################################.

if (F) {
  
  # Read in data files
  rv144_ank <- read.csv("RV144 data/Raw data/rv144_master_wk26.csv")
  rv144_ank_markers_1 <- read.csv(paste0("RV144 data/Raw data/rv144_data_wk26_",
                                         "correlates_long.csv"))
  rv144_ank_markers_2 <- read.csv(paste0("RV144 data/Raw data/Tomaras_IgG3_V2.",
                                         "csv"))
  rv144_ank_markers_2 %<>% dplyr::rename("pin"=ptid)
  rv144_ank_markers_3 <- read.csv(paste0("RV144 data/Raw data/rv144_ank_mdw_se",
                                         "lect.csv"))
  
  # Merge/process marker file 1
  rv144_ank_markers_1 %<>% dplyr::filter(
    outcome=="secondary_gt_iga_A1conenv03140CF"
  )
  rv144_ank_markers_1 %<>% subset(select=c(pin, y))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_1, by="pin")
  rv144_ank %<>% dplyr::rename("Day182iga_A1conenv03140CF"=y)
  
  # Merge/process marker file 2
  rv144_ank_markers_2a <- dplyr::filter(
    rv144_ank_markers_2, week==26 & outcome=="AEA244V1V2Tags293F"
  ) %>% subset(select=c(pin,y))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_2a, by="pin")
  rv144_ank %<>% dplyr::rename("Day182AEA244V1V2Tags293F"=y)
  rv144_ank_markers_2b <- dplyr::filter(
    rv144_ank_markers_2, week==26 & outcome=="C1086C_V1_V2Tags"
  ) %>% subset(select=c(pin,y))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_2b, by="pin")
  rv144_ank %<>% dplyr::rename("Day182C1086C_V1_V2Tags"=y)
  rv144_ank_markers_2c <- dplyr::filter(
    rv144_ank_markers_2, week==26 & outcome=="gp70_BCaseA_V1_V2"
  ) %>% subset(select=c(pin,y))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_2c, by="pin")
  rv144_ank %<>% dplyr::rename("Day182gp70_BCaseA_V1_V2"=y)
  rv144_ank_markers_2d <- dplyr::filter(
    rv144_ank_markers_2, week==26 & outcome=="gp70_C1086CV1V2293F"
  ) %>% subset(select=c(pin,y))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_2d, by="pin")
  rv144_ank %<>% dplyr::rename("Day182gp70_C1086CV1V2293F"=y)
  
  # Further data processing
  rv144_ank %<>% dplyr::rename("Trt"=trt)
  rv144_ank %<>% dplyr::filter(ittmod=="Yes" & cc_cohort==1)
  rv144_ank %<>% dplyr::mutate(
    Trt = ifelse(Trt=="VACCINE", 1, 0),
    ph1 = 1,
    ph2 = !is.na(Day182iga_A1conenv03140CF),
    infect = ifelse(infect=="Yes", 1, 0),
    dem_sex = ifelse(dem_sex=="Male", 1, 0),
    stratum = dplyr::case_when(
      Trt==1 & infect==1 ~ 1,
      Trt==1 & infect==0 & dem_sex==0 & vaccno==4 & perprot=="Yes" ~ 2,
      Trt==1 & infect==0 & dem_sex==0 & vaccno==4 & perprot=="No" ~ 3,
      Trt==1 & infect==0 & dem_sex==0 & vaccno==3 & perprot=="No" ~ 4,
      Trt==1 & infect==0 & dem_sex==0 & vaccno==2 & perprot=="No" ~ 5,
      Trt==1 & infect==0 & dem_sex==1 & vaccno==4 & perprot=="Yes" ~ 6,
      Trt==1 & infect==0 & dem_sex==1 & vaccno==4 & perprot=="No" ~ 7,
      TRUE ~ 8
    ),
    wt = dplyr::case_when(
      stratum==1 ~ 1,
      stratum==2 ~ 39.617,
      stratum==3 ~ 45.6,
      stratum==4 ~ 12.4,
      stratum==5 ~ 10.6,
      stratum==6 ~ 30.8,
      stratum==7 ~ 48.3,
      stratum==8 ~ 0,
      TRUE ~ 0
    )
  )
  
  # Filter out individuals in stratum 8
  rv144_ank %<>% dplyr::filter(!(stratum==8 & Trt==1))

  # # !!!!! strata not accounted for: male vaccno=1
  # nrow(dplyr::filter(rv144_ank, stratum==8 & Trt==1))
  # nrow(dplyr::filter(rv144_ank, stratum==8 & Trt==1 & dem_sex==0 & vaccno==1))
  # nrow(dplyr::filter(rv144_ank, stratum==8 & Trt==1 & dem_sex==1 & vaccno==3))
  # nrow(dplyr::filter(rv144_ank, stratum==8 & Trt==1 & dem_sex==1 & vaccno==2))
  # nrow(dplyr::filter(rv144_ank, stratum==8 & Trt==1 & dem_sex==1 & vaccno==1))
  
  # sum weights !!!!!!
  
  # Calculate weights
  # rv144_ank2 <- dplyr::filter(rv144_ank, ph2==1)
  # nrow(dplyr::filter(rv144_ank2, Trt==1 & dem_sex=="Female" & vaccno==3
  #                    & perprot=="Yes" & infect=="Yes"))
  #
  
  # DQA checks
  {
    chk_wt <- function(sex, vnum, pp) {
      d_case_ph1 <- dplyr::filter(rv144_ank, Trt==1 & dem_sex==sex & vaccno==vnum
                                  & perprot==pp & infect==1)
      d_ctrl_ph1 <- dplyr::filter(rv144_ank, Trt==1 & dem_sex==sex & vaccno==vnum
                                  & perprot==pp & infect==0)
      d_case_ph2 <- dplyr::filter(rv144_ank, Trt==1 & dem_sex==sex & vaccno==vnum
                                  & perprot==pp & infect==1 & ph2==T)
      d_ctrl_ph2 <- dplyr::filter(rv144_ank, Trt==1 & dem_sex==sex & vaccno==vnum
                                  & perprot==pp & infect==0 & ph2==T)
      message(paste0("# of cases    (PH-1): ", nrow(d_case_ph1)))
      message(paste0("# of controls (PH-1): ", nrow(d_ctrl_ph1)))
      message(paste0("# of cases    (PH-2): ", nrow(d_case_ph2)))
      message(paste0("# of controls (PH-2): ", nrow(d_ctrl_ph2)))
      message(paste0("IPS weight (cases): ", nrow(d_case_ph1)/nrow(d_case_ph2)))
      message(paste0("IPS weight 2 (cases): ", paste(unique(d_case_ph1$wt), collapse=",")))
      message(paste0("IPS weight (control): ", round(
        nrow(d_ctrl_ph1)/nrow(d_ctrl_ph2), 3
      )))
      message(paste0("IPS weight 2 (control): ", paste(unique(d_ctrl_ph1$wt), collapse=",")))
      invisible("")
    }
    chk_wt(sex=0, vnum=4, pp="Yes")
    chk_wt(sex=0, vnum=4, pp="No")
    chk_wt(sex=0, vnum=3, pp="No")
    chk_wt(sex=0, vnum=2, pp="No")
    chk_wt(sex=1, vnum=4, pp="Yes")
    chk_wt(sex=1, vnum=4, pp="No")
  }
  
  # These two should be equal
  nrow(dplyr::filter(rv144_ank, Trt==1))
  sum(dplyr::filter(rv144_ank, ph2==1)$wt)
  
  # Transform natural log values to log10 values
  cols <- c("Day182iga_A1conenv03140CF", "Day182AEA244V1V2Tags293F",
            "Day182C1086C_V1_V2Tags", "Day182gp70_BCaseA_V1_V2",
            "Day182gp70_C1086CV1V2293F")
  for(col in cols) {
    rv144_ank[[col]] <- log10(exp(rv144_ank[[col]]))
  }
  
  # Pull in breadth score
  rv144_ank_markers_3 %<>% subset(select=c("pin", "mdw.select.igg3.v1v2"))
  rv144_ank %<>% dplyr::left_join(rv144_ank_markers_3, by="pin")
  rv144_ank %<>% dplyr::relocate(mdw.select.igg3.v1v2, .before=ph1)
  
  # Export as csv
  write.table(rv144_ank, file="rv144_ank.csv", sep=",", row.names=F)
  
}



####################################.
##### RV144 vs. HVTN 705 plots #####
####################################.

if (F) {
  
  # New plotting function
  # !!!!! cfg2 and dat currently accessed globally
  create_plot2 <- function(plot_data, which, zoom_x=NA, zoom_y=NA,
                           zoom_y_max=NA, labs, hst_144, hst_705, rr_y_axis=F,
                           log10_x_axis=F, log10_y_axis=F, case_dots=F) {
    
    # Change curve labels to factors and set color scale
    curves <- c("RV 144", "HVTN 705")
    curve_colors <- c("darkgreen", "deepskyblue3")
    
    names(curve_colors) <- curves
    indices <- which(curves %in% unique(plot_data$curve))
    curve_colors <- curve_colors[indices]
    plot_data$curve <- factor(plot_data$curve, levels=curves[indices])
    
    # # Replace placeholder "Overall" X-values
    # plot_data[plot_data$overall=="Overall L","x"] <- min(hst$breaks)
    # plot_data[plot_data$overall=="Overall R","x"] <- max(hst$breaks)
    
    # Set default zoom levels
    if (is.na(zoom_x[1])) {
      z_x_L <- min(plot_data$x)
      z_x_R <- max(plot_data$x)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- min(zz, na.rm=T)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed llox") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- log10(cfg2$llox/2)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    }
    if (is.na(zoom_y[1])) {
      zoom_y <- c(0,1)
      zoom_y[2] <- zoom_y[2] + 0.05*(zoom_y[2]-zoom_y[1])
    } else if (zoom_y[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- min(zz$ci_lo, na.rm=T)
      z_y_U <- max(zz$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    } else if (zoom_y[1]=="zoomed (risk)") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- 0
      z_y_U <- max(plot_data$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    }
    if (!is.na(zoom_y_max)) { z_y_U <- min(zoom_y_max, z_y_U) }
    
    # Generate histogram data
    ymax_144 <- 0.6 * (hst_144$counts/max(hst_144$counts)) *
      (zoom_y[2]/1.05-zoom_y[1]) + zoom_y[1]
    ymax_705 <- 0.6 * (hst_705$counts/max(hst_705$counts)) *
      (zoom_y[2]/1.05-zoom_y[1]) + zoom_y[1]
    hist_data_144 <- data.frame(
      xmin = hst_144$breaks[-length(hst_144$breaks)],
      xmax = hst_144$breaks[-1],
      ymin = rep(zoom_y[1], length(hst_144$counts)),
      ymax = ymax_144, # max(ymax_144,ymax_705),
      trial = "RV 144"
    )
    hist_data_705 <- data.frame(
      xmin = hst_705$breaks[-length(hst_705$breaks)],
      xmax = hst_705$breaks[-1],
      ymin = rep(zoom_y[1], length(hst_705$counts)),
      ymax = ymax_705, # max(ymax_144,ymax_705),
      trial = "HVTN 705"
    )
    hist_data <- rbind(hist_data_144,hist_data_705)
    
    # Hack to get zoom_y value for HVTN 124 plots
    if (flags$hvtn124_plot) {
      zoom_x <<- zoom_x
      zoom_y <<- zoom_y
    }
    
    # Create and return ggplot2 object
    # Note: using geom_rect for the histogram so that it can be shifted up/down
    if (rr_y_axis) {
      syc_sec.axis <- sec_axis(~1-., breaks=seq(0,2,0.1),
                               name="Risk ratio (vaccine/placebo)")
    } else {
      syc_sec.axis <- waiver()
    }
    y_ticks <- ifelse(which=="CVE", 0.1, 0.01)
    syc_breaks <- seq(-1,1,y_ticks)
    if (log10_y_axis) { # !!!!! New section
      # zoom_y[2] <- max(plot_data$y, na.rm=T)
      # zoom_y[2] <- min(zoom_y[2], 0.985)
      zoom_y[2] <- 0.985
      plot_data %<>% mutate(
        y = ifelse(y==1, 0.999, y),
        ci_lo = ifelse(ci_lo==1, 0.999, ci_lo),
        ci_up = ifelse(ci_up==1, 0.999, ci_up)
      )
      syc_trans <- scales::trans_new(
        name = "log10_RR",
        transform = function(x) { -log10(1-x) },
        inverse = function(x) { 1 - 10^(-x) }
      )
      hist_data$ymax <- 1 - 10^(-(hist_data$ymax*(-log10(1-zoom_y[2]))))
      if (which=="CVE") {
        syc_breaks <- c(-1,0,0.5,0.75,0.9,0.95)
        if (rr_y_axis) {
          syc_sec.axis <- sec_axis(~1-., breaks=(1-syc_breaks),
                                   name="Risk ratio (vaccine/placebo)")
        }
      }
    } else {
      syc_trans <- "identity"
    }
    plot <- ggplot(plot_data, aes(x=x, y=y, color=curve)) +
      geom_ribbon(                               # Comment out to get histogram only
        aes(ymin=ci_lo, ymax=ci_up, fill=curve), # Comment out to get histogram only
        alpha = 0.05,                            # Comment out to get histogram only
        linetype = "dotted"                      # Comment out to get histogram only
      ) +                                        # Comment out to get histogram only
      geom_rect(
        aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, color=trial,
            fill=trial),
        hist_data,
        linewidth = 0.3, # asdf1
        alpha = 0.2,
        inherit.aes = F
      ) +
      geom_line(linewidth=0.5) + # 0.7           # Comment out to get histogram only
      scale_y_continuous(
        labels = scales::label_percent(accuracy=1),
        breaks = syc_breaks,
        minor_breaks = NULL,
        trans = syc_trans,
        sec.axis = syc_sec.axis
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
      if (is.null(cfg2$llox)) {
        x_axis <- draw.x.axis.cor(zoom_x, NA, cfg2$more_ticks)
      } else {
        x_axis <- draw.x.axis.cor(zoom_x, cfg2$llox, cfg2$more_ticks)
      }
      
      if (flags$janssen_id50_lloq) {
        x_axis$ticks[4] <- log10(2.7426)
        x_axis$labels[[4]] <- "LLOQ"
      }
      
      plot <- plot + scale_x_continuous(
        labels = do.call(expression,x_axis$labels),
        breaks = x_axis$ticks
      )
      
    }
    if (which=="Risk") {
      y_plac <- filter(plot_data, curve=="Placebo overall")[1,"y"]
      y_vacc <- filter(plot_data, curve=="Vaccine overall")[1,"y"]
      plot <- plot + annotate("text", label="Placebo overall", x=zoom_x[2],
                              y=y_plac, size=2.5, hjust=1.05, vjust=-0.5)
      plot <- plot + annotate("text", label="Vaccine overall", x=zoom_x[2],
                              y=y_vacc, size=2.5, hjust=1.05, vjust=-0.5)
    }
    
    if (case_dots) {
      s_cases <- dat$v$s[dat$v$delta==1]
      s_cases <- s_cases[!is.na(s_cases)]
      plot <- plot + geom_jitter(
        mapping = aes(x=x,y=y),
        data = data.frame(x=s_cases, y=0.05),
        width = 0.01,
        height = 0.01,
        inherit.aes = F,
        alpha = 0.9,
        color = "darkred",
        size = 1
      )
    }
    
    return(plot)
    
  }
  
  # New plotting function
  # !!!!! cfg2 and dat currently accessed globally
  create_plot3 <- function(plot_data, which, zoom_x=NA, zoom_y=NA,
                           zoom_y_max=NA, labs, hst_144, hst_705, rr_y_axis=F,
                           log10_x_axis=F, log10_y_axis=F, case_dots=F) {
    
    # Change curve labels to factors and set color scale
    curves <- c("RV 144", "HVTN 705")
    # curve_colors <- c("darkgreen", "deepskyblue3")
    curve_colors <- c("darkorange", "deepskyblue3")
    
    names(curve_colors) <- curves
    indices <- which(curves %in% unique(plot_data$curve))
    curve_colors <- curve_colors[indices]
    plot_data$curve <- factor(plot_data$curve, levels=curves[indices])
    
    # # Replace placeholder "Overall" X-values
    # plot_data[plot_data$overall=="Overall L","x"] <- min(hst$breaks)
    # plot_data[plot_data$overall=="Overall R","x"] <- max(hst$breaks)
    
    # Set default zoom levels
    if (is.na(zoom_x[1])) {
      z_x_L <- min(plot_data$x)
      z_x_R <- max(plot_data$x)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- min(zz, na.rm=T)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    } else if (zoom_x[1]=="zoomed llox") {
      zz <- dplyr::filter(plot_data, overall=="" & !is.na(y))$x
      z_x_L <- log10(cfg2$llox/2)
      z_x_R <- max(zz, na.rm=T)
      zoom_x <- c(z_x_L - 0.05*(z_x_R-z_x_L),
                  z_x_R + 0.05*(z_x_R-z_x_L))
    }
    if (is.na(zoom_y[1])) {
      zoom_y <- c(0,1)
      zoom_y[2] <- zoom_y[2] + 0.05*(zoom_y[2]-zoom_y[1])
    } else if (zoom_y[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- min(zz$ci_lo, na.rm=T)
      z_y_U <- max(zz$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    } else if (zoom_y[1]=="zoomed (risk)") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- 0
      z_y_U <- max(plot_data$ci_up, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    }
    if (!is.na(zoom_y_max)) { z_y_U <- min(zoom_y_max, z_y_U) }
    
    # Generate histogram/KDE data (RV144)
    max_y <- 0
    dens_height <- 0.6 * (zoom_y[2]/1.05-zoom_y[1])
    for (j in c(144,705)) {
      
      if (j==144) {
        dat_s <- dat_v_144$s
        dat_weights <- dat_v_144$weights
      } else {
        dat_s <- dat_v_705$s
        dat_weights <- dat_v_705$weights
      }
      
      density_type <- "kde edge"
      min_s <- min(dat_s, na.rm=T)
      p_edge <- mean(dat_s==min_s, na.rm=T) # !!!!! Make this weighted
      # print(paste0("P_edge, ", j, ": ",p_edge)) # !!!!!
      # if (p_edge<0.03 & density_type=="kde edge") { density_type <- "kde" } # !!!!!
      if (density_type=="kde") {
        
        df_dens <- data.frame(
          s = dat_s[!is.na(dat_s)],
          weights = dat_weights[!is.na(dat_s)]
        )
        df_dens$weights <- df_dens$weights / sum(df_dens$weights)
        dens <- stats::density(
          x = df_dens$s,
          bw = "ucv",
          # adjust = 2, # !!!!!
          weights = df_dens$weights
        )
        max_y <- max(max_y, dens$y)
        
      } else if (density_type=="kde edge") {
        
        df_dens <- data.frame(
          s = dat_s[!is.na(dat_s) & dat_s!=min_s],
          weights = dat_weights[!is.na(dat_s) & dat_s!=min_s]
        )
        df_dens$weights <- df_dens$weights / sum(df_dens$weights)
        dens <- stats::density(
          x = df_dens$s,
          bw = "ucv",
          # adjust = 2, # !!!!!
          weights = df_dens$weights
        )
        dens$y <- dens$y * (1-p_edge)
        
        plot_width <- zoom_x[2]-zoom_x[1]
        # rect_x <- c(zoom_x[1]+0.01*plot_width, zoom_x[1]+0.06*plot_width)
        rect_x <- c(min_s-0.025*plot_width, min_s+0.025*plot_width)
        rect_y <- p_edge / (rect_x[2]-rect_x[1])
        inds_to_remove <- dens$x>rect_x[2]
        dens$x <- dens$x[inds_to_remove]
        dens$y <- dens$y[inds_to_remove]
        dens$x[length(dens$x)+1] <- rect_x[1]
        dens$y[length(dens$y)+1] <- rect_y
        dens$x[length(dens$x)+1] <- rect_x[2]
        dens$y[length(dens$y)+1] <- rect_y
        dens$x[length(dens$x)+1] <- rect_x[2] + plot_width/10^5
        dens$y[length(dens$y)+1] <- zoom_y[1] # !!!!!
        max_y <- max(max_y, dens$y)
        
        assign(paste0("dens_",j), dens)
        rm(dens)
        
      }
      
    }
    
    for (j in c(144,705)) {
      dens <- get(paste0("dens_",j))
      if (j==144) {
        kde_data <- data.frame(
          x = dens$x,
          ymin = zoom_y[1],
          # ymax = dens_height * (dens$y/max(dens$y)) + zoom_y[1],
          ymax = dens_height * (dens$y/max_y) + zoom_y[1],
          trial = "RV 144"
        )
      } else {
        kde_data <- rbind(kde_data, data.frame(
          x = dens$x,
          ymin = zoom_y[1],
          # ymax = dens_height * (dens$y/max(dens$y)) + zoom_y[1],
          ymax = dens_height * (dens$y/max_y) + zoom_y[1],
          trial = "HVTN 705"
        ))
      }
    }
    
    # Hack to get zoom_y value for HVTN 124 plots
    if (flags$hvtn124_plot) {
      zoom_x <<- zoom_x
      zoom_y <<- zoom_y
    }
    
    # Create and return ggplot2 object
    # Note: using geom_rect for the histogram so that it can be shifted up/down
    if (rr_y_axis) {
      syc_sec.axis <- sec_axis(~1-., breaks=seq(0,2,0.1),
                               name="Risk ratio (vaccine/placebo)")
    } else {
      syc_sec.axis <- waiver()
    }
    y_ticks <- ifelse(which=="CVE", 0.1, 0.01)
    syc_breaks <- seq(-1,1,y_ticks)
    if (log10_y_axis) { # !!!!! New section
      # zoom_y[2] <- max(plot_data$y, na.rm=T)
      # zoom_y[2] <- min(zoom_y[2], 0.985)
      zoom_y[2] <- 0.985
      plot_data %<>% mutate(
        y = ifelse(y==1, 0.999, y),
        ci_lo = ifelse(ci_lo==1, 0.999, ci_lo),
        ci_up = ifelse(ci_up==1, 0.999, ci_up)
      )
      syc_trans <- scales::trans_new(
        name = "log10_RR",
        transform = function(x) { -log10(1-x) },
        inverse = function(x) { 1 - 10^(-x) }
      )
      hist_data$ymax <- 1 - 10^(-(hist_data$ymax*(-log10(1-zoom_y[2]))))
      if (which=="CVE") {
        syc_breaks <- c(-1,0,0.5,0.75,0.9,0.95)
        if (rr_y_axis) {
          syc_sec.axis <- sec_axis(~1-., breaks=(1-syc_breaks),
                                   name="Risk ratio (vaccine/placebo)")
        }
      }
    } else {
      syc_trans <- "identity"
    }
    plot <- ggplot(plot_data, aes(x=x, y=y, color=curve)) +
      geom_ribbon(                               # Comment out to get KDE only
        aes(ymin=ci_lo, ymax=ci_up, fill=curve), # Comment out to get KDE only
        alpha = 0.05,                            # Comment out to get KDE only
        linetype = "dotted",                     # Comment out to get KDE only
        linewidth = 0.2                          # Comment out to get KDE only
      ) +                                        # Comment out to get KDE only
      geom_ribbon(
        aes(x=x, ymin=ymin, ymax=ymax, fill=trial),
        data = kde_data,
        inherit.aes = F,
        color = "white",
        alpha = 0.2
      ) +
      geom_line(linewidth=0.5) + # 0.7           # Comment out to get KDE only
      scale_y_continuous(
        labels = scales::label_percent(accuracy=1),
        breaks = syc_breaks,
        minor_breaks = NULL,
        trans = syc_trans,
        sec.axis = syc_sec.axis
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
      if (is.null(cfg2$llox)) {
        x_axis <- draw.x.axis.cor(zoom_x, NA, cfg2$more_ticks)
      } else {
        x_axis <- draw.x.axis.cor(zoom_x, cfg2$llox, cfg2$more_ticks)
      }
      
      if (flags$janssen_id50_lloq) {
        x_axis$ticks[4] <- log10(2.7426)
        x_axis$labels[[4]] <- "LLOQ"
      }
      
      plot <- plot + scale_x_continuous(
        labels = do.call(expression,x_axis$labels),
        breaks = x_axis$ticks
      )
      
    }
    if (which=="Risk") {
      y_plac <- filter(plot_data, curve=="Placebo overall")[1,"y"]
      y_vacc <- filter(plot_data, curve=="Vaccine overall")[1,"y"]
      plot <- plot + annotate("text", label="Placebo overall", x=zoom_x[2],
                              y=y_plac, size=2.5, hjust=1.05, vjust=-0.5)
      plot <- plot + annotate("text", label="Vaccine overall", x=zoom_x[2],
                              y=y_vacc, size=2.5, hjust=1.05, vjust=-0.5)
    }
    
    return(plot)
    
  }
  
  # Set this manually
  i <- 1
  
  # Load data objects
  folder_144 <- "RV144 plots/"
  folder_705 <- "HVTN 705 (compare RV144) plots/"
  hst_144 <- readRDS(paste0(folder_144, "hst_", i, ".rds"))
  hst_705 <- readRDS(paste0(folder_705, "hst_", i, ".rds"))
  cutoffs_144 <- readRDS(paste0(folder_144, "cutoffs_", i, ".rds"))
  cutoffs_705 <- readRDS(paste0(folder_705, "cutoffs_", i, ".rds"))
  cfg2_144 <- readRDS(paste0(folder_144, "cfg2_", i, ".rds"))
  cfg2_705 <- readRDS(paste0(folder_705, "cfg2_", i, ".rds"))
  dat_v_144 <- readRDS(paste0(folder_144, "dat_v_", i, ".rds"))
  dat_v_705 <- readRDS(paste0(folder_705, "dat_v_", i, ".rds"))
  plot_data_cve_144 <- readRDS(paste0(folder_144,"plot_data_cve_",i,".rds")) %>%
    dplyr::filter(curve=="CVE, nonparametric") %>%
    trim_plot_data(cutoffs_144, cfg2) %>%
    mutate(curve="RV 144")
  plot_data_cve_705 <- readRDS(paste0(folder_705,"plot_data_cve_",i,".rds")) %>%
    dplyr::filter(curve=="CVE, nonparametric") %>%
    trim_plot_data(cutoffs_705, cfg2) %>%
    mutate(curve="HVTN 705")
  plot_data_merged <- rbind(plot_data_cve_144,plot_data_cve_705)
  

  # Hard-code labels
  cfg2 <- cfg2_705 # !!!!!
  cfg2$lab_y <- "Controlled VE against HIV"
  cfg2$lab_title <- strsplit(cfg2$lab_title, ":", fixed=T)[[1]][1]
  
  # Generate plots
  # plot <- create_plot2(
  plot <- create_plot3(
    plot_data = plot_data_merged,
    which = "CVE",
    zoom_x = cfg2$zoom_x,
    zoom_y = cfg2$zoom_y_cve,
    labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
    hst_144 = hst_144,
    hst_705 = hst_705,
    rr_y_axis = T,
    log10_x_axis = T,
    log10_y_axis = F,
    case_dots = F
  )
  
  ggsave(
    filename = paste0("plot_cve_", i, ".pdf"),
    plot=plot, device="pdf", width=6, height=4
  )
  
}



###############################################.
##### ARCHIVE: Summary stats / DQA checks #####
###############################################.

if (F) {
  
  library(ggfortify)
  
  # Alias vectors
  ind_tx <- df_tx[[cfg2$v$event]]
  ind_ct <- df_ct[[cfg2$v$event]]
  time_tx <- df_tx[[cfg2$v$time]]
  time_ct <- df_ct[[cfg2$v$time]]
  
  # Number of cases in each group
  num_case_tx <- sum(ind_tx)
  num_case_ct <- sum(ind_ct)
  num_case_tx_t_0 <- sum(ind_tx[time_tx<=cfg2$t_0])
  num_case_ct_t_0 <- sum(ind_ct[time_ct<=cfg2$t_0])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste0("Number of cases in vaccine group: ", num_case_tx))
  print(paste0("Number of cases in control group: ", num_case_ct))
  print(paste0("Number of cases by day ", cfg2$t_0, " in vaccine group: ",
               num_case_tx_t_0))
  print(paste0("Number of cases by day ", cfg2$t_0, " in control group: ",
               num_case_ct_t_0))
  print(paste0("Number at-risk in vaccine group: ", num_atrisk_tx))
  print(paste0("Number at-risk in control group: ", num_atrisk_ct))
  print(paste0("Naive P(COVID by day ", cfg2$t_0, ") in vaccine group: ",
               round(num_case_tx_t_0/num_atrisk_tx,3)))
  print(paste0("Naive P(COVID by day ", cfg2$t_0, ") in control group: ",
               round(num_case_ct_t_0/num_atrisk_ct,3)))
  print(paste0("Naive vaccine efficacy: ",
               round(1 - (num_case_tx_t_0/num_atrisk_tx) /
                       (num_case_ct_t_0/num_atrisk_ct),3)))
  
  # Fraction of point mass at edge
  s <- dat$v$s
  round(sum(s==min(s,na.rm=T),na.rm=T) / sum(!is.na(s)), 3)
  
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



########################################.
##### ARCHIVE: Grenander debugging #####
########################################.

if (F) {
  
  # Debugging: Examine intermediate objects
  
  x_vals <- ests$Phi_n(ests$grid)
  inds <- !base::duplicated(x_vals)
  x_vals <- x_vals[inds]
  y_vals <- -1 * ests$Gamma_os_n(ests$grid[inds])
  GCM_f <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots,
                     method="linear", rule=1)
  df_A <- data.frame(x=x_vals, y=y_vals, y2=GCM_f(x_vals))
  grid2 <- round(seq(0,1,0.02),2)
  df_B <- data.frame(x=grid2, y=ests$r_Mn_Gr(grid2))
  plot1 <- ggplot(df_A, aes(x=x,y=y)) + geom_point(alpha=0.4) + geom_line(aes(y=y2)) +
    labs(title="x=Phi(S), y=Gamma_n(S)")
  plot2 <- ggplot(df_B, aes(x=x,y=y)) + geom_line() + labs(title="r_Mn estimates")
  
  ggsave(
    filename = paste0(cfg2$analysis," plots/debug_A_",cfg2$tid,".pdf"),
    plot=plot1, device="pdf", width=6, height=4
  )
  ggsave(
    filename = paste0(cfg2$analysis," plots/debug_B_",cfg2$tid,".pdf"),
    plot=plot2, device="pdf", width=6, height=4
  )
  
  # grid <- round(seq(0,1,0.01),2)
  # gcm <- approxfun(x=ests$gcm$x.knots, y=ests$gcm$y.knots, ties="ordered")
  # omega_n <- Vectorize(function(s) {
  #   ests$omega_n(x=c(0,0),s,y=100,delta=0)
  # })
  # etastar_n <- Vectorize(function(s) {
  #   ests$etastar_n(s,x=c(0,0))
  # })
  # Q_n <- Vectorize(function(s) {
  #   ests$Q_n(t=cfg2$t_0, x=c(0,0), s)
  # })
  # 
  # int_data <- data.frame(
  #   x = rep(grid,3),
  #   y = c(ests$Psi_n(grid), gcm(grid), ests$dGCM(grid)),
  #   which = rep(c("Psi_n (-1*Theta_os_n)","gcm","dGCM (-1*r_Mn)"), each=101)
  # )
  # plot1 <- ggplot(int_data, aes(x=x, y=y, color=which)) +
  #   geom_line() +
  #   theme(legend.position="bottom")
  # 
  # ggsave(
  #   filename = paste0(cfg2$analysis," plots/debug_",cfg2$tid,".pdf"),
  #   plot=plot1, device="pdf", width=6, height=4
  # )
  
  if (F) {
    
    # Q_n: conditional survival function (as a function of S)
    # !!!!! Continue
    # as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
    int_data2 <- data.frame(
      x = grid,
      y = Q_n(grid)
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
    
    # # etastar_n
    # int_data4 <- data.frame(
    #   x = grid,
    #   y = etastar_n(grid)
    #   # which = rep(c("",""),each=101)
    # )
    # ggplot(int_data4, aes(x=x, y=y)) + # color=which
    #   geom_line() +
    #   theme(legend.position="bottom") +
    #   labs(title="etastar_n")
    
  }
  
}

if (F) {
  
  # Debugging: Grenander variance scale factor components
  
  print("Grenander variance scale factor components")
  print("deriv_r_Mn")
  print(ests$deriv_r_Mn(seq(0,1,0.05)))
  print("f_s_n")
  print(ests$f_s_n(seq(0,1,0.05)))
  # print("gamma_n")
  # print(ests$gamma_n(seq(0,1,0.05)))
  print("deriv_r_Mn*f_s_n")
  print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$f_s_n(seq(0,1,0.05)))
  # print("deriv_r_Mn*gamma_n")
  # print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  # print("f_s_n*gamma_n")
  # print(ests$f_s_n(seq(0,1,0.05))*ests$gamma_n(seq(0,1,0.05)))
  # print("deriv_r_Mn*f_s_n*gamma_n")
  # print(ests$deriv_r_Mn(seq(0,1,0.05))*ests$f_s_n(seq(0,1,0.05))*
  #         ests$gamma_n(seq(0,1,0.05)))
}



##########################################.
##### ARCHIVE: Data analysis (Qbins) #####
##########################################.

if (F) {
  
  if ("Qbins" %in% cfg2$estimators$cr) {
    
    # Obtain estimates
    s_orig <- dat$v$s[!is.na(dat$v$s)]
    s_grid <- seq(from=min(s_orig), to=max(s_orig), length.out=101)
    ests <- est_curve(
      dat_orig = dat$v,
      estimator = "Qbins",
      params = cfg2$params,
      points = s_grid,
      dir = "decr",
      return_extra = c()
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_q_",cfg2$tid,".rds"))
    }
    
    run_cve <- as.logical("CVE" %in% cfg2$plots)
    ests2 <- process_ests(ests, s_grid, run_cve=run_cve,
                          lab_risk="Risk, Qbins", lab_cve="CVE, Qbins")
    plot_data_risk <- rbind(plot_data_risk, ests2$risk)
    if (run_cve) { plot_data_cve <- rbind(plot_data_cve, ests2$cve) }
    
  }
  
}



##########################.
##### ARCHIVE: Other #####
##########################.

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
    rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-cfg2$t_0))]
    ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-cfg2$t_0))]
    ci_up_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-cfg2$t_0))]
    var_ct <- ((ci_up_ct-ci_lo_ct)/3.92)^2
    
    # Calculate treatment group survival (KM; with SE)
    srv_tx <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_tx
    )
    rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-cfg2$t_0))]
    ci_lo_tx <- 1 - srv_tx$upper[which.min(abs(srv_tx$time-cfg2$t_0))]
    ci_up_tx <- 1 - srv_tx$lower[which.min(abs(srv_tx$time-cfg2$t_0))]
    var_tx <- ((ci_up_tx-ci_lo_tx)/3.92)^2
    
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
    rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_tx_sub$time-cfg2$t_0))]
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
      x = c(x1,rep(s_grid,2),x1,x2),
      y = c(ests_cve[1],ests_cve,ests_gcomp,rep(ve_overall,2)),
      which = c("Controlled VE",rep(c("Controlled VE","Cox model"),
                                    each=length(p_grid)),
                rep("Overall VE",2)),
      ci_lo = c(ci_lo[1],ci_lo,ests_gcomp,rep(ve_ci[1],2)),
      ci_up = c(ci_up[1],ci_up,ests_gcomp,rep(ve_ci[2],2))
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
