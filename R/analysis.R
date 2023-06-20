# SAP: https://www.overleaf.com/project/604a54625b885d0da667de4b

#################.
##### Setup #####
#################.

{
  # Choose analysis
  which_analysis <- "Janssen (partA)" # "Janssen" "Moderna" "AMP" "AZD1222" "Janssen (partA)" "Profiscov" "HVTN 705 (primary)" "HVTN 705 (all)"
  
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
  
  # Uncomment this code to run multiple analyses (e.g. 1=4=Janssen, 5-14=Moderna)
  ..tid <- as.integer(Sys.getenv(.tid_var))
  if (..tid %in% c(1:4)) {
    which_analysis <- "Janssen"
    .tid_lst = list(as.character(round(..tid)))
  } else if (..tid %in% c(5:14)) {
    which_analysis <- "Moderna"
    .tid_lst = list(as.character(round(..tid-4)))
  } else if (..tid %in% c(15:26)) {
    which_analysis <- "HVTN 705 (ICS)"
    .tid_lst = list(as.character(round(..tid-14)))
  } else if (..tid %in% c(27:84)) {
    which_analysis <- "Janssen (partA)"
    .tid_lst = list(as.character(round(..tid-26)))
  }
  names(.tid_lst) = .tid_var
  do.call(Sys.setenv, .tid_lst)
  
  # Set seed
  set.seed(1)
  
  # Set configuration variables
  # Note: NA param values are set below based on the task ID
  cfg2 <- list(
    analysis = which_analysis,
    run_dqa = F,
    run_debug = list(gren_var=F, objs=F),
    run_hyptest = F
  )
  
  # Set analysis-specific flags
  # Note: some flags are set at the end of the "Setup" block because they are
  #       dependent on cfg2 variables
  flags <- list(
    hvtn705_abstract_fig = F,
    table_of_vals = F,
    save_data_objs = F,
    paper_npcve = F,
    paper_cox = F,
    hvtn124_plot = F
  )
  
  # Set up analysis-specific configuration variables. Each row in the cfg2$map
  #     dataframe represents the set of indices to use for a particular
  #     analysis.
  if (cfg2$analysis=="Janssen") {
    
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29ADCP")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to RBD: Day 29", "PsV Neutralization 50% Titer: Day 29", "Phagocytic Score: Day 29")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Phagocytic Score (=s)")
    cfg2$endpoint <- "COVID"
    cfg2$t_0 <- 54
    # Note: "janssen_pooled_real_..." changed to "janssen_pooled_EUA_"
    cfg2$dataset <- c("janssen_pooled_real_data_processed_with_riskscore.csv", "janssen_pooled_realADCP_data_processed_with_riskscore.csv")
    cfg2$txct <- T
    cfg2$cr2_trial <- c("janssen_pooled_real", "janssen_pooled_realADCP")
    cfg2$cr2_COR <- "D29IncludeNotMolecConfirmedstart1"
    cfg2$cr2_marker <- c(1,2,3)
    cfg2$edge_corr <- "min"
    cfg2$v <- list(
      id = "Ptid",
      time = "EventTimePrimaryIncludeNotMolecConfirmedD29",
      event = "EventIndPrimaryIncludeNotMolecConfirmedD29",
      wt = "wt.D29start1",
      ph1 = "ph1.D29start1",
      ph2 = "ph2.D29start1",
      covariates = "~. + risk_score + as.factor(Region)"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0,0.95),
      "CVE, nonparametric" = c(0,0.95),
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox model" = c(0,0.975),
      "CVE, Cox model" = c(0,0.975),
      "Risk, Cox (basic)" = c(0,0.975),
      "CVE, Cox (basic)" = c(0,0.975)
    )
    cfg2$zoom_x <- NA
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "Janssen data/"
    cfg2$folder_cluster <- "Z:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$params = list(g_n_type="binning", deriv_type="m-spline",
                       ci_type="logit", q_n_type="zero", # ci_type="trunc", # Used historically
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=1,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4),
      lab_x = c(1,2,3,4),
      lab_title = c(1,2,3,4),
      t_0 = c(1,1,1,1),
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
      v_covariates = c(1,1,1,1),
      zoom_x = c(1,1,1,1),
      zoom_y_cve = c(1,1,1,1),
      zoom_y_risk = c(1,1,1,1),
      more_ticks = c(1,1,1,1)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:4),
      map_row = c(1:4)
      # tid = c(1:8),
      # edge_ind = rep(c(F,T), 4),
      # map_row = rep(c(1:4), each=2)
    )
    
  }
  
  if (cfg2$analysis=="Moderna") {
    
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day29bindSpike", "Day57bindSpike", "Day29bindRBD", "Day57bindRBD", "Day29pseudoneutid50", "Day57pseudoneutid50", "Day29pseudoneutid80", "Day57pseudoneutid80", "Day29liveneutmn50", "Day57liveneutmn50")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to Spike: Day 57", "Binding Antibody to RBD: Day 29", "Binding Antibody to RBD: Day 57", "PsV Neutralization 50% Titer: Day 29", "PsV Neutralization 50% Titer: Day 57", "PsV Neutralization 80% Titer: Day 29", "PsV Neutralization 80% Titer: Day 57", "Live Virus Micro Neut 50% Titer: Day 29", "Live Virus Micro Neut 50% Titer: Day 57")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Pseudovirus-nAb ID80 (IU80/ml) (=s)", "Live Virus-mnAb ID50 (IU50/ml) (=s)")
    cfg2$endpoint <- "COVID"
    cfg2$t_0 <- c(126,100) # Try changing to 0
    cfg2$dataset <- "P3001ModernaCOVEimmunemarkerdata_correlates_processed_v1.1_lvmn_added_Jan14_2022.csv"
    cfg2$txct <- T
    cfg2$cr2_trial <- "moderna_real"
    cfg2$cr2_COR <- c("D29", "D57")
    cfg2$cr2_marker <- c(1,2,3,4,5)
    cfg2$edge_corr <- c("none", "min")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD29", "EventTimePrimaryD57"),
      event = c("EventIndPrimaryD29", "EventIndPrimaryD57"),
      wt = c("wt.D29", "wt.D57"),
      ph1 = c("ph1.D29", "ph1.D57"),
      ph2 = c("ph2.D29", "ph2.D57"),
      covariates = "~. + MinorityInd + HighRiskInd + risk_score"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.05,0.95),
      "CVE, nonparametric" = c(0.05,0.95),
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975),
      "Risk, Cox GAM" = c(0.025,0.975),
      "CVE, Cox GAM" = c(0.025,0.975),
      "Risk, Cox (analytic)" = c(0.025,0.975),
      "CVE, Cox (analytic)" = c(0.025,0.975),
      "Risk, Cox (basic)" = c(0.025,0.975),
      "CVE, Cox (basic)" = c(0.025,0.975),
      "Risk, Cox (spline 4 df)" = c(0.025,0.975),
      "CVE, Cox (spline 4 df)" = c(0.025,0.975),
      "Risk, Cox (edge)" = c(0.025,0.975),
      "CVE, Cox (edge)" = c(0.025,0.975)
    )
    # cfg2$zoom_x <- list(
    #   c(0.9,3.34),
    #   c(0.9,4.04),
    #   c(0.9,3.44),
    #   c(0.9,4.24),
    #   "zoomed llox"
    # )
    cfg2$zoom_x <- "zoomed"
    # cfg2$zoom_y_cve <- list(c(0.58,1.02)) # "zoomed"
    cfg2$zoom_y_risk <- list(c(-0.002,0.072))
    cfg2$zoom_y_cve <- "zoomed"
    # cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "Moderna data/"
    cfg2$folder_cluster <- "Z:/covpn/p3001/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$llox_label <- "LOD" # NEW
    cfg2$llox <- c(0.3076,1.594,2.42,15.02,22.66) # NEW
    cfg2$params = list(g_n_type="binning", deriv_type="m-spline",
                       ci_type="regular", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=1,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4,5,6,7,8,9,10),
      lab_x = c(1,1,2,2,3,3,4,4,5,5),
      lab_title = c(1,2,3,4,5,6,7,8,9,10),
      t_0 = c(1,2,1,2,1,2,1,2,1,2),
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
      v_covariates = c(1,1,1,1,1,1,1,1,1,1),
      # zoom_x = c(1,2,3,4,5,5,5,5,5,5),
      zoom_x = c(1,1,1,1,1,1,1,1,1,1),
      zoom_y_cve = c(1,1,1,1,1,1,1,1,1,1),
      zoom_y_risk = c(1,1,1,1,1,1,1,1,1,1),
      more_ticks = c(1,1,1,1,1,1,1,1,1,1),
      llox_label = c(1,1,1,1,1,1,1,1,1,1),
      llox = c(1,1,2,2,3,3,4,4,5,5)
    )
    
    # Secondary map for variations within a graph; map_row corresponds to which
    #     row of cfg2$map to use
    cfg2$map2 <- data.frame(
      tid = c(1:10),
      map_row = c(1:10)
      # tid = c(1:20),
      # # Q_n_type = rep(c("Cox", "Super Learner"),10),
      # # ci_type = rep(c("regular", "logit"),10),
      # edge_ind = rep(c(F,T), 10),
      # map_row = rep(c(1:10), each=2)
    )
    
    # Flag-specific operation
    if (flags$paper_npcve || flags$paper_cox) {
      
      cfg2$zoom_x <- c("zoomed", "zoomed llox")
      cfg2$map$zoom_x <- c(1,1,1,1,2,1,2,1,2,1)
      
    }
    
  }
  
  if (cfg2$analysis=="AMP") {
    
    cfg2$estimators <- list(overall="KM", cr=c("Grenander")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk")
    cfg2$marker <- "bweight"
    cfg2$lab_title <- c("HVTN703/HPTN081", "HVTN704/HPTN085", "Pooled AMP trials")
    cfg2$lab_x <- "Body Weight (kg)"
    cfg2$endpoint <- "HIV-1 infection"
    cfg2$amp_protocol <- c("HVTN 703", "HVTN 704", "Pooled")
    cfg2$amp_tx <- c("C3", "T1", "T2", "T1+T2")
    cfg2$t_0 <- c(595,595) # c(601,609)
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
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y_cve <- "zoomed"
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "AMP data/"
    cfg2$folder_cluster <- "Z:/vaccine/p704/analysis/datashare/avi_kenny/adata/"
    cfg2$params = list(g_n_type="binning", deriv_type="linear",
                       ci_type="regular", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=10,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      amp_protocol = c(1,1,1,1,2,2,2,2,3,3,3,3),
      amp_tx = c(1,2,3,4,1,2,3,4,1,2,3,4),
      marker = rep(1,12),
      lab_x = rep(1,12),
      lab_title = c(1,1,1,1,2,2,2,2,3,3,3,3),
      t_0 = c(1,1,1,1,2,2,2,2,1,1,1,1),
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
      v_covariates = c(1,1,1,1,1,1,1,1,2,2,2,2),
      zoom_x = rep(1,12),
      zoom_y_cve = rep(1,12),
      zoom_y_risk = rep(1,12),
      more_ticks = rep(1,12)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:12),
      map_row = c(1:12)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (all)") {
    
    # cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Cox gcomp", "Cox GAM", "Cox edge")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day210ELCZ", "Day210ELMo", "Day210ADCPgp140C97ZAfib", "Day210ADCPgp140Mos1fib", "Day210IgG3gp140C97ZAfibritin40delta", "Day210IgG3gp140Mos1fibritin40delta", "Day210IgG340mdw_gp120", "Day210IgG340mdw_gp140", "Day210IgG340mdw_V1V2", "Day210IgG3gp4140delta", "Day210IgG340mdw_multi", "Day210IgG340mdw_gp120_gp140_vm", "Day210IgG50mdw_V1V2", "Day210mdw_xassay", "Day210ADCCCAP8_pk", "Day210ADCCCH58_pk", "Day210ADCCWITO_pk", "Day210ADCCCAP8_pAUC", "Day210ADCCCH58_pAUC", "Day210ADCCWITO_pAUC", "Day210ICS4AnyEnvIFNg_OR_IL2", "Day210ICS8AnyEnvIFNg_OR_IL2", "Day210IgG3AE.A244.V1V2.Tags_293F40delta", "Day210IgG3C.1086C.V1.V2.Tags40delta", "Day210IgG3gp70.001428.2.42.V1V240delta", "Day210IgG3gp70.1012.11.TC21.3257.V1V240delta", "Day210IgG3gp70.1394C9G1.V1V240delta", "Day210IgG3gp70.BF1266.431a.V1V240delta", "Day210IgG3gp70.Ce1086.B2.V1V240delta", "Day210IgG3gp70.B.CaseA.V1.V240delta", "Day210IgGAE.A244.V1V2.Tags_293F50delta", "Day210IgGC.1086C.V1.V2.Tags50delta", "Day210IgGgp70_001428.2.42.V1V250delta", "Day210IgGgp70_1012.11.TC21.3257.V1V250delta", "Day210IgGgp70_1394C9G1.V1V250delta", "Day210IgGgp70_9004SS.A3.4.V1V250delta", "Day210IgGgp70_BF1266.431a.V1V250delta", "Day210IgGgp70_Ce1086.B2.V1V250delta", "Day210IgGgp70.B.CaseA.V1.V250delta")
    cfg2$lab_title <- c("IgG to VT-C (EU/ml): Day 210", "IgG to VT-M (EU/ml): Day 210", "Average phagocytosis score to gp140 C97ZA: Day 210", "Average phagocytosis score to gp140 Mos1: Day 210", "IgG3 Net MFI to gp140 C97ZA: Day 210", "IgG3 Net MFI to gp140 Mosaic: Day 210", "IgG3 gp120 breadth (Weighted avg log10 Net MFI): Day 210", "IgG3 gp140 breadth (Weighted avg log10 Net MFI): Day 210", "IgG3 V1V2 breadth (Weighted avg log10 Net MFI): Day 210", "IgG3 Net MFI to gp41: Day 210", "IgG3 multi-epitope breadth (Wt avg log10 Net MFI): Day 210", "IgG3 gp120 + gp140 breadth (Wt avg log10 Net MFI): Day 210", "IgG V1V2 breadth (Wt avg log10 Net MFI): Day 210", "Overall maximal diversity score: Day 210", "Peak baseline-subtracted pct loss luc activity to CAP8: Day 210", "Peak baseline-subtracted pct loss luc activity to CH58: Day 210", "Peak baseline-subtracted pct loss luc activity to WITO: Day 210", "AUC baseline-subtracted pct loss luc activity to CAP8: Day 210", "AUC baseline-subtracted pct loss luc activity to CH58: Day 210", "AUC baseline-subtracted pct loss luc activity to WITO: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2: Day 210", "IgG3 Net MFI to AE.A244 V1V2 Tags 293F: Day 210", "IgG3 Net MFI to C.1086C V1V2 Tags: Day 210", "IgG3 Net MFI to gp70-001428.2.42 V1V2: Day 210", "IgG3 Net MFI to gp70-1012.11.TC21.3257 V1V2: Day 210", "IgG3 Net MFI to gp70-1394C9G1 V1V2: Day 210", "IgG3 Net MFI to gp70-BF1266 431a V1V2: Day 210", "IgG3 Net MFI to gp70-Ce1086 B2 V1V2: Day 210", "IgG3 Net MFI to gp70-B.CaseAV1V2: Day 210", "IgG Net MFI to AE.A244 V1V2 Tags 293F: Day 210", "IgG Net MFI to C.1086C V1V2 Tags: Day 210", "IgG Net MFI to gp70-001428.2.42 V1V2: Day 210", "IgG Net MFI to gp70-1012.11.TC21.3257 V1V2: Day 210", "IgG Net MFI to gp70-1394C9G1 V1V2: Day 210", "IgG Net MFI to gp70-9004SS.A3.4 V1V2: Day 210", "IgG Net MFI to gp70-BF1266.431a V1V2: Day 210", "IgG Net MFI to gp70-Ce1086.B2 V1V2: Day 210", "IgG Net MFI to gp70.B.CaseA V1V2: Day 210")
    cfg2$lab_x <- c("IgG to VT-C (=s)", "IgG to VT-M (=s)", "ADCP gp140 C97ZA (=s)", "ADCP gp140 Mos1 (=s)", "IgG3 gp140 C97ZA (=s)", "IgG3 gp140 Mosaic (=s)", "IgG3 gp120 breadth (=s)", "IgG3 gp140 breadth (=s)", "IgG3 V1V2 breadth (=s)", "IgG3 gp41 (=s)", "IgG3 multi-epitope breadth (=s)", "IgG3 gp120+gp140 breadth (=s)", "IgG V1V2 breadth (=s)", "Overall max diversity score (=s)", "ADCC Peak CAP8 (=s)", "ADCC Peak CH58 (=s)", "ADCC Peak WITO (=s)", "ADCC  AUC CAP8 (=s)", "ADCC AUC CH58 (=s)", "ADCC AUC WITO (=s)", "CD4+ T-cells IFN-g/IL-2 (=s)", "CD8+ T-cells IFN-g/IL-2 (=s)", "IgG3 AE.A244 V1V2 Tags 293F (=s)", "IgG3 C.1086C V1V2 Tags (=s)", "IgG3 gp70-001428.2.42 V1V2 (=s)", "IgG3 gp70-1012.11.TC21.3257 V1V2 (=s)", "IgG3 gp70-1394C9G1 V1V2 (=s)", "IgG3 gp70-BF1266 431a V1V2 (=s)", "IgG3 gp70-Ce1086 B2 V1V2 (=s)", "IgG3 gp70-B.CaseAV1V2 (=s)", "IgG AE.A244 V1V2 Tags 293F (=s)", "IgG C.1086C V1V2 Tags (=s)", "IgG gp70-001428.2.42 V1V2 (=s)", "IgG gp70-1012.11.TC21.3257 V1V2 (=s)", "IgG gp70-1394C9G1 V1V2 (=s)", "IgG gp70-9004SS.A3.4 V1V2 (=s)", "IgG gp70-BF1266.431a V1V2 (=s)", "IgG gp70-Ce1086.B2 V1V2 (=s)", "IgG gp70.B.CaseA V1V2 (=s)")
    cfg2$endpoint <- "HIV"
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_excludeELISpotmarkers.csv"
    cfg2$txct <- T
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$cr2_marker <- c(1:39)
    cfg2$edge_corr <- c("none", "min")
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.05,0.95), # c(0.1,0.9)
      "CVE, nonparametric" = c(0.05,0.95), # c(0.1,0.9)
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox GAM" = c(0.025,0.975),
      "CVE, Cox GAM" = c(0.025,0.975),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- "zoomed" # !!!!! Changed from NA
    cfg2$zoom_y_cve <- "zoomed"
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "HVTN 705 (all) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$params = list(g_n_type="binning", deriv_type="line",
                       ci_type="regular", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=10,x_tol=15,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1:39),
      lab_x = c(1:39),
      lab_title = c(1:39),
      t_0 = rep(1,39),
      dataset = rep(1,39),
      cr2_trial = rep(1,39),
      cr2_COR = rep(1,39),
      cr2_marker = c(1:39),
      edge_corr = c(1,1,1,1,1,1,2,1,2,1,1,1,2,1,2,2,2,2,2,2,
                    2,2,1,2,2,2,2,2,2,2,1,1,2,2,2,2,2,1,2),
      v_id = rep(1,39),
      v_time = rep(1,39),
      v_event = rep(1,39),
      v_wt = rep(1,39),
      v_ph1 = rep(1,39),
      v_ph2 = rep(1,39),
      v_covariates = rep(1,39),
      zoom_x = rep(1,39),
      zoom_y_cve = rep(1,39),
      zoom_y_risk = rep(1,39),
      more_ticks = rep(1,39)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:39),
      map_row = c(1:39)
    )
    
  }
  
  if (cfg2$analysis=="HVTN 705 (ICS)") {
    
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day210ICS4JMos1gp120IFNg_OR_IL2", "Day210ICS4JMos1gp41IFNg_OR_IL2", "Day210ICS4JMos2GagIFNg_OR_IL2", "Day210ICS4JMos2RNAseIntIFNg_OR_IL2", "Day210ICS4JMos2Sgp120IFNg_OR_IL2", "Day210ICS4JMos2Sgp41IFNg_OR_IL2", "Day210ICS8JMos1gp120IFNg_OR_IL2", "Day210ICS8JMos1gp41IFNg_OR_IL2", "Day210ICS8JMos2GagIFNg_OR_IL2", "Day210ICS8JMos2RNAseIntIFNg_OR_IL2", "Day210ICS8JMos2Sgp120IFNg_OR_IL2", "Day210ICS8JMos2Sgp41IFNg_OR_IL2")
    cfg2$lab_title <- c("Pct CD4+ T-cells expressing IFN-g/IL-2 JMos1 gp120: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos1 gp41: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2 Gag: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2 RNAseInt: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2s gp120: Day 210", "Pct CD4+ T-cells expressing IFN-g/IL-2 JMos2s gp41: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos1 gp120: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos1 gp41: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2 Gag: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2 RNAseInt: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2s gp120: Day 210", "Pct CD8+ T-cells expressing IFN-g/IL-2 JMos2s gp41: Day 210")
    cfg2$lab_x <- c("CD4+ T-cells IFN-g/IL-2 JMos1 gp120 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos1 gp41 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2 Gag (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2 RNAseInt (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2s gp120 (=s)", "CD4+ T-cells IFN-g/IL-2 JMos2s gp41 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos1 gp120 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos1 gp41 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2 Gag (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2 RNAseInt (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2s gp120 (=s)", "CD8+ T-cells IFN-g/IL-2 JMos2s gp41 (=s)")
    cfg2$endpoint <- "HIV"
    cfg2$t_0 <- 550
    cfg2$dataset <- "HVTN705_secondcasecontrolprocesseddata_excludeELISpotmarkersaddICSmarkers.csv"
    cfg2$txct <- T
    cfg2$cr2_trial <- "hvtn705second"
    cfg2$cr2_COR <- "D210"
    cfg2$cr2_marker <- c(1:12)
    cfg2$edge_corr <- c("none", "min")
    cfg2$v <- list(
      id = "Subjectid",
      time = "Ttilde.D210",
      event = "Delta.D210",
      wt = "wt.D210",
      ph1 = "Ph1ptids.D210",
      ph2 = "Ph2ptids.D210",
      covariates = "~. + RSA + Age + BMI + Riskscore"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0,0.95),
      "CVE, nonparametric" = c(0,0.95),
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox GAM" = c(0,0.975),
      "CVE, Cox GAM" = c(0,0.975),
      "Risk, Cox model" = c(0,0.975),
      "CVE, Cox model" = c(0,0.975),
      "Risk, Cox (basic)" = c(0,0.975),
      "CVE, Cox (basic)" = c(0,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y_cve <- "zoomed"
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "HVTN 705 (ICS) data/"
    cfg2$folder_cluster <- "Z:/vaccine/p705/analysis/lab/cc/copcor/"
    cfg2$params = list(g_n_type="binning", deriv_type="line",
                       ci_type="logit", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=10,x_tol=15,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1:12),
      lab_x = c(1:12),
      lab_title = c(1:12),
      t_0 = rep(1,12),
      dataset = rep(1,12),
      cr2_trial = rep(1,12),
      cr2_COR = rep(1,12),
      cr2_marker = c(1:12),
      edge_corr = c(2,2,2,2,2,2,2,2,2,2,2,2),
      v_id = rep(1,12),
      v_time = rep(1,12),
      v_event = rep(1,12),
      v_wt = rep(1,12),
      v_ph1 = rep(1,12),
      v_ph2 = rep(1,12),
      v_covariates = rep(1,12),
      zoom_x = rep(1,12),
      zoom_y_cve = rep(1,12),
      zoom_y_risk = rep(1,12),
      more_ticks = rep(1,12)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:12),
      map_row = c(1:12)
    )
    
  }
  
  if (cfg2$analysis=="AZD1222") {
    
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day29pseudoneutid50", "Day57pseudoneutid50",
                     "Day29bindSpike", "Day57bindSpike")
    cfg2$lab_title <- c("PsV Neutralization 50% Titer: Day 29", "PsV Neutralization 50% Titer: Day 57", "Binding Antibody to Spike: Day 29", "Binding Antibody to Spike: Day 57")
    cfg2$lab_x <- c("Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Anti Spike IgG (BAU/ml) (=s)")
    cfg2$endpoint <- "COVID"
    cfg2$t_0 <- c(117,92)
    cfg2$dataset <- c("azd1222_data_processed_with_riskscore.csv", "azd1222_bAb_data_processed_with_riskscore.csv")
    cfg2$txct <- T
    cfg2$cr2_trial <- c("azd1222", "azd1222_bAb")
    cfg2$cr2_COR <- c("D29", "D57")
    cfg2$cr2_marker <- c(1,2)
    cfg2$edge_corr <- c("none", "min")
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD29", "EventTimePrimaryD57"),
      event = c("EventIndPrimaryD29", "EventIndPrimaryD57"),
      wt = c("wt.D29", "wt.D57"),
      ph1 = c("ph1.D29", "ph1.D57"),
      ph2 = c("ph2.D29", "ph2.D57"),
      covariates = "~. + Age + risk_score"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.1,0.9),
      "CVE, nonparametric" = c(0.1,0.9),
      "Risk, Qbins" = c(0,1),
      "CVE, Qbins" = c(0,1),
      "Risk, Cox GAM" = c(0.025,0.975),
      "CVE, Cox GAM" = c(0.025,0.975),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- NA
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- NA
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "AZD1222 data/"
    cfg2$folder_cluster <- "Z:/covpn/p3002/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$params = list(g_n_type="binning", deriv_type="line",
                       ci_type="regular", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=10,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1,2,3,4),
      lab_x = c(1,1,2,2),
      lab_title = c(1,2,3,4),
      t_0 = c(1,2,1,2),
      dataset = c(1,1,2,2),
      cr2_trial = c(1,1,2,2),
      cr2_COR = c(1,2,1,2),
      cr2_marker = c(1,1,1,1),
      edge_corr = c(2,1,1,1),
      v_id = c(1,1,1,1),
      v_time = c(1,2,1,2),
      v_event = c(1,2,1,2),
      v_wt = c(1,2,1,2),
      v_ph1 = c(1,2,1,2),
      v_ph2 = c(1,2,1,2),
      v_covariates = c(1,1,1,1),
      zoom_x = c(1,1,1,1),
      zoom_y_cve = c(1,1,1,1),
      zoom_y_risk = c(1,1,1,1),
      more_ticks = c(1,1,1,1)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:4),
      map_row = c(1:4)
    )
    
  }
  
  if (cfg2$analysis=="Janssen (partA)") {
    
    # Initial analysis: 1-58
    # Second manuscript: 1-3,13-15,25-27,37-39,45-47,49-51,59-64
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    # cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander")) # Janssen partA manuscript main figures
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day29bindSpike", "Day29bindRBD", "Day29pseudoneutid50", "Day29ADCP", "Day29pseudoneutid50la", "Day29pseudoneutid50sa")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 29", "Binding Antibody to RBD: Day 29", "PsV Neutralization 50% Titer: Day 29", "Phagocytic Score: Day 29", "PsV Neutralization 50% Titer (LA): Day 29", "PsV Neutralization 50% Titer (SA): Day 29")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Pseudovirus-nAb ID50 (IU50/ml) (=s)", "Phagocytic Score (=s)", "Pseudovirus-nAb ID50 LA (IU50/ml) (=s)", "Pseudovirus-nAb ID50 SA (IU50/ml) (=s)")
    cfg2$endpoint <- "COVID"
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
    cfg2$txct <- T
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
    cfg2$cr2_marker <- c(1,2,3,4,5)
    cfg2$edge_corr <- c("min")
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
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0,0.9),
      "CVE, nonparametric" = c(0,0.9),
      "Risk, Cox model" = c(0,0.975),
      "CVE, Cox model" = c(0,0.975),
      "Risk, Cox GAM" = c(0,0.975),
      "CVE, Cox GAM" = c(0,0.975),
      "Risk, Cox (analytic)" = c(0,0.975),
      "CVE, Cox (analytic)" = c(0,0.975),
      "Risk, Cox (basic)" = c(0,0.975),
      "CVE, Cox (basic)" = c(0,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- 0.15
    cfg2$more_ticks <- c(1,2)
    cfg2$folder_local <- "Janssen (partA) data/"
    cfg2$folder_cluster <- "Z:/covpn/p3003/analysis/correlates/Part_A_Blinded_Phase_Data/adata/"
    cfg2$llox_label <- c("LOD", "LLOQ") # !!!!!
    cfg2$llox <- c(NA, 4.8975) # !!!!!
    cfg2$params = list(g_n_type="parametric (edge) 2", deriv_type="m-spline",
                       ci_type="logit", q_n_type="zero",
                       Q_n_type="Super Learner")
    C <- list(appx=list(t_0=1,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      lab_x = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      lab_title = c(rep(c(1:4), 13), c(5,5,5,6,6,5), c(1,2,3,1,2,3)), # A (marker)
      t_0 = rep(1, 64), # C (same)
      dataset = c(rep(c(1:11), each=4), rep(1,4), rep(7,4), c(7:11), 7, c(1,1,1,7,7,7)), # D (region X age)
      cr2_trial = c(rep(c(1:11), each=4), rep(1,4), rep(7,4), c(7:11), 7, c(1,1,1,7,7,7)), # D (region X age)
      cr2_COR = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      cr2_marker = c(rep(c(1:4),6), rep(c(1,2,3,5),5), c(1:4), c(1,2,3,5), rep(4,6), c(1,2,3,1,2,3)),
      edge_corr = rep(1, 64), # C (same)
      v_id = rep(1, 64), # C (same)
      v_time = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      v_event = c(rep(1,44), rep(2,8), rep(1,5), 2, rep(3,6)), # E (regular/severe/moderate)
      v_wt = rep(1, 64), # C (same)
      v_ph1 = rep(1, 64), # C (same)
      v_ph2 = rep(1, 64), # C (same)
      v_covariates = c(rep(1,11), rep(2,33), rep(1,4), rep(2,4), rep(2,6), c(1,1,1,2,2,2)), # B (pooled or not)
      zoom_x = rep(1, 64), # C (same)
      zoom_y_cve = rep(1, 64), # C (same)
      zoom_y_risk = rep(1, 64), # C (same)
      # more_ticks = replace(rep(1,64),c(3,47,61),2), # G (custom)
      more_ticks = rep(2, 64), # G (custom)
      llox_label = c(rep(c(1,1,2,2), 13), rep(2,6), c(1,1,2,1,1,2)), # F (llox)
      llox = c(rep(c(1,1,2,2), 13), rep(2,6), c(1,1,2,1,1,2)) # F (llox)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:64),
      map_row = c(1:64)
    )
    
  }
  
  if (cfg2$analysis=="Profiscov") {
    
    cfg2$estimators <- list(overall="Cox gcomp", cr=c("Grenander", "Cox gcomp")) # Types: "Cox import", "Grenander", "Cox gcomp", "Qbins", "Cox GAM", "Cox edge"
    cfg2$plots <- c("Risk", "CVE")
    cfg2$marker <- c("Day43bindSpike", "Day43bindSpike_B.1.1.7", "Day43bindSpike_B.1.351", "Day43bindSpike_P.1", "Day43bindRBD", "Day43bindRBD_B.1.1.7", "Day43bindRBD_B.1.351", "Day43bindRBD_P.1", "Day43bindN",
                     "Day91bindSpike", "Day91bindSpike_B.1.1.7", "Day91bindSpike_B.1.351", "Day91bindSpike_P.1", "Day91bindRBD", "Day91bindRBD_B.1.1.7", "Day91bindRBD_B.1.351", "Day91bindRBD_P.1", "Day91bindN",
                     "Day43liveneutmn50")
    cfg2$lab_title <- c("Binding Antibody to Spike: Day 43", "Binding Antibody to Spike B.1.1.7: Day 43", "Binding Antibody to Spike B.1.351: Day 43", "Binding Antibody to Spike P.1: Day 43", "Binding Antibody to RBD: Day 43", "Binding Antibody to RBD B.1.1.7: Day 43", "Binding Antibody to RBD B.1.351: Day 43", "Binding Antibody to RBD P.1: Day 43", "Binding Antibody to Nucleocapsid: Day 43",
                        "Binding Antibody to Spike: Day 91", "Binding Antibody to Spike B.1.1.7: Day 91", "Binding Antibody to Spike B.1.351: Day 91", "Binding Antibody to Spike P.1: Day 91", "Binding Antibody to RBD: Day 91", "Binding Antibody to RBD B.1.1.7: Day 91", "Binding Antibody to RBD B.1.351: Day 91", "Binding Antibody to RBD P.1: Day 91", "Binding Antibody to Nucleocapsid: Day 91",
                        "Live Virus Micro Neut 50% Titer: Day 43")
    cfg2$lab_x <- c("Anti Spike IgG (BAU/ml) (=s)", "Anti Spike B.1.1.7 IgG (BAU/ml) (=s)", "Anti Spike B.1.351 IgG (BAU/ml) (=s)", "Anti Spike P.1 IgG (BAU/ml) (=s)", "Anti RBD IgG (BAU/ml) (=s)", "Anti RBD B.1.1.7 IgG (BAU/ml) (=s)", "Anti RBD B.1.351 IgG (BAU/ml) (=s)", "Anti RBD P.1 IgG (BAU/ml) (=s)", "Anti N IgG (BAU/ml) (=s)", "")
    cfg2$endpoint <- "COVID"
    cfg2$t_0 <- c(114,66)
    cfg2$dataset <- c("profiscov_data_processed_with_riskscore.csv",
                      "profiscov_lvmn_data_processed_with_riskscore.csv")
    cfg2$txct <- T
    cfg2$cr2_trial <- c("profiscov", "profiscov_lvmn") # dummy; not using Cox estimates
    cfg2$cr2_COR <- c("D43", "D91")
    cfg2$cr2_marker <- c(1:9)
    cfg2$edge_corr <- "none"
    cfg2$v <- list(
      id = "Ptid",
      time = c("EventTimePrimaryD43", "EventTimePrimaryD91"),
      event = c("EventIndPrimaryD43", "EventIndPrimaryD91"),
      wt = c("wt.D43", "wt.D91"),
      ph1 = c("ph1.D43", "ph1.D91"),
      ph2 = c("ph2.D43", "ph2.D91"),
      covariates = "~.+ HighRiskInd + Sex + Age + BMI"
    )
    cfg2$qnt <- list(
      "Risk, nonparametric" = c(0.05,0.95),
      "CVE, nonparametric" = c(0.05,0.95),
      "Risk, Cox model" = c(0.025,0.975),
      "CVE, Cox model" = c(0.025,0.975)
    )
    cfg2$zoom_x <- "zoomed"
    cfg2$zoom_y_cve <- NA
    cfg2$zoom_y_risk <- "zoomed (risk)"
    cfg2$zoom_y_max <- 0.21
    cfg2$more_ticks <- 1
    cfg2$folder_local <- "Profiscov data/"
    cfg2$folder_cluster <- "Y:/cavd/Objective 4/GH-VAP/ID127-Gast/correlates/adata/"
    cfg2$llox_label <- c("LLOQ", "LOD")
    cfg2$llox <- c(49*0.009, 70*0.009, 72*0.009, 32*0.009, 35*0.0272, 224*0.0272, 53*0.0272, 91*0.0272, 46*0.00236, 27.56)
    cfg2$params = list(
      g_n_type="parametric", deriv_type="m-spline", # "binning"
      ci_type="logit", q_n_type="zero",
      Q_n_type="Super Learner"
    )
    C <- list(appx=list(t_0=1,x_tol=25,s=0.01))
    
    # Variable map; one row corresponds to one CVE graph
    cfg2$map <- data.frame(
      marker = c(1:19), # A (same)
      lab_x = c(rep(c(1:9),2),10), # B (marker type)
      lab_title = c(1:19), # A (same)
      t_0 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      dataset = c(rep(1, 18), 2), # D (regular vs. LVNT)
      cr2_trial = c(rep(1, 18), 2), # D (regular vs. LVNT)
      cr2_COR = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      cr2_marker = c(rep(c(1:9),2),1), # E (marker order cr2)
      edge_corr = rep(1, 19), # A (same)
      v_id = rep(1, 19), # A (same)
      v_time = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_event = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_wt = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_ph1 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_ph2 = c(rep(c(1,2), each=9), 1), # C (day 43 vs 91)
      v_covariates = rep(1, 19), # A (same)
      zoom_x = rep(1, 19), # A (same)
      zoom_y_cve = rep(1, 19), # A (same)
      zoom_y_risk = rep(1, 19), # A (same)
      more_ticks = rep(1, 19), # A (same)
      llox_label = c(rep(1, 18), 2), # D (regular vs. LVNT)
      llox = c(rep(c(1:9),2),10) # B (marker type)
    )
    
    # Secondary map for variations within a plot
    cfg2$map2 <- data.frame(
      tid = c(1:19),
      map_row = c(1:19)
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
  
  # Set config based on cfg2$map and cfg2$map2
  i <- cfg2$map2[cfg2$tid,"map_row"]
  for (x in c("marker", "lab_x", "lab_title", "day", "dataset", "cr2_trial",
              "cr2_COR", "cr2_marker", "t_0", "zoom_x", "zoom_y_cve",
              "zoom_y_risk", "more_ticks", "llox_label", "llox")) {
    cfg2[[x]] <- cfg2[[x]][[cfg2$map[i,x]]]
  }
  for (x in c("id", "time", "event", "wt", "ph1", "ph2", "covariates")) {
    cfg2$v[[x]] <- cfg2$v[[x]][[cfg2$map[i,paste0("v_",x)]]]
  }
  cfg2$params$edge_corr <- cfg2$edge_corr[cfg2$map[i,"edge_corr"]]
  cfg2$v$covariates <- formula(cfg2$v$covariates)
  
  # Override certain values with map2
  # cfg2$params$ci_type <- cfg2$map2[cfg2$tid,"ci_type"] # !!!!!
  # cfg2$params$edge_ind <- cfg2$map2[cfg2$tid,"edge_ind"] # !!!!!
  
  # AMP-specific code
  if (cfg2$analysis=="AMP") {
    cfg2[["amp_protocol"]] <-
      cfg2[["amp_protocol"]][[cfg2$map[i,"amp_protocol"]]]
    cfg2[["amp_tx"]] <- cfg2[["amp_tx"]][[cfg2$map[i,"amp_tx"]]]
  }
  
  # Moderna-specific code
  if ((i %in% c(5,7,9)) && cfg2$analysis=="Moderna") {
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
  
  # Create data frames specific to treatment and control groups
  if (cfg2$txct) {
    
    tx_rows <- which(df_ph1_full$Trt==1)
    ct_rows <- which(df_ph1_full$Trt==0)
    df_tx <- df_ph1_full[tx_rows,]
    df_ct <- df_ph1_full[ct_rows,]

  }

  # Create data structures to hold results
  plot_data_risk <- data.frame(
    x = double(),
    y = double(),
    curve = character(),
    ci_lo = double(),
    ci_hi = double(),
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
  
  # Create C$t_0 variable if it is not provided
  if (length(df_tx[["SubcohortInd"]])!=length(dat$v$z)) {
    stop("Error; lengths differ between df_tx and dat$v.")
  }
  SubcohortInd <- df_tx[["SubcohortInd"]]
  indices_1 <- which(dat$v$z==1 & dat$v$delta==1)
  indices_2 <- which(dat$v$z==1 & SubcohortInd==1)
  time_1 <- max(dat$v$y[indices_1])
  time_2 <- sort(dat$v$y[indices_2], decreasing=T)[15] - 1
  # print(paste("# of PH2 events:", length(indices_1))) # QA
  s_num <- sum(dat$v$s==min(dat$v$s, na.rm=T), na.rm=T)
  s_den <- sum(!is.na(dat$v$s))
  # print(paste("Edge mass:", round(s_num/s_den,2))) # QA
  if (cfg2$analysis=="Janssen (partA)") {
    C$t_0 <- min(time_1,time_2)
  } else {
    if (cfg2$t_0==0) {
      C$t_0 <- max(dat$v$y[dat$v$z==1 & dat$v$delta==1])
    } else {
      C$t_0 <- cfg2$t_0
    }
  }
  
  # Generate grid of points
  s_grid <- seq(
    from = min(dat$v$s, na.rm=T),
    to = max(dat$v$s, na.rm=T),
    length.out = 101
  )
  
}



###################################.
##### Overall ests of risk/VE #####
###################################.

if (cfg2$estimators$overall %in% c("Cox gcomp", "KM")) {
  
  if (cfg2$estimators$overall=="Cox gcomp") {
    method <- "Cox"
  } else if (cfg2$estimators$overall=="KM") {
    method <- "KM"
  } else {
    stop("cfg2$estimators$overall incorrectly specified.")
  }
  
  ests_ov <- vaccine::overall(dat=dat, t_0=C$t_0, method=method)
  
  if ("CVE" %in% cfg2$plots) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = c(999,999),
      y = rep(ests_ov[ests_ov$stat=="ve","est"], 2),
      curve = rep("Overall VE", 2),
      ci_lo = rep(ests_ov[ests_ov$stat=="ve","ci_lo"], 2),
      ci_hi = rep(ests_ov[ests_ov$stat=="ve","ci_hi"], 2),
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
      ci_lo = c(rep(ests_ov[ests_ov$group=="placebo","ci_lo"], 2),
                rep(ests_ov[ests_ov$group=="vaccine","ci_lo"], 2)),
      ci_hi = c(rep(ests_ov[ests_ov$group=="placebo","ci_hi"], 2),
                rep(ests_ov[ests_ov$group=="vaccine","ci_hi"], 2)),
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
    # if (length(xx) %in% c(2,3)) {
    #   for (i in 2:length(xx)) {
    #     x=xx[i-1]
    #     if (x>=3) {
    #       label <- bquote(3%*%10^.(x))
    #     } else {
    #       label <- 3*10^x
    #     }
    #     x_axis$ticks[length(x_axis$ticks)+1] <- x+log10(3)
    #     x_axis$labels[[length(x_axis$labels)+1]] <- label
    #   }
    # }
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
      ci_hi = rep(overall.ve[[3]], 2),
      overall = c("Overall L", "Overall R")
    ))
  }
  
  if (cfg2$estimators$cr=="Cox import" && "CVE" %in% cfg2$plots) {
    plot_data_cve <- rbind(plot_data_cve, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(1-risks$prob/res.plac.cont["est"]),
      curve = rep("CVE, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_cve_cis[1,]),
      ci_hi = as.numeric(cox_cve_cis[2,]),
      overall = rep("", length(as.numeric(risks$marker)))
    ))
  }
  
  if (cfg2$estimators$overall=="Cox import" && "Risk" %in% cfg2$plots) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = rep(999,4),
      y = c(rep(prev.plac["est"], 2), rep(prev.vacc["est"], 2)),
      curve = c(rep("Placebo overall",2), rep("Vaccine overall",2)),
      ci_lo = c(rep(prev.plac["2.5%"], 2), rep(prev.vacc["2.5%"], 2)),
      ci_hi = c(rep(prev.plac["97.5%"], 2), rep(prev.vacc["97.5%"], 2)),
      overall = rep(c("Overall L", "Overall R"),2)
    ))
  }
  
  if (cfg2$estimators$cr=="Cox import" && "Risk" %in% cfg2$plots) {
    plot_data_risk <- rbind(plot_data_risk, data.frame(
      x = as.numeric(risks$marker),
      y = as.numeric(risks$prob),
      curve = rep("Risk, Cox model", length(as.numeric(risks$marker))),
      ci_lo = as.numeric(cox_risk_cis[1,]),
      ci_hi = as.numeric(cox_risk_cis[2,]),
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
    ci_lo_risk <- ests$cr$ci_lo %>% pmax(0) %>% pmin(1)
    ci_hi_risk <- ests$cr$ci_hi %>% pmax(0) %>% pmin(1)

    # Compute CVE estimates
    if (run_cve) {
      ests_cve <- ests$cve$est
      ci_lo_cve <- ests$cve$ci_lo %>% pmin(1) # Reversing is intentional
      ci_hi_cve <- ests$cve$ci_hi %>% pmin(1) # Reversing is intentional
    }
    
    plot_data_risk <- data.frame(
      x = s_grid,
      y = ests_risk,
      curve = rep(lab_risk, length(ests_risk)),
      ci_lo = ci_lo_risk,
      ci_hi = ci_hi_risk,
      overall = rep("", length(ests_risk))
    )
    if (run_cve) {
      plot_data_cve <- data.frame(
        x = s_grid,
        y = ests_cve,
        curve = rep(lab_cve, length(ests_cve)),
        ci_lo = ci_lo_cve,
        ci_hi = ci_hi_cve,
        overall = rep("", length(ests_cve))
      )
    } else {
      plot_data_cve <- NA
    }
    
    return(list(risk=plot_data_risk, cve=plot_data_cve))
    
  }
  
}



#################################.
##### Data analysis (NPCVE) #####
#################################.

if ("Grenander" %in% cfg2$estimators$cr) {
  
  # Obtain estimates
  return_extra <- c()
  if (cfg2$run_debug$objs) {
    return_extra <- c(return_extra, "omega_n", "f_sIx_n", "Q_n", "grid",
                      "Phi_n", "Gamma_os_n", "gcm", "dGCM", "r_Mn_Gr")
  }
  if (cfg2$run_debug$gren_var) {
    return_extra <- c(return_extra, "deriv_r_Mn", "f_s_n", "gamma_n")
  }

  calc_ests <- T
  if (calc_ests) {
    
    ests <- vaccine::est_np(
      dat = dat,
      t_0 = C$t_0,
      s_out = s_grid,
      edge_corr = as.logical(cfg2$params$edge_corr=="min"),
      ci_type = cfg2$params$ci_type,
      placebo_risk_method = "Cox",
      params = list(surv_type = cfg2$params$Q_n_type,
                    density_type = cfg2$params$g_n_type,
                    deriv_type = cfg2$params$deriv_type,
                    q_n_type = cfg2$params$q_n_type),
      grid_size = list(y=101, s=101, x=5)
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_g_",cfg2$tid,".rds"))
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



##################################################.
##### !!!!! Janssen mediation analysis !!!!! #####
##################################################.

if (F) {
  
  # NDE
  est_risk <- ests$r_Mn_edge_est
  var_risk <- ests$sigma2_edge_est
  n <- ests$n
  ci_lo_risk <- est_risk - 1.96*sqrt(var_risk/n)
  ci_hi_risk <- est_risk + 1.96*sqrt(var_risk/n)
  risk_ct <- ests_ov[ests_ov$group=="placebo","est"]
  cve <- Vectorize(function(x) { 1 - x/risk_ct })
  est_cve <- cve(est_risk)
  sd_cve <- sqrt(var_risk) / risk_ct
  ci_lo_cve <- cve(ci_hi_risk) %>% pmin(1)
  ci_hi_cve <- cve(ci_lo_risk) %>% pmin(1)
  
  # PM
  est_ov <- overall.ve[[1]]
  sd_ov <- sd(1-res.vacc.cont[2:1001]/res.plac.cont[2:1001])
  ci_lo_ov <- overall.ve[[2]]
  ci_hi_ov <- overall.ve[[3]]
  est_pm <- 1 - log(1-est_cve)/log(1-est_ov)
  var_pm <- sd_cve^2 / ( n * (1-est_cve)^2 * (log(1-est_ov))^2 ) +
    sd_ov * (log(1-est_cve))^2 / ((1-est_ov)^2+(log(1-est_ov))^4)
  sd_pm <- sqrt(var_pm)
  ci_lo_pm <- est_pm - 1.96*sd_pm
  ci_hi_pm <- est_pm + 1.96*sd_pm
  
  # Save results
  saveRDS(
    list(est_cve=est_cve, sd_cve=sd_cve, ci_lo_cve=ci_lo_cve,
         ci_hi_cve=ci_hi_cve, est_pm=est_pm, sd_pm=sd_pm, ci_lo_pm=ci_lo_pm,
         ci_hi_pm=ci_hi_pm),
    paste0(cfg2$analysis," plots/mediation_results_",cfg2$tid,".rds")
  )
  
  # Process results
  if (F) {
    
    df_med <- data.frame(
      "est_cve" = double(),
      "sd_cve" = double(),
      "ci_lo_cve" = double(),
      "ci_hi_cve" = double(),
      "est_pm" = double(),
      "sd_pm" = double(),
      "ci_lo_pm" = double(),
      "ci_hi_pm" = double()
    )
    for (i in c(1:58)) {
      file <- paste0("Janssen (partA) plots/Mediation/",
                     "/mediation_results_",i,".rds")
      med <- readRDS(file)
      df_med[round(nrow(df_med)+1),] <- med
    }
    write.table(df_med, file="mediation_results.csv", sep=",", row.names=F)
    
  }
  
}



###################################.
##### Data analysis (Cox GAM) #####
###################################.

if ("Cox GAM" %in% cfg2$estimators$cr) {
  
  calc_ests <- T
  if (calc_ests) {
    
    ests <- vaccine::est_cox(
      dat = dat,
      t_0 = C$t_0,
      s_out = s_grid,
      ci_type = "logit",
      placebo_risk_method = "Cox",
      spline_df = 4
    )
    
    if (flags$save_data_objs) {
      saveRDS(ests, paste0(cfg2$analysis," plots/ests_z_",cfg2$tid,".rds"))
    }
    
  } else {
    
    ests <- readRDS(paste0(cfg2$analysis," plots/ests_z_",cfg2$tid,".rds"))
    
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
    
    ests <- vaccine::est_cox(
      dat = dat,
      t_0 = C$t_0,
      s_out = s_grid,
      ci_type = "logit",
      placebo_risk_method = "Cox"
      # return_extras = T
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
    
    ests <- vaccine::est_cox(
      dat = dat,
      t_0 = C$t_0,
      s_out = s_grid,
      ci_type = "logit",
      placebo_risk_method = "Cox",
      edge_ind = T
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
##### Data analysis (hypothesis test) #####
###########################################.

if (cfg2$run_hyptest) {
  
  test_results <- test_2(
    dat_orig = dat$v,
    alt_type = "decr",
    # alt_type = "two-tailed",
    params = list(
      type = c("simple (with constant)", "edge", "combined", "combined 2"),
      q_n_type = "zero",
      Q_n_type = "Super Learner"
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
  #'   - `ci_hi`: Upper confidence bound for CVE
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
  #' @return A ggplot2 plot object
  #' @notes
  #'   - This function plots pointwise estimates and confidence intervals from
  #'     both Cox model and nonparametric approaches
  create_plot <- function(plot_data, which, zoom_x=NA, zoom_y=NA, zoom_y_max=NA,
                          labs, hst, rr_y_axis=F, log10_x_axis=F) {
    
    # Change curve labels to factors and set color scale
    curves <- c(
      "Placebo overall", "Vaccine overall", "Overall VE",
      "Risk, Cox model", "CVE, Cox model", "Risk, Qbins", "CVE, Qbins",
      "Risk, Cox GAM", "CVE, Cox GAM",
      "Risk, nonparametric", "CVE, nonparametric",
      "Risk, Cox (analytic)", "CVE, Cox (analytic)",
      "Control", "VRC01 10mg/kg", "VRC01 30mg/kg", "VRC01 Pooled",
      "Risk, Cox (basic)", "CVE, Cox (basic)",
      "Risk, Cox (spline 4 df)", "CVE, Cox (spline 4 df)",
      "Risk, Cox (edge)", "CVE, Cox (edge)"
    )
    curve_colors <- c(
      "darkgrey", "darkgrey", "darkgrey",
      "darkorchid3", "darkorchid3", "firebrick3", "firebrick3", # !!!!!
      # "firebrick3", "firebrick3", "firebrick3", "firebrick3", # !!!!!
      "darkgreen", "darkgreen",
      "deepskyblue3", "deepskyblue3",
      "deepskyblue3", "deepskyblue3",
      "deepskyblue3", "darkorchid3", "firebrick3", "darkolivegreen3",
      "darkorchid3", "darkorchid3",
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
    if (is.na(zoom_y[1])) {
      zoom_y <- c(0,1)
      zoom_y[2] <- zoom_y[2] + 0.05*(zoom_y[2]-zoom_y[1])
    } else if (zoom_y[1]=="zoomed") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- min(zz$ci_lo, na.rm=T)
      z_y_U <- max(zz$ci_hi, na.rm=T)
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
    } else if (zoom_y[1]=="zoomed (risk)") {
      zz <- dplyr::filter(plot_data, x>=zoom_x[1] & x<=zoom_x[2])
      z_y_L <- 0
      z_y_U <- max(plot_data$ci_hi, na.rm=T)
      if (!is.null(zoom_y_max)) { z_y_U <- min(zoom_y_max, z_y_U) }
      zoom_y <- c(z_y_L - 0.05*(z_y_U-z_y_L),
                  z_y_U + 0.05*(z_y_U-z_y_L))
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
    
    # Hack to get zoom_y value for HVTN 124 plots
    if (flags$hvtn124_plot) {
      zoom_x <<- zoom_x
      zoom_y <<- zoom_y
    }
    
    # Create and return ggplot2 object
    # Note: using geom_rect for the histogram so that it can be shifted up/down
    if (rr_y_axis) {
      sec.axis <- sec_axis(~1-., breaks=seq(0,2,0.1),
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
      geom_line(linewidth=0.7) +
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
      xlim <- c(min(dat$v$s, na.rm=T), max(dat$v$s, na.rm=T))
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
  trim_plot_data <- function(plot_data) {
    
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
    plot_data[rows, c("y", "ci_lo", "ci_hi")] <- NA
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
      ind_hi_cve <- max(which(plot_data_cve$ci_lo<=-2))
      plot_data_cve[c(ind_lo_cve:ind_hi_cve), c("y", "ci_lo", "ci_hi")] <- NA
      ind_lo_risk <- min(which(plot_data_risk$curve=="Risk, nonparametric"))
      ind_hi_risk <- max(which(plot_data_risk$ci_hi>=0.5))
      plot_data_risk[c(ind_lo_risk:ind_hi_risk), c("y", "ci_lo", "ci_hi")] <- NA
    } else {
      stop("HVTN 705 config changed")
    }
  }
  
  if (flags$paper_npcve || flags$paper_cox) { cfg2$lab_title <- NULL }
  
  if (nrow(plot_data_risk)>0) {
    
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", C$t_0)
    
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_risk),
      which = "Risk",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y_risk,
      zoom_y_max = cfg2$zoom_y_max,
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      log10_x_axis = T
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
      write.table(trim_plot_data(plot_data_risk),
                  file=paste0(cfg2$analysis," plots/risk_",cfg2$tid,".csv"),
                  sep=",",
                  row.names=FALSE)
    }
    
  }
  
  if (nrow(plot_data_cve)>0) {
    
    cfg2$lab_y <- paste0("Controlled VE against ", cfg2$endpoint,
                         " by day ", C$t_0)
    
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
    
    plot <- create_plot(
      plot_data = trim_plot_data(plot_data_cve),
      # plot_data = trim_plot_data(filter(plot_data_cve, curve!="CVE, Cox model")),
      which = "CVE",
      zoom_x = cfg2$zoom_x,
      zoom_y = cfg2$zoom_y_cve,
      labs = list(title=cfg2$lab_title, x=cfg2$lab_x, y=cfg2$lab_y),
      hst = hst,
      rr_y_axis = T,
      log10_x_axis = T
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
    
    if (flags$table_of_vals) {
      write.table(trim_plot_data(plot_data_cve),
                  file=paste0(cfg2$analysis," plots/cve_",cfg2$tid,".csv"),
                  sep=",",
                  row.names=FALSE)
    }
    
  }
  
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
    linewidth = 0.4,
    fill = "brown4",
    alpha = 0.1,
    inherit.aes = F
  ) + coord_cartesian(xlim=zoom_x, ylim=zoom_y, expand=F)
  
}



######################################################.
##### Generate and save AMP plots (run manually) #####
######################################################.

if (F) {
  
  # Temp
  C$t_0 <- 595
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
      risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_0))]
      ci_lo_ov <- 1 - srv_ov$upper[which.min(abs(srv_ov$time-C$t_0))]
      ci_hi_ov <- 1 - srv_ov$lower[which.min(abs(srv_ov$time-C$t_0))]
      km <- data.frame(
        x = c(999,999),
        y = rep(risk_ov, 2),
        curve = rep("Overall risk", 2),
        ci_lo = rep(ci_lo_ov, 2),
        ci_hi = rep(ci_hi_ov, 2),
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
    cfg2$lab_y <- paste0("Probability of ", cfg2$endpoint, " by day ", C$t_0)
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
  mrk_lab <- "Expanded Multi-epitope functions: Day 210"
  title_lab <- "Figure 1-5. HVTN 705 Expanded Multi-epitope functions at Month 7"
  
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



###############################################.
##### ARCHIVE: Summary stats / DQA checks #####
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
  num_case_tx_t_0 <- sum(ind_tx[time_tx<=C$t_0])
  num_case_ct_t_0 <- sum(ind_ct[time_ct<=C$t_0])
  num_atrisk_tx <- length(ind_tx)
  num_atrisk_ct <- length(ind_ct)
  print(paste0("Number of cases in vaccine group: ", num_case_tx))
  print(paste0("Number of cases in control group: ", num_case_ct))
  print(paste0("Number of cases by day ", C$t_0, " in vaccine group: ",
               num_case_tx_t_0))
  print(paste0("Number of cases by day ", C$t_0, " in control group: ",
               num_case_ct_t_0))
  print(paste0("Number at-risk in vaccine group: ", num_atrisk_tx))
  print(paste0("Number at-risk in control group: ", num_atrisk_ct))
  print(paste0("Naive P(COVID by day ", C$t_0, ") in vaccine group: ",
               round(num_case_tx_t_0/num_atrisk_tx,3)))
  print(paste0("Naive P(COVID by day ", C$t_0, ") in control group: ",
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

if (cfg2$run_debug$objs) {
  
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
  #   ests$Q_n(t=C$t_0, x=c(0,0), s)
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

if (cfg2$run_debug$gren_var) {
  
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
    rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_0))]
    ci_lo_ct <- 1 - srv_ct$upper[which.min(abs(srv_ct$time-C$t_0))]
    ci_hi_ct <- 1 - srv_ct$lower[which.min(abs(srv_ct$time-C$t_0))]
    var_ct <- ((ci_hi_ct-ci_lo_ct)/3.92)^2
    
    # Calculate treatment group survival (KM; with SE)
    srv_tx <- survfit(
      Surv(EventTimePrimaryIncludeNotMolecConfirmedD29,
           EventIndPrimaryIncludeNotMolecConfirmedD29)~1,
      data = df_tx
    )
    rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-C$t_0))]
    ci_lo_tx <- 1 - srv_tx$upper[which.min(abs(srv_tx$time-C$t_0))]
    ci_hi_tx <- 1 - srv_tx$lower[which.min(abs(srv_tx$time-C$t_0))]
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
    rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_tx_sub$time-C$t_0))]
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
