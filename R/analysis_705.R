
# !!!!! Need to update this file based on analysis_janssen.R

#################.
##### Setup #####
#################.

print(paste("Check 0:",Sys.time()))

# Vector of markers
markers <- c("Day210ELCZ", "Day210ELMo", "Day210ADCPgp140C97ZAfib",
             "Day210IgG3gp140C97ZAfibritin40delta",
             "Day210IgG3gp140Mos1fibritin40delta", "Day210IgG340mdw_gp120",
             "Day210IgG340mdw_gp140", "Day210IgG340mdw_V1V2",
             "Day210IgG3gp4140delta", "Day210IgG340mdw_multi",
             "Day210IgG340mdw_gp120_gp140_vm", "Day210mdw_xassay")

# Only focus on five "primary" markers
# marker_num will be in c(1:5), referring to this shortened list
markers <- markers[c(1,3,8,11,12)]

# Set configuration
cfg2 <- list(
  # tid = 1, # !!!!!
  tid = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")),
  S_n_type = "Cox PH", # !!!!! "Cox PH" "Super Learner"
  run_import_data = F,
  run_recomp_1 = T,
  run_recomp_2 = T,
  run_recomp_3 = T,
  run_recomp_4 = T,
  run_recomp_5 = T,
  save_fns = F,
  run_dqa = F,
  run_graphs = T
)
cfg2$folder = paste0("705 functions/705 (tid ",cfg2$tid,")")

# !!!!!
if (cfg2$tid %in% c(1:5)) {
  cfg2$marker_num <- cfg2$tid
  cfg2$S_n_type <- "Cox PH"
} else if (cfg2$tid %in% c(6:10)) {
  cfg2$marker_num <- cfg2$tid-5
  cfg2$S_n_type <- "Super Learner"
} else {
  stop("Error with SLURM_ARRAY_TASK_ID")
}
cfg2$edge_spread <- ifelse(cfg2$marker_num==3, TRUE, FALSE)
cfg2$which_marker <- markers[cfg2$marker_num]

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

  # Read in raw data; DO NOT STORE THIS LOCALLY !!!!!
  df_raw <- read.csv(paste0("Z:/vaccine/p705/analysis/lab/cc/copcor/HVTN705_fi",
                            "rstcasecontrolprocesseddata.csv"))
  
  # Subset data frames
  df_ctrl <- filter(df_raw, Trt==0)
  df_trt <- filter(df_raw, Trt==1)
  
  # Create data structure for analysis
  dat_orig <- list(
    "id" = df_trt[["Subjid"]],
    "y_star" = df_trt[["Ttilde.D210"]],
    "delta_star" = df_trt[["Delta.D210"]],
    "w" = data.frame(
      "w1" = df_trt[["RSA"]],
      "w2" = df_trt[["Riskscore"]],
      "w3" = df_trt[["BMI"]],
      "w4" = df_trt[["Age"]]
    ),
    "weights" = df_trt[["wt.D210"]],
    "a_list" = list(df_trt[[markers[1]]], df_trt[[markers[2]]],
                    df_trt[[markers[3]]], df_trt[[markers[4]]],
                    df_trt[[markers[5]]]),
    "delta" = as.integer(df_trt[["Ph2ptids.D210"]])
  )

  # Stabilize weights (rescale to sum to sample size)
  dat_orig$weights <- ifelse(dat_orig$delta==1, dat_orig$weights, 0)
  s <- sum(dat_orig$weights) / length(dat_orig$delta)
  dat_orig$weights <- dat_orig$weights / s

  saveRDS(dat_orig, file="705 data/dat_orig_705.rds")
  saveRDS(df_ctrl, file="705 data/df_ctrl_705.rds")

} else {

  dat_orig <- readRDS("705 data/dat_orig_705.rds")
  df_ctrl <- readRDS("705 data/df_ctrl_705.rds")

}



###############################.
##### Data quality checks #####
###############################.

if (cfg2$run_dqa) {
  # !!!!!
}



##############################.
##### Data analysis prep #####
##############################.

{
  
  # Set up end time of interest `t_e`
  C <- list(appx=cfg$appx, t_e=550)

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

  # Round A, W2
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  dat_orig$w$w2 <- round(dat_orig$w$w2, 0) # !!!!!
  dat_orig$w$w3 <- round(dat_orig$w$w3, 0) # !!!!!

  # Perform "spread" edge correction
  if (cfg2$edge_spread) {
    width <- 0.02
    noise_0 <- round(runif(length(dat_orig$a))*width, -log10(C$appx$a))
    noise_1 <- round(1 - runif(length(dat_orig$a))*width, -log10(C$appx$a))
    dat_orig$a <- ifelse(dat_orig$a==0, noise_0, dat_orig$a)
    dat_orig$a <- ifelse(dat_orig$a==1, noise_1, dat_orig$a)
    # dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  }

  # Set estimation tuning parameters
  params <- list(S_n_type=cfg2$S_n_type, g_n_type="binning",
                 ecdf_type="linear (mid)", deriv_type="m-spline",
                 gamma_type="kernel", ci_type="trunc")
  
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
  return_extra = c("gcomp")
)

# !!!!! continue (see Janssen)
asdf

# Calculate control group survival
srv_ct <- survfit(
  Surv(Ttilde.D210,Delta.D210)~1,
  data = df_ctrl
)
rate_ct <- 1 - srv_ct$surv[which.min(abs(srv_ct$time-C$t_e))]

# Calculate VE
calc_ve <- FALSE
if (calc_ve) {
  
  # VE, overall
  srv_tx <- survfit(
    Surv(Ttilde.D210,Delta.D210)~1,
    data = df_trt
  )
  rate_tx <- 1 - srv_tx$surv[which.min(abs(srv_tx$time-C$t_e))]
  ve_overall <- 1 - (rate_tx/rate_ct)
  print(ve_overall)
  
  # VE, within subcohort
  df_trt_sub <- filter(df_trt, Ph2ptids.D210==1)
  srv_tx_sub <- survfit(coxph(
    Surv(Ttilde.D210,Delta.D210)~1,
    data = df_trt_sub,
    weights = wt.D210
  ))
  rate_tx_sub <- 1 - srv_tx_sub$surv[which.min(abs(srv_ct$time-C$t_e))]
  ve_subcohort <- 1 - (rate_tx_sub/rate_ct)
  print(ve_subcohort)
  
}



##################.
##### Graphs #####
##################.

if (cfg2$run_graphs) {

  # CVE graph
  {
    
    # Extract results
    ests_gren <- rep(NA, length(p_grid))
    ci_lo <- rep(NA, length(p_grid))
    ci_hi <- rep(NA, length(p_grid))
    for (i in c(1:length(p_grid))) {
      ests_gren[i] <- ests[[i]]$est
      ci_lo[i] <- ests[[i]]$ci_lo
      ci_hi[i] <- ests[[i]]$ci_hi
    }
    ests_gcomp <- ests[["gcomp"]](p_grid)
    
    cve <- Vectorize(function(x) { 1 - x/rate_ct })
    
    # Generate CIs
    tau_ns <- tau_n(p_grid)
    qnt <- 1.00
    n_orig <- length(dat_orig$delta)
    if (params$ci_type=="regular") {
      ci_lo <- cve(ests_gren - qnt * (tau_ns/(n_orig^(1/3))))
      ci_hi <- cve(ests_gren + qnt * (tau_ns/(n_orig^(1/3))))
    } else if (params$ci_type=="logit") {
      ci_lo <- expit(
        logit(cve(ests_gren)) -
          qnt * (tau_ns*deriv_logit(cve(ests_gren)))/(rate_ct*n_orig^(1/3))
      )
      ci_hi <- expit(
        logit(cve(ests_gren)) +
          qnt * (tau_ns*deriv_logit(cve(ests_gren)))/(rate_ct*n_orig^(1/3))
      )
    } else if (params$ci_type=="trunc") {
      ci_lo <- cve(ests_gren - qnt * (tau_ns/(n_orig^(1/3)))) %>%
        # pmax(0) %>% pmin(1)
        pmin(1)
      ci_hi <- cve(ests_gren + qnt * (tau_ns/(n_orig^(1/3)))) %>%
        # pmax(0) %>% pmin(1)
        pmin(1)
    }
    ests_gcomp <- cve(ests_gcomp)
    ests_gren <- cve(ests_gren)

    # Marginal distribution of A
    df_marg <- data.frame(
      x = p_grid*a_scale-a_shift,
      ymin = 0,
      ymax = (1.6/max(f_a_n(p_grid))) * f_a_n(p_grid)
    )
    
    # Truncate at 5/95 or 10/90 quantiles
    # which <- p_grid>=Phi_n_inv(0.05) & p_grid<=Phi_n_inv(0.95)
    if (cfg2$marker_num==3) {
      which <- p_grid<=Phi_n_inv(0.9)
    } else {
      which <- p_grid>=Phi_n_inv(0.1) & p_grid<=Phi_n_inv(0.9)
    }
    ests_gren <- ifelse(which,ests_gren,NA)
    ci_lo <- ifelse(which,ci_lo,NA)
    ci_hi <- ifelse(which,ci_hi,NA)
    
    # Labels
    # x_labs <- c("IgG to VT-C (EU/ml)",
    #             "IgG to VT-M (EU/ml)",
    #             "Average phagocytosis score to gp140 C97ZA",
    #             "IgG3 Net MFI to gp140 C97ZA",
    #             "IgG3 Net MFI to gp140 Mosaic",
    #             "IgG3 gp120 breadth (Weighted average log10 Net MFI)",
    #             "IgG3 gp140 breadth (Weighted average log10 Net MFI)",
    #             "IgG3 V1V2 breadth (Weighted average log10 Net MFI)",
    #             "IgG3 Net MFI to gp41",
    #             "IgG3 multi-epitope breadth (Wt average log10 Net MFI)",
    #             "IgG3 gp120 + gp140 breadth (Wt average log10 Net MFI)",
    #             "Overall maximal diversity score")
    x_labs <- c("IgG Vx-VT-C (=s)", "ADCP Vx-C97ZA (=s)",
                "IgG3 V1V2 Breadth (=s)", "IgG3 Env Breadth (=s)",
                "Multi-Epitopes Functions (=s)")
    x_lab <- x_labs[cfg2$marker_num]
    y_lab <- "Controlled VE against HIV by day 550"
    
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
      geom_hline(yintercept=0.141, alpha=0.3, size=0.4) + # !!!!! Later calculate this manually
      geom_hline(yintercept=c(-0.193,0.401), alpha=0.3, size=0.4, linetype="dotted") + # !!!!! Later calculate this manually
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.1,
        linetype = "dotted",
        fill = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin-1, ymax=ymax-1), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_color_manual(values=c("purple", "darkblue")) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(-1,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      scale_x_continuous(label=math_format(10^.x)) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color="Estimator") +
      geom_line()
    
    plot_2 <- ggplot(
      data.frame(
        x = p_grid*a_scale-a_shift,
        y = ests_gren,
        ci_lo = ci_lo,
        ci_hi = ci_hi
      ),
      aes(x=x, y=y)
    ) +
      geom_hline(yintercept=0.141, alpha=0.3, size=0.4) + # !!!!! Later calculate this manually
      geom_hline(yintercept=c(-0.193,0.401), alpha=0.3, size=0.4, linetype="dotted") + # !!!!! Later calculate this manually
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.1,
        linetype = "dotted",
        fill = "darkblue",
        color = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin-1, ymax=ymax-1), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(-1,1),
                         breaks=seq(-1,1,0.1), minor_breaks=NULL) +
      theme(panel.grid.major=element_line(colour="white", size=0.3),
            panel.grid.minor=element_line(colour="white", size=0.3)) +
      scale_x_continuous(label=math_format(10^.x)) +
      theme(legend.position="bottom") +
      labs(x=x_lab, y=y_lab, color="Estimator") +
      geom_line(color="darkblue")
    
    # Save plots
    name_1 <- paste0("705 plots/plot1_",cfg2$tid,".pdf")
    name_2 <- paste0("705 plots/plot2_",cfg2$tid,".pdf")
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
      ymax = f_a_n(p_grid)*0.001
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
      ymax = f_a_n(p_grid)*0.01
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
