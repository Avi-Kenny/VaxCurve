
#################.
##### Setup #####
#################.

local <- FALSE
recompute_fns <- FALSE



###########################.
##### Data processing #####
###########################.

if (local) {
  
  # Read in raw data; DO NOT STORE THIS LOCALLY !!!!!
  df_raw <- read.csv(paste0("Z:/covpn/p3001/download_data/Moderna COVE mRNA 1273",
                            "P301_immune_20210915/moderna_real_data_processed_with_",
                            "riskscore.csv"))
  
  # df_raw_copy <- df_raw
  # df_raw <- df_raw_copy
  
  # Filter out placebo patients
  df_raw %<>% filter(Trt==1)
  
  # Filter out records with missing weights
  # !!!!! Check why some weights are missing
  df_raw %<>% filter(!is.na(df_raw$wt.D57))
  
  # Create data structure for analysis
  dat_orig <- list(
    "id" = df_raw$Ptid,
    "y_star" = df_raw$EventTimePrimaryD57,
    "delta_star" = df_raw$EventIndPrimaryD57,
    "a" = df_raw$Day57pseudoneutid50,
    "w" = data.frame(
      "w1" = df_raw$MinorityInd,
      "w2" = df_raw$standardized_risk_score,
      "w3" = df_raw$HighRiskInd
    ),
    "weights" = df_raw$wt.D57,
    "delta" = as.integer(!is.na(df_raw$Day57pseudoneutid50))
    # "delta" = as.integer(df_raw$SubcohortInd) # TwophasesampIndD57
  )
  
  # Stabilize weights (rescale to sum to sample size)
  dat_orig$weights <- ifelse(dat_orig$delta==1, dat_orig$weights, 0)
  s <- sum(dat_orig$weights) / length(dat_orig$delta)
  dat_orig$weights <- dat_orig$weights / s

  saveRDS(dat_orig, file="moderna/dat_orig_moderna.rds")
  
} else {
  
  dat_orig <- readRDS("moderna/dat_orig_moderna.rds")
  
}



###############################.
##### Data quality checks #####
###############################.

if (local) {
  
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

# !!!!! Eventually, replace this with a call to est_curve()

# Set up end time of interest
C <- list(t_e=100, appx=cfg$appx)

# Rescale A to lie in [0,1]
a_scale <- ceiling(10*max(dat_orig$a, na.rm=T))/10
dat_orig$a <- dat_orig$a / a_scale

# Round A, W2
# !!!!! Temporary; later use round(...,1)
dat_orig$a <- round(dat_orig$a, 2)
dat_orig$w$w2 <- round(dat_orig$w$w2, 0)

# Create truncated data object
dat <- ss(dat_orig, which(dat_orig$delta==1))

# Create val_list
vlist <- create_val_list(dat, C$appx)
vlist$AW_grid <- NA
vlist$omega <- NA
vlist$W_grid <- NA

# Set estimation tuning parameters
params <- list(S_n_type="Cox PH", g_n_type="binning",
               ecdf_type="linear (mid)", deriv_type="m-spline",
               gamma_type="kernel", ci_type="regular")

# Construct component functions; save individually
if (recompute_fns) {
  
  # Construct/save component functions related to point estimation
  Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
  Phi_n_inv <- construct_Phi_n(dat, which="inverse", type=params$ecdf_type)
  S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type) # vals=NA
  Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE) # vals=NA
  f_aIw_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type, k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  g_n <- construct_g_n(f_aIw_n, f_a_n)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  fns <- c("Phi_n","Phi_n_inv","S_n","Sc_n","f_aIw_n","f_a_n","g_n","omega_n")
  for (fn in fns) { saveRDS(eval(as.name(fn)), paste0("moderna/", fn, ".rds")) }
  
  # Construct/save Gamma_n (one-step and plug-in)
  Gamma_n <- construct_Gamma_n(dat, vlist$A_grid, omega_n, S_n, g_n)
  Gamma_n_pi <- construct_Gamma_n(dat, vlist$A_grid, omega_n, S_n, g_n,
                               type="plug-in")
  fns <- c("Gamma_n","Gamma_n_pi")
  for (fn in fns) { saveRDS(eval(as.name(fn)), paste0("moderna/", fn, ".rds")) }
  
  # Construct/save additional functions
  Psi_n <- Vectorize(function(x) {
    Gamma_n(round(Phi_n_inv(x), -log10(C$appx$a))) # Gamma_n(Phi_n_inv(x))
  })
  gcm <- gcmlcm(x=seq(0,1,0.0001), y=Psi_n(seq(0,1,0.0001)), type="lcm")
  dGCM <- Vectorize(function(x) {
    index <- which(round(x,5)<=gcm$x.knots)[1]-1
    if (index==0) { index <- 1 }
    return(gcm$slope.knots[index])
  })
  theta_n_Gr <- function(x) { dGCM(Phi_n(x)) }
  theta_n <- theta_n_Gr
  fns <- c("Psi_n","gcm","dGCM","theta_n_Gr","theta_n")
  for (fn in fns) { saveRDS(eval(as.name(fn)), paste0("moderna/", fn, ".rds")) }
  
  # Component functions related to variance estimation
  f_aIw_delta1_n <- construct_f_aIw_n(dat, vlist$AW_grid, type=params$g_n_type,
                                      k=15, delta1=TRUE)
  f_a_delta1_n <- construct_f_a_n(dat_orig, vlist$A_grid,
                                  f_aIw_delta1_n)
  gamma_n <- construct_gamma_n(dat_orig, dat, vlist$A_grid,
                               type=params$gamma_type, omega_n, f_aIw_n,
                               f_a_n, f_a_delta1_n)
  deriv_theta_n <- construct_deriv_theta_n(theta_n, type=params$deriv_type,
                                           dir="decr")
  tau_n <- construct_tau_n(deriv_theta_n, gamma_n, f_a_n)
  fns <- c("f_aIw_delta1_n","f_a_delta1_n","gamma_n","deriv_theta_n","tau_n")
  for (fn in fns) { saveRDS(eval(as.name(fn)), paste0("moderna/", fn, ".rds")) }
  
  # G-computation estimator for comparison
  gcomp <- construct_gcomp_n(dat_orig, vlist$A_grid, S_n=S_n)
  fns <- c("gcomp")
  for (fn in fns) { saveRDS(eval(as.name(fn)), paste0("moderna/", fn, ".rds")) }
  
} else {
  
  fns <- c("Phi_n","Phi_n_inv","S_n","Sc_n","f_aIw_n","f_a_n","g_n","omega_n",
           "Gamma_n","Gamma_n_pi","Psi_n","gcm","dGCM","theta_n_Gr","theta_n",
           "f_aIw_delta1_n","f_a_delta1_n","gamma_n","deriv_theta_n","tau_n",
           "gcomp")
  for (fn in fns) {
    assign(fn, readRDS(paste0("moderna/", fn, ".rds")))
  }
  
}

if (local) {
  
  # CVE graph
  # !!!!! plot 10% quantile cutoffs
  {
    # Generate point estimates
    grid <- seq(0,1,0.01)
    ests_gcomp <- 1 - gcomp(grid)/0.061
    ests_gren <- 1 - theta_n(grid)/0.061
    
    # Generate CIs
    tau_ns <- tau_n(grid)
    qnt <- 1.00
    n_orig <- length(dat_orig$delta)
    if (params$ci_type=="regular") {
      ci_lo <- ests_gren - (qnt*tau_ns)/(n_orig^(1/3))
      ci_hi <- ests_gren + (qnt*tau_ns)/(n_orig^(1/3))
    } else if (params$ci_type=="logit") {
      ci_lo <- expit(
        logit(ests_gren) - (qnt*tau_ns*deriv_logit(ests_gren))/(n_orig^(1/3))
      )
      ci_hi <- expit(
        logit(ests_gren) + (qnt*tau_ns*deriv_logit(ests_gren))/(n_orig^(1/3))
      )
    }
    
    # Marginal distribution of A
    df_marg <- data.frame(
      x = grid*a_scale,
      ymin = 0,
      ymax = f_a_n(grid)*0.2
    )
    
    # Plot
    # Export: 6" x 5"
    ggplot(
      data.frame(
        x = rep(grid*a_scale,2),
        y = c(ests_gcomp,ests_gren),
        which = rep(c("G-comp", "Grenander"), each=length(grid)),
        ci_lo = c(ests_gcomp,ci_lo),
        ci_hi = c(ests_gcomp,ci_hi)
      ),
      aes(x=x, y=y, color=which)
    ) +
      geom_ribbon(
        aes(ymin=ci_lo, ymax=ci_hi),
        alpha = 0.2,
        linetype = "dotted",
        fill = "darkblue"
      ) +
      geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                  data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
      scale_color_manual(values=c("purple", "darkblue")) +
      scale_y_continuous(labels=label_percent(accuracy=1), limits=c(0,1.05),
                         breaks=seq(0,1,0.1)) +
      theme(legend.position="bottom") +
      labs(x="Day 57 pseudoneut-ID50", y="Controlled vaccine efficacy",
           color="Estimator") +
      geom_line()
      
  }
  
  # Gamma_n graph
  # Export: 6" x 5"
  grid <- seq(0,1,0.01)
  df_marg <- data.frame(
    x = grid*a_scale,
    ymin = 0,
    ymax = f_a_n(grid)*0.001
  )
  ggplot(
    data.frame(
      x = grid*a_scale,
      y = Gamma_n(grid)
    ),
    aes(x=x, y=y)
  ) +
    geom_ribbon(aes(x=x, ymin=ymin, ymax=ymax), inherit.aes=F,
                data=df_marg, fill="forestgreen", color=NA, alpha=0.3) +
    theme(legend.position="bottom") +
    labs(x="Day 57 pseudoneut-ID50", y="Gamma_n",
         color="Estimator") +
    geom_line()
  
  # theta graph
  # Export: 6" x 5"
  grid <- seq(0,1,0.01)
  df_marg <- data.frame(
    x = grid*a_scale,
    ymin = 0,
    ymax = f_a_n(grid)*0.01
  )
  ests_gcomp <- gcomp(grid)
  ests_gren <- theta_n(grid)
  ggplot(
    data.frame(
      x = rep(grid*a_scale,2),
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
