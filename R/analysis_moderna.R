
#################.
##### Setup #####
#################.

library(dplyr)
library(ggplot2)
library(magrittr)



#########################.
##### Data cleaning #####
#########################.


# Variables used:
#     Trt: randomization (1=vax, 0=placebo)
#     risk_score: (covariate)
#     standardized_risk_score: (covariate)
#     MinorityInd: (covariate)
#     EventTimePrimaryD1:
#     x:

# dat1_raw <- read.csv(paste0("Z:/covpn/p3001/download_data/Moderna COVE mRNA 1273P30",
#                        "1_immune_20210915/immune1.csv"))

dat_raw <- read.csv(paste0("Z:/covpn/p3001/download_data/Moderna COVE mRNA 1273P3",
                       "01_immune_20210915/moderna_real_data_processed_with_ri",
                       "skscore.csv"))

# Quality checks
xtabs(~dat$URMforsubcohortsampling, addNA=T)
xtabs(~dat$MinorityInd, addNA=T)

# Rename columns
dat <- dat_raw %>% rename(
  "y_star" = EventTimePrimaryD57,
  "delta_star" = EventIndPrimaryD57,
  "a" = Day57pseudoneutid50,
  "w1" = MinorityInd,
  "w2" = standardized_risk_score,
  "weight" = wt.D57,
  "delta" = TwophasesampIndD57
)

# Filter out placebo patients
dat %<>% filter(Trt==1)

# Select needed columns
dat %<>% subset(select=c(y_star,delta_star,a,w1,w2,weight,delta))

# Filter data to subcohort
dat %<>% filter(delta==1)
dat %<>% filter(!is.na(w2))     # !!!!!
dat %<>% filter(!is.na(weight)) # !!!!!

# Rescale A to lie in [0,1]
a_scale <- ceiling(10*max(dat$a))/10
dat$a <- dat$a / a_scale

# Quality checks
c(min(dat$a),max(dat$a))
hist <- function(x){ ggplot(data.frame(x=x),aes(x=x))+geom_histogram() }
sum(is.na(dat$y_star))
sum(is.na(dat$delta_star))
sum(is.na(dat$a))
sum(is.na(dat$w1))
sum(is.na(dat$w2))
dat_raw$Day29pseudoneutid50 %>% hist()
dat_raw$Day57pseudoneutid50 %>% hist()

# Delta (selected in second-stage sample)
# Weight (sampling weight)



#########################.
##### Data analysis #####
#########################.

# Construct component functions
Phi_n <- construct_Phi_n(d$a, d$weights, type=params$ecdf_type) # !!!!!
Phi_n_inv <- construct_Phi_n(d$a, d$weights, which="inverse",
                             type=params$ecdf_type) # !!!!!
S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                      csf=TRUE)
f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                             type=params$g_n_type, k=15)
f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
g_n <- construct_g_n(vlist$AW_grid, f_aIw_n, f_a_n)
omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)

# Construct Gamma_n
Gamma_n <- construct_Gamma_n(dat_orig, vlist$A_grid, omega_n, S_n, g_n)

# Construct theta_n
Psi_n <- Vectorize(function(x) { Gamma_n(Phi_n_inv(x)) })
gcm <- gcmlcm(x=seq(0,1,0.0001), y=Psi_n(seq(0,1,0.0001)), type="gcm")
dGCM <- Vectorize(function(x) {
  # The round deals with a floating point issue
  index <- which(round(x,5)<=gcm$x.knots)[1]-1
  if (index==0) { index <- 1 }
  return(gcm$slope.knots[index])
})

# Construct Grenander-based theta_n
theta_n_Gr <- function(x) { dGCM(Phi_n(x)) }
theta_n <- theta_n_Gr

grid <- seq(0,1,0.01)
ggplot(
  data.frame(x=grid, y=theta_n(grid)),
  aes(x=x, y=y)
) +
  geom_line()
