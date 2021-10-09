
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

df_raw <- read.csv(paste0("Z:/covpn/p3001/download_data/Moderna COVE mRNA 1273P3",
                       "01_immune_20210915/moderna_real_data_processed_with_ri",
                       "skscore.csv"))

# Rename columns
df <- df_raw %>% rename(
  "y_star" = EventTimePrimaryD57,
  "delta_star" = EventIndPrimaryD57,
  "a" = Day57pseudoneutid50,
  "w1" = MinorityInd,
  "w2" = standardized_risk_score,
  "weight" = wt.D57,
  "delta" = TwophasesampIndD57
)

# Quality checks
xtabs(~df$URMforsubcohortsampling, addNA=T)
xtabs(~df$MinorityInd, addNA=T)

# Filter out placebo patients
df %<>% filter(Trt==1)

# Select needed columns
df %<>% subset(select=c(y_star,delta_star,a,w1,w2,weight,delta))

# Filter data to subcohort
df %<>% filter(delta==1)
df %<>% filter(!is.na(w2))     # !!!!!
df %<>% filter(!is.na(weight)) # !!!!!

# Rescale A to lie in [0,1]
a_scale <- ceiling(10*max(df$a))/10
df$a <- df$a / a_scale

# Quality checks
c(min(df$a),max(df$a))
hist <- function(x){ ggplot(data.frame(x=x),aes(x=x))+geom_histogram() }
sum(is.na(df$y_star))
sum(is.na(df$delta_star))
sum(is.na(df$a))
sum(is.na(df$w1))
sum(is.na(df$w2))
df_raw$Day29pseudoneutid50 %>% hist()
df_raw$Day57pseudoneutid50 %>% hist()

# Delta (selected in second-stage sample)
# Weight (sampling weight)



#########################.
##### Data analysis #####
#########################.

# !!!!! TO DO; copy from est_curve()
