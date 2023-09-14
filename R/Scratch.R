
# Alternative to log(1-CVE) transformation
if (F) {
  
  # log transform
  
  # logRRs <- seq(-5,5,0.1)
  logRRs <- seq(-5,0,0.1)
  RRs <- exp(logRRs)
  log2 <- function(x) { log(x+0.01) }
  exp2 <- function(x) { exp(x) - 0.01 }
  logRRs2 <- log2(RRs)
  df_plot <- data.frame(
    x = rep(RRs,2),
    y = c(logRRs,logRRs2),
    which = rep(c("log(x)", "log(x+0.01)"), each=length(RRs))
  )
  ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_line() +
    labs(color="Transformation") +
    theme(legend.position="bottom")
  
  # logit transform
  
  expit <- function(x) { exp(x) / ( 1 + exp(x) ) }
  logit <- function(x) { log(x/(1-x)) }
  logitRs <- seq(-5,5,0.1)
  # logitRs <- seq(-5,0,0.1)
  Rs <- expit(logitRs)
  logitRs2 <- logit(Rs+0.02*(0.5-Rs))
  df_plot <- data.frame(
    x = rep(Rs,2),
    y = c(logitRs,logitRs2),
    which = rep(c("logit(x)", "logit*(x)"), each=length(Rs))
  )
  ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_line() +
    labs(color="Transformation") +
    theme(legend.position="bottom")
  

}

# Distribution of nondifferentiable function of normal variable
if (F) {
  
  hist <- function(x, bins=100, xlim=NA){
    ggplot(data.frame(x=x),aes(x=x)) + geom_histogram(bins=bins)
  }
  
  n_reps <- 10000
  var <- rep(NA, n_reps)
  x <- rnorm(n_reps)
  g <- function(x) { max(x, -2) }
  for (i in c(1:n_reps)) {
    # var[i] <- x[i]
    var[i] <- g(x[i])
  }
  hist(var)
  
}

# Debugging
if (F) {
  
  # data <- dat.vac.seroneg <- subset(dat.mock, Trt==1 & ph1)
  data <- data.frame(
    EventTimePrimary = dat$v$y,
    EventIndPrimary = dat$v$delta,
    s = dat$v$s,
    risk_score = dat$v$x$x1,
    ph2 = dat$v$z,
    Wstratum = dat$v$strata
  )
  marker.name <- "s"
  tfinal.tpeak <- 181
  # data.ph2=subset(data, ph2==1)
  data.ph2 <- dplyr::filter(data, ph2==1)
  form.s <- Surv(EventTimePrimary, EventIndPrimary) ~ 1 # !!!!!
  form.0 <- update(form.s, as.formula("~. + risk_score")) # !!!!!
  f1=update(form.0, as.formula(paste0("~.+",marker.name)))
  fit.risk.1 <- svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=data))
  # out=marginalized.risk(fit.risk.1, marker.name, data.ph2, t=t, ss=ss, weights=data.ph2$wt, categorical.s=categorical.s)
  

  marginalized.risk.svycoxph.boot=function(marker.name, type, data, t, B, ci.type="quantile", numCores=1) {  

    
    # used in both point est and bootstrap
    # many variables are not passed but defined in the scope of marginalized.risk.svycoxph.boot
    fc.1=function(data.ph2, data, categorical.s, n.dean=FALSE){
      # non-competing risk implementation
      # inline design object b/c it may also throw an error
      if ( !inherits(fit.risk.1, "try-error" )) {
        if (n.dean) c(n.dean= last(coef(fit.risk.1)/sqrt(diag(fit.risk.1$var))) * sqrt(1/fit.risk.1$n + 1/fit.risk.1$nevent), out) else out
      } else {
        rep(NA, ifelse(n.dean,1,0)+length(ss))
      }
    }
    
    if (type==1) {
      # conditional on S=s (quantitative)
      # don't sort ss or do ss=ss[!duplicated(ss)] because e.g. 15% will be lost and later code depends on that
      ss=sort(c(
        # Lars quantiles so that to be consistent with his analyses, also add every 5% to include s1 and s2 for sensitivity analyses
        report.assay.values(data[[marker.name]][data$EventIndPrimary==1], marker.name.to.assay(marker.name)), 
        # 2.5% and 97.5% as the leftmost and rightmost points 
        wtd.quantile(data[[marker.name]], data$wt, c(0.025,0.05,0.95,0.975)),
        # equally spaced values so that the curves look good  
        seq(min(data[[marker.name]], na.rm=TRUE), max(data[[marker.name]], na.rm=TRUE), length=100)[-c(1,100)],
        # useful for reports
        if (log10(100)>min(data[[marker.name]], na.rm=TRUE) & log10(100)<max(data[[marker.name]], na.rm=TRUE)) log10(100)
      ))
      
      prob = fc.1(data.ph2, data, n.dean=TRUE, categorical.s=F)
      if (!comp.risk) {
        n.dean=prob[1]
        prob=prob[-1]
      } 
      
    }
    
    # bootstrap
    if(config$case_cohort) ptids.by.stratum=get.ptids.by.stratum.for.bootstrap (data)     
    seeds=1:B; names(seeds)=seeds
    out=mclapply(seeds, mc.cores = numCores, FUN=function(seed) {   
      seed=seed+560
      if (verbose>=2) myprint(seed)
      
      if(config$case_cohort) {
        dat.b = get.bootstrap.data.cor (data, ptids.by.stratum, seed) 
      } else {
        dat.b = bootstrap.case.control.samples(data, seed, delta.name="EventIndPrimary", strata.name="tps.stratum", ph2.name="ph2") 
      }        
      dat.b.ph2=subset(dat.b, ph2==1)     
      
      if(type==1) {
        # conditional on s
        fc.1(dat.b.ph2, dat.b, categorical.s=F, n.dean=T)
        
      } else if (type==2) {
        # conditional on S>=s
        fc.2(dat.b.ph2)        
        
      } else if (type==3) {
        # conditional on a categorical S
        fc.1(dat.b.ph2, dat.b, n.dean=F, categorical.s=T)
        
      } else if (type==4) {
        # conditional on S=s (quantitative)
        fit.risk.b=try(svycoxph(f1, design=twophase(id=list(~1,~1), strata=list(NULL,~Wstratum), subset=~ph2, data=dat.b)))
        if ( class (fit.risk.b)[1] != "try-error" ) {
        } else {
          NA
        }
        
      } else stop("wrong type")
      
    })
    res=do.call(cbind, out)
    if (type==1 & !comp.risk) {
      # the first row is n.dean
      boot.n.dean=res[1,]
      res=res[-1,]
    }
    res=res[,!is.na(res[1,])] # remove NA's
    if (verbose) str(res)
    
    # restore rng state 
    assign(".Random.seed", save.seed, .GlobalEnv)    
    
    if (ci.type=="quantile") {
      ci.band=t(apply(res, 1, function(x) quantile(x, c(.025,.975), na.rm=T)))
    } else {
      stop("only quantile bootstrap CI supported for now")
    }
    
    ret = list(marker=if(type==3) names(prob) else ss, prob=prob, boot=res, lb=ci.band[,1], ub=ci.band[,2], if(type==1 & !comp.risk) n.dean=c(n.dean, boot.n.dean))   
    if (type==1 & !comp.risk) names(ret)[length(ret)]="n.dean" # this is necessary because when using if, that element won't have a name
    ret  
  }
  
}

# Truncation issue
if (F) {
  
  library(survival)
  library(magrittr)
  library(dplyr)
  
  n <- 10000
  time_birth <- 10*runif(n)
  lmbd <- 1/5
  # lmbd <- 1/20
  H_0_inv <- function(t) { t/lmbd }
  beta_1 <- 0.5
  tx <- as.integer(time_birth>7)
  # z <- rbinom(n, size=1, prob=0.3)
  # beta_2 <- 5
  # lin <- beta_1*tx + beta_2*z
  lin <- beta_1*tx
  U <- runif(n)
  age_death <- H_0_inv(-1*log(U)*exp(-1*lin))
  time_study_start <- 5
  time_death <- time_birth+age_death
  comm_id <- sample(c(1:20), n, replace=T)
  
  # # Option 1: include all children born after study start
  # age_start <- rep(0, length(U))
  # incl <- as.integer(time_birth>time_study_start)
  
  # Option 2: include all children alive at study start
  age_start <- pmax(time_study_start-time_birth,0)
  incl <- as.integer(time_death>time_study_start)
  
  time_cens <- 10
  age_cens <- time_cens-time_birth
  age_end <- pmin(age_death,age_cens)
  delta <- as.integer(age_end==age_death)
  
  # Run model
  df_cox <- filter(data.frame(age_start=age_start, age_end=age_end, delta=delta,
                              tx=tx, incl=incl, comm_id=comm_id), incl==1)
  model <- coxph(Surv(age_start, age_end, delta)~tx, data=df_cox)
  summary(model)
  print(exp(beta_1))
  
  model2 <- coxme(
    Surv(age_start, age_end, delta) ~ tx + (1|comm_id),
    data = df_cox
  )
  summary(model2)
  
  
}

# Plotting spline
if (F) {
  
  # x <- c(0, rep(c(60,365,730,1095), each=2), 1500)
  # y <- rep(c(1.412983,.9819709,.8390149,.9493304,1.012977), each=2)
  # ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) +
  #   geom_line()
  
  # s_coeffs <- log(c(1.002517,.7429874,1.488451,.8494693,1.112713)) # Only cal time trend
  # s_coeffs <- log(c(1.000547,.9162426,1.14562,.9338781))
  # s_coeffs <- log(c(1.00107,.8885394,1.197337,.9192548))
  x <- seq(0,1500, length.out=100)
  knts <- c(0,90,365,730,1460)
  K <- length(knts)
  dk <- function(x,i) {
    if (i==1) {
      return(x)
    } else {
      i <- round(i-1)
      num <- max((x-knts[i])^3,0) - (1/(knts[K]-knts[K-1])) * (
        max((x-knts[K-1])^3,0)*(knts[K]-knts[i]) -
          max((x-knts[K])^3,0)*(knts[K-1]-knts[i])
      )
      den <- (knts[K]-knts[1])^2
      return(num/den)
    }
  }
  n1 <- sapply(x, function(x) {dk(x,1)})
  n2 <- sapply(x, function(x) {dk(x,2)})
  n3 <- sapply(x, function(x) {dk(x,3)})
  n4 <- sapply(x, function(x) {dk(x,4)})
  # spl_vals2 <- s_coeffs[1]*n1 + s_coeffs[2]*n2 + s_coeffs[3]*n3 + s_coeffs[4]*n4
  
  ggplot(data.frame(x=x, y=exp(spl_vals2)), aes(x=x, y=y)) +
    geom_line()
  
  # Calculating lincom multipliers: 0-1 years
  print(mean(sapply(c(1:(365*1)), function(x) {dk(x,1)}))) # intgrl_1
  print(mean(sapply(c(1:(365*1)), function(x) {dk(x,2)}))) # intgrl_2
  print(mean(sapply(c(1:(365*1)), function(x) {dk(x,3)}))) # intgrl_3
  print(mean(sapply(c(1:(365*1)), function(x) {dk(x,4)}))) # intgrl_4
  
  # Calculating lincom multipliers: 0-2 years
  print(mean(sapply(c(1:(365*2)), function(x) {dk(x,1)}))) # intgrl_1
  print(mean(sapply(c(1:(365*2)), function(x) {dk(x,2)}))) # intgrl_2
  print(mean(sapply(c(1:(365*2)), function(x) {dk(x,3)}))) # intgrl_3
  print(mean(sapply(c(1:(365*2)), function(x) {dk(x,4)}))) # intgrl_4
  
  # Calculating lincom multipliers: 0-3 years
  print(mean(sapply(c(1:(365*3)), function(x) {dk(x,1)}))) # intgrl_1
  print(mean(sapply(c(1:(365*3)), function(x) {dk(x,2)}))) # intgrl_2
  print(mean(sapply(c(1:(365*3)), function(x) {dk(x,3)}))) # intgrl_3
  print(mean(sapply(c(1:(365*3)), function(x) {dk(x,4)}))) # intgrl_4
  
  
  knts <- c(0,90,365,730,1460)
  K <- length(knts)
  dk <- function(x,i) {
    if (i==1) { return(x) } else {
      i <- round(i-1)
      num <- max((x-knts[i])^3,0) - (1/(knts[K]-knts[K-1])) * (
        max((x-knts[K-1])^3,0)*(knts[K]-knts[i]) -
          max((x-knts[K])^3,0)*(knts[K-1]-knts[i])
      )
      den <- (knts[K]-knts[1])^2
      return(num/den)
    }
  }
  n1 <- sapply(x, function(x) {dk(x,1)})
  n2 <- sapply(x, function(x) {dk(x,2)})
  n3 <- sapply(x, function(x) {dk(x,3)})
  n4 <- sapply(x, function(x) {dk(x,4)})
  
  model <- coxme(Surv(`_t0`, `_t`, `_d`) ~ intervention + s_int1 + s_int2 + s_int3 + s_int4 + s_cal1 + ...)
  c <- as.numeric(model$coefficients)
  
  
}

# DWD spline basis 2
if (F) {
  
  x <- c(0, 180, 266, 13)
  K <- 4
  knts <- c(0,100,200,300)
  dk <- function (x,k) {
    num <- max((x-knts[k])^3,0) - max((x-knts[K])^3,0)
    den <- knts[K] - knts[k]
    return(num/den)
  }
  n2 <- x
  n3 <- sapply(x, function(x) {dk(x,1)}) - sapply(x, function(x) {dk(x,K-1)})
  n4 <- sapply(x, function(x) {dk(x,2)}) - sapply(x, function(x) {dk(x,K-1)})
  print(n2); print(n3); print(n4);

}

# DWD spline basis
if (F) {
  
  library(ggplot2)
  library(splines)
  
  n <- 20
  x <- runif(n)
  y <- 1 + 2*x + sin(5*x) + rnorm(n, sd=0.3)
  # basis <- ns(x, df=4, intercept=F, Boundary.knots=quantile(x, c(0.05,0.95)))
  basis <- ns(x, knots=c(0.25,0.5,0.75), intercept=F, Boundary.knots=c(0,1))
  
  b1 <- basis[,1]
  b2 <- basis[,2]
  b3 <- basis[,3]
  b4 <- basis[,4]
  # b5 <- basis[,5] # !!!!!
  
  # !!!!! Construct basis manually
  {
    K <- 5
    knts <- c(0,0.25,0.5,0.75,1)
    dk <- function (x,k) {
      num <- max((x-knts[k])^3,0) - max((x-knts[K])^3,0)
      den <- knts[K] - knts[k]
      return(num/den)
    }
    n2 <- x
    n3 <- sapply(x, function(x) {dk(x,1)}) - sapply(x, function(x) {dk(x,K-1)})
    n4 <- sapply(x, function(x) {dk(x,2)}) - sapply(x, function(x) {dk(x,K-1)})
    n5 <- sapply(x, function(x) {dk(x,3)}) - sapply(x, function(x) {dk(x,K-1)})
  }
  
  model1 <- lm(y~b1+b2+b3+b4-1)
  model2 <- lm(y~n2+n3+n4+n5-1)
  y_pred1 <- predict(model1)
  y_pred2 <- predict(model2)
  
  df_plot <- data.frame(
    x = rep(x,2),
    y = c(y_pred1,y_pred2),
    which = rep(c("ns","manual"), each=length(x))
  )
  ggplot(data.frame(x=x, y=y), aes(x=x, y=y)) +
    geom_point(alpha=0.3) +
    geom_line(data=df_plot, color="darkgreen") +
    facet_wrap(~which)
  
}

# Comparing new package to old estimates
if (F) {
  
  var <- "omega_4"
  x1 <- filter(sim$results, use_package==T)
  x2 <- filter(sim$results, use_package==F)
  df_plot <- data.frame(
    x = c(x1[,var], x2[,var]),
    grp = c(rep("vaccine", length(x1[,var])),
            rep("VaxCurve", length(x2[,var])))
  )
  ggplot(df_plot, aes(x=x, group=grp, fill=factor(grp))) +
    geom_histogram(color="white") +
    facet_wrap(~grp, ncol=1)

}

# Illustration for Nima
if (F) {
  
  # f <- Vectorize(function(x) {
  #   k <- 5
  #   dbeta(x, shape1=1, shape2=k)
  # })
  # f <- function(x) { 1 }
  
  library(truncnorm)
  library(ggplot2)
  
  grid <- seq(0.01,1,0.01)
  vec_x <- vec_y <- vec_d <- vec_p <- c()
  for (distr in c("Uniform", "Normal")) {
    for (k in c(1,2,3,4,5)) {
      
      if (distr=="Uniform") {
        f <- function(x) { 1 }
      } else if (distr=="Normal") {
        f <- Vectorize(function(x) {
          dtruncnorm(x, a=0, b=1, mean=0.5, sd=0.2)
        })
      }
      
      if (k==1) { p <- c(5,1) }
      if (k==2) { p <- c(2,1) }
      if (k==3) { p <- c(1,1) }
      if (k==4) { p <- c(1,2) }
      if (k==5) { p <- c(1,5) }
      
      f_mod <- function(x) { f(x) * dbeta(x, shape1=p[1], shape2=p[2]) }
      int_f <- Vectorize(function(x) { integrate(f_mod, lower=0, upper=x)$value })
      f2 <- function(x) { grad(int_f,x)/int_f(1) }
      # f2 <- function(x) { int_f(x)/int_f(1) }
      
      vec_x <- c(vec_x, grid)
      vec_y <- c(vec_y, f2(grid))
      vec_d <- c(vec_d, rep(distr, length(grid)))
      vec_p <- c(vec_p, rep(paste0(k,". Beta(",p[1],",",p[2],")"), length(grid)))
    }
  }
  
  df_plot <- data.frame(x=vec_x, y=vec_y, d=vec_d, p=vec_p)
  ggplot(df_plot, aes(x=x, y=y)) + geom_line() +
    facet_grid(rows=dplyr::vars(d), cols=dplyr::vars(p),
               scales="free_y")
  
}

# Loop through files and extract edge mass and # events
if (F) {
  
  files <- dir("Debugging/TEMP")
  files <- files[-c(1,2)]
  file_stub <- "slurm-4559435_"
  
  numPH2 <- rep(NA,44)
  edgeMass <- rep(NA,44)
  
  for (i in c(1:44)) {
    
    file <- paste0("Debugging/TEMP/",file_stub,i,".out")
    fileString <- readChar(file, file.info(file)$size)
    numPH2[i] <- as.numeric(sub("\".*", "", sub(".*# of PH2 events: ", "", fileString)))
    edgeMass[i] <- as.numeric(sub("\".*", "", sub(".*Edge mass: ", "", fileString)))
    
  }
  
  print(data.frame(i=c(1:44), numPH2=numPH2, edgeMass=edgeMass))
  
}

# DWD
if (F) {
  
  library(microbenchmark)
  
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir)
  dat <- ss(dat_orig, which(dat_orig$z==1))
  vlist <- create_val_list(dat_orig)
  
  fit_CoxPH <- construct_Q_n(dat, vlist$Q_n, type="Cox PH")
  fit_CoxPH2 <- construct_Q_n(dat, vlist$Q_n, type="Cox PH2")
  fit_CoxPH3 <- construct_Q_n(dat, vlist$Q_n, type="Cox PH3")
  
  Q_CoxPH <- fit_CoxPH$srv
  Q_CoxPH2 <- fit_CoxPH2$srv
  Q_CoxPH3 <- fit_CoxPH3$srv
  
  times <- round(seq(0,200,10))
  n <- length(times)
  x_a <- as.data.frame(cbind(x1=rep(0.2,n), x2=rep(1,n)))
  
  microbenchmark({
    Q_CoxPH(t=150, x=c(1,1), s=0.5)
  }, times=100L)
  microbenchmark({
    Q_CoxPH2(t=150, x=c(1,1), s=0.5)
  }, times=100L)
  microbenchmark({
    Q_CoxPH3(t=150, x=c(1,1), s=0.5)
  }, times=100L)
  
  microbenchmark({
    Q_CoxPH(t=times, x=x_a, s=rep(0.2,n))
  }, times=100L)
  microbenchmark({
    Q_CoxPH2(t=times, x=x_a, s=rep(0.2,n))
  }, times=100L)
  microbenchmark({
    sapply(times, function(t) {
      Q_CoxPH3(t=t, x=c(0.2,1), s=0.2)
    })
  }, times=100L)
  microbenchmark({
    unlist(lapply(times, function(t) {
      Q_CoxPH3(t=t, x=c(0.2,1), s=0.2)
    }))
  }, times=100L)
  
}

# Profiling construct_S_n
if (F) {
  
  C <- list(
    points = round(seq(0,1,0.02),2), # round(seq(0,1,0.1),2)
    alpha_1 = 0.5,
    alpha_2 = 0.7,
    t_0 = 200,
    appx = list(t_0=1, x_tol=25, s=0.01) # !!!!!
  )
  L <- list(
    n=300, alpha_3=-1, dir="decr",
    sc_params=list(lmbd=2e-4, v=1.5, lmbd2=5e-5, v2=1.5),
    distr_S="N(0.3+0.4x2,0.09)", edge="none", surv_true="Cox PH",
    sampling="two-phase (50%)", wts_type="true"
  )
  
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_S, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir, L$wts_type)
  
  # Rescale S to lie in [0,1] and round values
  s_min <- min(dat_orig$s,na.rm=T)
  s_max <- max(dat_orig$s,na.rm=T)
  s_shift <- -1 * s_min
  s_scale <- 1/(s_max-s_min)
  dat_orig$s <- (dat_orig$s+s_shift)*s_scale
  dat_orig <- round_dat(dat_orig)
  
  # Setup
  n_orig <- length(dat_orig$z)
  dat <- ss(dat_orig, which(dat_orig$z==1))
  vlist <- create_val_list(dat_orig)
  
  # Construct Q_n
  srvSL <- construct_Q_n(dat, vlist$Q_n, type="Cox PH", print_coeffs=T)
  
}

# Histograms
if (F) {
  
  # [1] "sim_uid"     "level_id"    "rep_id"      "n"           "alpha_3"    
  # [6] "dir"         "sc_params"   "distr_S"     "edge"        "surv_true"  
  # [11] "sampling"    "wts_type"    "test"        "runtime"     "type_1"     
  # [16] "reject_1"    "p_val_1"     "beta_n_1"    "var_n_1"     "type_2"     
  # [21] "reject_2"    "p_val_2"     "beta_n_2"    "var_n_2"     "Theta_0.1"  
  # [26] "Theta_0.4"   "Theta_0.8"   "etastar_0.1" "etastar_0.4" "etastar_0.8"  
  
  r <- sim$results
  hist <- function(x, bins=30, xlim=NA){
    ggplot(data.frame(x=x),aes(x=x)) + geom_histogram(bins=bins)
    # xlim(-0.01,0.01)
  }
  hist(r$Theta_0.1, bins=50)
  hist(r$Theta_0.4, bins=50)
  hist(r$Theta_0.8, bins=50)
  
}

# Testing coefficients
if (F) {
  
  n <- 10000
  x <- runif(n)
  # x <- rbeta(n, shape1=10, shape2=0.1)
  y <- 0.8*x^2 + rnorm(n, sd=0.1)
  # y <- 0.8*x^2 + 0.6*x + rnorm(n, sd=0.1)
  # y <- 0.8*x^2 + 0.6*x + 0.5 + rnorm(n, sd=0.1)
  
  # Intermediates
  s_x1 <- mean(x)
  s_x2 <- mean(x^2)
  s_x3 <- mean(x^3)
  s_x4 <- mean(x^4)
  s_y1 <- mean(y)
  s_xy <- mean(x*y)
  s_x2y <- mean(x^2*y)
  
  coeff_noint <- (s_xy*s_x3-s_x2y*s_x2) / (s_x3^2-s_x2*s_x4)
  num <- (s_x3-s_x1*s_x2)*(s_xy-s_x1*s_y1) + (s_x2y-s_x2*s_y1)*(s_x1^2-s_x2)
  den <- (s_x1*s_x2-s_x3)^2 + (s_x2^2-s_x4)*(s_x2-s_x1^2)
  coeff_withint <- num/den
  print(den) # !!!!!
  
  model_1 <- lm(y~0+x+I(x^2))
  print("Without intercept")
  print(coeff_noint)
  print(model_1$coefficients[[2]])
  
  model_2 <- lm(y~x+I(x^2))
  print("With intercept")
  print(coeff_withint)
  print(model_2$coefficients[[3]])
  
}

# DWD
if (F) {
  
  # One
  r <- sim$results
  df_plot_1 <- data.frame(
    x = c(r$Theta_0.0,r$Gamma_0.0,
          r$Theta_0.2,r$Gamma_0.2,
          r$Theta_0.5,r$Gamma_0.5,
          r$Theta_0.7,r$Gamma_0.7,
          r$Theta_1.0,r$Gamma_1.0),
    which = rep(c(rep("Theta",500),rep("Gamma",500)),5),
    point = c(rep("0.0",1000), rep("0.2",1000), rep("0.5",1000),
              rep("0.7",1000), rep("1.0",1000))
  )
  df_plot_2 <- data.frame(
    x = c(r$Theta_0.0, r$Theta_0.2, r$Theta_0.5, r$Theta_0.7, r$Theta_1.0),
    point = c(rep("0.0",500), rep("0.2",500), rep("0.5",500),
              rep("0.7",500), rep("1.0",500))
  )
  ggplot(df_plot_2, aes(x=x)) +
    facet_wrap(~point, scales="free") +
    geom_histogram(bins=50)
  
  # Two
  r <- sim$results
  nr <- nrow(r)
  df_plot <- data.frame(
    x = c(r$Theta_0.5, r$Theta2_0.5, r$Th_cmp_1,
          r$Th_cmp_3, r$etastar1, r$etastar2),
    which = c(rep("Theta (0.5)",nr),
              rep("Theta2 (0.5)",nr),
              rep("Theta (cmp 1)",nr),
              rep("Theta (cmp 3)",nr),
              rep("etastar1",nr),
              rep("etastar2",nr))
  )
  ggplot(df_plot, aes(x=x)) +
    facet_wrap(~which, scales="free") +
    geom_histogram(bins=50)
    # labs(title=t)
  
  # Three
  r1 <- filter(sim$results, tmp=="old etastar")
  r2 <- filter(sim$results, tmp=="new etastar")
  mean(r1$reject_1)
  mean(r2$reject_1)
  nr1 <- nrow(r1)
  nr2 <- nrow(r2)
  df_plot <- data.frame(
    x = c(r1$Theta_0.2, r2$Theta_0.2, r1$Theta_0.5, r2$Theta_0.5,
          r1$etastar_0.2, r2$etastar_0.2, r1$etastar_0.5, r2$etastar_0.5,
          r1$p_val_1, r2$p_val_1, r1$beta_n_1, r2$beta_n_1,
          r1$var_n_1, r2$var_n_1, sqrt(r1$var_n_1), sqrt(r2$var_n_1)),
    which = c(rep("Theta_0.2 (old eta)",nr1),
              rep("Theta_0.2 (new eta)",nr2),
              rep("Theta_0.5 (old eta)",nr1),
              rep("Theta_0.5 (new eta)",nr2),
              rep("etastar_0.2 (old eta)",nr1),
              rep("etastar_0.2 (new eta)",nr2),
              rep("etastar_0.5 (old eta)",nr1),
              rep("etastar_0.5 (new eta)",nr2),
              rep("P-val (old eta)",nr1),
              rep("P-val (new eta)",nr2),
              rep("beta_n (old eta)",nr1),
              rep("beta_n (new eta)",nr2),
              rep("var_n (old eta)",nr1),
              rep("var_n (new eta)",nr2),
              rep("sd_n (old eta)",nr1),
              rep("sd_n (new eta)",nr2))
  )
  ggplot(df_plot, aes(x=x)) +
    facet_wrap(~which, scales="free", ncol=4) +
    geom_histogram(bins=30)
  
}

# Basic simulation of hypothesis test
if (F) {
  
  sim <- new_sim()
  
  create_data <- function(n, mu_x, mu_y, sigma_x, sigma_y, rho) {
    mu <- c(mu_x,mu_y)
    Sigma <- rbind(
      c(sigma_x^2,rho*sigma_x*sigma_y),
      c(rho*sigma_x*sigma_y,sigma_y^2)
    )
    xy <- mvrnorm(n=n, mu=mu, Sigma=Sigma)
    return(list(x=xy[,1], y=xy[,2]))
  }
  
  sim %<>% set_levels(
    n = c(500), # 100
    mu_x = c(0,1),
    mu_y = c(0,1),
    sigma_x = 8,
    sigma_y = 20,
    rho = c(0.1,0.5,0.9)
  )
  
  sim %<>% set_config(num_sim=10000)
  
  sim %<>% set_script(function() {
    
    dat <- create_data(n=L$n, mu_x=L$mu_x, mu_y=L$mu_y, sigma_x=L$sigma_x,
                       sigma_y=L$sigma_y, rho=L$rho)
    x <- dat$x
    y <- dat$y
    rho_n <- cor(x,y)
    
    stat_x <- (sqrt(L$n)*mean(x))/sd(x)
    stat_y <- (sqrt(L$n)*mean(y))/sd(y)
    stat_xy_unscaled <- ( sqrt(L$n)*(mean(x)+mean(y)) ) /
      sqrt( var(x) + 2*cov(x,y) + var(y) )
    stat_xy_scaled <- sqrt(L$n/(2+2*rho_n)) * (mean(x)/sd(x) + mean(y)/sd(y))
    
    stat_xy_l2norm <- L$n*(mean(x)/sd(x))^2 +
      L$n*((sd(x)*mean(y)-rho_n*sd(y)*mean(x))/(sd(x)*sd(y)*sqrt(1-rho_n^2)))^2
    
    p_x <- pchisq(stat_x^2, df=1, lower.tail=FALSE)
    p_y <- pchisq(stat_y^2, df=1, lower.tail=FALSE)
    p_xy_unscaled <- pchisq(stat_xy_unscaled^2, df=1, lower.tail=FALSE)
    p_xy_scaled <- pchisq(stat_xy_scaled^2, df=1, lower.tail=FALSE)
    p_xy_l2norm <- pchisq(stat_xy_l2norm, df=2, lower.tail=FALSE)
    
    return (list(
      "rej_x" = ifelse(p_x<0.05,1,0),
      "rej_y" = ifelse(p_y<0.05,1,0),
      "rej_xy_unscaled" = ifelse(p_xy_unscaled<0.05,1,0),
      "rej_xy_scaled" = ifelse(p_xy_scaled<0.05,1,0),
      "rej_xy_l2norm" = ifelse(p_xy_l2norm<0.05,1,0),
      "rej_bonf" = ifelse(p_x<0.025||p_y<0.025,1,0),
      "rej_holm" = ifelse(min(p_x,p_y)<0.025||max(p_x,p_y)<0.05,1,0)
    ))
    
  })
  
  sim %<>% run()
  
  SimEngine::summarize(sim) %>% dplyr::rename(
    "rej_x" = mean_rej_x,
    "rej_y" = mean_rej_y,
    "rej_xy_unscaled" = mean_rej_xy_unscaled,
    "rej_xy_scaled" = mean_rej_xy_scaled,
    "rej_xy_l2norm" = mean_rej_xy_l2norm,
    "rej_bonf" = mean_rej_bonf,
    "rej_holm" = mean_rej_holm
  )
  
  # ggplot(data.frame(x=test_stat), aes(x=x)) + geom_histogram(bins=100)
  # mean(test_stat)
  # var(test_stat)
  
}

# Hyp test generalized primitive
if (F) {
  
  In <- as.integer
  grid <- seq(0,1,0.001)
  theta_0 <- Vectorize(function(x) {
    q <- 100
    In(x<=(1/q))*(q*x) + In(x>(1/q))
  })
  k <- 100
  G_0 <- function(x) {
    k_x <- 0.5/k
    m <- (0.5*k)/(k-0.5)
    In(x<=k_x)*(k*x) + In(x>k_x)*(m*x+(1-m))
  }
  G_0_inv <- function(x) {
    m <- (k-0.5)/(0.5*k)
    In(x<=0.5)*(x/k) + In(x>0.5)*(m*x+(1-m))
  }
  g_0 <- function(x) { grad(G_0,x) }
  Theta_0 <- Vectorize(function(x) {
    integrate(theta_0, lower=0, upper=x)$value
  })
  Gamma_0 <- Vectorize(function(x) {
    integrate(function(y) {
      theta_0(y) * g_0(y)
    }, lower=0, upper=x)$value
  })
  Gamma_0s <- Vectorize(function(x) { Gamma_0(G_0_inv(x)) })
  fns <- c("theta_0", "G_0", "G_0_inv", "g_0", "Gamma_0", "Gamma_0s", "Theta_0")
  
  df_plot <- data.frame(
    x = rep(grid, 7),
    y = c(theta_0(grid), G_0(grid), G_0_inv(grid),
          g_0(grid), Gamma_0(grid), Gamma_0s(grid),
          Theta_0(grid)),
    which = rep(factor(fns, levels=fns), each=length(grid))
  )
  ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_line() +
    facet_wrap(~which, ncol=3, scales="free") +
    theme(legend.position="none")
  
}

# Linear spline
if (F) {
  
  # !!!!!
  grid <- seq(0,1,0.001)
  m <- 0.25
  f <- function(x) {
    In()
  }
  
  In <- as.integer
  a <- c(-1,1,2,3)
  f1 <- function(x) {
    a1 <- a[1]; a2 <- a[2]; a3 <- a[3]; a4 <- a[4];
    a1 + In(x<=a3)*a2*x + In(x>a3)*(a3*(a2-a4)+a4*x)
  }
  f2 <- function(x) {
    a1 <- a[1]; a2 <- a[2]; a3 <- a[3]; a4 <- a[4];
    a1 + a2*x + (a4-a2)*pmax(0,x-a3)
  }
  grid <- seq(0,3,0.01)
  ggplot(
    data.frame(
      x = rep(grid,2),
      y = c(f1(grid),f2(grid)),
      which = rep(c("f1","f2"), each=length(grid))),
    aes(x=x, y=y, color=which)
  ) +
    geom_line() +
    facet_wrap(~which)
    

}

# New superfunc pattern (v2)
if (F) {
  
  # Benchmarking
  i <- round(seq(10,200,10),-1)
  w_long <- as.data.frame(
    matrix(rep(c(0.5,1),length(i)), ncol=length(c(0.5,1)), byrow=T)
  )
  a_long <- rep(0.8,length(i))
  microbenchmark({
    S_n(i,w_long,a_long) # Mean: 1178 microsec
  }, times=1000L)
  microbenchmark({
    S_n2(i,w_long,a_long) # Mean: 1178 microsec
  }, times=1000L)
  
  for (j in c(1:1000)) {
    S_n(i,w_long,a_long)
  }
  
  S_n2 <- construct_superfunc2(S_n, vec=c(1,2,1))
  
  construct_superfunc2 <- function(fnc, aux=NA, vec=TRUE, vals=NA, rnd=NA) {
    
    htab <- new.env()
    ..new_fnc <- function() {
      
      ..e <- parent.env(environment())
      ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame()) # Simplify?
      
      for (j in 1:length(..e$arg_names)) {
        if (j==2) { # ..e$vec[j]==2
          ..mc[[..e$arg_names[j]]] <- as.list(
            as.data.frame(t(..mc[[..e$arg_names[j]]]), row.names=NA)
          )
        }
      }
      FUN <- function(...) {
        keylist <- lapply(..e$arg_names, function(arg_name) {
          as.numeric(list(...)[[arg_name]])
        })
        key <- paste(keylist, collapse=";")
        val <- ..e$htab[[key]]
        if (is.null(val)) { # !!!!! We can remove this by pre-calculating
          val <- do.call(..e$fnc, list(...))
          ..e$htab[[key]] <- val
        }
        return(val)
      }
      
      # Return value
      return(do.call(mapply, c(FUN=FUN, ..mc, MoreArgs=NULL, USE.NAMES=F)))
      
    }
    
    # Transform `vec` (if needed) and validate
    if (identical(vec[1],T)) { vec <- rep(1, length(names(formals(fnc)))) }
    if (!is.numeric(vec) || any(!(vec %in% c(0,1,2)))) {
      stop("`vec` must be a vector of 0s, 1s, and 2s")
    }
    
    # Set formals and set up environment
    formals(..new_fnc) <- formals(fnc)
    f_env <- new.env(parent=environment(fnc))
    f_env$arg_names <- names(formals(fnc))
    f_env$vec <- vec
    f_env$rnd <- rnd
    f_env$htab <- htab
    f_env$aux <- aux
    f_env$fnc <- fnc
    environment(..new_fnc) <- f_env
    
    # Run function on vals list
    if (is.list(vals)) { do.call(..new_fnc, vals) }
    
    return(..new_fnc)
    
  }
  
}

# New superfunc pattern
if (F) {
  
  construct_superfunc2 <- function(fnc, aux=NA, vec=TRUE, vals=NA, rnd=NA) {
    
    htab <- new.env()
    ..new_fnc <- function() {
      
      ..e <- parent.env(environment())
      ..mc <- lapply(as.list(match.call())[-1L], eval, parent.frame())
      
      keylist <- lapply(..e$arg_names, function(arg_name) {
        as.numeric(..mc[[arg_name]])
      })
      key <- paste(keylist, collapse=";")
      val <- ..e$htab[[key]]
      if (is.null(val)) {
        val <- do.call(..e$fnc, ..mc)
        ..e$htab[[key]] <- val
      }
      res <- val
      
      # Return value
      return(res)
      
    }
    
    # Set formals and set up environment
    formals(..new_fnc) <- formals(fnc)
    f_env <- new.env(parent=environment(fnc))
    f_env$arg_names <- names(formals(fnc))
    f_env$vec <- vec
    f_env$rnd <- rnd
    f_env$htab <- htab
    f_env$aux <- aux
    f_env$fnc <- fnc
    environment(..new_fnc) <- f_env
    
    # Run function on vals list
    if (is.list(vals)) { do.call(..new_fnc, vals) }
    
    return(..new_fnc)
    
  }
  

  S_n2 <- construct_superfunc2(S_n, vec=F)
  
  S_n3 <- (function() {
    .cache <- new.env()
    function(t,w,a) {
      key <- paste(c(t,w,a), collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(t,w,a) {
          S_n(t,w,a)
        })(t,w,a)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  

  
  
  
  n <- length(dat$a)
  microbenchmark({
    S_n(rep(C$t_e, length(dat$a)),dat$w,dat$a)
  }, times=50L)
  microbenchmark({
    unlist(lapply(c(1:n), function(i) { S_n(C$t_e, dat$w[i,], dat$a[i]) }))
  }, times=50L)
  microbenchmark({
    unlist(lapply(c(1:n), function(i) { S_n2(C$t_e, dat$w[i,], dat$a[i]) }))
  }, times=50L)
  microbenchmark({
    unlist(lapply(c(1:n), function(i) { S_n3(C$t_e, dat$w[i,], dat$a[i]) }))
  }, times=50L)

  w_lst <- list()
  dat_lst <- list()
  for (i in c(1:n)) {
    dat_lst[[i]] <- list(
      t = C$t_e,
      w = as.numeric(dat$w[i,]),
      a = dat$a[i]
    )
    w_lst[[i]] <- as.numeric(dat$w[i,])
  }
  microbenchmark({
    unlist(lapply(dat_lst, function(l) { S_n3(l$t, l$w, l$a) }))
  }, times=50L)
  
  microbenchmark({
    do.call(mapply, c(
      FUN = S_n3,
      list(t=rep(C$t_e,n), w=w_lst, a=dat$a),
      USE.NAMES = F
    ))
  }, times=50L)
  
  
  
  
  
  dat2 <- list(
    w = data.frame(w1=c(1,2,3), w2=c(11,22,33)),
    w_alt = data.frame(w1=c(1,2,3), w2=c(11,22,33)),
    a = c(0.1,0.2)
  )
  n <- length(dat2$a)
  
  fnc_plain <- function(w,a) { sum(w*a) }
  fnc_plain(as.numeric(dat2$w[1,]), dat2$a[1])
  fnc_plain(as.numeric(dat2$w[2,]), dat2$a[2])
  
  # Current
  fnc_plain(dat2$w,dat2$a)
  
  # New
  sapply(c(1:n), function(i) { fnc_plain(dat2$w[i,],dat2$a[i]) })
  
  fnc_vec <- function(fnc) {
    
  }
  
}

# Misc
if (F) {
  
  dat_orig <- generate_data(300, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir) # wts_type="estimated"
  if (T) {
    a_min <- min(dat_orig$a,na.rm=T)
    a_max <- max(dat_orig$a,na.rm=T)
    a_shift <- -1 * a_min
    a_scale <- 1/(a_max-a_min)
    dat_orig$a <- (dat_orig$a+a_shift)*a_scale
    dat_orig <- round_dat(dat_orig)
  }
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  vlist <- create_val_list(dat_orig)
  srvSL <- construct_S_n(dat, vlist$S_n, type="Cox PH", print_coeffs=T)
  S_n <- srvSL$srv
  Sc_n <- srvSL$cens
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n, type="estimated")
  
  # Profile this line
  # omega_n(dat$w,dat$a,dat$y_star,dat$delta_star)
  sapply(c(1:length(dat$a)), function(i) {
    omega_n(as.numeric(dat$w[i,]),dat$a[i],dat$y_star[i],dat$delta_star[i])
  })
  omega_n(dat$w[1,],dat$a[1],dat$y_star[1],dat$delta_star[1]) # -0.4827077
  omega_n(dat$w[2,],dat$a[2],dat$y_star[2],dat$delta_star[2]) # 0.3532617
  
}

# Coming up with a DGM where the COx model fails
if (F) {
  
  # Setup
  C <- list(points=round(seq(0,1,0.02),2), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10, w_tol=25, a=0.01))
  L <- list(
    n=5000, alpha_3=-2, dir="decr",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
    distr_A="Unif(0,1)", edge="none", surv_true="Complex", # "Complex" "Non PH"
    sampling="iid", wts_type="true"
  )
  
  # Generate dataset
  d <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir, L$wts_type)
  
  # Fit a Cox model and extract coefficients
  df <- data.frame(time=d$y_star, ev=d$delta_star, w1=d$w$w1, w2=d$w$w2, a=d$a)
  model <- coxph(Surv(time, ev)~w1+w2+a, data=df)
  beta_n <- as.numeric(model$coefficients)
  
  # Survival estimator (at a point)
  bh <- basehaz(model, centered=F)
  Lambda_n <- bh$hazard[which.min(abs(C$t_e-bh$time))]
  S_n <- (function() {
    .cache <- new.env()
    function(z) {
      key <- paste(z, collapse=" ")
      val <- .cache[[key]]
      if (is.null(val)) {
        val <- (function(z) {
          # exp(-exp(sum(z*beta_n))*Lambda_n(C$t_e))
          exp(-exp(sum(z*beta_n))*Lambda_n)
        })(z)
        .cache[[key]] <- val
      }
      return(val)
    }
  })()
  
  theta_hat <- Vectorize(function(a) {
    1 - mean(sapply(c(1:length(d$a)), function(i) {
      S_n(c(df$w1[i],df$w2[i],a))
    }))
  })
  
  # Plot results
  grid <- round(seq(0,1,0.02),2)
  plot_data <- data.frame(x = rep(grid,2),
                          y = c(attr(d, "theta_true"), theta_hat(grid)),
                          which = rep(c("true", "estimated"), each=51))
  ggplot(plot_data, aes(x=x, y=y, color=which)) +
    geom_line() +
    ylim(c(0,1))
  
}

# Graphs of Phi_n
if (F) {
  
  # Summarize results
  summ_biasP <- list()
  for (i in c(1:51)) {
    m <- format(round(i/50-0.02,2), nsmall=1)
    summ_biasP[[i]] <- list(
      name = paste0("biasP_",m),
      estimate = paste0("estP_",m),
      truth = paste0("Phi_",m)
    )
  }
  summ <- summarize(sim, bias=summ_biasP)
  p_data <- pivot_longer(
    data = summ,
    cols = -c(level_id,n,alpha_3,sc_params,distr_A,edge,
              surv_true,sampling,estimator,dir),
    names_to = c("stat","point"),
    names_sep = "_"
  )
  p_data %<>% mutate(point = as.numeric(point))
  df_vlines <- data.frame(
    x = c(qunif(0.1,0,1), qunif(0.9,0,1)),
    distr_A = rep("Unif(0,1)",2)
  )
  
  # Bias plot
  # Export: 10" x 6"
  # Note: change "bias" to "biasG" for Gamma and "biasP" for Phi
  ggplot(
    filter(p_data, stat=="biasP"),
    aes(x=point, y=value, color=factor(sampling), group=factor(sampling))
  ) +
    geom_vline(aes(xintercept=x), data=df_vlines, color="orange",
               linetype="dashed") +
    geom_line() +
    facet_grid(rows=dplyr::vars(surv_true), cols=dplyr::vars(estimator)) +
    scale_y_continuous(limits=plot_lims$b) + # labels=percent
    # scale_color_manual(values=m_colors) +
    theme(legend.position="bottom") +
    labs(title="Bias", x="A", y=NULL, color="n")
  
}

# Checking Gamma_os_n_star vs Gamma_os_n_star2
if (F) {
  
  # Check 1
  unique_rows <- function(df) { df[!duplicated(df), ] }
  weights_i <- dat$weights
  weights_f <- dat_orig$weights
  n_orig <- round(sum(weights_i))
  (1/n_orig) * sum((1-weights_f))
  
  # Check 2
  w_f <- dat_orig$w
  y_star_f <- dat_orig$y_star
  delta_star_f <- dat_orig$delta_star
  (1/n_orig) * sum((1-weights_f)*q_n(w_f, y_star_f, delta_star_f, 0.7))
  (1/n_orig) * sum(sapply(c(1:n_orig), function(i) {
    (1-weights_f[i])*q_n(as.numeric(w_f[i,]), y_star_f[i], delta_star_f[i], 0.7)
  }))
  
  # Check 3
  (1/n_orig) * sum((1-weights_f)*eta_n(0.7,w_f))
  (1/n_orig) * sum(sapply(c(1:n_orig), function(i) {
    (1-weights_f[i])*eta_n(0.7,as.numeric(w_f[i,]))
  }))
  
  # Check 4
  str_f <- dat_orig$strata
  (1/n_orig) * sum((1-weights_f)*str_f)
  (1/n_orig) * sum((1-weights_f)*(str_f^2))
  
  # Check 5
  # dat_orig2 <- dat_orig # estimated weights
  # dat_orig1 <- dat_orig # true weights
  df_test <- cbind(dat_orig$w,
                   delta = dat_orig$delta,
                   delta_star = dat_orig$delta_star,
                   y_star = dat_orig$y_star,
                   weights = dat_orig$weights,
                   strata = dat_orig$strata)
  head(df_test)
  df_unique <- unique_rows(df_test)
  df_unique %>% arrange(w1,w2,delta_star)
  
}

# Convex least squares
if (F) {
  
  n <- 100
  x <- runif(n)
  # y <- x^2 + rnorm(n, sd=0.1)
  y <- x + rnorm(n, sd=0.1)
  fit <- cvx.lse.reg(t=x, z=y)
  pred_x <- seq(0,1,0.01)
  pred_y <- predict(fit, newdata=pred_x)
  ggplot(data.frame(x=x,y=y), aes(x=x, y=y)) +
    geom_point(alpha=0.5) +
    geom_line(data=data.frame(x=pred_x,y=pred_y), color="red")
  
}

# Comparing plots
if (F) {
  
  C <- list(points=round(seq(0,1,0.02),2), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10, w_tol=25, a=0.01))
  L <- list(
    n=500, alpha_3=-2, dir="decr",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-4, v2=1.5),
    distr_A="N(0.5,0.04)", edge="none", surv_true="exp",
    sampling="iid", estimator=list(est="Grenander",params=params)
  )
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir)
  
  
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  na_head <- sum(round(points,-log10(C$appx$a))<round(a_min,-log10(C$appx$a)))
  points <- round((points+a_shift)*a_scale, -log10(C$appx$a))
  na_tail <- sum(points>1)
  if (na_head>0) {
    points <- points[-c(1:na_head)]
  }
  if (na_tail>0) {
    points <- points[-c((length(points)-na_tail+1):length(points))]
  }
  
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  vlist <- create_val_list(dat_orig)
  Phi_n1 <- construct_Phi_n(dat, type="linear (mid)")
  Phi_n2 <- construct_Phi_n(dat, type="step")
  Phi_n3 <- function(x) { ptruncnorm(x, a=0, b=1, mean=0.5, sd=0.2) }
  
  m <- 10^5
  a <- rtruncnorm(m, a=0, b=1, mean=0.5, sd=0.2)
  Gamma_os_n_star <- Vectorize(function(x) {
    mean( as.integer(a<=x) * (1-exp(-1*L$sc_params$lmbd*C$t_e)) )
  })
  
  grid <- sort(unique(dat$a))
  x_vals1 <- Phi_n1(grid)
  x_vals2 <- Phi_n2(grid)
  x_vals3 <- Phi_n3(grid)
  inds1 <- !base::duplicated(x_vals1)
  inds2 <- !base::duplicated(x_vals2)
  inds3 <- !base::duplicated(x_vals3)
  x_vals1 <- x_vals1[inds1]
  x_vals2 <- x_vals2[inds2]
  x_vals3 <- x_vals3[inds3]
  y_vals1 <- -1 * Gamma_os_n_star(grid[inds1])
  y_vals2 <- -1 * Gamma_os_n_star(grid[inds2])
  y_vals3 <- -1 * Gamma_os_n_star(grid[inds3])
  gcm1 <- gcmlcm(x=x_vals1, y=y_vals1, type="gcm")
  gcm2 <- gcmlcm(x=x_vals2, y=y_vals2, type="gcm")
  gcm3 <- gcmlcm(x=x_vals3, y=y_vals3, type="gcm")
  
  GCM_f1 <- approxfun(x=gcm1$x.knots, y=gcm1$y.knots, method="linear", rule=1)
  GCM_f2 <- approxfun(x=gcm2$x.knots, y=gcm2$y.knots, method="linear", rule=1)
  GCM_f3 <- approxfun(x=gcm3$x.knots, y=gcm3$y.knots, method="linear", rule=1)
  
  dGCM1 <- approxfun(x=gcm1$x.knots[-length(gcm1$x.knots)],
                     y=gcm1$slope.knots, method="constant",rule=2,f=0)
  dGCM2 <- approxfun(x=gcm2$x.knots[-length(gcm2$x.knots)],
                     y=gcm2$slope.knots, method="constant",rule=2,f=0)
  dGCM3 <- approxfun(x=gcm3$x.knots[-length(gcm3$x.knots)],
                     y=gcm3$slope.knots, method="constant",rule=2,f=0)
  theta_n_Gr1 <- Vectorize(function(x) { -1 * dGCM1(Phi_n1(x)) })
  theta_n_Gr2 <- Vectorize(function(x) { -1 * dGCM2(Phi_n2(x)) })
  theta_n_Gr3 <- Vectorize(function(x) { -1 * dGCM3(Phi_n3(x)) })
  df_A <- data.frame(
    x = c(x_vals1,x_vals2,x_vals3),
    y = c(y_vals1,y_vals2,y_vals3),
    y2 = c(GCM_f1(x_vals1), GCM_f2(x_vals2), GCM_f3(x_vals3)),
    which = factor(c(rep("Phi_n (lin mid)",length(x_vals1)),
                     rep("Phi_n (step)",length(x_vals2)),
                     rep("Phi_0",length(x_vals3))))
  )
  grid2 <- round(seq(0,1,0.02),2)
  df_B <- data.frame(
    x = rep(grid2, 3),
    y = c(theta_n_Gr1(grid2), theta_n_Gr2(grid2), theta_n_Gr3(grid2)),
    which = factor(c(rep("Phi_n (lin mid)",51), rep("Phi_n (step)",51),
                     rep("Phi_0",51)))
  )
  ggplot(df_A, aes(x=x,y=y,color=which)) + geom_point(alpha=0.4) + 
    geom_line(aes(y=y2)) +
    facet_wrap(~which,ncol=3) +
    labs(title="x=Phi(A), y=Gamma_0(A)")
  ggplot(df_B, aes(x=x,y=y,color=which)) + geom_line() + 
    facet_wrap(~which,ncol=3) +
    geom_hline(yintercept=0.1813, color="grey") +
    labs(title="theta_n estimates")
  
}

# Testing spline models
if (F) {
  
  # Existing model
  {
    # Fit Cox model
    fml <- "Surv(y_star,delta_star)~a"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star, "a"=dat$a,
                dat$w, "weights"=dat$weights)
    model <- coxph(fml, data=df, weights=dat$weights)
    beta_n <- coefficients(model)
    
    # Get Breslow estimator
    bh <- basehaz(model, centered=FALSE)
    index <- max(which((bh$time<C$t_e)==T))
    est_bshz <- bh$hazard[index]
    
    # Construct conditional survival function
    S_n <- function(w,a) {
      exp(-1*est_bshz*exp(sum(as.numeric(beta_n)*c(a,w))))
    }
    
    # Construct marginalized survival function
    r_M <- Vectorize(function(a) {
      1 - mean(sapply(c(1:length(dat_orig$a)), function(i) {
        S_n(as.numeric(dat_orig$w[i,]),a)
      }))
    })
  }
  
  # Square model
  {
    # Fit Cox model
    fml <- "Surv(y_star,delta_star)~a+a2"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star, "a"=dat$a,
                "a2"=(dat$a)^2,
                dat$w, "weights"=dat$weights)
    model2 <- coxph(fml, data=df, weights=dat$weights)
    beta_n2 <- coefficients(model2)
    
    # Get Breslow estimator
    bh2 <- basehaz(model2, centered=FALSE)
    index2 <- max(which((bh2$time<C$t_e)==T))
    est_bshz2 <- bh2$hazard[index2]
    
    # Construct conditional survival function
    S_n2 <- function(w,a) {
      exp(-1*est_bshz2*exp(sum(as.numeric(beta_n2)*c(a,a^2,w))))
    }
    
    # Construct marginalized survival function
    r_M2 <- Vectorize(function(a) {
      1 - mean(sapply(c(1:length(dat_orig$a)), function(i) {
        S_n2(as.numeric(dat_orig$w[i,]),a)
      }))
    })
  }
  
  # Cube model
  {
    # Fit Cox model
    fml <- "Surv(y_star,delta_star)~a+a2+a3"
    for (i in 1:length(dat$w)) {
      fml <- paste0(fml, "+w",i)
    }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star, "a"=dat$a,
                "a2"=(dat$a)^2, "a3"=(dat$a)^3,
                dat$w, "weights"=dat$weights)
    model3 <- coxph(fml, data=df, weights=dat$weights)
    beta_n3 <- coefficients(model3)
    
    # Get Breslow estimator
    bh3 <- basehaz(model3, centered=FALSE)
    index3 <- max(which((bh3$time<C$t_e)==T))
    est_bshz3 <- bh3$hazard[index3]
    
    # Construct conditional survival function
    S_n3 <- function(w,a) {
      exp(-1*est_bshz3*exp(sum(as.numeric(beta_n3)*c(a,a^2,a^3,w))))
    }
    
    # Construct marginalized survival function
    r_M3 <- Vectorize(function(a) {
      1 - mean(sapply(c(1:length(dat_orig$a)), function(i) {
        S_n3(as.numeric(dat_orig$w[i,]),a)
      }))
    })
  }
  
  # NCS model
  {
    
    # Generate spline basis (up to 6 degrees of freedom)
    degf <- 4
    qnt <- as.numeric(quantile(dat$a, seq(0,1,1/degf)))
    qnt <- unique(qnt)
    spl <- list(
      K = qnt[-c(1,length(qnt))],
      B = c(qnt[1],qnt[length(qnt)]),
      L = length(qnt)-1
    )
    ns_basis <- ns(dat$a, knots=spl$K, Boundary.knots=spl$B)
    
    # Fit Cox model
    fml <- "Surv(y_star,delta_star)~b1"
    for (i in 2:spl$L) { fml <- paste0(fml, "+b",i) }
    for (i in 1:length(dat$w)) { fml <- paste0(fml, "+w",i) }
    fml <- formula(fml)
    df <- cbind("y_star"=dat$y_star, "delta_star"=dat$delta_star,
                dat$w, "weights"=dat$weights)
    for (i in 1:spl$L) { df[paste0("b",i)] <- ns_basis[,i] }
    model4 <- coxph(fml, data=df, weights=dat$weights)
    beta_n4 <- coefficients(model4)
    
    # Get Breslow estimator
    bh4 <- basehaz(model4, centered=FALSE)
    index4 <- max(which((bh4$time<C$t_e)==T))
    est_bshz4 <- bh4$hazard[index4]
    
    # Construct conditional survival function
    S_n4 <- function(w,a) {
      lin <- c(as.numeric(ns(a, knots=spl$K, Boundary.knots=spl$B)),w)
      exp(-1*est_bshz4*exp(sum(as.numeric(beta_n4)*lin)))
    }
    
    # Construct marginalized survival function
    r_M4 <- Vectorize(function(a) {
      1 - mean(sapply(c(1:length(dat_orig$a)), function(i) {
        S_n4(as.numeric(dat_orig$w[i,]),a)
      }))
    })
  }
  
  # Plot results
  grid <- seq(log10(100), log10(3000), length.out=10)
  df_plot <- data.frame(
    x = rep(grid,4),
    y = c(r_M(grid), r_M2(grid), r_M3(grid), r_M4(grid)),
    which = c(rep("linear",10), rep("square",10),
              rep("cube",10), rep("NCS",10))
  )
  ggplot(df_plot, aes(x=x, y=y, color=which)) +
    geom_line() +
    xlim(c(log10(100),log10(30000))) +
    ylim(c(0,0.08))
  
}

# Kaplan-Meier estimates
if (F) {
  
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  
  dat1 <- ss(dat, which(dat$a==2))
  dat2 <- ss(dat, which(dat$a>2 & dat$a<=log10(140)))
  dat3 <- ss(dat, which(dat$a>log10(140)))
  # cve <- function(r_v) { 1-r_v/0.06 }
  
  # srv_ov <- survfit(Surv(dat_orig$y_star,dat_orig$delta_star)~1)
  # risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
  # # print(paste0("Risk: ",round(risk_ov,3),", CVE: ",round(cve(risk_ov),3)))
  # print(paste0("Risk: ",round(risk_ov,3),", CVE: ",round(cve(risk_ov),3)))
  
  srv_ov <- survfit(Surv(dat1$y_star,dat1$delta_star)~1, weights=dat1$weights)
  risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
  # print(paste0("Risk: ",round(risk_ov,3),", CVE: ",round(cve(risk_ov),3)))
  print(paste0("Risk: ",round(risk_ov,3)))
  
  srv_ov <- survfit(Surv(dat2$y_star,dat2$delta_star)~1, weights=dat2$weights)
  risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
  # print(paste0("Risk: ",round(risk_ov,3),", CVE: ",round(cve(risk_ov),3)))
  print(paste0("Risk: ",round(risk_ov,3)))
  
  srv_ov <- survfit(Surv(dat3$y_star,dat3$delta_star)~1, weights=dat3$weights)
  risk_ov <- 1 - srv_ov$surv[which.min(abs(srv_ov$time-C$t_e))]
  # print(paste0("Risk: ",round(risk_ov,3),", CVE: ",round(cve(risk_ov),3)))
  print(paste0("Risk: ",round(risk_ov,3)))
  
}

# Testing whether calculation of Gamma_true is correct
if (F) {
  
  C <- list(points=round(seq(0,1,0.02),2), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10, w_tol=25, a=0.01))
  L <- list(
    n=500, alpha_3=-2, dir="decr",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=4e-5, v2=1.5),
    # distr_A="Unif(0,1)", edge="none", surv_true="exp",
    distr_A="N(0.5,0.01)", edge="none", surv_true="exp",
    sampling="iid", estimator=list(est="Grenander",params=params)
  )
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir)
  attr(dat_orig, "Gamma_true")
  sapply(seq(0,1,0.02), function(x) {
    # F_0 <- function(x) { x }
    F_0 <- function(x) { pnorm(x, mean=0.5, sd=0.1) }
    (1-exp(-1*1e-3*C$t_e)) * F_0(x)
  })
  
}

# Debugging standard errors
if (F) {
  
  # Adapted from coxsurv.fit()
  
  # I_tilde_inv
  mf <- stats::model.frame(model)
  offset <- rep(0,n)
  Y <- model[["y"]]
  n <- nrow(Y)
  risk <- rep(exp(offset-mean(offset)), length=n)
  Terms <- terms(model)
  Terms2 <- Terms
  Terms2 <- delete.response(Terms)
  Terms2 <- Terms2[-ss$terms]
  
  Call <- call("survfit", model, newdata=newdata, stype=2, ctype=1)
  tcall <- Call[c(1, match(c("id", "na.action"), names(Call),
                           nomatch = 0))]
  tcall$data <- newdata
  tcall$formula <- Terms2
  tcall$xlev <- NULL # model$xlevels[match(attr(Terms2, "term.labels"), names(model$xlevels), nomatch = 0)]
  tcall[[1L]] <- quote(stats::model.frame)
  mf2 <- eval(tcall)
  y2 <- model.extract(mf2, "response")
  y2 <- y2[, 1:2, drop = F]
  offset2 <- 0
  x2 <- model.matrix(Terms2, mf2)[, -1, drop = FALSE]
  
  beta <- model$coefficients
  xcenter <- 0
  risk2 <- exp(c(x2 %*% beta) + offset2 - xcenter)
  strata2 <- factor(rep(0, nrow(mf2)))
  id2 <- rep(1, nrow(mf2))
  unlist <- TRUE
  ctype <- 1
  stype <- 2
  se.fit <- T
  varmat <- model$var
  cluster <- NULL
  y <- Y
  x <- model.matrix(model, data=mf)
  wt <- model.weights(mf)
  position <- NULL
  strata <- NULL
  oldid <- NULL
  
  ank(ctype=ctype, stype=stype, se.fit=se.fit, varmat=varmat, cluster=cluster,
      y=y, x=x, wt=wt, risk=risk, position=position, strata=strata, oldid=oldid,
      y2=y2, x2=x2, risk2=risk2, strata2=strata2, id2=id2, unlist=unlist)
  
  
  ank <- function (ctype, stype, se.fit, varmat, cluster, y, x, wt, risk, 
            position, strata, oldid, y2, x2, risk2, strata2, id2, unlist = TRUE) 
  {
    if (missing(strata) || length(strata) == 0) {
      strata <- rep(0L, nrow(y))
    }
    if (is.factor(strata)) {
      ustrata <- levels(strata)
    } else {
      ustrata <- sort(unique(strata))
    }
    nstrata <- length(ustrata)
    survlist <- vector("list", nstrata)
    names(survlist) <- ustrata
    survtype <- if (stype == 1) {
      1
    } else {
      ctype + 1
    }
    vartype <- survtype
    if (is.null(wt))  wt <- rep(1, nrow(y))
    if (is.null(strata)) strata <- rep(1L, nrow(y))
    for (i in 1:nstrata) {
      indx <- which(strata == ustrata[i])
      survlist[[i]] <- agsurv(y[indx, , drop = F], x[indx, 
                                                     , drop = F], wt[indx], risk[indx], survtype, vartype)
    }
    expand <- function(fit, x2, varmat, se.fit) {
      # if (survtype == 1)  {
      #   surv <- cumprod(fit$surv)
      # } else {
      #   surv <- exp(-fit$cumhaz)
      # }
      # if (is.matrix(x2) && nrow(x2) > 1) {
      #   fit$surv <- outer(surv, risk2, "^")
      #   dimnames(fit$surv) <- list(NULL, row.names(x2))
      #   if (se.fit) {
      #     varh <- matrix(0, nrow = length(fit$varhaz), 
      #                    ncol = nrow(x2))
      #     for (i in 1:nrow(x2)) {
      #       dt <- outer(fit$cumhaz, x2[i, ], "*") - fit$xbar
      #       varh[, i] <- (cumsum(fit$varhaz) + rowSums((dt %*% 
      #                                                     varmat) * dt)) * risk2[i]^2
      #     }
      #     fit$std.err <- sqrt(varh)
      #   }
      #   fit$cumhaz <- outer(fit$cumhaz, risk2, "*")
      # }
      # else {
      #   fit$surv <- surv^risk2
      #   if (se.fit) {
      #     dt <- outer(fit$cumhaz, c(x2)) - fit$xbar
      #     varh <- (cumsum(fit$varhaz) + rowSums((dt %*% 
      #                                              varmat) * dt)) * risk2^2
      #     fit$std.err <- sqrt(varh)
      #   }
      #   fit$cumhaz <- fit$cumhaz * risk2
      # }
      # fit
    }
    if (missing(id2) || is.null(id2)) {
      result <- lapply(survlist, expand, x2, varmat, se.fit)
    } else {
      onecurve <- function(slist, x2, y2, strata2, risk2, se.fit) {
        ntarget <- nrow(x2)
        surv <- vector("list", ntarget)
        n.event <- n.risk <- n.censor <- varh1 <- varh2 <- time <- surv
        hazard <- vector("list", ntarget)
        stemp <- as.integer(strata2)
        timeforward <- 0
        for (i in 1:ntarget) {
          slist <- survlist[[stemp[i]]]
          # indx <- which(slist$time > y2[i, 1] & slist$time <=
          #                 y2[i, 2])
          indx <- 18 # !!!!!
          if (length(indx) == 0) {
            timeforward <- timeforward + y2[i, 2] - y2[i,
                                                       1]
          }
          else {
            time[[i]] <- diff(c(y2[i, 1], slist$time[indx]))
            time[[i]][1] <- time[[i]][1] + timeforward
            timeforward <- y2[i, 2] - max(slist$time[indx])
            hazard[[i]] <- slist$hazard[indx] * risk2[i]
            if (survtype == 1)
              surv[[i]] <- slist$surv[indx]^risk2[i]
            n.event[[i]] <- slist$n.event[indx]
            n.risk[[i]] <- slist$n.risk[indx]
            n.censor[[i]] <- slist$n.censor[indx]
            dt <- outer(slist$cumhaz[indx], x2[i, ]) -
              slist$xbar[indx, , drop = F]
            varh1[[i]] <- slist$varhaz[indx] * risk2[i]^2
            varh2[[i]] <- rowSums((dt %*% varmat) * dt) *
              risk2[i]^2
          }
        }
        cumhaz <- cumsum(unlist(hazard))
        if (survtype == 1)
          surv <- cumprod(unlist(surv))
        else surv <- exp(-cumhaz)
        if (se.fit)
          list(n = as.vector(table(strata)[stemp[1]]),
               time = cumsum(unlist(time)), n.risk = unlist(n.risk),
               n.event = unlist(n.event), n.censor = unlist(n.censor),
               surv = surv, cumhaz = cumhaz, std.err = sqrt(cumsum(unlist(varh1)) +
                                                              unlist(varh2)))
        else list(n = as.vector(table(strata)[stemp[1]]),
                  time = cumsum(unlist(time)), n.risk = unlist(n.risk),
                  n.event = unlist(n.event), n.censor = unlist(n.censor),
                  surv = surv, cumhaz = cumhaz)
      }
      if (all(id2 == id2[1])) {
        result <- list(onecurve(survlist, x2, y2, strata2, 
                                risk2, se.fit))
      }
      else {
        uid <- unique(id2)
        result <- vector("list", length = length(uid))
        for (i in 1:length(uid)) {
          indx <- which(id2 == uid[i])
          result[[i]] <- onecurve(survlist, x2[indx, , 
                                               drop = FALSE], y2[indx, , drop = FALSE], strata2[indx], 
                                  risk2[indx], se.fit)
        }
        names(result) <- uid
      }
    }
    if (unlist) {
      if (length(result) == 1) {
        if (se.fit) 
          result[[1]][c("n", "time", "n.risk", "n.event", 
                        "n.censor", "surv", "cumhaz", "std.err")]
        else result[[1]][c("n", "time", "n.risk", "n.event", 
                           "n.censor", "surv", "cumhaz")]
      }
      else {
        temp <- list(n = unlist(lapply(result, function(x) x$n), 
                                use.names = FALSE), time = unlist(lapply(result, 
                                                                         function(x) x$time), use.names = FALSE), n.risk = unlist(lapply(result, 
                                                                                                                                         function(x) x$n.risk), use.names = FALSE), n.event = unlist(lapply(result, 
                                                                                                                                                                                                            function(x) x$n.event), use.names = FALSE), 
                     n.censor = unlist(lapply(result, function(x) x$n.censor), 
                                       use.names = FALSE), strata = sapply(result, 
                                                                           function(x) length(x$time)))
        names(temp$strata) <- names(result)
        if ((missing(id2) || is.null(id2)) && nrow(x2) > 
            1) {
          temp$surv <- t(matrix(unlist(lapply(result, 
                                              function(x) t(x$surv)), use.names = FALSE), 
                                nrow = nrow(x2)))
          dimnames(temp$surv) <- list(NULL, row.names(x2))
          temp$cumhaz <- t(matrix(unlist(lapply(result, 
                                                function(x) t(x$cumhaz)), use.names = FALSE), 
                                  nrow = nrow(x2)))
          if (se.fit) 
            temp$std.err <- t(matrix(unlist(lapply(result, 
                                                   function(x) t(x$std.err)), use.names = FALSE), 
                                     nrow = nrow(x2)))
        }
        else {
          temp$surv <- unlist(lapply(result, function(x) x$surv), 
                              use.names = FALSE)
          temp$cumhaz <- unlist(lapply(result, function(x) x$cumhaz), 
                                use.names = FALSE)
          if (se.fit) 
            temp$std.err <- unlist(lapply(result, function(x) x$std.err), 
                                   use.names = FALSE)
        }
        temp
      }
    }
    else {
      names(result) <- ustrata
      result
    }
  }
  
}

# Combine hypothesis test results
if (F) {
  
  res <- data.frame(
    "analysis" = character(),
    "plot" = character(),
    "reject" = integer(),
    "p_val" = double(),
    "beta_n" = double()
  )
  for (i in c(1:10)) {
    d <- read.csv(paste0("Moderna plots/hyptest_",i,".csv"))
    res[i,] <- list(analysis="Moderna", plot=i, reject=d$reject,
                    p_val=round(d$p_val,5), beta_n=round(d$beta_n,5))
  }
  for (i in c(1:4)) {
    d <- read.csv(paste0("Janssen plots/hyptest_",i,".csv"))
    res[round(i+10),] <- list(analysis="Janssen", plot=i, reject=d$reject,
                    p_val=round(d$p_val,5), beta_n=round(d$beta_n,5))
  }
  res
  
  grid <- round(seq(0,1,0.01),2)
  pdata <- data.frame(
    x = rep(grid,2),
    y = c(Theta_os_n_1(grid), Theta_os_n_2(grid)),
    which = rep(c("1","2"), each=101)
  )
  ggplot(pdata, aes(x=x, y=y, color=which)) + geom_line()
  
  lm1 <- lm(y~x+I(x^2), data=filter(pdata, which==1))
  lm2 <- lm(y~x+I(x^2), data=filter(pdata, which==2))
  summary(lm1)
  summary(lm2)
  pred1 <- as.numeric(predict(lm1, newdata=data.frame(x=grid)))
  pred2 <- as.numeric(predict(lm2, newdata=data.frame(x=grid)))
  ggplot(data.frame(x=grid, y=pred1), aes(x=x, y=y)) + geom_line()
  ggplot(data.frame(x=grid, y=pred2), aes(x=x, y=y)) + geom_line()
  
}

# Comparing isotonic regression PAVA CIs to mean CI
if (F) {
  
  boot_reps <- 100
  v_x <- c()
  v_iso <- c()
  v_grp <- c()
  v_clr <- c()
  for (i in c(1:boot_reps)) {
    n <- 100
    x <- runif(n)
    y <- x + rnorm(n, sd=0.2)
    # x <- 0.3 + 0.4*rbinom(n, size=1, prob=0.5) + 0.01*runif(n)
    # y <- x + rnorm(n, sd=0.01)
    df <- data.frame(x=x, y=y)
    # df_boot <- df[sample(c(1:n), replace=T),]
    df %<>% arrange(x)
    iso <- pava(y=df$y)
    v_x <- c(v_x, df$x, c(0,1))
    v_iso <- c(v_iso, iso, rep(mean(df$y), 2))
    v_grp <- c(v_grp, rep(i,n), rep(i+0.5,2))
    v_clr <- c(v_clr, rep("Iso",n), rep("Mean",2))
  }
  
  res <- data.frame(
    x = v_x,
    y = v_iso,
    grp = v_grp,
    clr = v_clr
  )
  
  ggplot(data=res, aes(x=x, y=y, group=grp, color=clr)) +
    geom_abline(slope=1, intercept=0, alpha=0.3) +
    geom_abline(slope=0, intercept=0.5, alpha=0.3) +
    geom_point(data=df, aes(group=NA, color=NA), alpha=0.5) +
    geom_line(alpha=0.3) +
    theme(legend.position="none") +
    scale_color_manual(values=viridis::viridis(4)[2:3])
  
}

# Qbins cutoffs
if (F) {
  
  # Generate cutoffs
  n_bins <- 4
  n <- 1000
  x <- runif(n)
  pi <- rbinom(n, size=1, prob=0.4)
  # pi <- rep(0,n)
  a <- x * (1-pi)
  wt <- rep(1,n)
  dat <- list(a=a, weights=wt)
  Phi_n_inv <- construct_Phi_n(dat, which="inverse", type="step")
  cutoffs <- round(Phi_n_inv(seq(0,1,length.out=n_bins+1)),3)
  
  num_0 <- sum(cutoffs==0)
  if (num_0==1) {
    cutoffs2 <- cutoffs
  } else {
    a_sub <- a[a!=0]
    wt_sub <- wt[a!=0]
    dat_sub <- list(a=a_sub, weights=wt_sub)
    Phi_n_inv_sub <- construct_Phi_n(dat_sub, which="inverse", type="step")
    cutoffs2 <- c(0,round(Phi_n_inv_sub(seq(0,1,length.out=n_bins)),3))
  }
  cutoffs2
  
  # View percent in each bin
  a_binned <- cut(a, breaks=cutoffs2, right=F, include.lowest=T)
  xtabs(~a_binned)/length(a_binned)
  
}

# Cox models
if (F) {
  
  # Subset datasets
  dat_J <- ss(dat_orig_J, which(dat_orig_J$delta==1))
  dat_M <- ss(dat_orig_M, which(dat_orig_M$delta==1))
  
  # Right-censor data
  t_e_J <- 66
  indices_J <- which(dat_J$y_star>t_e_J)
  dat_J$y_star[indices_J] <- t_e_J
  dat_J$delta_star[indices_J] <- 0
  t_e_M <- 126
  indices_M <- which(dat_M$y_star>t_e_M)
  dat_M$y_star[indices_M] <- t_e_M
  dat_M$delta_star[indices_M] <- 0
  
  # Create analysis dataframes
  df_J <- data.frame(time=dat_J$y_star, ev=dat_J$delta_star, a=dat_J$a,
                     w1=dat_J$w$w1, w2=dat_J$w$w2, w3=dat_J$w$w3,
                     wts1=dat_J$weights)
  df_M <- data.frame(time=dat_M$y_star, ev=dat_M$delta_star, a=dat_M$a,
                     w1=dat_M$w$w1, w2=dat_M$w$w2, w3=dat_M$w$w3,
                     wts1=dat_M$weights)
  
  # Run Cox models
  coxph(Surv(time,ev)~a+w1+w2+w3, data=df_J, weights=wts1) %>% summary()
  coxph(Surv(time,ev)~a+w1+w2+w3, data=df_M, weights=wts1) %>% summary()
  
}

# Approximating Gamma_0
if (F) {
  
  Theta_true <- attr(dat_orig,"Theta_true")
  Gamma_0 <- Vectorize(function(x) {
    Theta_true[which.min(abs(x-seq(0,1,0.02)))[1]]
  })
  
  Gamma_appx <- function(x) { sqrt(x)/3 }
  
  grid <- seq(0,1,0.02)
  ggplot(data.frame(
    x = rep(grid,2),
    y = c(Gamma_0(grid), Gamma_appx(grid)),
    which = rep(c("Gamma_0","Gamma_appx"), each=51)
  ), aes(x=x, y=y, color=which)) + geom_line()  
  
}

# Super Learner for estimating regressions
if (F) {
  
  # Generate data
  n <- 100
  w1 <- sample(seq(0,1,0.1), size=n, replace=T)
  w2 <- rbinom(n, size=1, prob=0.5)
  y <- 1 + 2*w1 + 3*w2 + rnorm(n)
  df <- data.frame(w1=w1,w2=w2,y=y)
  newX <- data.frame(w1=c(0.1,0.3,0.8), w2=c(1,0,1))
  
  # Construct true regression function
  reg_true <- function(w1,w2) { 1 + 2*w1 + 3*w2 }
  
  # True regression values
  apply(X=newX, MARGIN=1, FUN = function(r) { reg_true(r[["w1"]], r[["w2"]]) })
  
  # Fit superleaner model
  model_sl <- SuperLearner(
    Y = df$y,
    X = subset(df, select=-y),
    newX = newX,
    family = "gaussian",
    SL.library = c("SL.mean", "SL.glm", "SL.gam", "SL.randomForest"),
    verbose = FALSE
  )
  as.numeric(model_sl$SL.predict)
  
  # Fit a linear regression
  model <- lm(y~w1+w2, data=df)
  reg_lin <- Vectorize(function(w1,w2) {
    cf <- as.numeric(coef(model))
    return( cf[1] + cf[2]*w1 + cf[3]*w2 )
  })
  apply(X=newX, MARGIN=1, FUN = function(r) { reg_lin(r[["w1"]], r[["w2"]]) })
  
}

# Debugging influence functions
if (F) {
  
  # Construct dataframe
  df <- subset(sim$results, select=c(
    n, sampling, if1_mean, if2_mean, r_1n, r_2n, Psi_1_est, Psi_2_est
  ))
  
  # Scatterplots (Psi_1)
  ggplot(df, aes(x=Psi_1_est, y=if1_mean)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  ggplot(df, aes(x=Psi_1_est, y=r_1n)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  ggplot(df, aes(x=if1_mean, y=r_1n)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  
  # Scatterplots (Psi_2)
  ggplot(df, aes(x=Psi_2_est, y=if2_mean)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  ggplot(df, aes(x=Psi_2_est, y=r_2n)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  ggplot(df, aes(x=if2_mean, y=r_2n)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  
  # Scatterplots (both)
  ggplot(df, aes(x=Psi_1_est, y=Psi_2_est)) + geom_point(alpha=0.3) +
    facet_grid(rows=dplyr::vars(sampling), cols=dplyr::vars(n)) +
    xlim(c(-2e-3,2e-3)) + ylim(c(-2e-3,2e-3))
  
  # Remainder plots
  ggplot(df, aes(x=n, y=sqrt(n)*r_1n)) + geom_point(alpha=0.3)
  ggplot(df, aes(x=n, y=sqrt(n)*r_2n)) + geom_point(alpha=0.3)
  
  # Toy example: sample variance
  var_ests <- c()
  if_means <- c()
  r_ns <- c()
  for (n in c(100,200,400,800)) {
    for (i in c(1:100)) {
      x <- rnorm(n)
      m <- mean(x)
      var_est <- mean((x-m)^2) - 1
      if_mean <- mean(
        # (x-0)^2 - 1
        (x-m)^2 - (var_est+1)
      )
      r_n <- var_est - if_mean
      var_ests <- c(var_ests, var_est)
      if_means <- c(if_means, if_mean)
      r_ns <- c(r_ns, r_n)
    }
  }
  ns <- rep(c(100,200,400,800), each=100)
  df2 <- data.frame(var_est=var_ests, if_mean=if_means,
                    r_n=r_ns, n=ns)
  
  # Toy example: scatterplots
  ggplot(df2, aes(x=var_est, y=if_mean)) + geom_point(alpha=0.3) +
    facet_grid(cols=dplyr::vars(n)) + ylim(c(-1e-4,1e-4))
  ggplot(df2, aes(x=var_est, y=r_n)) + geom_point(alpha=0.3) +
    facet_grid(cols=dplyr::vars(n))
  ggplot(df2, aes(x=if_mean, y=r_n)) + geom_point(alpha=0.3) +
    facet_grid(cols=dplyr::vars(n))
  ggplot(df2, aes(x=n, y=sqrt(n)*r_n)) + geom_point(alpha=0.3)
  
}

# Plotting VE
if (F) {
  
  x <- seq(0,1,0.1)
  ve <- (function(x) {
    0.8 + 0.15*x
  })(x)
  ve_to_risk <- Vectorize(function(ve) {
    risk_p <- 0.065 # Moderna
    # risk_p <- 0.048 # Janssen
    risk_p*(1-ve)
  })
  risk <- ve_to_risk(ve)
  
  ggplot(data.frame(x=x,y=ve), aes(x=x,y=y)) +
    geom_line() +
    ylim(c(0,1)) +
    labs(y="VE")
  ggplot(data.frame(x=x,y=(1-ve)), aes(x=x,y=y)) +
    geom_line() +
    ylim(c(0,1)) +
    labs(y="1-VE")
  ggplot(data.frame(x=x,y=(1-ve)), aes(x=x,y=y)) +
    geom_line() +
    scale_y_continuous(trans='log10') +
    # ylim(c(0,1)) +
    labs(y="1-VE (log scale))")
  ggplot(data.frame(x=x,y=risk), aes(x=x,y=y)) +
    geom_line() +
    ylim(c(0,0.1)) +
    labs(y="Risk")
  
}

# Debugging eta_n runtime
if (F) {
  
  params = list(g_n_type="true", S_n_type="true", omega_n_type="true")
  C <- list(
    points=seq(0,1,0.02), alpha_1=0.5, alpha_2=0.7, t_e=200,
    appx=list(t_e=10, w_tol=25, y_star=1, a=0.01)
  )
  L <- list(
    n=200, alpha_3=-2, dir="decr", sampling="two-phase (72%)",
    sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
    distr_A="Unif(0,1)", edge="expit 0.2", surv_true="Cox PH"
  )
  dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                            L$sc_params, L$sampling, L$dir)
  
  # Setup
  .default_params <- list(
    var="asymptotic", ecdf_type="step", g_n_type="binning",
    S_n_type="Super Learner", omega_n_type="estimated", cf_folds=1
  )
  for (i in c(1:length(.default_params))) {
    if (is.null(params[[names(.default_params)[i]]])) {
      params[[names(.default_params)[i]]] <- .default_params[[i]]
    }
  }
  a_lims <- c(min(dat_orig$a,na.rm=T),max(dat_orig$a,na.rm=T))
  a_shift <- -1 * a_lims[1]
  a_scale <- 1/(a_lims[2]-a_lims[1])
  dat_orig$a <- (dat_orig$a+a_shift)*a_scale
  dat_orig$a <- round(dat_orig$a, -log10(C$appx$a))
  for (i in c(1:length(dat_orig$w))) {
    rnd <- 8
    tol <- C$appx$w_tol
    n_unique <- tol + 1
    while(n_unique>tol) {
      rnd <- rnd - 1
      n_unique <- length(unique(round(dat_orig$w[,i],rnd)))
    }
    dat_orig$w[,i] <- round(dat_orig$w[,i],rnd)
  }
  n_orig <- length(dat_orig$delta)
  dat <- ss(dat_orig, which(dat_orig$delta==1))
  weights <- dat$weights
  vlist <- create_val_list(dat_orig)
  S_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type)
  Sc_n <- construct_S_n(dat, vlist$S_n, type=params$S_n_type, csf=TRUE)
  
  # !!!!! New constructor: omega_n
  construct_omega_n <- function(vals=NA, S_n, Sc_n, type="estimated") {
    H_n <- function(t,w,a) { -1 * log(S_n(t,w,a)) }
    fnc <- function(w,a,y_star,delta_star) {
      k <- round(min(y_star,C$t_e))
      if (k==0) { integral <- 0 } else {
        i <- c(1:k)
        w_long <- as.data.frame(matrix(rep(w,length(i)), ncol=length(w), byrow=T))
        a_long <- rep(a,length(i))
        integral <- 0.5 * sum(
          ( H_n(i,w_long,a_long) - H_n(i-1,w_long,a_long) ) * (
            ( S_n(i,w_long,a_long) * Sc_n(i,w_long,a_long) )^-1 +
              ( S_n(i-1,w_long,a_long) * Sc_n(i-1,w_long,a_long))^-1
          )
        )
      }

      return(S_n(C$t_e,w,a) * (
        (delta_star*as.integer(y_star<=C$t_e))/(S_n(k,w,a)*Sc_n(k,w,a))-integral
      ))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(2,1,1,1), vals=vals))
  }

  # # !!!!! New constructor: eta_n
  # construct_eta_n <- function(dat, vals=NA, S_n) {
  #   n_orig <- sum(dat$weights)
  #   fnc <- function(x,w) {
  #     w_long <- as.data.frame(matrix(rep(w,length(dat$a)), ncol=length(w), byrow=T))
  #     return(
  #       (1/n_orig) * sum(
  #         dat$weights * as.integer(dat$a<=x)*(1-S_n(rep(C$t_e,length(dat$a)),w_long,dat$a))
  #       )
  #     )
  #   }
  #   return(construct_superfunc(fnc, aux=NA, vec=c(1,2), vals=vals))
  # }
  
  # !!!!! New constructor calls
  # eta_n <- construct_eta_n(dat, vlist$AW_grid, S_n)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n,
                               type=params$omega_n_type)
  
  time_omega_n <- system.time({
    for (x in seq(0,1,0.01)) {
      for (w1 in seq(0,1,0.1)) {
        z <- omega_n(c(w1,0),x,10,0) # C$t_e
      }
    }
  })
  time_eta_n <- system.time({
    for (x in seq(0,1,0.01)) {
      for (w1 in seq(0,1,0.1)) {
        z <- eta_n(x,c(w1,0))
      }
    }
  })
  print(time_omega_n)
  print(time_eta_n)
  
}

# Cross-validation wrapper
if (F) {
  
  # CV prep
  n_folds <- 5
  folds <- sample(cut(c(1:nrow(df)), breaks=n_folds, labels=FALSE))
  bws <- seq(2*C$appx$a, 0.3, length.out=30)
  
  best <- list(bw=999, sum_sse=999)
  for (bw in bws) {
    
    sum_sse <- 0
    for (i in c(1:n_folds)) {
      
      df_train <- df[-which(folds==i),]
      df_test <- df[which(folds==i),]
      
      ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal", bandwidth=bw)
      reg <- Vectorize(function(a) {
        index <- which.min(abs(a-ks$x))
        return(ks$y[index])
      })
      
      sum_sse <- sum_sse + sum((reg(df_test$a)-df_test$po)^2, na.rm=T)
      
    }
    
    # print(paste0("bw: ",bw,", sum_sse: ",sum_sse))
    if (sum_sse<best$sum_sse || best$sum_sse==999) {
      best$bw <- bw
      best$sum_sse <- sum_sse
    }
    
    
  }
  
  # Construct optimal function from true data
  ks <- ksmooth(x=df_train$a, y=df_train$po, kernel="normal",
                bandwidth=best$bw)
  reg <- Vectorize(function(a) {
    index <- which.min(abs(a-ks$x))
    return(ks$y[index])
  })
  
  # !!!!!
  # x <- runif(100)
  # y <- sin(x*10)+rnorm(100,sd=0.1)+1
  # df=data.frame(a=x,po=y)
  print(best)
  grid <- seq(0,1,0.01)
  y_range <- c(0,1)
  ggplot(df, aes(x=a, y=po)) +
    geom_point() +
    geom_line(
      aes(x=x, y=y),
      data.frame(x=grid, y=reg(grid))
    ) +
    ylim(y_range) +
    labs(title=paste0("Regression (CV; ",y_range[1],"--",y_range[2],")"))
  
}

# Checking variance components
if (F) {
  
  sim %>% summarize(
    mean = list(
      list(x="psi1psi2_var", name="var_est_psi1psi2"),
      list(x="Psi_1_var_est", name="var_est_psi1"),
      list(x="Psi_2_var_est", name="var_est_psi2")
    ),
    var = list(
      list(x="psi1psi2", name="var_emp_psi1psi2"),
      list(x="Psi_1_est", name="var_emp_psi1"),
      list(x="Psi_2_est", name="var_emp_psi2")
    )
  )
  
}

# Debugging infl_fn_1
if (F) {
  
  ######################################################.
  ##### Old plot of estimated Var against true Var #####
  ######################################################.
  
  if (F) {
    # Var estimates (True Var is 0.00025)
    ggplot(data.frame(x=Psi_6_var_est), aes(x=x)) +
      geom_histogram(alpha=0.7) +
      geom_vline(xintercept=var(Psi_6_est), color="forestgreen", linetype="dashed")
    
    # Var estimates (True Var is ???)
    print(paste("Observed variance:", round(var(Psi_4_est),7)))
    print(paste("Estimated variance:", round(mean(Psi_4_var_est),7)))
  }
  
  
  
  #################.
  ##### Psi_6 #####
  #################.
  
  # IF constructor
  construct_infl_fn_6 <- function(Phi_n) {
    fnc <- function(a_i) {
      as.integer(a_i<=0.5) - Phi_n(0.5)
    }
    return(fnc) # return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_6_est <- rep(NA, n_reps)
  Psi_6_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF
    n_orig <- 1000
    a <- round(runif(n_orig), 3)
    Phi_n <- ecdf(a)
    
    # Construct estimator, IF, and variance estimator
    Psi_6_est[i] <- Phi_n(0.5)
    infl_fn_Psi_6 <- construct_infl_fn_6(Phi_n)
    Psi_6_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_6(a))^2)
    
    # Hypothesis test (Null: Phi_0(0.5)=0.5)
    Phi_0 <- 0.5
    test_stat <- (Psi_6_est[i]-Phi_0)^2/Psi_6_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_4 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_4 <- function(a,Phi_n,Gamma_0) {
    a_j <- a
    piece_1 <- Gamma_0(a_j)
    piece_2 <- 2*mean(Phi_n(a_j)*Gamma_0(a_j))
    fnc <- function(a_i) {
      Phi_n(a_i)*Gamma_0(a_i) + mean(as.integer(a_i<=a_j)*piece_1) - piece_2
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_4_est <- rep(NA, n_reps)
  Psi_4_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF
    n_orig <- 1000
    a <- round(runif(n_orig), 3)
    Phi_n <- ecdf(a)
    
    # Construct estimator, IF, and variance estimator
    Psi_4_est[i] <- mean(Phi_n(a)*Gamma_0(a))
    infl_fn_Psi_4 <- construct_infl_fn_4(a,Phi_n,Gamma_0)
    Psi_4_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_4(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.13335
    test_stat <- (Psi_4_est[i]-Phi_0)^2/Psi_4_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_7 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_7 <- function(a,Phi_n,Gamma_0,lambda_3n) {
    a_j <- a
    piece_1 <- (Phi_n(a_j))^2
    piece_2 <- mean(Phi_n(a_j)*Gamma_0(a_j))
    piece_3 <- Gamma_0(a_j)
    fnc <- function(a_i) {
      (3*mean(as.integer(a_i<=a_j)*piece_1) + (Phi_n(a_i))^3 - 6*lambda_3n) *
        piece_2 +
        lambda_3n * (
          mean(as.integer(a_i<=a_j)*piece_3) + Phi_n(a_i)*Gamma_0(a_i)
        )
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_7_est <- rep(NA, n_reps)
  Psi_7_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_3n <- mean((Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_7_est[i] <- mean(lambda_3n*Phi_n(a)*Gamma_0(a))
    infl_fn_Psi_7 <- construct_infl_fn_7(a,Phi_n,Gamma_0,lambda_3n)
    Psi_7_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_7(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.03333
    test_stat <- (Psi_7_est[i]-Phi_0)^2/Psi_7_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_8 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  construct_infl_fn_8 <- function(a,Phi_n,Gamma_0,lambda_2n) {
    a_j <- a
    piece_1 <- Phi_n(a_j)
    piece_2 <- mean((Phi_n(a_j))^2*Gamma_0(a_j))
    piece_3 <- Phi_n(a_j)*Gamma_0(a_j)
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_1) + (Phi_n(a_i))^2 - 6*lambda_2n) * 
        piece_2 +
        lambda_2n *
        (2*mean(as.integer(a_i<=a_j)*piece_3) + (Phi_n(a_i))^2*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 100
  Psi_8_est <- rep(NA, n_reps)
  Psi_8_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_2n <- mean((Phi_n(a))^2)
    
    # Construct estimator, IF, and variance estimator
    Psi_8_est[i] <- mean(lambda_2n*(Phi_n(a))^2*Gamma_0(a))
    infl_fn_Psi_8 <- construct_infl_fn_8(a,Phi_n,Gamma_0,lambda_2n)
    Psi_8_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_8(a))^2)
    
    # Hypothesis test
    Phi_0 <- 0.031746
    test_stat <- (Psi_8_est[i]-Phi_0)^2/Psi_8_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  #################.
  ##### Psi_1 #####
  #################.
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  rho_n <- function(a,Phi_n,Gamma_0,x) {
    mean( (Phi_n(a))^x * Gamma_0(a) )
  }
  construct_infl_fn_1 <- function(a,Phi_n,Gamma_0,lambda_2n,lambda_3n) {
    a_j <- a
    rho_1 <- rho_n(a,Phi_n,Gamma_0,1)
    rho_2 <- rho_n(a,Phi_n,Gamma_0,2)
    # xi_01 <- construct_xi_n(a,Phi_n,Gamma_0,0,1)
    # xi_20 <- construct_xi_n(a,Phi_n,Gamma_0,2,0)
    # xi_10 <- construct_xi_n(a,Phi_n,Gamma_0,1,0)
    # xi_11 <- construct_xi_n(a,Phi_n,Gamma_0,1,1)
    piece_01 <- Gamma_0(a_j)
    piece_20 <- (Phi_n(a_j))^2
    piece_10 <- Phi_n(a_j)
    piece_11 <- Phi_n(a_j) * Gamma_0(a_j)
    
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_10)+(Phi_n(a_i))^2-6*lambda_2n)*rho_2 +
        lambda_2n*(2*mean(as.integer(a_i<=a_j)*piece_11)+(Phi_n(a_i))^2*Gamma_0(a_i)) -
        (3*mean(as.integer(a_i<=a_j)*piece_20)+(Phi_n(a_i))^3-6*lambda_3n)*rho_1 -
        lambda_3n*(mean(as.integer(a_i<=a_j)*piece_01)+Phi_n(a_i)*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 1000
  Psi_1_est <- rep(NA, n_reps)
  Psi_1_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    # a <- round(runif(n_orig), 3)
    a <- runif(n_orig)
    Phi_n <- ecdf(a)
    lambda_2n <- mean((Phi_n(a))^2)
    lambda_3n <- mean((Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_1_est[i] <- mean(
      lambda_2n*(Phi_n(a))^2*Gamma_0(a) - lambda_3n*Phi_n(a)*Gamma_0(a)
    )
    infl_fn_Psi_1 <- construct_infl_fn_1(a,Phi_n,Gamma_0,lambda_2n,lambda_3n)
    Psi_1_var_est[i] <- (1/n_orig^2) * sum((infl_fn_Psi_1(a))^2)
    
    # Hypothesis test
    Phi_0 <- -0.0015878
    test_stat <- (Psi_1_est[i]-Phi_0)^2/Psi_1_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
  
  
  ################################.
  ##### Psi_1 (with weights) #####
  ################################.
  
  # Setup
  
  C <- list(points=seq(0,1,0.02), alpha_1=0.5, alpha_2=0.7, t_e=200,
            appx=list(t_e=10,w1=0.1,y_star=1,a=0.01))
  
  # Approximate Gamma_0
  Gamma_0 <- function(x) { sqrt(x)/3 }
  
  # IF constructor
  rho_n <- function(a,Phi_n,Gamma_0,x) {
    mean( (Phi_n(a))^x * Gamma_0(a) )
  }
  construct_infl_fn_1 <- function(a,Phi_n,Gamma_0,lambda_2n,lambda_3n) {
    a_j <- a
    rho_1 <- rho_n(a,Phi_n,Gamma_0,1)
    rho_2 <- rho_n(a,Phi_n,Gamma_0,2)
    # xi_01 <- construct_xi_n(a,Phi_n,Gamma_0,0,1)
    # xi_20 <- construct_xi_n(a,Phi_n,Gamma_0,2,0)
    # xi_10 <- construct_xi_n(a,Phi_n,Gamma_0,1,0)
    # xi_11 <- construct_xi_n(a,Phi_n,Gamma_0,1,1)
    piece_01 <- Gamma_0(a_j)
    piece_20 <- (Phi_n(a_j))^2
    piece_10 <- Phi_n(a_j)
    piece_11 <- Phi_n(a_j) * Gamma_0(a_j)
    
    fnc <- function(a_i) {
      (2*mean(as.integer(a_i<=a_j)*piece_10)+(Phi_n(a_i))^2-6*lambda_2n)*rho_2 +
        lambda_2n*(2*mean(as.integer(a_i<=a_j)*piece_11)+(Phi_n(a_i))^2*Gamma_0(a_i)) -
        (3*mean(as.integer(a_i<=a_j)*piece_20)+(Phi_n(a_i))^3-6*lambda_3n)*rho_1 -
        lambda_3n*(mean(as.integer(a_i<=a_j)*piece_01)+Phi_n(a_i)*Gamma_0(a_i))
    }
    return(construct_superfunc(fnc, aux=NA, vec=c(1), vals=NA))
  }
  
  n_reps <- 10
  Psi_1_est <- rep(NA, n_reps)
  Psi_1_var_est <- rep(NA, n_reps)
  reject <- rep(NA, n_reps)
  for (i in c(1:n_reps)) {
    
    # Generate data, construct ECDF and lambda_3
    n_orig <- 1000
    dat_orig <- generate_data(1000, 0, "Unif(0,1)", "none", "Cox PH",
                              list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
                              "iid", "decr") # "two-phase (72%)"
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    a <- dat$a
    weights <- dat$weights
    
    Phi_n <- construct_Phi_n(dat, type="step")
    lambda_2n <- (1/n_orig) * sum(weights*(Phi_n(a))^2)
    lambda_3n <- (1/n_orig) * sum(weights*(Phi_n(a))^3)
    
    # Construct estimator, IF, and variance estimator
    Psi_1_est[i] <- (1/n_orig) * sum(weights*(
      lambda_2n*(Phi_n(a))^2*Gamma_0(a) - lambda_3n*Phi_n(a)*Gamma_0(a)
    ))
    infl_fn_Psi_1 <- construct_infl_fn_1(a,Phi_n,Gamma_0,lambda_2n,lambda_3n)
    Psi_1_var_est[i] <- (1/n_orig^2) * sum(weights*(infl_fn_Psi_1(a))^2)
    
    # Hypothesis test
    Phi_0 <- -0.0015878
    test_stat <- (Psi_1_est[i]-Phi_0)^2/Psi_1_var_est[i]
    p_val <- pchisq(test_stat, df=1, lower.tail=FALSE)
    reject[i] <- as.integer(p_val<0.05)
    
  }
  
  # False positive rate
  mean(reject)
  
}



# Testing infl_fn_1
if (F) {
  
  L <- list(n=1000, alpha_3=-2, dir="decr",
            sc_params=list(lmbd=1e-3, v=1.5, lmbd2=5e-5, v2=1.5),
            distr_A="Unif(0,1)", edge="none", surv_true="Cox PH", # Unif(0,1) N(0.5,0.04)
            ecdf_type="true", sampling="iid", # two-phase (72%)
            estimator=list(est="Grenander",params=params)
  )
  
  n_reps <- 50
  ests <- rep(NA,n_reps)
  vars <- rep(NA,n_reps)
  for (i in c(1:n_reps)) {
    
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    n_orig <- length(dat_orig$delta)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    weights <- dat$weights
    vlist <- create_val_list(dat_orig)
    
    Phi_n <- construct_Phi_n(dat, type=params$ecdf_type)
    lambda_2 <- lambda(dat,2,Phi_n)
    lambda_3 <- lambda(dat,3,Phi_n)
    xi_n <- construct_xi_n(Phi_n, lambda_2, lambda_3)
    rho_n <- function(x) { 0 }
    
    # Partial stat
    Gamma_0 <- Vectorize(function(a) {
      Theta_true <- attr(dat_orig,"Theta_true")
      grid <- seq(0,1,0.02)
      index <- which.min(abs(a-seq(0,1,0.02)))
      return(Theta_true[index])
    })
    partial_est <- (1/n_orig) * sum( weights * (
      (lambda_2*(Phi_n(dat$a))^2 - lambda_3*Phi_n(dat$a)) * Gamma_0(dat$a)
    ))
    
    # Variance estimate of partial stat
    infl_fn_1b <- construct_infl_fn_1(dat, Gamma_0, Phi_n, xi_n, rho_n,
                                      lambda_2, lambda_3)
    partial_var <- (1/n_orig^2) * sum((weights * (infl_fn_1b(dat$a)))^2)
    
    ests[i] <- partial_est
    vars[i] <- partial_var
    
    print(paste0("Rep ", i, " of ", n_reps))
    
  }
  
  # Process results
  print(paste0("sd(ests), iid:", sd(ests)))
  print(paste0("mean(sqrt(vars)), iid:", mean(sqrt(vars))))
  ggplot(data.frame(x=ests), aes(x=x)) + geom_histogram() + labs(title="ests")
  ggplot(data.frame(x=vars), aes(x=x)) + geom_histogram() + labs(title="var")
  
}

# Figure out why returned functions are so large
if (F) {
  
  S_n <- readRDS("705 (SL, marker 8)/S_n.rds")
  
  objs <- ls(environment(get("fnc",envir=environment(gamma_n))))
  for (obj in objs) {
    print(obj)
    print(object.size(get(
      obj,
      envir = environment(get("fnc",envir=environment(gamma_n)))
    )))
  }

}

# Debugging
if (F) {
  
  # Set up containers
  df <- data.frame(
    "rep" = integer(),
    "point" = double(),
    "theta_n" = double(),
    "Gamma_n" = double(),
    "Psi_n" = double()
  )
  df_list <- list()
  
  # Run for-loop
  for (i in 1:2) {
    
    print(paste0("Rep ",i,": ",Sys.time()))
    set.seed(i)
    dat_orig <- generate_data(L$n, L$alpha_3, L$distr_A, L$edge, L$surv_true,
                              L$sc_params, L$sampling, L$dir)
    dat <- ss(dat_orig, which(dat_orig$delta==1))
    vlist <- create_val_list(dat_orig)
    vlist$AW_grid <- NA; vlist$omega <- NA; vlist$W_grid <- NA;
    Gamma_os_n <- construct_Gamma_os_n(dat, vlist$A_grid, omega_n, S_n, g_n)
    Psi_n <- Vectorize(function(x) {
      -1 * Gamma_os_n(round(Phi_n_inv(x), -log10(C$appx$a)))
    })
    gcm <- gcmlcm(
      x = seq(0,1,C$appx$a),
      y = Psi_n(seq(0,1,C$appx$a)), # rev(Psi_n(seq(0,1,C$appx$a)))
      type = "gcm"
    )
    dGCM <- Vectorize(function(x) {
      # The round deals with a floating point issue
      index <- which(round(x,5)<=gcm$x.knots)[1]-1
      if (index==0) { index <- 1 }
      return(gcm$slope.knots[index])
    })
    theta_n_Gr <- Vectorize(function(x) { min(max(-1 * dGCM(Phi_n(x)),0),1) })
    
    for (j in 1:51) {
      df[nrow(df)+1,] <- c(
        i, points[j], theta_n_Gr(points[j]), Gamma_os_n(points[j]),
        Psi_n(points[j])
      )
    }
    
    df_list[[i]] <- list(
      dat_orig = dat_orig,
      Phi_n = Phi_n,
      Phi_n_inv = Phi_n_inv,
      S_n = S_n,
      Sc_n = Sc_n,
      g_n = g_n,
      omega_n = omega_n,
      Gamma_os_n = Gamma_os_n
    )
    
  }
  
  # Plot results
  df$Psi_n <- -1 * df$Psi_n
  ggplot(df, aes(x=point, y=Psi_n, group=rep)) + # Gamma_n Psi_n
    geom_line(alpha=0.4) +
    ylim(c(0,0.4)) # 0.4
  
  # !!!!!
  ggplot(
    data.frame(
      x = c(gcm$x.knots, seq(0,1,0.01)),
      y = c(gcm$y.knots, rev(Psi_n(seq(0,1,0.01)))),
      grp = c(rep("gcm",length(gcm$x.knots)),rep("Psi_n",101))
    ),
    aes(x=x, y=y, color=factor(grp))) +
    geom_line()
  
  # !!!!!
  identical(
    sim$results_complex$sim_uid_1$dat_orig,
    df_list[[1]]$dat_orig
  )
  sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.1))
  df_list[[1]]$Phi_n(seq(0,1,0.1))
  
  # !!!!!
  sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.1))
  ggplot(
    data.frame(
      x = seq(0,1,0.01),
      y = sim$results_complex$sim_uid_1$Phi_n(seq(0,1,0.01))
    ),
    aes(x=x, y=y)
  ) + geom_line()
  df_list[[1]]$Phi_n(seq(0,1,0.1))
  ggplot(
    data.frame(
      x = seq(0,1,0.01),
      y = df_list[[1]]$Phi_n(seq(0,1,0.01))
    ),
    aes(x=x, y=y)
  ) + geom_line()
  
  # return(function(a) { ptruncnorm(a, a=0, b=1, mean=0.5, sd=0.2) })
  
  
}

# Unit tests for superfunc
if (F) {
  
  f <- function(a,b,c) {
    if (T) { Sys.sleep(1) }
    (a-b)+sum(c)/5
  }
  fs <- construct_superfunc(f, vec=c(0,1,2))
  # f(9,2,c(10,40))
  # f(9,3,c(20,50))
  # f(9,4,c(30,60))
  # system.time({x<-fs(9,2,c(10,40)); print(x);})
  # system.time({x<-fs(9,3,c(20,50)); print(x);})
  # system.time({x<-fs(9,4,c(30,60)); print(x);})
  # system.time({x<-fs(9,c(2,3,4),data.frame(c(10,20,30),c(40,50,60))); print(x)})
  
  # Test rounding
  fs(9,2,c(10,40))
  fs(9,c(2,3,4),data.frame(c(10,20,30),c(40,50,60)))
  
}

# Deconstruct vectorize
if (F) {
  
  # # !!!!!
  # (function(){
  #   (function(){
  #     parent.frame()
  #   })()
  #   # parent.frame()
  # })()
  
  Vectorize2 <- function (FUN, vectorize.args = arg.names, SIMPLIFY = TRUE, 
            USE.NAMES = TRUE) 
  {
    arg.names <- as.list(formals(FUN))
    arg.names[["..."]] <- NULL
    arg.names <- names(arg.names)
    vectorize.args <- as.character(vectorize.args)
    if (!length(vectorize.args)) 
      return(FUN)
    if (!all(vectorize.args %in% arg.names)) 
      stop("must specify names of formal arguments for 'vectorize'")
    collisions <- arg.names %in% c("FUN", "SIMPLIFY", "USE.NAMES", 
                                   "vectorize.args")
    if (any(collisions)) 
      stop(sQuote("FUN"), " may not have argument(s) named ", 
           paste(sQuote(arg.names[collisions]), collapse = ", "))
    rm(arg.names, collisions)
    (function() {
      FUNV <- function() {
        print("parent.frame()")
        print(parent.frame())
        args <- lapply(as.list(match.call())[-1L], eval, 
                       parent.frame())
        print("args")
        print(args)
        names <- if (is.null(names(args))) 
          character(length(args))
        else names(args)
        print("names")
        print(names)
        dovec <- names %in% vectorize.args
        print("dovec")
        print(dovec)
        print("args[dovec]")
        print(args[dovec])
        print("args[!dovec]")
        print(args[!dovec])
        do.call("mapply", c(FUN = FUN, args[dovec], MoreArgs = list(args[!dovec]), 
                            SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))
      }
      formals(FUNV) <- formals(FUN)
      environment(FUNV) <- parent.env(environment())
      FUNV
    })()
  }
  
  x <- Vectorize2(function(a,b,c,d) {a+b+c+d}, vectorize.args=c("a","b"))
  x(1,2,3,4)
  
}

# Deconstruct memoise
if (F) {
  
  mem2 <- function (f) {
    memo_f <- function(...) {
      encl <- parent.env(environment())
      v <- encl$`_val`
    }
    memo_f_env <- new.env(parent=environment(f))
    memo_f_env$`_val` <- 3
    environment(memo_f) <- memo_f_env
    memo_f
  }
  
  f1 <- function() { print("there") }
  f2 <- mem2(f1)
  f2()
  
}

# construct_superfunc: fixing environments issue
if (F) {
  
  x <- 1
  f0 <- function() { return(x) }
  f1 <- memoise(f)
  f2 <- construct_superfunc(f)
  f3 <- Vectorize(f)
  
  f2 <- function() {
    x <- 99
    f()
    # f1()
  }
  f2()
  
  f2(x=4) # This should give same answer as f2b
  f2b <- function(x) {
    f1b <- function() { return(x) }
    f1b(a=3)
  }
  f2b(x=4)
  
  # !!!!! Try both passing function anonymously and named
  # !!!!! Test use of aux
  
}

# Smoothed ecdf
if (F) {
  
  # Create dataset
  n_orig <- 100
  dat <- data.frame(weights=rep(1,n_orig), a=runif(n_orig, min=0, max=1))
  
  # Calculate estimators
  dat %<>% arrange(a)
  vals_x <- unique(dat$a)
  vals_y <- c()
  for (j in 1:length(vals_x)) {
    indices <- which(dat$a==vals_x[j])
    weights_j <- dat$weights[indices]
    new_y_val <- (1/n_orig) * sum(weights_j)
    vals_y <- c(vals_y, new_y_val)
  }
  vals_y <- cumsum(vals_y)
  
  vals_x_ <- c(vals_x[1], vals_x[1:(length(vals_x)-1)]+diff(vals_x)/2,
               vals_x[length(vals_x)])
  vals_y_ <- c(0, vals_y[1:(length(vals_y)-1)], 1)
  
  rval <- approxfun(vals_x_, vals_y_, method="linear", yleft=0,
                    yright=1, ties="ordered")
  
  rval_inv <- approxfun(vals_y_, vals_x_, method="linear", yleft=min(vals_x),
                        yright=max(vals_x), ties="ordered")
  
  grid <- seq(0,1,0.01)
  jitter <- 0
  # jitter <- rnorm(2*length(grid), sd=0.003)
  df <- data.frame(
    x = rep(grid,2),
    y = c(rval(grid), rval_inv(grid)) + jitter,
    which = rep(c("ecdf", "smoothed"), each=length(grid))
  )
  ggplot(df, aes(x=x, y=y, color=which)) + geom_line() +
    geom_abline(slope=1, color="grey")
  
  rval(rval_inv(seq(0,1,0.1)))
  rval_inv(rval(seq(0,1,0.1)))
  
}

# Histogram for edge values
if (F) {
  
  r_cox <- filter(sim$results, surv_true=="Cox PH")
  r_com <- filter(sim$results, surv_true=="Complex")
  cox <- r_cox$est_0.0
  com <- r_com$est_0.0
  tr_cox <- r_cox[1,"theta_0.0"]
  tr_com <- r_com[1,"theta_0.0"]
  
  (mean(cox)-tr_cox)/mean(cox)
  (mean(com)-tr_com)/mean(com)
  
  df <- data.frame(
    x = c(cox,com),
    grp = c(rep("Cox PH", length(cox)),
            rep("Complex", length(com)))
  )
  
  ggplot(df, aes(x=x, group=grp, fill=factor(grp))) +
    geom_histogram(color="white") +
    facet_wrap(~grp, ncol=2) +
    geom_vline(
      aes(xintercept=tr),
      data = data.frame(tr=c(tr_cox,tr_com), grp=c("Cox PH","Complex"))
    )
  
}

# Monotone spline for derivative
if (F) {
  
  x <- sort(runif(5))
  y <- sort(runif(5))
  
  theta_n_smoothed <- splinefun(x=x, y=y, method="monoH.FC")
  
  ggplot(
    data.frame(x=seq(0,1,0.01), y=theta_n_smoothed(seq(0,1,0.01))),
    aes(x=x,y=y)) +
    geom_line() +
    labs(title=mmm) +
    geom_point(data=data.frame(x=x, y=y))
  
}

# Use softer cutoffs instead of step functions
if (F) {
  
  ff <- function(w1) { as.integer(abs(w1-0.5)<0.11) }
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=ff)
  
  ff2 <- function(w1) { pmax(0,1-4*abs(w1-0.5)) }
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) + stat_function(fun=ff2)
  
}

#
if (F) {
  
  # Generate data
  dat_orig <- generate_data(
    n = 2000,
    alpha_3 = -3,
    distr_A = "N(1.5+w1,1.5+w2)", # "Unif(0,1)"
    edge = "none",
    surv_true = "Complex",
    sc_params = L$sc_params,
    sampling = "two-phase (72%)",
    dir = "decr"
  )
  
  # Prep
  n_orig <- nrow(dat_orig)

  # Construct dataframes of values to pre-compute functions on
  vlist <- create_val_list(dat_orig)
  
  # Construct component functions
  S_n <- construct_S_n(dat_orig, vlist$S_n, type="true")
  Sc_n <- construct_S_n(dat_orig, vlist$S_n, type="true",
                        csf=TRUE)
  f_aIw_n <- construct_f_aIw_n(dat_orig, vlist$AW_grid,
                               type="true", k=15)
  f_a_n <- construct_f_a_n(dat_orig, vlist$A_grid, f_aIw_n)
  g_n <- construct_g_n(f_aIw_n, f_a_n)
  omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
  
  offset_25 <- (dat$weights *
     as.integer(dat$a<=0.25) *
     omega_n(dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star)
  ) / g_n(dat$a,dat$w1,dat$w2)
  offset_75 <- (dat$weights *
                  as.integer(dat$a<=0.75) *
                  omega_n(dat$w1,dat$w2,dat$a,dat$y_star,dat$delta_star)
  ) / g_n(dat$a,dat$w1,dat$w2)
  print(mean(offset_25))
  print(mean(offset_75))
  # print(head(sort(offset)))
  # print(tail(sort(offset)))
  # print(mean(sort(offset)[2:length(offset)]))
    
}

# Variance of one-step edge estimator
if (F) {
  
  ##### Edge mass #####
  
  n_reps <- 30
  edge <- "expit"
  ests <- c()
  sigma2s <- c()
  
  for (i in 1:n_reps) {
    
    # Generate data
    dat_orig <- generate_data(
      n = 1000,
      alpha_3 = -3,
      distr_A = "N(0.5,0.04)",
      edge = edge,
      surv_true = "Cox PH",
      sc_params = L$sc_params,
      sampling = "two-phase (72%)",
      dir = "decr"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    g_sn <- construct_g_sn(dat, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat, g_sn, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat, g_sn, S_n, omega_n, theta_os_n_est)
    
    ests <- c(ests, theta_os_n_est)
    sigma2s <- c(sigma2s, sigma2_os_n_est/n_orig)
    
  }
  
  v1 <- var(ests)
  m1 <- mean(sigma2s)
  
  ##### Edge mass #####
  
  n_reps <- 30
  edge <- "none"
  ests <- c()
  sigma2s <- c()
  
  for (i in 1:n_reps) {
    
    # Generate data
    dat_orig <- generate_data(
      n = 1000,
      alpha_3 = -3,
      distr_A = "N(0.5,0.04)",
      edge = edge,
      surv_true = "Cox PH",
      sc_params = L$sc_params,
      sampling = "two-phase (72%)",
      dir = "decr"
    )
    
    n_orig <- nrow(dat_orig)
    vlist <- create_val_list(dat_orig)
    
    S_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type)
    Sc_n <- construct_S_n(dat_orig, vlist$S_n, type=params$S_n_type,
                          csf=TRUE)
    omega_n <- construct_omega_n(vlist$omega, S_n, Sc_n)
    g_sn <- construct_g_sn(dat, vlist$W_grid, type="logistic")
    theta_os_n_est <- theta_os_n(dat, g_sn, S_n, omega_n)
    sigma2_os_n_est <- sigma2_os_n(dat, g_sn, S_n, omega_n, theta_os_n_est)
    
    ests <- c(ests, theta_os_n_est)
    sigma2s <- c(sigma2s, sigma2_os_n_est/n_orig)
    
  }
  
  v2 <- var(ests)
  m2 <- mean(sigma2s)
  
  print("Edge mass")
  print(v1)
  print(m1)
  print("No edge mass")
  print(v2)
  print(m2)

}

# Testing cross-fitting
if (F) {
  
  system.time({
    ests_reg <- est_curve(
      dat_orig = dat_tp,
      estimator = "Grenander",
      params = list(S_n_type=params$S_n_type, g_n_type=params$g_n_type, ci_type="logit",
                    cf_folds=1),
      points = C$points
    )
  })
  
  system.time({
    ests_cf <- est_curve(
      dat_orig = dat_tp,
      estimator = "Grenander",
      params = list(S_n_type=params$S_n_type, g_n_type=params$g_n_type, ci_type="logit",
                    cf_folds=3),
      points = C$points
    )
  })
  
  system.time({
    reject_reg <- test_2(
      dat_orig = dat_tp,
      alt_type = "incr",
      params = list(
        var = "asymptotic",
        S_n_type = params$S_n_type,
        g_n_type = params$g_n_type,
        est_known_nuis = FALSE,
        cf_folds = 1
      )
    )
  })
  
  system.time({
    reject_cf <- test_2(
      dat_orig = dat_tp,
      alt_type = "incr",
      params = list(
        var = "asymptotic",
        S_n_type = params$S_n_type,
        g_n_type = params$g_n_type,
        est_known_nuis = FALSE,
        cf_folds = 3
      )
    )
  })
  
  
}

# Does memoising work if function is passed?
if (F) {
  
  f <- function() {
    print(get("hey", envir=parent.frame()))
  }
  
  (function(hey) {
    f()
  })("there")
  
  fn_mem <- memoise(function(x) {
    Sys.sleep(2)
    return(9)
  })
  fn_mem(1)
  
  fn2 <- function(fn) {
    fn(1)
  }
  fn3 <- function(fn) {
    fn(1)
  }
  fn2(fn_mem)
  fn3(fn_mem)
  
}

# Fastest way of getting a vectorized/memoised function (part 2)
if (F) {
  
  construct_S_n_mem <- function(dat, vals) {
      model <- coxph(Surv(y_star, delta_star)~w1+w2+a, data=dat)
      c_1 <- model$coefficients[[1]]
      c_2 <- model$coefficients[[2]]
      c_3 <- model$coefficients[[3]]
      bh <- basehaz(model, centered=FALSE)
      H_0 <- c()
      for (t in 0:C$t_e) {
        index <- which.min(abs(bh$time-t))
        H_0[t+1] <- bh$hazard[index]
      }
      # return(Vectorize(memoise(function(t, w1, w2, a){
      #   return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
      # })))
      return(memoise(Vectorize(function(t, w1, w2, a){
        return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
      })))
  }
  
  construct_S_n_htab <- function(dat, vals) {
    model <- coxph(Surv(y_star, delta_star)~w1+w2+a, data=dat)
    c_1 <- model$coefficients[[1]]
    c_2 <- model$coefficients[[2]]
    c_3 <- model$coefficients[[3]]
    bh <- basehaz(model, centered=FALSE)
    H_0 <- c()
    for (t in 0:C$t_e) {
      index <- which.min(abs(bh$time-t))
      H_0[t+1] <- bh$hazard[index]
    }
    fn <- function(t, w1, w2, a){
      return(exp(-1*H_0[t+1]*exp(c_1*w1+c_2*w2+c_3*a)))
    }
    return (create_htab(fn, vals))
  }
  
  vals_S_n <- expand.grid(t=c(0:C$t_e),w1=seq(0,1,0.1),w2=c(0,1),a=seq(0,1,0.1))
  S_n_mem <- construct_S_n_mem(dat, vals_S_n)
  for (i in 1:nrow(vals_S_n)) {
    row <- vals_S_n[i,]
    x <- do.call(S_n_mem, as.list(as.numeric(row)))
  }
  S_n_htab <- construct_S_n_htab(dat, vals_S_n)
  
  v1 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        # key <- "200;0;0;0.1"
        key <- paste(x, collapse=";")
        val <- fn[[key]]
        return(val)
      }
    )
  }
  v2 <- function(fn, ...) {
    do.call("mapply", c(
      FUN = function(...) {
        # key <- "200;0;0;0.1"
        key <- paste(c(...), collapse=";")
        return(fn[[key]])
      },
      list(...)
    ))
  }
  v3 <- memoise(function(...) {
    do.call("mapply", c(
      FUN = function(...) {
        key <- paste(c(...), collapse=";")
        return(S_n_htab[[key]])
      },
      list(...)
    ))
  })
  
  S_n_mem(C$t_e,0,0,round(dat$a,1))[1:10]
  v1(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  v2(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  v3(C$t_e,0,0,round(dat$a,1))[1:10]
  v4("S_n_htab",C$t_e,0,0,round(dat$a,1))[1:10]
  
  microbenchmark({
    # S_n_mem(C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    S_n_mem(C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v1(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v1(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v2(S_n_htab,C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v3(C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  microbenchmark({
    # v2(S_n_htab,C$t_e,0,0,c(0.1,0.2,0.3))[1:10]
    v4("S_n_htab",C$t_e,0,0,round(dat$a,1))[1:10]
  }, times=100L)
  
}

# Fastest way of getting a vectorized/memoised function (part 1)
if (F) {
  
  n <- 4*10^4
  fn <- function(x,y) { x^2+y }
  
  # Method 1: memoise/Vectorize
  # mean: 838 microsec
  fn_1 <- Vectorize(memoise(fn))
  for(i in 1:n) {
    x <- fn_1(i,10)
  }
  microbenchmark({
    fn_1(c(100:200),10)
  }, times=100L)
  
  # Method 2: hash table
  # mean: 205 microsec
  htab <- new.env()
  for(i in 1:n) {
    htab[[paste(c(i,10), collapse=";")]] <- fn_1(i,10)
  }
  fn_2 <- function(...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        key <- paste(c(x), collapse=";")
        val <- htab[[key]]
        return(val)
      }
    )
  }
  microbenchmark({
    fn_2(c(100:200),10)
  }, times=100L)
  
  # Create, populate, and return hash table
  for (i in 1:n) {
    row <- vals[i,]
    key <- paste(row, collapse=";")
    htab[[key]] <- do.call(fn, as.list(as.numeric(row)))
    # v <- do.call(fn, as.list(as.numeric(row)))
    # if (is.null(v)) {
    #   stop(paste0("key: ",key))
    # }
    # htab[[key]] <- v
  }
  return (htab)
  
}

# !!!! TEMP
if (F) {
  
  # fn_mem <- memoise(Vectorize(function(x) {x^2}))
  # htab <- new.env()
  # 
  # for(i in 1:10000) {
  #   x <- fn_mem(i)
  #   row <- vals[i,]
  #   key <- paste(row, collapse=";")
  #   htab[[key]] <- do.call(fn, as.list(as.numeric(row)))
  # }
  
  S_n <- construct_S_n(dat, vals_S_n, type=params$S_n_type)
  
  
  
  
  S_n_htab_vc <- Vectorize(function(w,x,y,z) {
    S_n_htab[["200;0;0;0.1"]]
    # S_n_htab[[paste(c(w,x,y,z), collapse=";")]]
  })
  
  
  microbenchmark({
    S_n_htab_vc(200,0,0,rep(0.1,100))
    # S_n_htab[[paste0("200;0;0;",0.1)]]
    # v(S_n_htab,C$t_e,0,0,round(dat$a,1))
  }, times=100L)
  
  microbenchmark({
    S_n_mem(200,0,0,rep(0.1,100))
    # S_n_mem(C$t_e,0,0,round(dat$a,1))
  }, times=100L)
  
  (function(a,b,c) {
    mc <- match.call()
    called_args <- as.list(mc)[-1]
    print(called_args)
  })(1,2,3)

  function (t, w1, w2, a) 
  {
    mc <- match.call()
    encl <- parent.env(environment())
    called_args <- as.list(mc)[-1]
    default_args <- encl$`_default_args`
    default_args <- default_args[setdiff(names(default_args), 
                                         names(called_args))]
    called_args[encl$`_omit_args`] <- NULL
    args <- c(lapply(called_args, eval, parent.frame()), lapply(default_args, 
                                                                eval, envir = environment()))
    key <- encl$`_hash`(c(encl$`_f_hash`, args, lapply(encl$`_additional`, 
                                                       function(x) eval(x[[2L]], environment(x)))))
    res <- encl$`_cache`$get(key)
    if (inherits(res, "key_missing")) {
      mc[[1L]] <- encl$`_f`
      res <- withVisible(eval(mc, parent.frame()))
      encl$`_cache`$set(key, res)
    }
    if (res$visible) {
      res$value
    }
    else {
      invisible(res$value)
    }
  }
  
  
  
  
  
  library(stringi)
  
  v1 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[paste(c(x), collapse=";")]]
      }
    )
  }
  
  v2 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[paste(x, collapse=";")]]
      }
    )
  }
  
  v3 <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {
        fn[[stri_join(x,collapse=";")]]
      }
    )
  }
  
  microbenchmark({
    v1(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
  microbenchmark({
    v2(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
  microbenchmark({
    v3(S_n,C$t_e,round(w1,1),w2,round(dat_orig$a,1))
  }, times=100L)
  
}

# New hashing/memoising structure for creators
if (F) {
  
  # Creator option #1
  creator_1 <- function(vals) {
    
    # Declare function
    fn <- function(x,y) { x*y }
    
    # Set environment
    htab <- new.env()
    
    # Run function on vals
    for (i in 1:nrow(vals)) {
      row <- vals[i,]
      key <- paste(row, collapse=";")
      htab[[key]] <- do.call("fn", as.list(as.numeric(row)))
    }
    
    return (htab)
    
  }
  
  # Usage
  vals <- expand.grid(a=c(1:n), w=c(1:n))
  fn_1 <- creator_1(vals)
  z <- function(...) {
    paste(c(...), collapse=";")
  }
  v <- function(fn, ...) {
    apply(
      X = cbind(...),
      MARGIN = 1,
      FUN = function(x) {fn[[z(x)]]}
    )
  }
  
  fn_1[[z(1,2)]]
  fn_1[[z(4,5)]]
  fn_1[[z(7,8)]]
  v(fn_1,c(1,4,7),c(2,5,8))

  fn_1[[z(2,2)]]
  fn_1[[z(2,5)]]
  fn_1[[z(2,8)]]
  v(fn_1,3,c(2,5,8))
  
  n <- 10000
  vals <- data.frame(a=c(1:n), w=c(2:(n+1)))
  
  microbenchmark({
    fn_1[["77;98"]]
  }, times=1000L)
  microbenchmark({
    fn_1[[paste(c(77,98), collapse=";")]]
  }, times=1000L)
  microbenchmark({
    fn_1[[z(77,98)]]
  }, times=1000L)

  # memoise() vs. new.env()
  mem_fn <- memoise(fn)
  htab_old <- new.env()
  
  # Initialize hash structures
  # n= 200:  0.4 sec
  # n= 400:  1.8 sec
  # n= 800: 10.4 sec
  # n=1600: 10.4 sec
  n <- 1600
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        htab_old[[paste0(i,";",j)]] <- fn(i,j)
      }
    }
  })
  
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        x <- mem_fn(i,j)
      }
    }
  })
  
  # Run benchmarks
  microbenchmark({x <- mem_fn(5,6)}, times=1000L)
  microbenchmark({x <- htab[[paste0(5,";",6)]]}, times=1000L)
  
}

# Checking lambda factors with mixture distribution
if (F) {
  
  n <- 10^6
  prob <- 0.9 # 0.9
  mix <- rbinom(n, size=1, prob=prob)
  unif <- runif(n)
  a <- mix*unif
  F_a <- ecdf(a)
  mean((F_a(a))^2)
  mean((F_a(a))^3)
  
}

# Probability integral transform
if (F) {
  
  x <- rbeta(5,0.1,0.1)
  F_n <- ecdf(x)
  print(x)
  print(F_n(x))
  
  n <- 100
  i <- c(1:n)
  sum(i^2)
  n^3/3 + n^2/2 + n/6
  sum(i^3)
  n^4/4 + n^3/2 + n^2/4
  
}

# Benchmarking hashing/memoising structures
if (F) {
  
  library(memoise)
  library(microbenchmark)
  
  # Unmemoised function
  fn <- function(x,y) { x*y*1000 }
  
  # memoise() vs. new.env()
  mem_fn <- memoise(fn)
  htab <- new.env()
  htab_old <- new.env()
  
  # Initialize hash structures
  # n= 200:  0.4 sec
  # n= 400:  1.8 sec
  # n= 800: 10.4 sec
  # n=1600: 10.4 sec
  n <- 1600
  htab_old <- new.env()
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        htab_old[[paste0(i,";",j)]] <- fn(i,j)
      }
    }
  })
  
  system.time({
    for (i in c(1:n)) {
      for (j in c(1:n)) {
        x <- mem_fn(i,j)
      }
    }
  })
  
  # Run benchmarks
  microbenchmark({x <- mem_fn(5,6)}, times=1000L)
  microbenchmark({x <- htab[[paste0(5,";",6)]]}, times=1000L)
  
}

# TEMP
if (F) {
  
  int <- (3/5)*3^5 + (2/3)*3^3
  print(int)
  
  G <- function(x) { x^3 + 2*x }
  h <- function(x) { x^2 }
  m <- 1000000
  i <- c(1:m)
  int_appx <- sum(
    (G((3*i)/m)-G((3*i-3)/m)) * (h((3*i)/m))
  )
  print(int_appx)
  
}

# TEMP
if (F) {
  
  dat <- c(0,0,0,0,1,1,1,2,2,3)+0.5
  n <- length(dat)
  I <- as.integer
  
  # Version 1
  {
    dens <- Vectorize(function(x, phi1, phi2, phi3) {
      phi4 <- 4-(phi1+phi2+phi3)
      bin <- which.min(x>c(0,1,2,3,4))-1
      lik <- phi1*I(bin==1)+phi2*I(bin==2)+phi3*I(bin==3)+phi4*I(bin==4)
      return(lik)
    }, vectorize.args=c("x"))
    
    wlik <- function(par) {
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        lik <- dens(dat[i], par[1], par[2], par[3])
        return(log(lik))
      }))
      return(-1*sum_loglik)
    }
    
    opt <- optim(par=c(1,1,1), fn=wlik, method="CG") # BFGS CG SANN
    if (opt$convergence!=0) { warning("nonconvergence of optim()") }
    dens(c(0,1,2,3)+0.5, phi1=opt$par[1], phi2=opt$par[2], phi3=opt$par[3])
    # sum(d)
  }
  
  # Version 2 (only dens)
  {
    dens <- Vectorize(function(x, par) {
      bin <- which.min(x>c(0,1,2,3,4))-1
      k <- 4 # num bins
      hz <- sapply(c(1:(ifelse(bin==k,k-1,bin))), function(j) {
        expit(par[j])
      })
      p1 <- ifelse(bin==k, 1, hz[bin])
      p2 <- ifelse(bin==1, 1, prod(1-hz[1:(bin-1)]))
      lik <- k*p1*p2
      return(lik)
    }, vectorize.args=c("x"))

    wlik <- function(par) {
      sum_loglik <- sum(sapply(c(1:n), function(i) {
        lik <- dens(dat[i], par)
        return(log(lik))
      }))
      return(-1*sum_loglik)
    }
    
    opt <- optim(par=c(0.5,0.5,0.5), fn=wlik, method="CG") # BFGS CG SANN
    if (opt$convergence!=0) { warning("nonconvergence of optim()") }
    dens(c(0,1,2,3)+0.5, par=c(opt$par[1],opt$par[2],opt$par[3]))
    sum(d)
    
  }
  
}

# Testing Ted's density estimator function
if (F) {
  
  library(ctsCausal)
  
  # Sample data
  n <- 1000
  WW <- data.frame(W1 = runif(n))
  Z <- rbinom(n, size = 1, prob = 1/(1 + exp(2-WW$W1)))
  AA <- (1-Z) * rnorm(n, mean = WW$W1, sd = abs(1 + WW$W1))
  fit <- cmdSuperLearner(
    A = AA,
    W = WW,
    control = list(
      SL.library = c("SL.mean","SL.glm","SL.gam","SL.earth"),
      verbose = FALSE,
      n.bins = c(2:10)
    )
  )  
  
}


# Testing G-comp variance estimator
if (F) {
  
  # Generate data
  ests <- c()
  for (i in 1:100) {
    n <- 100
    alpha_0 <- -1.5
    alpha_1 <- 0.3
    alpha_2 <- 0.7
    alpha_3 <- -3
    w1 <- runif(n)
    w2 <- rbinom(n, size=1, prob=0.5)
    a <- rbeta(n, shape1=0.9, shape2=(1.1 + 0.4*w2))
    probs <- expit(alpha_0 + alpha_1*w1 + alpha_2*w2 + alpha_3*a)
    y <- rbinom(n, size=1, prob=probs)
    dat <- data.frame(w1=w1, w2=w2, a=a, y=y)
    
    # Fit logistic regression model
    model <- glm(y~w1+w2+a, data=dat, family="binomial")
    coefs <- as.numeric(summary(model)$coefficients[,1])
    
    # Run gcomp
    x <- 0.5
    gcomp_i <- apply(
      X = dat,
      MARGIN = 1,
      FUN = function(r) {
        M_i <- c(1,r[["w1"]],r[["w2"]],x)
        return(expit(as.numeric(M_i %*% coefs)))
      }
    )
    
    ests <- c(ests, mean(gcomp_i))
  }
  sd(ests)
  
  # Run histogram
  (function(x){ggplot(data.frame(x=x),aes(x=x))+geom_histogram()})(ests)
  
  # Calculate variance estimate
  I_inv <- n*vcov(model)
  mu_n <- function(w1,w2,a) {
    expit(coefs[1] + coefs[2]*w1 + coefs[3]*w2 + coefs[4]*a)
  }
  exp_lin <- function(w1,w2,a) {
    exp(coefs[1] + coefs[2]*w1 + coefs[3]*w2 + coefs[4]*a)
  }
  v1 <- function(w1,w2,a) {
    t(as.matrix(c(1,w1,w2,a)))
  }
  v2 <- function(y_i,w1_i,w2_i,a_i,w1_j,w2_j,a) {
    y_mins_mu <- y_i - mu_n(w1_j,w2_j,a)
    as.matrix(c(
      y_mins_mu,
      w1_i*y_mins_mu,
      w2_i*y_mins_mu,
      a_i*y_mins_mu
    ))
  }
  
  w1 <- dat$w1
  w2 <- dat$w2
  y <- dat$y
  a <- dat$a
  var_est <- (1/(n^4)) * sum(sapply(c(1:n), function (i) {
    (sum(sapply(c(1:n), function (j) {
      
      mu_i <- mu_n(w1[i],w2[i],x)
      mu_j <- mu_n(w1[j],w2[j],x)
      mu_i - mu_j + (
        ( mu_j / (1+exp_lin(w1[j],w2[j],x)) ) *
        ( v1(w1[j],w2[j],x) %*%
          I_inv %*%
          v2(y[i],w1[i],w2[i],a[i],w1[j],w2[j],x) )
      )
      
    })))^2
  }))
  print(var_est)
  
  # var_est <- (1/(n^4)) * sum(
  #   
  #   mu_n(w1_i_long,w2_i_long,x) - mu_n(w1_j_long,w2_j_long,x) + (
  #     ( mu_n(w1_j_long,w2_j_long,x) / (1 + exp_lin(w1_j_long,w2_j_long,x)) ) *
  #     ( v1(w1_j_long,w2_j_long,x) %*%
  #       I_inv %*%
  #       v2(y_i_long,w1_i_long,w2_i_long,a_i_long,w1_j_long,w2_j_long,x) )
  #   )
  #   
  # )
  
}

# Testing conditional distribution estimators
if (F) {
  
  w <- 10
  gamma <- c(0,-3,-3)
  hz1 <- expit(gamma[1]+0.2*w)
  hz2 <- expit(gamma[2]+0.2*w)
  hz3 <- expit(gamma[3]+0.2*w)
  hz4 <- expit(gamma[4]+0.2*w)
  p1 <- hz1
  p2 <- hz2 * (1-hz1)
  p3 <- hz3 * (1-hz2)*(1-hz1)
  p4 <- (1-hz3)*(1-hz2)*(1-hz1) # Modification of Diaz and VDL
  ggplot(data.frame(
    x = c(1:4),
    y = c(p1,p2,p3,p4)
  ), aes(x=x, y=y)) + geom_bar(stat="identity") + ylim(c(0,1))
  #
  
  
  n <- 10000
  w1 <- runif(n)
  w2 <- runif(n)
  a <- rbeta(n, shape1=0.3+w1, shape2=0.3+w2)
  ggplot(data.frame(a=a), aes(x=a)) + geom_histogram(bins=50)
  
}

# Experimenting with the profiler
if (F) {
  
  y <- function() {
    rnorm(10^7)
  }
  z <- function() {
    y()
  }
  x1 <- rnorm(10^7)
  x2 <- 999
  x3 <- y()
  x4 <- 999
  x5 <- z()
  x6 <- 999
  
}

# Testing whether memoising works with constructors
if (F) {
  
  construct_test <- function() {
    
    return(memoise(Vectorize(function(a) {
      mean(mu_n(a,dat$w1,dat$w2))
    })))
    
  }
  t1 <- construct_test()
  t2 <- construct_test()
  environment(t1)
  environment(t2)
  
  construct_test2 <- function(test_func) {
    
    return(function(a) {
      print(environment(test_func))
    })
    
  }
  
  tt1 <- construct_test2(t1)
  tt2 <- construct_test2(t2)
  tt1()
  tt2()
  
}

# Adapting ECDF and inverse ECDF to two-phase sampling
if (F) {
  
  dat <- data.frame(
    a = c(0.5,0.75,0.75,NA,NA),
    wts = c(6,1,1,NA,NA)
  )
  
  # Old functions
  Phi_n <- ecdf(dat$a)
  Phi_n_inv <- function(x) {
    dat %<>% filter(!is.na(a))
    qemp(p=x, obs=dat$a, discrete=T)
  }
  
  # New functions
  construct_Phi_n2 <- function(dat, type=params$ecdf_type) {
    dat <- cbind(dat, wts=wts(dat, scale="sum 1"))
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    s <- sum(dat$wts) / n
    Vectorize(function(x) { (1/(n*s))*sum(dat$wts*as.integer(dat$a<=x)) })
  }
  construct_Phi_n3 <- function (dat, which="ecdf", type=params$ecdf_type) {
    # Adaptation of ecdf() source code
    dat <- cbind(dat, wts=wts(dat, scale="sum 1"))
    dat %<>% arrange(a)
    n <- nrow(dat)
    dat %<>% filter(!is.na(a))
    vals_x <- unique(dat$a)
    vals_y <- c()
    s <- sum(dat$wts) / n
    for (j in 1:length(vals_x)) {
      indices <- which(dat$a==vals_x[j])
      weights_j <- dat$wts[indices]
      new_y_val <- sum(weights_j) / (n*s)
      vals_y <- c(vals_y, new_y_val)
    }
    vals_y <- cumsum(vals_y)
    if (type=="ecdf") {
      rval <- approxfun(vals_x, vals_y, method="constant", yleft=0,
                        yright=1, f=0, ties="ordered")
    } else if (type=="inverse") {
      rval <- approxfun(vals_y, vals_x, method="constant", yleft=min(vals_x),
                        yright=max(vals_x), f=1, ties="ordered")
    }
    return(rval)
  }
  Phi_n2 <- construct_Phi_n2(dat, type=params$ecdf_type)
  Phi_n3 <- construct_Phi_n3(dat, type=params$ecdf_type)
  Phi_n_inv3 <- construct_Phi_n3(dat, which="inverse", type=params$ecdf_type)
  
  
  # # Run benchmarks
  # microbenchmark(Phi_n(0.5), times=100)
  # microbenchmark(Phi_n2(0.5), times=100)
  # microbenchmark(Phi_n3(0.5), times=100)
  
  # ECDF
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n2)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n3)
  
  # Inverse ECDF
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n_inv)
  ggplot(data.frame(x=seq(0,1,0.01)), aes(x=x)) +
    stat_function(fun = Phi_n_inv3)
  
  # Correct behavior at knots
  Phi_n(0.5)
  Phi_n2(0.5)
  Phi_n3(0.5)
  Phi_n_inv(1/3)
  Phi_n_inv3(1/3)
  
}
