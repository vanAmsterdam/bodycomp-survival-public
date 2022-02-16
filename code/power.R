library(ggplot2); theme_set(theme_bw()) # plotting
library(purrr) # mapping and walking 
library(survival) # survival modeling
library(data.table) # data handling
source(here::here('code', "pgwsource.R")) # source for power-generalized weibull
source(here::here('code', '00_functions.R')) # helper functions
library(copula)

# loghrs on original variable scale
# from sjoblom 2016; 
loghrs <- c(
  age = log(0.99), # -0.01
  sex_male = log(0.77), # -0.26
  histoother = 0.2,    # these ones are not given so these are assumed
  histosquamous = 0.3, # these ones are not given so these are assumed
  bmi = -0.01,
  ecog_bin1 = log(1.24), # 0.215
  ecog_bin2 = log(1.89), # 0.637
  smd    = log(0.84) / 10, # take univariable hr of 10 SMDs (take 10 for higher precision)  
  smi    = log(0.99), 
  smismd = 0.5 * (log(0.84) / 10 + log(0.99)) # take midpoint of smd and smi; will be scaled by sd of smismd inside of data simulation
)
fwrite(data.table(term=names(loghrs), hazard_ratio=exp(loghrs), log_hazard_ratio=loghrs),
       here::here('resources', 'samplesize_hrs.csv'),
       row.names=F)

# stats for smi (calculate from stats for age <75 and age>=75)
p1 <- (1 - 0.18) 
p2 <- 0.18 # age >=75
mu1 <- 48.1; mu2 <- 45.6
var1 <- 7.6^2; var2 <- 7.4^2
mu_smi <- mu1 * p1 + mu2 * p2
sd_smi <- sqrt(p1*var1 + p2*var2 + p1*p2*(mu1 - mu2)^2)
# stats for smd
mu1 <- 38.2; mu2 <- 33.0
var1 <- 8.3^2; var2 <- 8.1^2
mu_smd <- mu1 * p1 + mu2 * p2
sd_smd <- sqrt(p1*var1 + p2*var2 + p1*p2*(mu1 - mu2)^2)

marginstats <- list(
  age=list(mean=65.4, sd=9.4),
  bmi=list(mean=23.4, sd=3.9),
  smi=list(mean= mu_smi, sd = sd_smi),
  smd=list(mean= mu_smd, sd = sd_smd) 
)


# get freqs of stages based on dolan 2020 (rt only, stage I, II and III);
# and abbass 2020 (rt only, stages III and IV)
p1_123_d = .353
p2_123_d = .185
p3_123_d = .345

p3_34_a = .373
p4_34_a = .627

# make system of equations and solve for unknowns
Amat <- matrix(c(
  1 - p1_123_d, - p1_123_d, - p1_123_d, 0,
  p2_123_d, 1 - p2_123_d, - p2_123_d, 0,
  0, 0, 1 - p3_34_a, - p3_34_a, 
  1, 1, 1, 1
), nrow=4, byrow = T)
b <- c(0, 0, 0, 1)
p_c_stage_int <- solve(Amat, b)

catfreqs <- list(
  histo=c(adeno= 1 - (0.063 + .176 + 0.233), squamous= 0.233, other= 0.063+.176), # adeno, squamous, other
  ecog_bin=c(`0` = 0.225, `1` = 0.55, `2` = 0.225), # 0, 1, >=2
  c_stage_int=purrr::set_names(p_c_stage_int, 1:4)
)

p_sex_male=c(0.572)

# initiate copula
dimtot <- 8L # age, bmi, ecog_bin, histo, sex_male, c_stage_int, smi, bmi
margindists <- c(
  age='norm', bmi='norm', ecog_bin='norm',
  histo='norm', sex_male='binom',
  c_stage_int='norm',
  smi='norm', smd='norm'
)

# copula will be defined on scaled variables
marginprms <- list(
  age=list(mean=0, sd=1),
  bmi=list(mean=0, sd=1),
  ecog_bin=list(mean=0, sd=1),
  histo=list(mean=0, sd=1),
  sex_male=list(size=1, p=p_sex_male),
  c_stage_int=list(mean=0, sd=1),
  smi=list(mean=0, sd=1),
  smd=list(mean=0, sd=1)
)

# define function to turn matrix of copula samples into data frame,
#' handling the discretation of histology and ecog
copula_smps_to_df <- function(copsmps, col.names, catfreqs, contstats) {
  colnames(copsmps) <- col.names
  copdf <- as.data.table(copsmps)
  
  # scale continuous vars to raw scale
  contvars <- names(contstats)
  walk2(contvars, contstats, function(varname, varstats) {
    setnames(copdf, varname, paste0(varname, '_raw'))
    set(copdf, j=varname, value=copdf[[paste0(varname, '_raw')]] * varstats$sd + varstats$mean)
  })
  
  # determine categorical vars
  # use marginal probabilities to get cut-offs from normal distribution
  catvars <- names(catfreqs)
  walk2(catvars, catfreqs, function(varname, varfreqs) {
    stopifnot(sum(varfreqs) == 1)
    varcutprobs <- cumsum(varfreqs)
    varcutpoints <- c(-Inf, qnorm(varcutprobs))
    varcatnames <- names(varfreqs)
    setnames(copdf, varname, paste0(varname, '_raw'))
    set(copdf, j=varname, value=cut(copdf[[paste0(varname, '_raw')]],
                                    breaks=varcutpoints,
                                    labels=varcatnames,
                                    include.lowest=T))
  })
  return(copdf)
}

# create sampling function for covariate data
sample_X <- function(n=100, copula_prm=1) {
  # create copula with right parameter for covariance
  copulai <- claytonCopula(copula_prm, 8L)
  # create multivariate distribution to sample from
  mvdi <- mvdc(copula = copulai, 
               margins=margindists,
               paramMargins = marginprms)
  # sample raw values
  copsmps <- rMvdc(n, mvdi)
  # turn in to pretty data.table
  copdf <- copula_smps_to_df(copsmps, names(marginprms), catfreqs, marginstats)
  return(copdf)
}

# function for adding extra columns to xmat in place
augment_X <- function(xdf) {
  # assign interaction term between smi and smd
  set(xdf, j='smismd', value=(xdf[['smi']] - mean(xdf[['smi']])) * (xdf[['smd']] - mean(xdf[['smd']])))
  
  # create dummy variables of factor variables
  colclasses <- map_chr(xdf, class)
  facvars <- names(colclasses)[colclasses == 'factor']
  walk(facvars, function(varname) {
    faclevels <- levels(xdf[[varname]])
    walk(tail(faclevels, -1), function(levelname) {
      set(xdf, j=paste0(varname, levelname), value=as.integer(xdf[[varname]] == levelname))
    })
  })
}

# make function to add beta column based on log hazard ratios
calculate_lhr <- function(data, loghrs) {
  beta_ = numeric(length=nrow(data))
  Xmat = data[, as.matrix(.SD), .SDcols=names(loghrs)]
  return(Xmat %*% loghrs)
}

#' get marginal statistics of survival and censoring from our own data
#' including intercepts
#' if the data is not available, read the marginal parameters from a file
if (file.exists(here::here('data', 'clinical_and_measurements_complete.csv'))) {
  dfo <- fread(here::here('data', 'clinical_and_measurements_complete.csv'))
  maxtime <- max(dfo$time)
  margfits <- dfo[!is.na(c_stage_int), list(fit=list(pgwfit(Surv(time,deceased)~1))), by=c_stage_int]
  margstats <- copy(margfits)
  margstats[, `:=`(
    c_stage_int = factor(c_stage_int, levels=1:4),
    alpha0 = map(fit, 'mle') %>% map_dbl(function(mat) mat['a_(Intercept)', 'mle']),
    beta0  = map(fit, 'mle') %>% map_dbl(function(mat) mat['b_(Intercept)', 'mle']),
    nu0    = map(fit, 'mle') %>% map_dbl(function(mat) mat['n_(Intercept)', 'mle']),
    fit = NULL
  )]
  setkey(margstats, c_stage_int)
  fwrite(margstats, here::here('resources', 'margstats.csv'), row.names=F)
  
  # calculate marginal censoring distribution
  censfit <- pgwfit(Surv(time, 1-deceased)~1, data=dfo)
  censstats <- c(
    'alpha0' = censfit$mle['a_(Intercept)', 'mle'],
    'beta0'  = censfit$mle['b_(Intercept)', 'mle'],
    'nu0'    = censfit$mle['n_(Intercept)', 'mle']
  )
  
  dump(c('censstats'), here::here('resources', 'censstats.txt'))
  
  tmin <- min(dfo$time)
  tmax <- max(dfo$time)
  tvec <- seq(tmin, tmax, len=100)
  
  #' also create kaplan-meier estimates to compare the parametric fit with
  #' the non-parametric estimator
  kmfits <- dfo[!is.na(c_stage_int), list(fit=list(survfit(Surv(time,deceased)~1))), by=c_stage_int]
  margfits[kmfits, kmfit:=i.fit, on='c_stage_int']
  png(here::here('results', 'marginal_survivals.png'), width=1200, height=1200, pointsize=24)
  par(mfrow=c(2,2))
  margfits[order(c_stage_int), walk2(fit, kmfit, function(pgwfit, kmfit) {
    plot(kmfit, main=paste0('stage: ', c_stage_int), xlim=c(0, tmax), xlab='years', ylab='cumulative survival probability')
    lines(tvec, pgwS(rbind(fitvec(pgwfit, grp=1)), tvec), col=2)
  }), by=c_stage_int]
  dev.off()
} else {
  margstats <- fread(here::here('resources', 'margstats.csv'))
  setkey(margstats, c_stage_int)
  source(here::here('resources', 'censstats.txt'))
}

simulate_data <- function(nsim=100, copula_prm=1.0, loghrs) {
  # sample Xs using copula
  simdata <- sample_X(nsim, copula_prm)
  # augment deterministic variables
  augment_X(simdata)
  
  # add in intercepts of survival and censoring
  simdata <- simdata[margstats, on='c_stage_int']
  walk2(names(censstats), censstats, ~set(simdata, j=paste0(.x, 'c'), value=.y))
  
  # update log hazard ratios based on standard deviation of smismd
  smd_sd <- sd(simdata$smd)
  smismd_sd <- sd(simdata$smismd)
  loghrs[['smismd']] <- loghrs[['smismd']] * smd_sd / smismd_sd
  
  # calculate hazard ratios
  simdata[, lhr:=calculate_lhr(.SD, loghrs)]
  simdata[, beta_:=beta0 + lhr - mean(lhr)] # keep the center right
  
  # now sample noise terms
  simdata[, u:=runif(.N)]
  simdata[, time:=pgwsim(cbind(0, beta_, alpha0, nu0), u)]
  
  # sample censoring times
  simdata[, uc:=runif(.N)]
  simdata[, timec:=pgwsim(cbind(0, beta0c, alpha0c, nu0c), uc)]
  
  simdata[, deceased:=time < timec]
  
  # fix extremes
  simdata[time>=maxtime | !is.finite(time), `:=`(time=maxtime, deceased=0)]
  
  return(simdata)
}

set.seed(1)

# generate datasets
## check reasonable values for the copula parameter
copulaprms <- seq(0.1, 1.0, by=0.1)
copulas <- map(copulaprms, claytonCopula)
mvdcs <- map(copulas, mvdc, 
             margins=c('norm', 'norm'),
             paramMargins=list(list(mean=0, sd=1), list(mean=0, sd=1)))
copsmps <- map(mvdcs, rMvdc, n=1000)
corrs <- map_dbl(copsmps, ~cor.test(.x[, 1], .x[, 2])$estimate)
plot(copulaprms, corrs)
fwrite(data.frame(copula_prm=copulaprms, pearson_correlation=corrs), 
       here::here('resources', 'copula_prm_vs_pearson_cor.csv'), row.names=F)
if (FALSE) {cvsp <- fread(here::here('resources', 'copula_prm_vs_pearson_cor.csv'))}

# lay out simulation settings
nreps <- 1000
set.seed(123454321)
simsettings <- CJ(
  simsize = c(500, 1000, 2000),
  copula_prm = c(0.1, 0.25, 1.0),
  dk = 1L
)
simsettings[, simsettingidx:=.I]
simgrid <- simsettings[data.table(simrep=1:nreps, dk=1L), on='dk', allow.cartesian=T]
simgrid[, simkey:=paste0('simsettingidx', simsettingidx, 'rep', simrep)]
setkey(simgrid, 'simkey')

# get simulated data
simdats <- simgrid[, list(simdat=list(simulate_data(
  copula_prm = copula_prm, n=simsize, loghrs = loghrs))), by='simkey']
simdf <- simdats[, rbindlist(set_names(simdat, simkey), idcol='simkey')]
nevents <- simdf[, list(n_events=sum(deceased)), by='simkey']
dim(simdf) / 1e6

# fit outcome models
outcomeformula <- Surv(time,deceased)~age+sex_male+histoother+histosquamous+bmi+ecog_bin1+ecog_bin2+smd+smi+smismd+strata(c_stage_int)
simfits <- simdf[, list(fit=list(coxph(outcomeformula, data=.SD))), by='simkey']
simfitcoefs <- simfits[, rbindlist(set_names(map(fit, broom::tidy), simkey), idcol='simkey')]
setkey(simfitcoefs)
simfitcoefs[, n_terms:=.N, by='simkey']
simfitcoefs <- simfitcoefs[simgrid]
simfitcoefs <- simfitcoefs[nevents]
## make one-sided significance test
simfitcoefs[, p.value.leq0:=pt(statistic, df=n_events-n_terms, lower.tail=T)]

## calculate powers
powers <- simfitcoefs[, 
   list(power_two_sided = mean(p.value < 0.05),
        power_leq0 = mean(p.value.leq0 < 0.05 & estimate < 0)),
   by=c('simsettingidx', 'term')]
powers <- powers[simsettings, on='simsettingidx']
powers[, copula_prm:=factor(copula_prm)]
powers[term=='smismd']

## plot powers
p_smismd <- powers %>% 
  .[term=='smismd'] %>%
  ggplot(aes(x=simsize, y=power_leq0)) + 
  geom_line(aes(linetype=copula_prm)) + 
  geom_hline(aes(yintercept=0.8), linetype=2)

p <- powers %>% 
  ggplot(aes(x=simsize, y=power_two_sided, col=term)) + 
  geom_line(aes(linetype=copula_prm)) + 
  geom_hline(aes(yintercept=0.8), linetype=2) + 
  facet_wrap(~term)

ggsave(here::here('results', 'power_onesided.png'), p_smismd, width = 15, height = 15)
ggsave(here::here('results', 'power_twosided.png'), p, width=15, height=15)

## save things that take long to compute
fwrite(simdf, here::here('resources', 'powersimdf.csv'), row.names=F)
saveRDS(simfits, here::here('resources', 'powerfits.rds'))
fwrite(simfitcoefs, here::here('resources', 'powercoefs.csv'), row.names=F)
fwrite(powers, here::here('results', 'powerresults.csv'), row.names=F)

