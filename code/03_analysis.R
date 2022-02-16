library(stringr) # string operations
library(magrittr) # piping
library(ggplot2); theme_set(theme_bw()) # plotting
library(purrr) # mapping
library(survival) # survival modeling
library(data.table) # handling data
library(mice) # D1 function for model comparison and pool
library(mitools) # imputation utils
source(here::here('code', '00_functions.R')) # some helper functions
library(Hmisc) # splines
library(directlabels) # for plots
require(ggh4x) # for the multi-facetted plot
library(dplyr)

datadir <- here::here('data')

# read in and prepare data
dfo <- fread(here::here(datadir, 'clinical_and_measurements_complete.csv'))
dfo[, c_stage_int:=ordered(c_stage_int, levels=1:4)]
dfo[histo=='', histo:=NA]
dfo[, bmi_raw := weight / length^2]
dfo[, total_muscle:=copy(psoas_muscle)]
dfo[, pmd_raw := copy(mean_HU_long_spine_muscle)]
dfo[, pmi_raw := total_muscle / length^2]
npts <- nrow(dfo)

means_and_sds <- fread(here::here('resources', 'means_and_sds_cont.csv'))
means <- as.list(set_names(means_and_sds$mean, means_and_sds$variable))
sds <- as.list(set_names(means_and_sds$sd, means_and_sds$variable))

#### scale continuous variables
dfo[, age_raw:=copy(age)]
age_mean <- dfo[, mean(age_raw)]
age_sd <- dfo[, sd(age_raw)]
dfo[, age:=(age_raw - age_mean) / age_sd]
bmi_mean <- dfo[, mean(bmi_raw, na.rm=T)]
bmi_sd <- dfo[, sd(bmi_raw, na.rm=T)]
pmi_mean <- dfo[, mean(pmi_raw, na.rm=T)]
pmi_sd <- dfo[, sd(pmi_raw, na.rm=T)]
pmd_mean <- dfo[, mean(pmd_raw, na.rm=T)]
pmd_sd <- dfo[, sd(pmd_raw, na.rm=T)]
dfo[, bmi:=(bmi_raw - bmi_mean) / bmi_sd] 
dfo[, pmi:=(pmi_raw - pmi_mean) / pmi_sd] 
dfo[, pmd:=(pmd_raw - pmd_mean) / pmd_sd]
dfo[, pmipmd := pmi * pmd]
pmipmd_sd <- dfo[, sd((pmi_raw - pmi_mean) * (pmd_raw - pmd_mean), na.rm=T)]

# load imputed data
impfname <- 'imp_full_smcfcs_strat2_rcs5_i250_'

## potentially, multiple sets of imputations can exist (and be added to get more imputed datasets)
impmatches <- str_subset(list.files(datadir),
                         paste0('^', impfname, "m(\\d+)_s(\\d+)p.rds"))
impfname1 <- impmatches[1]
# impraws <- map(impmatches, ~readRDS(here::here('data', .x)))
impraws <- map(impmatches, ~readRDS(here::here('tmp', .x)))
## combine imputations
impsets <- do.call(c, map(impraws, 'impDatasets'))

# set ecog to a factor for easier interpretation of coefficients
impsets <- map(impsets, mutate, ecog_bin=factor(ecog_bin, ordered=F))
impsets <- map(impsets, mutate, ecog_bin1=ecog_bin=='1', ecog_bin2=ecog_bin=='2')

#' for some reason, adding ecog_bin:c_stage early adds three terms instead of 2,
#' also for the reference category of ecog_bin = 0
#' this is non-identifyable as it gets conflicted with the baseline hazard
#' add a dummy var
impsets <- map(impsets, mutate,
               ecog_bin_c_stage_early1 = as.integer((ecog_bin == '1') * c_stage_early),
               ecog_bin_c_stage_early2 = as.integer((ecog_bin == '2') * c_stage_early),
               histono_pa_c_stage_earlyTRUE = as.integer((histo=='no_pa') * c_stage_early),
               histoother_c_stage_earlyTRUE = as.integer((histo=='other') * c_stage_early),
               histosquamous_c_stage_earlyTRUE = as.integer((histo=='squamous') * c_stage_early)
               )

imp <- imputationList(impsets)
imp20 <- imputationList(impsets20)
nimp <- length(imp$imputations)
print(paste0("using ", nimp, " imputed datasets"))

# make formulas
lin_base <- 'age + sex_male + histo + ecog_bin + bmi + pmi + pmd + pmipmd'
rcs_base <- paste(lin_base, "+ age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 + pmipmd1 + pmipmd2 + pmi1pmd + pmi2pmd + age3 + bmi3 + pmi3 + pmd3 + pmipmd3 + pmi3pmd")

f_rcs_strat <- paste0(
  "Surv(time,deceased) ~ ",
  rcs_base,
  " + (", rcs_base, "):c_stage_early",
  " + strata(c_stage_int)"
)                  
f_rcs_pool <- paste0(
  "Surv(time,deceased) ~ ",
  rcs_base,
  " + strata(c_stage_int)"
)

f_lin_strat <- paste0(
  "Surv(time,deceased) ~ ",
  lin_base,
  " + (", lin_base, "):c_stage_early",
  " + strata(c_stage_int)"
)                  
f_lin_pool <- paste0(
  "Surv(time,deceased) ~ ",
  lin_base,
  " + strata(c_stage_int)"
)

## baseline formulas
f_lin_pool1 <- "Surv(time,deceased)~age + sex_male + histo + ecog_bin + bmi + pmi + pmd + strata(c_stage_int)"
f_lin_pool0 <- "Surv(time,deceased)~age + sex_male + histo + ecog_bin + bmi + strata(c_stage_int)"

# main hypothesis
fs_2 <- with(imp, coxph(as.formula(f_rcs_strat))) # full version with stratification and splines
fspmipmd_2l <- with(imp, coxph( # linear version without stratification of pmipmd
  Surv(time,deceased) ~ age + sex_male + histo + ecog_bin + bmi + pmi + pmd +
    age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 +
    age3 + bmi3 + pmi3 + pmd3 + pmipmd +
    (age + sex_male + histo + ecog_bin + bmi + pmi + pmd + 
    age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 + 
    age3 + bmi3 + pmi3 + pmd3):c_stage_early + strata(c_stage_int)
))

## null-models for comparison
fs_1 <- with(imp, coxph(
  Surv(time,deceased) ~ age + sex_male + histo + ecog_bin + bmi + pmi + pmd +
    age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 +
    age3 + bmi3 + pmi3 + pmd3 + 
    (age + sex_male + histo + ecog_bin + bmi + pmi + pmd + 
    age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 + 
    age3 + bmi3 + pmi3 + pmd3):c_stage_early + strata(c_stage_int)
))

test_interaction_nonlinear <- D1(fs_2, fs_1) # non-linear interaction and stratification
test_interaction_linear <- D1(fspmipmd_2l, fs_1) # linear interaction
test_interaction_nonlinear
test_interaction_linear

prms_spmipmd_2l <- summary(pool(fspmipmd_2l))
fwrite(prms_spmipmd_2l, here::here('results', 'allprms.csv'), row.names=F)

## make plot of effect of pmi per pmd
# generate data where the continuous variable with spline terms varies between min and max
testdata <- make_testdat(impfname1)
testdatpmipmd <- testdata[term=='pmipmd']
  
## plot effect of pmi for values of pmd
### to get the confidence intervals we will use bootstrap
#' method 2, "MI Boot" (see Bootstrap inference with multiple imputation, Schomaker and Heumann 2018)
#' 
#' get M imputed datasets
#' make test data (range of pmd values)
#' for each imputed dataset i
#'   fit survival model
#'   predict on test data
#'   store mean prediction
#'   
#'   get B bootstrap samples of imputed dataset i
#'   for each bootstrap sample j
#'     fit survival model
#'     predict on test data
#'     
#'   calculate s.e. of predictions using bootstrap samples (usual sd(...) function)
#'   store s.e.
#' use Rubin's rules on means and s.e.s for final prediction

nboots = 100

get_boot_preds <- function(dfimp, testdata) {
  # get nboots bootstrap samples of imputed dataset
  bootiis <- map(1:nboots, ~sample(npts, replace=T))
  bootdfs <- map(bootiis, function(iis) dfimp[iis,])
  bootdf <- rbindlist(bootdfs, idcol='bootidx')
  
  # fit and predict on 'full' imputed dataset
  fitfull <- coxph(Surv(time,deceased)~
       age+age1+age2+age3+sex_male+histo+ecog_bin+bmi+bmi1+bmi2+bmi3+pmi+pmd+pmd1+pmd2+pmd3+pmipmd+pmipmd1+pmipmd2+pmipmd3+
       strata(c_stage_int), data=dfimp)
  predfull <- predict(fitfull, newdata=testdata, se.fit=F)
  predfulldf <- data.table(pmd_raw=testdata$pmd_raw, estimate=predfull, pmi=testdata$pmi)
  predfulldfw <- dcast(predfulldf, pmd_raw ~ pmi, value.var='estimate')
  predfulldfw[, estimate:=`1` - `0`]
  
  # fit and predict for each imputed dataset
  bootfits <- bootdf[, list(
    fit = list(coxph(Surv(time,deceased)~
       age+age1+age2+age3+sex_male+histo+ecog_bin+bmi+bmi1+bmi2+bmi3+pmi+pmd+pmd1+pmd2+pmd3+pmipmd+pmipmd1+pmipmd2+pmipmd3+
       strata(c_stage_int), data=.SD))),
    by='bootidx']
  bootpreds <- map(bootfits$fit, predict, newdata=testdata, se.fit=F)
  bootpreddf <- rbindlist(map(bootpreds, ~data.table(pred=.x, pmd_raw=testdata$pmd_raw, pmi=testdata$pmi)), idcol='bootidx')
  
  # calculate difference between pmi = 1 and pmi = 0
  bootpreddfw <- dcast(bootpreddf, ...~pmi, value.var = 'pred')
  bootpreddfw[, lhr:=`1` - `0`]
  bootpredses <- bootpreddfw[, list(std.error=sd(lhr)), by=c('pmd_raw')]
  
  bootpredsout <- bootpredses[predfulldfw[, list(pmd_raw, estimate)], on='pmd_raw']
}
st <- system.time({
  bp1 <- get_boot_preds(imp$imputations[[1]], testdatpmipmd)
})
stopifnot(nrow(bp1) == nrow(testdatpmipmd[pmi==0]))

# this may take ~5 minutes per 100 imputations
print(paste0('running bootstrap fits to get std.error with ', nboots, ' bootstrap samples for ', nimp, ' imputed datasets.'))
print(paste0('will take approx ', round(st[3] * nimp / 60, 1), ' minutes'))
bootpreds <- map(imp$imputations, get_boot_preds, testdata=testdatpmipmd)
bootpreddf <- rbindlist(bootpreds, idcol='impidx')
## save this intermediate result to be able to pick up later if needed
fwrite(bootpreddf, here::here('resources', 'pmdpmibootpreddf.csv'), row.names=F)
if (FALSE) {
  bootpreddf <- fread(here::here('resources', 'pmdpmibootpreddf.csv'))
}

# summarize and plot
bootpreddf[, u:=std.error^2]
bootpredsum <- bootpreddf[, list(pooled=list(pool.scalar(estimate, u))), by='pmd_raw']
bootpredsum[, `:=`(
  estimate=map_dbl(pooled, 'qbar'),
  std.error=map(pooled, 't') %>% map_dbl(~sqrt(.x)) 
)]
bootpredsum[, `:=`(estimate = estimate/sds$pmi, std.error=std.error/sds$pmi)]
bootpredsum[, `:=`(
  ci_lo=estimate - 1.96 * std.error,
  ci_hi=estimate + 1.96 * std.error
)]

# define breaks on y-axis
ylabels_hr <- c(0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25)
ybreaks <- log(ylabels_hr)

# check where 90% of observations lie
pmd_range90 <- dfo[, quantile(pmd_raw, probs=c(0.05, 0.95), na.rm=T)]
pmd_range95 <- dfo[, quantile(pmd_raw, probs=c(0.025, 0.975), na.rm=T)]
pmd_range98 <- dfo[, quantile(pmd_raw, probs=c(0.01, 0.99), na.rm=T)]
pmd_raw_min <- pmd_range98[1]
pmd_raw_max <- pmd_range98[2]

# combine data 
fig1data <- rbindlist(list(
  # histogram = dfo[pmd_raw>=pmd_raw_min & pmd_raw <= pmd_raw_max, list(pmd=pmd_raw)],
  histogram = dfo[, list(pmd=pmd_raw)],
  # pmdpmi = bootpredsum[pmd_raw>=pmd_raw_min & pmd_raw <= pmd_raw_max]
  pmdpmi = bootpredsum
), fill=T, idcol='facet')
fig1data[, facet:=factor(facet,
                         levels=c('pmdpmi', 'histogram'),
                         labels=c('Association between PMI and overall survival per PMD',
                                  'histogram')
                         )]  

p1 <- fig1data %>%
  ggplot(aes(x=pmd_raw)) +
  geom_line(aes(y=estimate)) +
  geom_ribbon(aes(ymin=ci_lo, ymax=ci_hi), linetype=2, alpha=0.1) + 
  geom_hline(aes(yintercept=0), linetype=2, alpha=0.5) +
  geom_vline(aes(xintercept=pmd_raw_min), alpha=0.5) + 
  geom_vline(aes(xintercept=pmd_raw_max), alpha=0.5) + 
  labs(x='PMD (Hounsfield Units)', y=bquote(atop(
         "hazard ratio for 1 standard deviation increase in PMI",
         "Improved overall survival" %<->% "Worse overall survival")))+
  scale_y_continuous(breaks=ybreaks, labels=ylabels_hr) +
  geom_histogram(aes(x=pmd)) + 
  facet_grid(facet~., scales='free') + 
  force_panelsizes(row=c(8,2)) +
  facetted_pos_scales(
    y=list(
      # scale_y_continuous(breaks=ybreaks, labels=ylabels_hr),
      scale_y_continuous(),
      scale_y_continuous()
    ))
p1
ggsave(here::here('results', 'hr_for_pmi_per_pmd.png'), p1,
       width = 10, height=10)  

ggsave(here::here('results', 'hr_for_pmi_per_pmd.pdf'), p1,
       width = 10, height=10, dpi=1200)  

