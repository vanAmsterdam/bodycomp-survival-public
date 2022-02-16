library(stringr) # string operations
library(magrittr) # piping
library(purrr) # mapping
library(survival) # survival modeling
library(data.table) # data handling
library(mice) # setting up imputation methods
library(smcfcs) # actual imputation package
library(Hmisc) # splines
library(furrr) # parallelization
library(tictoc) # timing

## set flags
N_WORKERS = 15 # number of cores to use in parallel
N_IMPS = 180 # number of imputations (should be a multitude of N_WORKERS)
NUMIT = 250 # number of iterations per imputation
IMP_SEED = 3 # seed to use for random imputations

## for main analysis
N_KNOTS = 5 # number of knots for restricted cubic splines
CCSMI = F # use complete cases wrt to body composition
STRAT2 = T # stratify by early vs late stage

#' for sensitivity analysis:
#' run once with:
# N_KNOTS = 0; CCSMI = F; STRAT2 = F

# read data 
dfo <- fread(here::here('data', 'clinical_and_measurements_complete.csv'))
dfo[histo=='', histo:=NA]

# select the right columns
colnamedf <- fread(here::here('resources', 'L3_colnames.csv'))
all_measurement_columns <- colnamedf[type=='measurement', unique(newname)]
measurement_columns <- all_measurement_columns %>%
  str_subset('ratio', negate=T) %>%
  str_subset('muscle')

allvars <- c(
  'episode',
  'time', 'deceased',
  'age','sex_male',
  'length', 'weight',
  'ecog_bin', 'histo',
  'c_stage_int',
  measurement_columns
)
facvars <- c('c_stage_int', 'ecog_bin', 'histo')
contvars <- c('age', 'length', 'weight', measurement_columns)

df <- dfo[, .SD, .SDcols=allvars]

walk(facvars, ~set(df, j=.x, value=factor(df[[.x]])))

df[, ecog_bin:=as.ordered(ecog_bin)]
df[, c_stage_int:=as.ordered(c_stage_int)]
df[, c_stage_early:=as.integer(c_stage_int < '3')]

### add in deterministic variables
df[, total_muscle:=copy(psoas_muscle)]
df[, bmi_raw := weight / length^2]
df[, pmi_raw := total_muscle / length^2]
df[, pmd_raw := mean_HU_psoas_muscle]
df[, pmipmd_raw := (pmi_raw - mean(pmi_raw, na.rm=T))*(pmd_raw - mean(pmd_raw, na.rm=T))]

## write out means and sds for continuous variables
contxs <- c('age', 'bmi', 'pmi', 'pmd')
contvars <- c('age', 'bmi_raw', 'pmd_raw', 'pmi_raw', 'pmipmd_raw', 'weight', 'length', 'total_muscle')
names(contxs) <- contxs
dfm <- melt(df, id.col='episode', measure.vars=contvars)
meanssds <- dfm[, list(mean=mean(value, na.rm=T), sd=sd(value, na.rm=T)), by='variable']
meanssds[, variable:=str_replace(variable, '_raw$', '')]
fwrite(meanssds, here::here('resources', 'means_and_sds_cont.csv'), row.names=F)
df[, pmipmd_raw:=NULL]

#### scale continuous variables
df[, age_raw:=copy(age)]
age_mean <- df[, mean(age)]
age_sd <- df[, sd(age)]
df[, age:=(age_raw - age_mean) / age_sd]
bmi_mean <- df[, mean(bmi_raw, na.rm=T)]
bmi_sd <- df[, sd(bmi_raw, na.rm=T)]
pmi_mean <- df[, mean(pmi_raw, na.rm=T)]
pmi_sd <- df[, sd(pmi_raw, na.rm=T)]
pmd_mean <- df[, mean(pmd_raw, na.rm=T)]
pmd_sd <- df[, sd(pmd_raw, na.rm=T)]
df[, bmi:=(bmi_raw - bmi_mean) / bmi_sd] 
df[, pmi:=(pmi_raw - pmi_mean) / pmi_sd] 
df[, pmd:=(pmd_raw - pmd_mean) / pmd_sd]
df[, pmipmd := pmi * pmd]
pmipmd_sd <- df[, sd((pmi_raw - pmi_mean) * (pmd_raw - pmd_mean), na.rm=T)]

if (N_KNOTS > 0) {
  ### define spline basis
  ### find placement of knots (taken from Harrell's RMS book, chapter 2.4.)
  spline_probs <- c(0.05, 0.275, 0.5, 0.725, 0.95)
  ### define knot positions for all continuous variables
  rcsknots <- map(contxs, function(x) quantile(df[[x]], probs=spline_probs, na.rm=T))
    
  ## create rcs variables
  walk2(names(rcsknots), rcsknots, function(varname, knotspos) {
    for (i in 1:(length(knotspos) - 2)) { # omit the last knot as it will be forced to linearity
      set(df, j=paste0(varname, i), value=rcspline.eval(df[[varname]], knots = knotspos)[,i])
    }
  })
  knotnames <- do.call(c, map(contxs, paste0, 1:(N_KNOTS - 2)))
    
  ## create product terms between pmi and pmd
  df[, pmipmd1:=pmi*pmd1]
  df[, pmipmd2:=pmi*pmd2]
  df[, pmi1pmd:=pmi1*pmd]
  df[, pmi2pmd:=pmi2*pmd]
  df[, pmipmd3:=pmi*pmd3]
  df[, pmi3pmd:=pmi3*pmd]
    
  knotnames <- c(knotnames,
                 'pmipmd1', 'pmipmd2', 'pmipmd3',
                 'pmi1pmd', 'pmi2pmd', 'pmi3pmd')
}  
### optionally remove rows with missing values for pmd
if (CCSMI) df <- df[!is.na(pmi)]

### remove raw measurement vars
### we are using total muscle and pmd
keepcols <- setdiff(colnames(df), measurement_columns)
  
dfi <- as.data.frame(df[, .SD, .SDcols=keepcols])

### first initialize to make predictor matrix
ini <- mice(dfi, maxit=0, print=F, defaultMethod = c('norm', 'logreg', 'mlogit', 'podds'))
meth <- ini$method

### assign the non-default imputation methods
if (!CCSMI) meth['pmd'] = "(pmd_raw - pmd_mean) / pmd_sd"
meth['bmi_raw'] = "weight / length^2"
meth['bmi'] = "(bmi_raw - bmi_mean) / bmi_sd"
if (!CCSMI) {
  meth['pmi_raw'] = "total_muscle / length^2"
  meth['pmi'] = "(pmi_raw - pmi_mean) / pmi_sd"
  meth['pmipmd'] = "pmi * pmd"
}
meth['c_stage_early'] = "c_stage_int < '3'"

if (N_KNOTS > 0) {
  ## add splines
  meth['bmi1'] = "rcspline.eval(bmi, knots=rcsknots$bmi)[, 1]"
  meth['bmi2'] = "rcspline.eval(bmi, knots=rcsknots$bmi)[, 2]"
  meth['bmi3'] = "rcspline.eval(bmi, knots=rcsknots$bmi)[, 3]"
  if (!CCSMI) {
    meth['pmi1'] = "rcspline.eval(pmi, knots=rcsknots$pmi)[, 1]"
    meth['pmi2'] = "rcspline.eval(pmi, knots=rcsknots$pmi)[, 2]"
    meth['pmi3'] = "rcspline.eval(pmi, knots=rcsknots$pmi)[, 3]"
    meth['pmd1'] = "rcspline.eval(pmd, knots=rcsknots$pmd)[, 1]"
    meth['pmd2'] = "rcspline.eval(pmd, knots=rcsknots$pmd)[, 2]"
    meth['pmd3'] = "rcspline.eval(pmd, knots=rcsknots$pmd)[, 3]"
    meth['pmipmd1'] = "pmi*pmd1"
    meth['pmipmd2'] = "pmi*pmd2"
    meth['pmipmd3'] = "pmi*pmd3"
    meth['pmi1pmd'] = "pmi1*pmd"
    meth['pmi2pmd'] = "pmi2*pmd"
    meth['pmi3pmd'] = "pmi3*pmd"
  }
}
  
# define predictor matrix
# (columns predict the rows)
predmat <- ini$predictorMatrix
predmat[, 'episode'] <- 0 # dont predict with episode
predmat[, c('pmi', 'pmipmd')] <- 0 # dont predict with deterministic variables
if (N_KNOTS > 0) {
  predmat[, knotnames] <- 0
}
predmat[, c('pmi_raw', 'pmd')] <- 0
predmat[, c('bmi', 'bmi_raw')] <- 0
predmat[, c('c_stage_early')] <- 0
predmat['total_muscle', 'pmd_raw'] <- 0 # if pmd_raw is known, total muscle is always known
# however, in scans with IV contrast, only total_muscle is known and not pmd_raw (not reliable when there is contrast in the muscle)

# define the outcome formula used during imputation
lin_base <- 'age + sex_male + histo + ecog_bin + bmi + pmi + pmd + pmipmd'
rcs_base <- paste(lin_base, "+ age1 + age2 + bmi1 + bmi2 + pmi1 + pmi2 + pmd1 + pmd2 + pmipmd1 + pmipmd2 + pmi1pmd + pmi2pmd + age3 + bmi3 + pmi3 + pmd3 + pmipmd3 + pmi3pmd")

formstring_rcs_strat <- paste0(
  "Surv(time,deceased) ~ ",
  rcs_base,
  " + (", rcs_base, "):c_stage_early",
  " + strata(c_stage_int)"
)

formstring_lin_pool <- paste0(
  "Surv(time,deceased) ~ ",
  lin_base,
  " + strata(c_stage_int)"
)

if (N_KNOTS > 0) {
  formstring <- formstring_rcs_strat
} else {
  formstring <- formstring_lin_pool
}

print("imputing clinical data")

# make name for saving
outname <- 'imp'
if (CCSMI) {
  outname <- paste0(outname, '_ccpmi')
} else {
  outname <- paste0(outname, '_full')
}
outname <- paste0(outname, '_smcfcs')
if (STRAT2) {
  outname <- paste0(outname, '_strat2')
} else {
  outname <- paste0(outname, '_pool')
}
if (N_KNOTS > 0) {
  outname <- paste0(outname, '_rcs', N_KNOTS)
} else {
  outname <- paste0(outname, '_lin')
}
outname <- paste0(outname,
                  '_i', NUMIT,
		  '_m', N_WORKERS
		  )
		  #)
                  #'_s', IMP_SEED,
                  #'p.rds')
print(paste0("outname: ", outname))

print(formstring)

## test if saving works to not waste long computation
saveRDS(c('teststring'), here::here('tmp', paste0('test-', outname)))
saveRDS(c('teststring'), here::here('data', paste0('test-', outname)))

## save knots
if (N_KNOTS > 0) saveRDS(rcsknots, here::here('resources', 'rcsknots', outname))

## start multi-core processing
plan(multicore, workers=N_WORKERS)
tic()
n_per_worker = floor(N_IMPS / N_WORKERS)
seedi = IMP_SEED
for (i in 1:n_per_worker) {
  seedi = seedi + 1
  outnamei = paste0(outname, '_s', seedi, 'p.rds')
  imps <- future_map(1:N_WORKERS, function(j) {
    smcfcs(as.data.frame(dfi),
  	 smtype='coxph',
  	 smformula = formstring,
  	 predictorMatrix = predmat, method=meth, 
  	 numit=NUMIT,
  	 rjlimit=5000,
  	 m=1)},
    .options=furrr_options(seed=seedi)
  )
  toc()
  
  ## combine all imputations created in parallel together
  impsets <- do.call(c, map(imps, "impDatasets"))
  imp <- imps[[1]]
  imp$impDatasets <- impsets
  
  warnings()
  
  saveRDS(imp, here::here('tmp', outnamei))
  #saveRDS(imp, here::here('data', outnamei))
}
