# helper function definitions
require(data.table)
require(Hmisc)

# function to make dummy coxph object with pooled parameters
make_pooled_cph <- function(fits) {
  # get one fit
  fit1 <- fits[[1]]
  coef1 <- fit1$coefficients
  vcov1 <- vcov(fit1)
  
  # get pooled parameters
  pooled <- summary(pool(fits))
  
  # update parameters
  coefp <- pooled[, 'estimate']
  names(coefp) <- pooled[, 'term']
  stopifnot(all(pooled[, 'term'] == names(coef1)))
  fit1$coefficients <- coefp
  
  # update covariance matrix
  diag(fit1$var) <- pooled[, 'std.error'] ^ 2
  return(fit1)
}

make_testdat <- function(impfname='imp_full_smcfcs_strat2_rcs5_i250_', stratecog=FALSE) {
  # define ranges of variables on normalized scale
  agerange <- sort(unique(c(seq(-3.6, 2.6, length.out=100), 0)))
  bmirange <- sort(unique(c(seq(-3.6, 2.6, length.out=100), 0)))
  smirange <- sort(unique(c(seq(-4.1, 4.4, length.out=100), 0)))
  smdrange <- sort(unique(c(seq(-8.597, 4.951, length.out=150), 0)))
  
  # retrieve rcs knots and means / sds of raw variables
  if (file.exists(here::here('resources', 'rcsknots', impfname))) {
    rcsknots <- readRDS(here::here('resources', 'rcsknots', impfname))
  } else {
    rcsknots <- readRDS(here::here('resources', 'rcsknots', str_replace(impfname, '_s(\\d+)p.rds', '')))
  }
  means_and_sds <- fread(here::here('resources', 'means_and_sds_cont.csv'))
  means <- set_names(means_and_sds$mean, means_and_sds$variable)
  sds <- set_names(means_and_sds$sd, means_and_sds$variable)
  
  # create datasets where one variable is changing
  testdatage <- CJ(
    age = agerange,
    sex_male = F,
    histo = 'adeno',
    ecog_bin = '1',
    ecog_bin_c_stage_early1 = 0,
    ecog_bin_c_stage_early2 = 0,
    histono_pa_c_stage_earlyTRUE = 0,
    histoother_c_stage_earlyTRUE = 0,
    histosquamous_c_stage_earlyTRUE = 0,
    bmi = 0.0,
    c_stage_int = ordered(1),
    smi=0,
    smd=0
  )
  testdatage[, rangeval:=copy(age)]
  
  testdatbmi <- CJ(
    age = 0.0,
    sex_male = F,
    histo = 'adeno',
    ecog_bin = '1',
    ecog_bin_c_stage_early1 = 0,
    ecog_bin_c_stage_early2 = 0,
    histono_pa_c_stage_earlyTRUE = 0,
    histoother_c_stage_earlyTRUE = 0,
    histosquamous_c_stage_earlyTRUE = 0,
    bmi = bmirange,
    # c_stage_int = ordered(1:4),
    c_stage_int = ordered(1),
    smi=0,
    smd=0
  )
  testdatbmi[, rangeval:=copy(bmi)]
  
  testdatsmi <- CJ(
    age = 0.0,
    sex_male = F,
    histo = 'adeno',
    ecog_bin = '1',
    ecog_bin_c_stage_early1 = 0,
    ecog_bin_c_stage_early2 = 0,
    histono_pa_c_stage_earlyTRUE = 0,
    histoother_c_stage_earlyTRUE = 0,
    histosquamous_c_stage_earlyTRUE = 0,
    bmi = 0.0,
    smi = smirange,
    # c_stage_int = ordered(1:4),
    c_stage_int = ordered(1),
    smd=0.0
  )
  testdatsmi[, rangeval:=copy(smi)]
  
  testdatsmd <- CJ(
    age = 0.0,
    sex_male = F,
    histo = 'adeno',
    ecog_bin = '1',
    ecog_bin_c_stage_early1 = 0,
    ecog_bin_c_stage_early2 = 0,
    histono_pa_c_stage_earlyTRUE = 0,
    histoother_c_stage_earlyTRUE = 0,
    histosquamous_c_stage_earlyTRUE = 0,
    smd = smdrange,
    bmi = 0.0,
    # c_stage_int = ordered(1:4),
    c_stage_int = ordered(1),
    smi=0.0
  )
  testdatsmd[, rangeval:=copy(smd)]
  
  testdatsmismd <- CJ(
    age = 0.0,
    sex_male = F,
    histo = 'adeno',
    ecog_bin = '1',
    ecog_bin_c_stage_early1 = 0,
    ecog_bin_c_stage_early2 = 0,
    histono_pa_c_stage_earlyTRUE = 0,
    histoother_c_stage_earlyTRUE = 0,
    histosquamous_c_stage_earlyTRUE = 0,
    bmi = 0.0,
    # c_stage_int = ordered(1:4),
    c_stage_int = ordered(1),
    smi=0:1,
    smd=smdrange
  )
  testdatsmismd[, rangeval:=copy(paste0('smd', smd, 'smi', smi))]
  
  testdat <- rbindlist(list(
    age=testdatage,
    bmi=testdatbmi,
    smi=testdatsmi,
    smd=testdatsmd,
    smismd=testdatsmismd), idcol='term', use.names=T)
  
  # add deterministically determined variables
  testdat[, smismd:=smi*smd]
  testdat[, age_raw:=(age * sds['age']) + means['age']]
  testdat[, bmi_raw:=(bmi * sds['bmi']) + means['bmi']]
  testdat[, smi_raw:=(smi * sds['smi']) + means['smi']]
  testdat[, smd_raw:=(smd * sds['smd']) + means['smd']]
  testdat[, smismd_raw:=(smismd * sds['smismd']) + means['smismd']]
  testdat[term=='age', rangeval_raw:=age_raw]
  testdat[term=='bmi', rangeval_raw:=bmi_raw]
  testdat[term=='smd', rangeval_raw:=smd_raw]
  testdat[term=='smi', rangeval_raw:=smi_raw]
  
  
  # make rcs variables
  testdat[, `:=`(
    c_stage_early=c_stage_int < 3,
    age1 = rcspline.eval(age, knots=rcsknots$age)[, 1],
    age2 = rcspline.eval(age, knots=rcsknots$age)[, 2],
    age3 = rcspline.eval(age, knots=rcsknots$age)[, 3],
    bmi1 = rcspline.eval(bmi, knots=rcsknots$bmi)[, 1],
    bmi2 = rcspline.eval(bmi, knots=rcsknots$bmi)[, 2],
    bmi3 = rcspline.eval(bmi, knots=rcsknots$bmi)[, 3],
    smi1 = rcspline.eval(smi, knots=rcsknots$smi)[, 1],
    smi2 = rcspline.eval(smi, knots=rcsknots$smi)[, 2],
    smi3 = rcspline.eval(smi, knots=rcsknots$smi)[, 3],
    smd1 = rcspline.eval(smd, knots=rcsknots$smd)[, 1],
    smd2 = rcspline.eval(smd, knots=rcsknots$smd)[, 2],
    smd3 = rcspline.eval(smd, knots=rcsknots$smd)[, 3]
  )]
  
  # make interaction terms
  testdat[, `:=`(
    smismd1 = smi*smd1,
    smismd2 = smi*smd2,
    smismd3 = smi*smd3,
    smi1smd = smi1*smd,
    smi2smd = smi2*smd,
    smi3smd = smi3*smd
  )]
  return(testdat)
}  

#' calculate std.error for covariate + interaction term
#' @param fit result of a lm, glm, coxph
#' @param var1 termname of variable1
#' @param var2 termname of variable2
add_se <- function(fit , var1, var2){
  vcovmat <- vcov(fit)  #grab the standard error of the coefficients
  ses <- sqrt(diag(vcovmat))
  cov12 <- vcovmat[var1, var2]
  var_out <- ses[var1]^2 + ses[var2]^2 + 2 * cov12
  return(var_out^(0.5))
}

#' for making better variable names
rename_terms <- function(term) {
  term %>% 
  str_replace_all("sex_maleTRUE", "male sex") %>%
  str_replace_all("bmi", "BMI") %>%
  str_replace_all("histo", "") %>%
  str_replace_all("no_pa", "no PA") %>%
  str_replace_all("other", "PA: other") %>%
  str_replace_all("squamous", "PA: squamous") %>%
  str_replace_all("ecog_bin2", "PS 2 or higher") %>%
  str_replace_all("ecog_bin", "PS ") %>%
  str_replace_all("pmd", "PMD") %>%
  str_replace_all("pmi", "PMI") %>%
  str_replace_all("smd", "PMD") %>%
  str_replace_all("smi", "PMI")
}

