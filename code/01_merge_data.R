### merge data from different sources

library(data.table)

datadir <- here::here('data')

## make names based on csv file with colnames
colnamedf <- fread(here::here('resources', 'L3_colnames.csv'))
measurementvars <- colnamedf[type=='measurement', unique(newname)]

## load measurements
l3 <- readRDS(here::here(datadir, 'l3.rds'))
setDT(l3)

## import clinical data
lc <- fread(here::here(datadir, '20210608_in_clinicaldata.csv'))
lc[, pid:=as.character(pid)]
lc[histo=='', histo:=NA]

## merge dates
l3dates <- l3[, list(pid, scan_date)]
datedf <- l3dates[lc, on='pid', nomatch=0]
datedf[, days_before_fup_start:=as.numeric(fup_start - scan_date)]
datedf[order(abs(days_before_fup_start)), timediff_scan_fupstart_order:=1:.N, by='episode']
closestdates <- datedf[timediff_scan_fupstart_order==1]

## pick only dates within 90 days before or after fup start
chosendates <- closestdates[days_before_fup_start > -30 & days_before_fup_start < 90]

l3[, .N, by=c('pid', 'scan_date')][, table(N)]
l3e <- l3[chosendates[, list(pid, scan_date, days_before_fup_start)], on=c('pid', 'scan_date'), nomatch=0]
mvars <- c('scan_date', 'days_before_fup_start')
lc[chosendates, (mvars):=mget(paste0('i.', mvars)), on='episode']

## remove duplicated patients by taking only the first followupdate
lc[order(fup_start), episode_order:=1:.N, by='pid']
lc <- lc[episode_order==1]
lc[, episode_order:=NULL]

## get measurement columns
mvars <- c(measurementvars, 'l3_slice_index', 'dicom_serie_pks', 'series_uid')
# lc[l3, (mvars):=mget(paste0('i.', mvars)), on=c('pid', 'scan_date')]
lc[l3e, (mvars):=mget(paste0('i.', mvars)), on=c('pid', 'scan_date')]

## remove density measurements in scans with iv contrast
contrastdf <- fread(here::here(datadir, '20220121_ivcontrast_check.csv'))
contrastdf[, `:=`(
  pid=str_extract(record_id, '^\\d+(?=s)'),
  dicom_serie_pks=as.integer(str_extract(record_id, '(?<=s)\\d+(?=.png)'))
  )]
lc[contrastdf, `:=`(contrasttype=i.label_name), on=c('pid', 'dicom_serie_pks')]
hu_vars <- str_subset(measurementvars, '^mean_HU')
lc[contrasttype=='iv', (hu_vars):=NA]

## export
fwrite(lc, here::here(datadir, 'clinical_and_measurements_complete.csv'), row.names=F)
