

# load graphics libraries
library(plyr)
library(ggplot2)
library(reshape2)
library(tidyverse)
library(knitr)
library(kableExtra)
library(plotly)
library(stringr)
library(bsselectR)
library(lubridate)
library(textclean)

# returns the previous month if today's date is earlier than the 15th
get_month <- function(){
  if (day(Sys.Date()) < 15)
    return (Sys.Date() %m-% months(1))
  else
    return (Sys.Date())
}

prepare_init_df <- function(df, mri, sites){
  
  #for the purposes of this report - only consider the first session of the baseline scan - might need to modify later if want to count followup
  mri <- subset(mri, mri$redcap_repeat_instrument =="optimum_neuro_mri")
  mri <- subset(mri, mri$redcap_event_name =="baseline_arm_6"|mri$redcap_event_name =="6_mth_fu_arm_6")
  
  #merge mri and df dataframes
  df <- merge(df, mri, by=c('record_id', 'redcap_event_name'), all.x=TRUE)
  df$redcap_event_name <- as.factor(df$redcap_event_name)
  #pull out site code from subject IDs
  df$site <- substring(df$record_id, 1, 2)
  
  #make sure that site data is a factor
  df$site <- as.factor(df$site)
  
  
  #make sure all site codes are uppercase
  df$site <- toupper(df$site)
  
  #remove any sites that don't match expected site codes  
  df <- df[(df$site == 'CU' | df$site == 'LA' | df$site == 'UP' | df$site == 'UT' | df$site == 'WU'),]
  
  #rename timepoint variable
  names(df)[names(df)=="redcap_event_name"] <- "timepoint"
  
  #rename factors in timepoint variable
  levels(df$timepoint) <- c(levels(df$timepoint), "baseline", "6 mth FU", "24 mth FU") 
  df$timepoint[df$timepoint == "baseline_arm_6"]  <- "baseline"
  df$timepoint[df$timepoint == "6_mth_fu_arm_6"]  <- "6 mth FU" 
  df$timepoint[df$timepoint == "24_mth_fu_arm_6"]  <- "24 mth FU" 
  df$timepoint <- droplevels(df$timepoint, exclude = c("baseline_arm_6", "6_mth_fu_arm_6", "24_mth_fu_arm_6"))
  
  return(df)
}

prepare_recruit_df <- function(recruit_df){
  
  #make sure important dates are in date format
  recruit_df$dkefs_date <- as.Date(recruit_df$dkefs_date_admin, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to date format
  recruit_df$rbans_date <- as.Date(recruit_df$rbans_date_admin, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to date format
  
  #get current date and add to dataframe
  recruit_df$currentYm <- as.character(Sys.Date()) #add current sysdate to df
  recruit_df$currentYm_str <- as.character(substr(recruit_df$currentYm, 1, 7)) #store month and year, as character
  
  #add column that indicates if consented, consent date, and if this month
  recruit_df$consented[is.na(recruit_df$meta_consent_date)] <- 0
  recruit_df$consented[is.na(recruit_df$consented)] <- 1
  recruit_df$meta_consent_date <- as.Date(recruit_df$meta_consent_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to date format
  recruit_df$meta_consent_mth <- as.character(substr(recruit_df$meta_consent_date, 1, 7)) #store month and year, as character
  recruit_df$meta_consent_mthT <- recruit_df$meta_consent_mth %in% recruit_df$currentYm_str #logical T or F- was consent this month?
  
  #add column that indicates if completed, complete date, and if this month
  recruit_df$completed <- ifelse(recruit_df$meta_terminate_reason == 3, 1, 0)
  recruit_df$meta_terminate_date <- as.Date(recruit_df$meta_terminate_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to date format
  recruit_df$meta_terminate_mth <- as.character(substr(recruit_df$meta_terminate_date, 1, 7)) #store month and year, as character
  recruit_df$meta_terminate_mthT <- recruit_df$meta_terminate_mth %in% recruit_df$currentYm_str #logical T or F- was complete this month?
  recruit_df$meta_complete_mthT <- ifelse(recruit_df$completed == 1 & recruit_df$meta_terminate_mthT == TRUE, TRUE, FALSE)
  
  #add column that indicates if terminated
  recruit_df$terminated <- ifelse(recruit_df$meta_terminate_reason == 3, 0, 1)
  
  #add columns that indicate timepoint at which terminated (for reasons other than completion)
  #GABI note to self: I don't think this would work...
  recruit_df$terminate_baseline <- ifelse(recruit_df$meta_terminate_reason != 3 & recruit_df$timepoint == 'baseline', 1, NA)
  recruit_df$terminate_FU_06_MTH <- ifelse(recruit_df$meta_terminate_reason != 3 & recruit_df$timepoint == '6 mth FU', 1, NA)
  recruit_df$terminate_FU_24_MTH <- ifelse(recruit_df$meta_terminate_reason != 3 & recruit_df$timepoint == '24 mth FU', 1, NA)
  
  return(recruit_df)
}

prepare_enroll_df <- function(enroll_df, enroll_criteria=2) {
  #NOTE TO SELF FOR GABI: It seems like this is where enroll df starts only counting baseline. I think that I could just filter out 6 month follow up from baseline and calculate that way???  
  #NOTE TO SELF FOR GABI: I JUST COMMENTED THIS OUT 
  #take only baseline
  enroll_df <- enroll_df[(enroll_df$timepoint == 'baseline'),] 
  
  today_str <- as.character(Sys.Date())
  
  #first, check if T1 complete, and if date is this month
  #note - here, 'complete' means just T1 complete
  names(enroll_df)[names(enroll_df) == "mr_t1"] <- "enroll_mri" #change name of status column
  names(enroll_df)[names(enroll_df) == "mr_date"] <- "enroll_mri_date" #change name of date column
  
  enroll_df$enroll_mri <- ifelse(enroll_df$enroll_mri == 1, 1, NA) #turn no to NA
  enroll_df$enroll_mri_init <- ifelse(enroll_df$enroll_mri == 1, 1, NA) #turn no to NA
  enroll_df$enroll_mri_date <- as.Date(enroll_df$enroll_mri_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  enroll_df$enroll_mri_mth <- as.character(substr(enroll_df$enroll_mri_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_mri_mthT <- enroll_df$enroll_mri_mth %in% enroll_df$currentYm_str #logical T or F- was MRI this month?
  enroll_df$mri_fu_due <-  ifelse(enroll_df$enroll_mri_date %m+% months(6) < today_str, 1, 0)
  enroll_df$mri_fu_7_mos <-  ifelse(enroll_df$enroll_mri_date %m+% months(7) < today_str, 1, 0)
  enroll_df$mri_fu_mth <- as.character(substr(enroll_df$enroll_mri_date %m+% months(6), 1, 7)) 
  enroll_df$mri_fu_mthT <- enroll_df$mri_fu_mth %in% enroll_df$currentYm_str
  #note: one thing to consider, do we want to count 6mos from that specific assessment, or anchor everything to their consent date? leaning toward the latter tbh.
  
  #second, check if blood complete, and if date is this month
  #note - here, 'complete' means that just the plasma is complete (don't care about serum)
  enroll_df$enroll_bld <- ifelse(enroll_df$plasma_blood_status == 1, 1, NA) #turn no and partial to NA
  enroll_df$enroll_bld_date <- substr(enroll_df$plasma_blood_date, 1, 11)#cut out the time contained in datetime
  enroll_df$enroll_bld_date <- as.Date(enroll_df$enroll_bld_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert to datetime 
  enroll_df$enroll_bld_mth <- as.character(substr(enroll_df$enroll_bld_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_bld_mthT <- enroll_df$enroll_bld_mth %in% enroll_df$currentYm_str #logical T or F - was blood this month?    
  enroll_df$bld_fu_due <-  ifelse(enroll_df$enroll_bld_date %m+% months(6) < today_str, 1, 0)
  enroll_df$bld_fu_7_mos <-  ifelse(enroll_df$enroll_bld_date %m+% months(7) < today_str, 1, 0) 
  
  #third, check if neuropsych complete, and if data is this month
  #note - here, 'complete' means rbans and dkefs done in entirity
  #make new columns with completion and date data for 2 assessments we care about - dkefs and rbans
  enroll_df$dkefs <- ifelse(enroll_df$dkefs_complete == 1, 1, NA) #turn no and partial to NA
  enroll_df$enroll_dkefs_date <- as.Date(enroll_df$dkefs_date_admin, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  enroll_df$enroll_dkefs_mth <- as.character(substr(enroll_df$enroll_dkefs_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_dkefs_mthT <- enroll_df$enroll_dkefs_mth %in% enroll_df$currentYm_str #logical T or F- was dkefs this month?
  
  enroll_df$rbans <- ifelse(enroll_df$rbans_complete == 1, 1, NA) #turn no and partial to NA
  enroll_df$enroll_rbans_date <- as.Date(enroll_df$rbans_date_admin, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  enroll_df$enroll_rbans_mth <- as.character(substr(enroll_df$enroll_rbans_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_rbans_mthT <- enroll_df$enroll_rbans_mth %in% enroll_df$currentYm_str #logical T or F- was rbans this month?
  
  #take the sum of columns - want both of the central neuropsych assessments (dkefs, rbans) to be complete
  enroll_df[, c('dkefs', 'rbans')] <- sapply(enroll_df[, c('dkefs', 'rbans')], as.numeric) #change to numeric
  enroll_df$np_count <- rowSums(enroll_df[,c('dkefs', 'rbans')]) #sum
  
  #make an enrolled column for neuropsych - if sum is 2 
  enroll_df$enroll_np <- ifelse(enroll_df$np_count == 2, 1, 0)
  
  #identify the second date in ordered dates - i.e., when enrollment for np criteria is met
  np_dates <- c('dkefs_date', 'rbans_date')
  enroll_df$enroll_np_date <- NA
  for (n in 1:nrow(enroll_df)){
    dates <- enroll_df[n, np_dates] # pull out NP assessment dates
    dates <- dates[order(dates)] # order dates
    if (length(dates[!is.na(dates)]) == 2) { # if there are 2 assesments
      enroll_df$enroll_np_date[n] <- dates[[2]]   # pull out the second assessment date
    } else {
      enroll_df$enroll_np_date[n] <- NA         # otherwise print 'NA'
    }
  }
  
  #determine if np criteria was met this month
  enroll_df$enroll_np_date <- as.Date(enroll_df$enroll_np_date, origin = "1970-01-01") #turn back into date format
  enroll_df$enroll_np_mth <- as.character(substr(enroll_df$enroll_np_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_np_mthT <- enroll_df$enroll_np_mth %in% enroll_df$currentYm_str #logical T or F- was criteria met this month?
  enroll_df$np_fu_due <-  ifelse(enroll_df$enroll_np_date %m+% months(6) < today_str, 1, 0)
  enroll_df$np_fu_7_mos <-  ifelse(enroll_df$enroll_np_date %m+% months(7) < today_str, 1, 0)
  
  
  #determine if enrollment criteria is met, i.e., at least 2/3 of mri, blood, and np completed
  enroll_df$enroll <- rowSums(enroll_df[,c('enroll_mri', 'enroll_bld', 'enroll_np')], na.rm = T) #calculate sum
  enroll_df$enroll <- ifelse(enroll_df$enroll >= enroll_criteria, 1, 0)
  
  #determine date of 2/3 of mri, blood, and psy completed
  enroll_dates <- c('enroll_mri_date','enroll_bld_date', 'enroll_np_date')
  enroll_df$enroll_date <- NA
  for (n in 1:nrow(enroll_df)){
    dates <- enroll_df[n, enroll_dates] # pull out component dates
    dates <- dates[order(dates)] # order dates
    if (length(dates[!is.na(dates)]) >= 2) { # if there are 2 or more components
      enroll_df$enroll_date[n] <- dates[[2]]   # pull out the second component dates
    } else {
      enroll_df$enroll_date[n] <- NA         # otherwise print 'NA'
    }
  }
  enroll_df$enroll_date <- as.Date(enroll_df$enroll_date, origin = "1970-01-01") #turn back into date format
  
  #determine if enrollment was met this month
  enroll_df$enroll_mth <- as.character(substr(enroll_df$enroll_date, 1, 7)) #store month and year, as character
  enroll_df$enroll_mthT <- enroll_df$enroll_mth %in% enroll_df$currentYm_str #logical T or F- was criteria met this month?
  
  return(enroll_df)
  
}

prep_fu_df <- function(fu_df){
  
  fu_df$currentYm <- as.character(get_month()) #add current sysdate to df
  fu_df$currentYm_str <- as.character(substr(fu_df$currentYm, 1, 7)) #store month and year, as character
  
  
  #first, check if T1 complete, and if date is this month
  #note - here, 'complete' means just T1 complete
  names(fu_df)[names(fu_df) == "mr_t1"] <- "enroll_mri" #change name of status column
  names(fu_df)[names(fu_df) == "mr_date"] <- "enroll_mri_date" #change name of date column
  
  
  
  fu_df$enroll_mri <- ifelse(fu_df$enroll_mri == 1, 1, NA) #turn no to NA
  fu_df$enroll_mri_date <- as.Date(fu_df$enroll_mri_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  fu_df$enroll_mri_mth <- as.character(substr(fu_df$enroll_mri_date, 1, 7)) #store month and year, as character
  fu_df$enroll_mri_mthT <- fu_df$enroll_mri_mth %in% fu_df$currentYm_str #logical T or F- was MRI this month?  #OKAY SO they must also have T1 completed
  
  
  # if they got an MRI this month then: fu_df %>% filter(enroll_mri==1, enroll_mri_mthT==TRUE)
  
  #second, check if blood complete, and if date is this month
  #note - here, 'complete' means that just the plasma is complete (don't care about serum)
  fu_df$enroll_bld <- ifelse(fu_df$plasma_blood_status.x == 1, 1, NA) #turn no and partial to NA
  fu_df$enroll_bld_date <- substr(fu_df$plasma_blood_date.x, 1, 11)#cut out the time contained in datetime
  fu_df$enroll_bld_date <- as.Date(fu_df$enroll_bld_date, format = "%Y-%m-%d", origin = "1970-01-01") #convert to datetime 
  fu_df$enroll_bld_mth <- as.character(substr(fu_df$enroll_bld_date, 1, 7)) #store month and year, as character
  fu_df$enroll_bld_mthT <- fu_df$enroll_bld_mth %in% fu_df$currentYm_str #logical T or F - was blood this month?    
  
  #third, check if neuropsych complete, and if data is this month
  #note - here, 'complete' means rbans and dkefs done in entirity
  #make new columns with completion and date data for 2 assessments we care about - dkefs and rbans
  fu_df$dkefs <- ifelse(fu_df$dkefs_complete.x == 1, 1, NA) #turn no and partial to NA
  fu_df$enroll_dkefs_date <- as.Date(fu_df$dkefs_date_admin.x, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  fu_df$enroll_dkefs_mth <- as.character(substr(fu_df$enroll_dkefs_date, 1, 7)) #store month and year, as character
  fu_df$enroll_dkefs_mthT <- fu_df$enroll_dkefs_mth %in% fu_df$currentYm_str #logical T or F- was dkefs this month?
  
  fu_df$rbans <- ifelse(fu_df$rbans_complete.x == 1, 1, NA) #turn no and partial to NA
  fu_df$enroll_rbans_date <- as.Date(fu_df$rbans_date_admin.x, format = "%Y-%m-%d", origin = "1970-01-01") #convert date to datetime format
  fu_df$enroll_rbans_mth <- as.character(substr(fu_df$enroll_rbans_date, 1, 7)) #store month and year, as character
  fu_df$enroll_rbans_mthT <- fu_df$enroll_rbans_mth %in% fu_df$currentYm_str #logical T or F- was rbans this month?
  
  #take the sum of columns - want both of the central neuropsych assessments (dkefs, rbans) to be complete
  fu_df[, c('dkefs', 'rbans')] <- sapply(fu_df[, c('dkefs', 'rbans')], as.numeric) #change to numeric
  fu_df$np_count <- rowSums(fu_df[,c('dkefs', 'rbans')]) #sum
  
  #make an enrolled column for neuropsych - if sum is 2 
  fu_df$enroll_np <- ifelse(fu_df$np_count == 2, 1, 0)
  
  #identify the second date in ordered dates - i.e., when enrollment for np criteria is met
  np_dates <- c('dkefs_date.x', 'rbans_date.x')
  fu_df$enroll_np_date <- NA
  for (n in 1:nrow(fu_df)){
    dates <- fu_df[n, np_dates] # pull out NP assessment dates
    dates <- dates[order(dates)] # order dates
    if (length(dates[!is.na(dates)]) == 2) { # if there are 2 assesments
      fu_df$enroll_np_date[n] <- dates[[2]]   # pull out the second assessment date
    } else {
      fu_df$enroll_np_date[n] <- NA         # otherwise print 'NA'
    }
  }
  
  #determine if np criteria was met this month
  fu_df$enroll_np_date <- as.Date(fu_df$enroll_np_date, origin = "1970-01-01") #turn back into date format
  fu_df$enroll_np_mth <- as.character(substr(fu_df$enroll_np_date, 1, 7)) #store month and year, as character
  fu_df$enroll_np_mthT <- fu_df$enroll_np_mth %in% fu_df$currentYm_str #logical T or F- was criteria met this month?
  
  #determine if enrollment criteria is met, i.e., at least 2/3 of mri, blood, and np completed
  fu_df$enroll <- rowSums(fu_df[,c('enroll_mri', 'enroll_bld', 'enroll_np')], na.rm = T) #calculate sum
  fu_df$enroll <- ifelse(fu_df$enroll >= 2, 1, 0)
  
  #determine date of 2/3 of mri, blood, and psy completed
  enroll_dates <- c('enroll_mri_date','enroll_bld_date', 'enroll_np_date')
  fu_df$enroll_date <- NA
  for (n in 1:nrow(fu_df)){
    dates <- fu_df[n, enroll_dates] # pull out component dates
    dates <- dates[order(dates)] # order dates
    if (length(dates[!is.na(dates)]) >= 2) { # if there are 2 or more components
      fu_df$enroll_date[n] <- dates[[2]]   # pull out the second component dates
    } else {
      fu_df$enroll_date[n] <- NA         # otherwise print 'NA'
    }
  }
  fu_df$enroll_date <- as.Date(fu_df$enroll_date, origin = "1970-01-01") #turn back into date format
  
  #determine if enrollment was met this month
  fu_df$enroll_mth <- as.character(substr(fu_df$enroll_date, 1, 7)) #store month and year, as character
  fu_df$enroll_mthT <- fu_df$enroll_mth %in% fu_df$currentYm_str #logical T or F- was criteria met this month?
  
  return(fu_df)
  
}

make_enroll_table <- function(df, sites, targets){
  #initialize dataframe and names columns and rows
  out_table <- data.frame(matrix(ncol=11, nrow=6))
  
  target_month <- filter(targets, month==format(get_month(), '%Y-%m'))$per_month
  target_total <- filter(targets, month==format(get_month(), '%Y-%m'))$target
  
  #names of columns and rows on table
  names(out_table) <- c(paste0(c('n', '%', 'n', '%', 'n', '%', 'n', '%', 'n', '%', 'total')))
  row.names(out_table) <- c(
    paste0("completed neuropsych - ", as.yearmon(get_month())),
    "completed neuropsych to date", 
    paste0("completed blood - ", as.yearmon(get_month())),
    "completed blood to date",
    paste0("completed MRI - ", as.yearmon(get_month())),
    "completed MRI to date")
  
  df_24 <- filter(df, timepoint=='baseline')
  by_site <- group_by(df_24, site, .drop=FALSE)
  
  # TODO: fix this to use dplyr/summarise() instead of loops
  
  # Neuropsych
  neuro_24 <- filter(by_site, dkefs_complete==1, rbans_complete==1)
  i <- 1
  for (site in sites){
    # this month
    out_table[1,i] = sum( filter( neuro_24, as.yearmon(dkefs_date_admin) == as.yearmon(get_month()) )$site == site )
    # all time
    out_table[2,i] = sum(neuro_24$site == site)
    i <- i+2
  }
  
  # Blood
  blood_24 <- filter(by_site, plasma_blood_status==1)
  i <- 1
  for (site in sites){
    # this month
    out_table[3,i] = sum( filter( blood_24, as.yearmon(plasma_blood_date) == as.yearmon(get_month()) )$site == site )
    # all time
    out_table[4,i] = sum(blood_24$site == site)
    i <- i+2
  }
  
  # MRI 
  mri_24 <- filter(by_site, mr_t1==1)
  
  if (dim(mri_24)[1] !=0){
    i <- 1
    for (site in sites){
      # this month
      out_table[5,i] = sum( filter( mri_24, as.yearmon(mr_date) == as.yearmon(get_month()) )$site == site )
      # all time
      out_table[6,i] = sum(mri_24$site == site)
      i <- i+2
    }
  } 
  else {
    i <- 1
    for (site in sites){
      out_table[5,i] = 0
      out_table[6,i] = 0
      i <- i+2
    }
  }
  
  # add totals
  for (i in 1:nrow(out_table)){
    out_table[i,ncol(out_table)] = sum(out_table[i,], na.rm=TRUE)
  }
  
  # add percentages of targets
  i <- 2
  for (site in sites){
    out_table[1,i] = percent(out_table[1,i-1]/target_month)
    out_table[2,i] = percent(out_table[2,i-1]/target_total)
    out_table[3,i] = percent(out_table[3,i-1]/target_month)
    out_table[4,i] = percent(out_table[4,i-1]/target_total)
    out_table[5,i] = percent(out_table[5,i-1]/target_month)
    out_table[6,i] = percent(out_table[6,i-1]/target_total)
    i <- i+2
  }
  
  return(out_table)
}

make_recruit_table <- function(recruit_df, sites, targets) {
  
  # create vector of recruitment variables
  recruit_vars <- c(
    "meta_consent_mthT.x",       #consented current month
    "consented.x",               #consented to date
    "enroll_mthT",               #enrolled current month        no
    "enroll",                    #enrolled to date              no
    "meta_complete_mthT.x",      #completed current month
    "completed.x",               #completed to date
    "terminate_baseline.x",      #terminated baseline
    "terminate_FU_06_MTH.x",     #terminated 6 mth FU
    "terminate_FU_24_MTH.x")     #terminated 24 mth FU
  
  #turn all NAs into 0
  recruit_df[, recruit_vars] <- apply(recruit_df[, recruit_vars], 2, function(x){replace(x, is.na(x), 0)})
  
  #initialize dataframe and names columns and rows
  recruit_table <- data.frame(matrix(ncol=11, nrow=length(recruit_vars)))
  
  #names of columns and rows on demo_table
  names(recruit_table) <- c(paste0(c('n', '%', 'n', '%', 'n', '%', 'n', '%', 'n', '%', 'total')))
  row.names(recruit_table) <- c(
    paste0("consented - ", as.yearmon(get_month())),
    "consented to date", 
    paste0("enrolled - ", as.yearmon(get_month())),
    "enrolled to date",
    paste0("completed - ", as.yearmon(get_month())),
    "completed to date",
    "terminated during baseline",
    "terminated during 6 month follow-up",
    "terminated during 24 month follow-up")
  
  consented <- filter(recruit_df, !is.na(meta_consent_date.x)) 
  consented_this_month <- filter(consented, 
                                 as.yearmon(meta_consent_date.x) == as.yearmon(get_month()))
  
  # pt is considered enrolled if they have at least 1 of 3 measures at baseline
  enrolled <- filter(recruit_df, timepoint=='baseline') %>%
              mutate(mr_complete = mr_t1==1, 
                     blood_complete = plasma_blood_status.x==1) %>%
              mutate(np_complete = dkefs_complete.x == 1 & rbans_complete.x == 1) %>%
              mutate(n_measures = rowSums(select(., np_complete, mr_complete, blood_complete), na.rm=TRUE)
                     ) %>%
              filter(n_measures >= 1)
  
  enrolled_this_month <- filter(enrolled, 
                                 as.yearmon(enroll_date) == as.yearmon(get_month()))
  
  
  # Completed could be considered as having 1/3 measures completed at 24 mo:
  # completed <- filter(recruit_df, timepoint=='24 mth FU') %>%
  #               mutate(mr_complete = mr_t1==1,
  #                      blood_complete = plasma_blood_status.x==1) %>%
  #               mutate(np_complete = dkefs_complete.x == 1 & rbans_complete.x == 1) %>%
  #               mutate(n_measures = rowSums(select(., np_complete, mr_complete, blood_complete), na.rm=TRUE)
  #                       ) %>%
  #               filter(n_measures >= 1)
  
  # Pts are counted as completed when their baseline tracking form in REDCap is
  # marked as 'completed'
  completed <- filter(recruit_df, completed.x==1)
  
  # TODO: figure out what date is best to use to determine when completed
  # completed_this_month <- filter(completed, 
  #                               as.yearmon(mr_date) == as.yearmon(get_month()))
  
  # terminated <- filter(recruit_df, !is.na(meta_terminate_date.x)) 
  terminated = filter(recruit_df, terminated.x==1)
  terminated_bl <- filter(terminated, timepoint=='baseline')
  terminated_6m <- filter(terminated, timepoint=='6 mth FU')
  terminated_24 <- filter(terminated, timepoint=='24 mth FU')
  
  target_month <- filter(targets, month==format(get_month(), '%Y-%m'))$per_month
  target_total <- filter(targets, month==format(get_month(), '%Y-%m'))$target
  
  i <- 1
  for (site in sites){
    # Consent 
    recruit_table[1,i] = nrow(filter(consented_this_month, site.x==site))
    recruit_table[1,i+1] = percent( recruit_table[1,i]/target_month )
    
    recruit_table[2,i] = nrow(filter(consented, site.x==site))
    recruit_table[2,i+1] = percent( recruit_table[[2,i]]/target_total )
    
    # Enrollment
    recruit_table[3,i] = nrow(filter(enrolled_this_month, site.x==site))
    recruit_table[3,i+1] = percent( recruit_table[3,i]/target_month )
    
    recruit_table[4,i] = nrow(filter(enrolled, site.x==site))
    recruit_table[4,i+1] = percent( recruit_table[[4,i]]/target_total )
    
    # Completion
    # recruit_table[5,i] = nrow(filter(completed_this_month, site.x==site))
    # recruit_table[5,i+1] = percent( recruit_table[5,i]/target_month )
    
    recruit_table[6,i] = nrow(filter(completed, site.x==site))
    # TODO: figure out how many people should be completed based on 2 years prior
    # recruit_table[6,i+1] = percent( recruit_table[[6,i]]/target_total )
    
    # Termination
    recruit_table[7,i] = nrow(filter(terminated_bl, site.x==site))
    recruit_table[8,i] = nrow(filter(terminated_6m, site.x==site))
    recruit_table[9,i] = nrow(filter(terminated_24, site.x==site))
    
    i <- i+2
  }
  
  # add totals to last column
  recruit_table[1,11] = nrow(consented_this_month)
  recruit_table[2,11] = nrow(consented)
  
  recruit_table[3,11] = nrow(enrolled_this_month)
  recruit_table[4,11] = nrow(enrolled)
  
  recruit_table[6,11] = nrow(completed)
  
  recruit_table[7,11] = nrow(terminated_bl)
  recruit_table[8,11] = nrow(terminated_6m)
  recruit_table[9,11] = nrow(terminated_24)
  
  
  #add row names to make like demo table
  recruit_table <- cbind(metric = rownames(recruit_table), recruit_table)
  rownames(recruit_table) <- NULL
  
  
  return(recruit_table)
}

make_demo_table <- function(df, sites){
  
  # create a separate dataframe of data from just baseline timepoint
  demo_df <- df[(df$timepoint == 'baseline'),]
  
  # modify race variable to have 2 levels
  demo_df$demo_race <- ifelse(demo_df$demo_race != 2, 1, 2)
  
  # create vector of demographic variables
  demo_vars <- c(
    'demo_sex',
    'demo_ethnicity',
    'demo_race',
    'demo_age',
    'demo_edu')
  
  # initialize dataframe 
  demo_table <- data.frame(matrix(ncol=6, nrow=length(demo_vars)))
  
  # names of columns and rows on demo_table
  names(demo_table) <- c(paste0(c('CU', 'LA', 'UP', 'UT', 'WU', 'p')))
  row.names(demo_table) <- c('sex (M:F:O)', 'ethnicity (H:NH)', 'race (minority:non-minority)', 'age', 'education')
  
  # initialize counters (j = row, k = column)
  j <- 1
  k <- 1
  
  # for loop for categorical
  for (var in demo_vars[1:3]) {
    
    # count observations for each cluster
    for (site in sites) {
      out <- table(demo_df$site, demo_df[[var]])
      demo_table[j,k] <- paste(out[site,], collapse = ':')
      k <- k + 1}
    
    # run chi-square tests 
    chi_sq <- chisq.test(out)
    demo_table[j,k] <- sub("^(-?)0.", "\\1.", sprintf("%.3f", chi_sq$p.value))
    
    k <- 1
    j <- j + 1} 
  
  # initialize counters (j = row, k = column)
  j <- 4
  k <- 1
  
  # for loop for continuous variables
  for (var in demo_vars[4:5]) {
    
    # calculate means and SDs for each site
    for (site in sites) {
      M <- sprintf('%.02f', mean(demo_df[demo_df$site == site, var], na.rm = TRUE))
      SD <- sprintf('%.02f', sd(demo_df[demo_df$site == site, var], na.rm = TRUE))
      demo_table[j,k] <- paste( M,' (',SD,')', sep='') 
      
      k <- k + 1}
    
    # run one-way ANOVA with cluster as between-subjects variable
    F_test <- aov(demo_df[[var]] ~ demo_df$site, na.action=na.omit)
    # extract p-value
    F_test.p.value <- summary(F_test)[[1]][["Pr(>F)"]][[1]]
    # rounded p-value to 3 decimals and without leading zero
    F_test.p <- sub("^(-?)0.", "\\1.", sprintf("%.3f", F_test.p.value))
    demo_table[j,k] <- F_test.p 
    
    k <- 1 
    j <- j + 1}
  
  
  return(demo_table)
  
}

make_fu_table <- function(df, sites){
  
  #initialize dataframe and names columns and rows
  fu_table <- data.frame(matrix(ncol=11, nrow=6))
  
  #names of columns and rows on table
  names(fu_table) <- c(paste0(c('n', '%', 'n', '%', 'n', '%', 'n', '%', 'n', '%', 'total')))
  row.names(fu_table) <- c(
    paste0("completed neuropsych - ", as.yearmon(get_month())),
    "completed neuropsych to date", 
    paste0("completed blood - ", as.yearmon(get_month())),
    "completed blood to date",
    paste0("completed MRI - ", as.yearmon(get_month())),
    "completed MRI to date")
  
  by_site <- filter(df, timepoint=='6 mth FU') %>%
              group_by(site, .drop=FALSE)
  
  # TODO: fix this to use dplyr/summarise() instead of loops
  
  # Neuropsych
  neuro_df <- filter(by_site, dkefs_complete==1, rbans_complete==1)
  i <- 1
  for (site in sites){
    # this month
    fu_table[1,i] = sum( filter( neuro_df, as.yearmon(dkefs_date_admin) == as.yearmon(get_month()) )$site == site )
    # all time
    fu_table[2,i] = sum(neuro_df$site == site)
    i <- i+2
  }
  
  # Blood
  blood_df <- filter(by_site, plasma_blood_status==1)
  i <- 1
  for (site in sites){
    # this month
    fu_table[3,i] = sum( filter( blood_df, as.yearmon(plasma_blood_date) == as.yearmon(get_month()) )$site == site )
    # all time
    fu_table[4,i] = sum(blood_df$site == site)
    i <- i+2
  }
  
  # MRI 
  mri_df <- filter(by_site, mr_t1==1)
  
  if (dim(mri_df)[1] !=0){
    i <- 1
    for (site in sites){
      # this month
      fu_table[5,i] = sum( filter( mri_df, as.yearmon(mr_date) == as.yearmon(get_month()) )$site == site )
      # all time
      fu_table[6,i] = sum(mri_df$site == site)
      i <- i+2
    }
  } 
  else {
    i <- 1
    for (site in sites){
      fu_table[5,i] = 0
      fu_table[6,i] = 0
      i <- i+2
    }
  }
  
  # add totals
  for (i in 1:nrow(fu_table)){
    fu_table[i,ncol(fu_table)] = sum(fu_table[i,], na.rm=TRUE)
  }
  
  # add percentages of targets
  neuro_due <- filter(df, timepoint=='baseline',
                      dkefs_complete==1,
                      as.Date(dkefs_date_admin) %m+% months(6) < get_month()) %>%
              count(site)
  blood_due <- filter(df, timepoint=='baseline',
                      plasma_blood_status==1,
                      as.Date(plasma_blood_date) %m+% months(6) < get_month()) %>%
              count(site)
  mri_due   <- filter(df, timepoint=='baseline',
                      mr_t1==1,
                      as.Date(mr_date) %m+% months(6) < get_month()) %>%
              count(site)
  
  i <- 2
  for (s in sites){
    fu_table[2,i] = percent(fu_table[2,i-1]/(neuro_due %>% filter(site==s))$n)
    fu_table[4,i] = percent(fu_table[4,i-1]/(blood_due %>% filter(site==s))$n)
    fu_table[6,i] = percent(fu_table[6,i-1]/(mri_due %>% filter(site==s))$n)
    i <- i+2
  }
  return(fu_table)
}

make_fu_24_table <- function(df, sites) {
  #initialize dataframe and names columns and rows
  fu_table <- data.frame(matrix(ncol=11, nrow=6))
  
  #names of columns and rows on table
  names(fu_table) <- c(paste0(c('n', '%', 'n', '%', 'n', '%', 'n', '%', 'n', '%', 'total')))
  row.names(fu_table) <- c(
    paste0("completed neuropsych - ", as.yearmon(get_month())),
    "completed neuropsych to date", 
    paste0("completed blood - ", as.yearmon(get_month())),
    "completed blood to date",
    paste0("completed MRI - ", as.yearmon(get_month())),
    "completed MRI to date")

  df_24 <- filter(df, timepoint=='24 mth FU')
  by_site <- group_by(df_24, site, .drop=FALSE)
  
  # TODO: fix this to use dplyr/summarise() instead of loops
  
  # Neuropsych
  neuro_24 <- filter(by_site, dkefs_complete==1, rbans_complete==1)
  i <- 1
  for (site in sites){
    # this month
    fu_table[1,i] = sum( filter( neuro_24, as.yearmon(dkefs_date_admin) == as.yearmon(get_month()) )$site == site )
    # all time
    fu_table[2,i] = sum(neuro_24$site == site)
    i <- i+2
  }
  
  # Blood
  blood_24 <- filter(by_site, plasma_blood_status==1)
  i <- 1
  for (site in sites){
    # this month
    fu_table[3,i] = sum( filter( blood_24, as.yearmon(plasma_blood_date) == as.yearmon(get_month()) )$site == site )
    # all time
    fu_table[4,i] = sum(blood_24$site == site)
    i <- i+2
  }
  
  # MRI 
  mri_24 <- filter(by_site, mr_t1==1)
  
  if (dim(mri_24)[1] !=0){
    i <- 1
    for (site in sites){
      # this month
      fu_table[5,i] = sum( filter( mri_24, as.yearmon(mr_date) == as.yearmon(get_month()) )$site == site )
      # all time
      fu_table[6,i] = sum(mri_24$site == site)
      i <- i+2
    }
  } 
  else {
    i <- 1
    for (site in sites){
      fu_table[5,i] = 0
      fu_table[6,i] = 0
      i <- i+2
    }
  }
  
  # add totals
  for (i in 1:nrow(fu_table)){
    fu_table[i,ncol(fu_table)] = sum(fu_table[i,], na.rm=TRUE)
  }
  
  # add percentages of targets
  neuro_due <- filter(df, timepoint=='baseline',
                      dkefs_complete==1,
                      as.Date(dkefs_date_admin) %m+% months(24) < get_month()) %>%
    count(site)
  blood_due <- filter(df, timepoint=='baseline',
                      plasma_blood_status==1,
                      as.Date(plasma_blood_date) %m+% months(24) < get_month()) %>%
    count(site)
  mri_due   <- filter(df, timepoint=='baseline',
                      mr_t1==1,
                      as.Date(mr_date) %m+% months(24) < get_month()) %>%
    count(site)
  
  i <- 2
  for (s in sites){
    fu_table[2,i] = percent(fu_table[2,i-1]/(neuro_due %>% filter(site==s))$n)
    fu_table[4,i] = percent(fu_table[4,i-1]/(blood_due %>% filter(site==s))$n)
    fu_table[6,i] = percent(fu_table[6,i-1]/(mri_due %>% filter(site==s))$n)
    i <- i+2
  }
  
  return(fu_table)
  
}


