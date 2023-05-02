################################################################################
###Pulling together a cohort for analysis#######################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#First use command prompt to log in to Slade using SSH

#To install package
#library(devtools)
#install_github("drkgyoung/Exeter_Diabetes_aurum_package")


#Load tidyverse and lubridate
library(tidyverse)
library(lubridate)

#Load aurum package
library(aurum)
library(EHRBiomarkr)

################################################################################
###SETUP########################################################################

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Setting up/loading analysis test
analysis = cprd$analysis("Rhian_covid")

################################################################################
#Setting cohort information
#Rest of script should run using these inputs and output cohort table
index.date <- as.Date("2020-02-01")
cohort.name <- "feb2020" #this should be pasted into all relevant table names when caching
################################################################################

##Load table already generated (in t1t2_setup) or copy from another analysis (how to do this in MySQL is on the github, alternatively connect to other analysis, then reconnect and cache)
#Includes patid, diabetes diagnosis date, diabetes diagnosis flag, diabetes diagnosis age, dob, diabetes type
all_t1t2 <- all_t1t2 %>% analysis$cached("all_t1t2_cohort", unique_indexes="patid",indexes="dm_diag_date")

#ONS death dates
ons_deathdate <- cprd$tables$onsDeath %>%
  select(patid, ons_ddate = dod) %>%
  analysis$cached("all_ons_ddates", unique_indexes = "patid", indexes= "ons_ddate")

#Getting a cohort of people actively registered and alive on index date
active_cohort <- all_t1t2 %>% 
  #Don't need to drop those diagnoses after index?
  #filter(dm_diag_date<=index.date) %>% 
  #Categorise diagnosis age
  mutate(diag_age_cat = ifelse(dm_diag_age<18, "0-17", ifelse(dm_diag_age<30 & dm_diag_age>=18, "18-29", ifelse(dm_diag_age<50 & dm_diag_age>=30, "30-49", ifelse(dm_diag_age<65 & dm_diag_age>=50, "50-64", ifelse(dm_diag_age<75 & dm_diag_age>=65, "65-74", ifelse(dm_diag_age>=75, "75+", NA))))))) %>%
  #Join valid date lookup and filter people actively registered
  inner_join(cprd$tables$validDateLookup) %>% 
  filter(gp_ons_end_date>=index.date) %>% 
  inner_join(cprd$tables$patient) %>% 
  filter(regstartdate<=index.date) %>% 
  select(patid, pracid, gender, dob, dm_diag_date, dm_diag_flag, dm_diag_age, diag_age_cat, diabetes_type, regstartdate, regenddate, cprd_ddate) %>%
  #Join ONS death dates, take ONS death date as death date, if missing use CPRD death date if have one
  left_join (ons_deathdate, by = "patid") %>%
  mutate(dod = ifelse(!is.na(ons_ddate), ons_ddate, cprd_ddate)) %>%
  #Work out diabetes duration and age at the index date and categorise
  mutate(dm_duration_at_index=datediff(index.date,dm_diag_date)/365.25,age_at_index=datediff(index.date,dob)/365.25) %>% 
  mutate(age_cat = ifelse(age_at_index<40, "<40", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index< 90 & age_at_index>=80, "80-89", ifelse(age_at_index >= 90, "90+", NA)))))))) %>%
  mutate(duration_cat = ifelse(dm_duration_at_index<1, "<1", ifelse(dm_duration_at_index<3 & dm_duration_at_index>=1, "1-2", ifelse(dm_duration_at_index<6 & dm_duration_at_index>=3, "3-5", ifelse(dm_duration_at_index<10 & dm_duration_at_index>=6, "6-9", ifelse(dm_duration_at_index<15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index<20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index>=20, "20+", NA)))))))) %>%
  analysis$cached(paste0(cohort.name,"_cohort"), unique_indexes="patid")

active_cohort %>% count()
active_cohort_patids <- active_cohort %>% select(patid) %>% analysis$cached(paste0(cohort.name,"_cohort_patids"), unique_indexes="patid")

################################################################################
###SOCIOECONOMIC FACTORS########################################################

###Deprivation##################################################################
deprivation <- cprd$tables$patientImd2015 %>%
  mutate(imd_quintile = ifelse(imd2015_10 == 1 | imd2015_10 == 2, "1", ifelse(imd2015_10 == 3 | imd2015_10 == 4, "2", ifelse(imd2015_10 == 5 | imd2015_10 == 6, "3", ifelse(imd2015_10 == 7 | imd2015_10 == 8, "4", ifelse(imd2015_10 == 9 | imd2015_10 == 10, "5", "Missing")))))) %>%
  select(patid, imd_quintile) %>%
  analysis$cached("all_deprivation", unique_indexes= "patid")

###Ethnicity####################################################################

#Cache raw ethnicity table
ethnicity_raw <- cprd$tables$observation %>%
  inner_join(codes$ethnicity_5cat, by = "medcodeid") %>%
  left_join (codes$ethnicity_16cat, by = "medcodeid") %>%
  analysis$cached("all_ethnicity_raw", indexes= c("patid", "ethnicity_5cat_cat", "ethnicity_16cat_cat"))

#Drop ethnicity unknown categories (eth16 =17/ eth5 =5)
ethnicity_all <- ethnicity_raw %>% select(patid, obsdate, medcodeid, eth16= ethnicity_16cat_cat, eth5 = ethnicity_5cat_cat) %>% filter(eth5 != 5) %>% distinct() 

#Find most commonly recorded ethnicity for each person with an ethnicity observation
eth_most_common <- ethnicity_all %>% group_by(patid, eth5, eth16) %>% summarise(n = n()) %>% ungroup() %>% group_by(patid) %>% filter(n == max(n, na.rm =TRUE)) %>%
  select(patid, eth16,eth5) %>% distinct() %>% filter(n()==1) %>% ungroup(patid)

#For those with more than one equally most commonly recorded ethnicities, find the most recently recorded from these categories
eth_most_recent <- ethnicity_all %>% anti_join(eth_most_common, by = "patid") %>% group_by(patid, eth5, eth16) %>% summarise(n = n()) %>% ungroup() %>%
  group_by(patid) %>% filter(n == max(n, na.rm =TRUE)) %>% left_join(ethnicity_all) %>% filter(obsdate == max(obsdate)) %>% select(patid,eth16, eth5) %>% filter(n()==1) %>% ungroup(patid)

#Pull together
ethnicity <- eth_most_common %>% union(eth_most_recent) %>% analysis$cached("all_ethnicity", unique_indexes= "patid")


##Adding in HES ethnicity
#Categorising: 1 = White - > 1 = White, 5 = Indian + 6 = Pakistani + 7 = Bangladeshi + 8 = Other Asian -> 1 = South Asian,
#2 = Black Caribbean + 3 = Black African + 4 = Black Other -> 2 = Black, 9 = Chinese + 11 = Other -> 3 = Other, 10 = Mixed -> 4 = Mixed
hes_ethnicity <- cprd$tables$hesPatient %>% select(patid, gen_ethnicity) %>% distinct() %>% 
  mutate(hes_eth5_cat = ifelse(gen_ethnicity == 1, "0", ifelse(gen_ethnicity == 2 | gen_ethnicity == 3 | gen_ethnicity == 4, "2", ifelse(gen_ethnicity == 5 | gen_ethnicity == 6 | gen_ethnicity == 7 | gen_ethnicity == 8, "1", ifelse(gen_ethnicity == 9 | gen_ethnicity == 11, "3", ifelse(gen_ethnicity == 10, "4", NA )))))) %>%
  select(patid, hes_eth5_cat) %>% analysis$cached("all_ethnicity_hes_categorised", unique_indexes= "patid")


###Smoking######################################################################
smoking_raw <- cprd$tables$observation %>% inner_join(codes$smoking) %>% analysis$cached("all_smoking_raw", indexes= c("patid", "obsdate", "smoking_cat"))

#Get all valid smoking codes and flag is ever smoked
smoking_all <- active_cohort_patids %>% inner_join(smoking_raw) %>% filter(!is.na(obsdate) & obsdate < as.Date(index.date)) %>% select(patid, obsdate, smoking_cat) %>% distinct()
smoked_ever <- smoking_all %>% filter(smoking_cat == "Active smoker") %>% distinct(patid) %>% mutate(smoked_ever_flag = 1)

#Take most recently recorded code category, if this is non-smoker but have previous codes of active smoking then code as ex-smoker
smoking_recent <- smoking_all %>% group_by(patid) %>% filter(obsdate == max(obsdate)) %>% ungroup() %>% left_join(smoked_ever, by = "patid") %>% mutate(smoked_ever_flag = ifelse(is.na(smoked_ever_flag), 0, smoked_ever_flag)) %>%
  mutate (smoking_cat = ifelse((smoking_cat == "Non-smoker" & smoked_ever_flag == 1), "Ex-smoker", smoking_cat)) %>% select(patid, smoking_status = smoking_cat) %>% distinct() %>%
  group_by(patid) %>% filter(n()==1) %>% ungroup()

#For those with two different category codes on most recent day, look at next date before this, if still multiple categories then will be coded as unknown
smoking_next_recent <- smoking_all %>% anti_join(smoking_recent, by = "patid") %>% group_by(patid) %>% filter(obsdate != max(obsdate)) %>% filter(obsdate == max(obsdate)) %>% ungroup() %>% 
  left_join(smoked_ever, by = "patid") %>% mutate(smoked_ever_flag = ifelse(is.na(smoked_ever_flag), 0, smoked_ever_flag)) %>% mutate (smoking_cat = ifelse(smoking_cat == "Non-smoker" & smoked_ever_flag == 1, "Ex-smoker", smoking_cat)) %>% 
  select(patid, smoking_status = smoking_cat) %>% distinct() %>% group_by(patid) %>% filter(n()==1) %>% ungroup()

#Pull together
smoking <- smoking_recent %>% union(smoking_next_recent) %>% analysis$cached(paste0(cohort.name,"_smoking"), unique_indexes = "patid")

###Alcohol consumption##########################################################
alcohol_raw <- cprd$tables$observation %>% inner_join(codes$alcohol) %>% analysis$cached("all_alcohol_raw", indexes= c("patid", "obsdate", "alcohol_cat"))

#Get valid alcohol codes and flag is ever had a consumption level 3 code recorded 
alcohol_all <- active_cohort_patids %>% inner_join(alcohol_raw) %>% filter(!is.na(obsdate) & obsdate < as.Date(index.date)) %>% select(patid, obsdate, alcohol_cat) %>% distinct()
level3_ever <- alcohol_all %>% filter(alcohol_cat == "AlcoholConsumptionLevel3") %>% distinct(patid) %>% mutate(level3_ever_flag = 1)

#Take most recently recorded code category, if ever category 3 then code as 3
alcohol <- alcohol_all %>% group_by(patid) %>% filter(obsdate == max(obsdate)) %>% ungroup() %>% left_join(level3_ever, by = "patid") %>% mutate(level3_ever_flag = ifelse(is.na(level3_ever_flag), 0, level3_ever_flag)) %>%
  mutate(alcohol_consumption = ifelse(level3_ever_flag == 1, "AlcoholConsumptionLevel3", alcohol_cat)) %>% 
  mutate(level_number = ifelse(alcohol_consumption == "AlcoholConsumptionLevel0", 0, ifelse(alcohol_consumption == "AlcoholConsumptionLevel1", 1, ifelse(alcohol_consumption == "AlcoholConsumptionLevel2", 2, ifelse(alcohol_consumption == "AlcoholConsumptionLevel3", 3, NA))))) %>%
  group_by(patid) %>% filter(level_number == max(level_number)) %>% ungroup() %>% select(patid, alcohol_consumption) %>% distinct() %>%
  analysis$cached(paste0(cohort.name,"_alcohol"), unique_indexes = "patid")


###Region#######################################################################
region <- active_cohort %>% left_join(cprd$tables$practice) %>% left_join(cprd$tables$region, by = c("region" = "regionid")) %>% select(patid, region= description) %>%
  analysis$cached(paste0(cohort.name,"_region"), unique_indexes = "patid")

################################################################################
###BIOMARKERS###################################################################

starting_biomarkers <- c("bmi", "height", "weight", "creatinine_blood", "sbp", "totalcholesterol", "acr")

##Join codelist to obs table and cache raw biomarker tables
for (i in starting_biomarkers) {
  print(i)
  raw_tablename <- paste0("all_",i,"_raw")
  tabledata <- cprd$tables$observation %>% inner_join(codes[[i]]) %>% analysis$cached(raw_tablename,indexes=c("patid","obsdate", "testvalue","numunitid"))
  assign(raw_tablename,tabledata)
  print(tabledata %>% count())
}

###HbA1c
all_hba1c_raw <- cprd$tables$observation %>% inner_join(codes$hba1c) %>% analysis$cached("all_hba1c_raw",indexes=c("patid","obsdate","testvalue","numunitid"))
##Take raw table and clean dates, values, units, and take mean value if more than one observation on same date
#Note: there is already an all_hba1c_cleaned table generated in setup so this should just be called
#If there isn't already a cleaned HbA1c table, this will need to be done separately/ add hba1c unit conversion beforehand as below:
all_hba1c_cleaned <- all_hba1c_raw %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>%
  filter(year(obsdate)>=1990) %>% mutate(testvalue=ifelse(testvalue<=20,((testvalue-2.152)/0.09148),testvalue)) %>%  
  clean_biomarker_units(numunitid, "hba1c") %>% clean_biomarker_values(testvalue, "hba1c") %>% group_by(patid,obsdate) %>% summarise(testvalue=mean(testvalue)) %>% ungroup() %>%
  select(patid, obsdate, testvalue) %>% distinct() %>% analysis$cached("all_hba1c_cleaned",indexes=c("patid","obsdate","testvalue"))

for (i in starting_biomarkers) {
  print(i)
  raw_tablename <- paste0("all_",i,"_raw")
  clean_tablename <- paste0("all_",i,"_cleaned")
  tabledata <- get(raw_tablename) %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% filter(year(obsdate)>=1990) %>% 
  clean_biomarker_units(numunitid, i) %>% clean_biomarker_values(testvalue, i) %>% group_by(patid,obsdate) %>% summarise(testvalue=mean(testvalue)) %>% ungroup() %>%
  select(patid, obsdate, testvalue) %>% distinct() %>% analysis$cached(clean_tablename,indexes=c("patid","obsdate","testvalue"))
  assign(clean_tablename,tabledata)
  print(tabledata %>% count())
}

#Calculating BMI from height and weight
#Drop heights/ weights under age 18
adult_height <- all_height_cleaned %>% inner_join(all_t1t2) %>% mutate(obs_age=datediff(obsdate,dob)/365.25) %>% filter(obs_age >=18) %>% select(patid, obsdate, testvalue) %>% analysis$cached("adult_heights",indexes=c("patid","obsdate","testvalue"))
adult_weight <- all_weight_cleaned %>% inner_join(all_t1t2) %>% mutate(obs_age=datediff(obsdate,dob)/365.25) %>% filter(obs_age >=18) %>% select(patid, obsdate, testvalue) %>% analysis$cached("adult_weights",indexes=c("patid","obsdate","testvalue"))

#Calculating median height (no median SQL median function so have to work out the long way, believe below should match what comes out of R function normally)
#May give error, but creates table as should anyway?
median_height <- adult_height %>% arrange(patid,testvalue) %>% group_by(patid) %>% add_count() %>% mutate(medianrow_down = floor((n+1)/2), medianrow_up = ceiling((n+1)/2), rownumber = row_number()) %>% 
  filter(rownumber == medianrow_down | rownumber == medianrow_up) %>% select(patid, testvalue) %>% summarise(median_height = mean(testvalue)) %>% analysis$cached("adult_median_height", unique_indexes = "patid", indexes= "median_height")

#BMI calculation
bmi_calc <- adult_weight %>% inner_join(median_height, by = "patid") %>% mutate (testvalue = testvalue / ((median_height/100)^2)) %>% clean_biomarker_values(testvalue, "bmi") %>% 
  select(patid, obsdate, testvalue) %>% analysis$cached("all_bmi_calc",indexes=c("patid","obsdate","testvalue"))
#Combine with BMI values from code, where no coded value on obsdate take calculated value
all_bmi_cleaned <- all_bmi_cleaned %>% rename(bmi_code_value = testvalue) %>% union(bmi_calc) %>% rename(bmi_calc_value = testvalue) %>% group_by(patid,obsdate) %>%
  summarise(bmi_code_value = max(bmi_code_value), bmi_calc_value = max(bmi_calc_value)) %>% mutate(bmi_value = ifelse(is.na(bmi_code_value), bmi_calc_value, bmi_code_value)) %>%
  ungroup() %>% select(patid, obsdate, testvalue = bmi_value) %>% analysis$cached("all_bmi_combined", indexes=c("patid","obsdate","testvalue"))

#eGFR calculation
all_egfr_cleaned <- all_creatinine_blood_cleaned %>% left_join(all_t1t2) %>% mutate(age_egfr = datediff(obsdate,dob)/365.25) %>%
  mutate (creatinine_mgdl = testvalue * 0.0113) %>% left_join(cprd$tables$patient) %>% select(patid, obsdate, gender, age_egfr, creatinine_mgdl) %>%
  mutate(egfr_value= ifelse((creatinine_mgdl<=0.7 & gender ==2),(142*((creatinine_mgdl/0.7)^-0.241)*(0.9938^age_egfr) *1.012),
  ifelse((creatinine_mgdl>0.7 & gender==2),(142*((creatinine_mgdl/0.7)^-1.2)*(0.9938^age_egfr)*1.012),
  ifelse((creatinine_mgdl<=0.9 & gender==1),(142*((creatinine_mgdl/0.9)^-0.302)*(0.9938^age_egfr)),
  ifelse((creatinine_mgdl>0.9 & gender==1),(142*((creatinine_mgdl/0.9)^-1.2)*(0.9938^age_egfr)), NA))))) %>% filter(!is.na(egfr_value)) %>%
  dplyr::select(patid, obsdate, testvalue = egfr_value) %>% distinct() %>% analysis$cached("all_egfr_cleaned", indexes=c("patid","obsdate","testvalue"))
#will want to use this for CKD


##Adding biomarker values to cohort
final_biomarkers <- c("hba1c", "bmi", "egfr", "sbp", "totalcholesterol")

#Set biomarkers to cohort patids and create date 2 years than index date
biomarkers <- active_cohort_patids
two.years.earlier <- index.date - years(2)
#Find biomarker measurement closest to index date in the 2 years prior to and including index date
for (i in final_biomarkers) {
  cleaned_tablename <- paste0("all_",i,"_cleaned")
  new_date_variablename <- paste0(i,"_measure_date")
  new_value_variablename <- paste0(i,"_value")
  nearest_biomarker <- get(cleaned_tablename) %>% filter(obsdate <= index.date & obsdate >= two.years.earlier) %>% group_by(patid) %>% filter(obsdate == max(obsdate,na.rm=TRUE)) %>% ungroup() %>% rename({{new_date_variablename}}:= obsdate, {{new_value_variablename}}:= testvalue)
  biomarkers <- biomarkers %>% left_join(nearest_biomarker)
}
biomarkers <- biomarkers %>% analysis$cached(paste0(cohort.name,"_biomarkers"),unique_indexes="patid")


################################################################################
###COMORBIDITIES AND COMPLICATIONS##############################################

conditions <- c("heartfailure", "myocardialinfarction", "stroke", "retinopathy", "neuropathy", "diabeticnephropathy", "af", "angina", "asthma", "cld", "copd", "dementia",
                "hypertension", "ihd", "pad", "revasc", "solidorgantransplant", "tia", "haem_cancer", "solid_cancer", "otherneuroconditions", "pulmonaryfibrosis", "cysticfibrosis",
                "pulmonaryhypertension", "bronchiectasis")

#Join codelist to obs table and cache raw observation tables
for (i in conditions) {
  print(i)
  tablename <- paste0("all_",i,"_raw")
  tabledata <- cprd$tables$observation %>% inner_join(codes[[i]]) %>% analysis$cached(tablename,indexes=c("patid","obsdate"))
  assign(tablename,tabledata)
  print(tabledata %>% count())
}

#
for (i in conditions) {
  print(i)
  raw_tablename <- paste0("all_",i,"_raw")
  date_tablename <- paste0("first_",i,"_gp")
  tabledata <- get(raw_tablename) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% group_by(patid) %>% summarise(date =min(obsdate,na.rm=TRUE)) %>% ungroup() %>% analysis$cached(date_tablename,indexes=c("patid","date"))
  assign(date_tablename,tabledata)
  print(tabledata %>% count())
}

################################################################################
##HES comorbidities

#
for (i in conditions) {
  if(length(codes[[paste0("icd10_",i)]]) > 0) {
  print(i)
  tablename <- paste0("all_",i,"_HES")
  tabledata <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes[[paste0("icd10_",i)]], sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% 
    analysis$cached(tablename,indexes=c("patid","epistart"))
  assign(tablename,tabledata)
  print(tabledata %>% count())
  }
}

#
for (i in conditions) {
  raw_tablename <- paste0("all_",i,"_HES")
  if(exists(raw_tablename)) {
  print(i)
  date_tablename <- paste0("first_",i,"_hes")
  tabledata <- get(raw_tablename) %>% inner_join(cprd$tables$validDateLookup) %>% filter(epistart>=min_dob) %>% group_by(patid) %>% summarise(date =min(epistart,na.rm=TRUE)) %>% ungroup() %>% analysis$cached(date_tablename,indexes=c("patid","date"))
  assign(date_tablename,tabledata)
  print(tabledata %>% count())
  }
}


##Same for procedures tables

#
for (i in conditions) {
  if(length(codes[[paste0("opcs4_",i)]]) > 0) {
  print(i)
  tablename <- paste0("all_",i,"_op")
  tabledata <- cprd$tables$hesProceduresEpi %>% inner_join(codes[[paste0("opcs4_",i)]],by=c("OPCS"="opcs4")) %>% 
    analysis$cached(tablename,indexes=c("patid","evdate"))
  assign(tablename,tabledata)
  print(tabledata %>% count())
  }
}

#
for (i in conditions) {
  raw_tablename <- paste0("all_",i,"_op")
  if(exists(raw_tablename)) {
  print(i)
  date_tablename <- paste0("first_",i,"_op")
  tabledata <- get(raw_tablename) %>% inner_join(cprd$tables$validDateLookup) %>% filter(evdate>=min_dob) %>% group_by(patid) %>% summarise(date =min(evdate,na.rm=TRUE)) %>% ungroup() %>% analysis$cached(date_tablename,indexes=c("patid","date"))
  assign(date_tablename,tabledata)
  print(tabledata %>% count())
  }
}

################################################################################
###Combining
for(i in conditions) {
  print(i)
  gp_tablename <- paste0("first_",i,"_gp")
  hes_tablename <- paste0("first_",i,"_hes")
  op_tablename <- paste0("first_",i,"_op")
  combined_tablename <- paste0("first_",i, "_all")
  tabledata <- get(gp_tablename)
  if(exists(hes_tablename)) {
    hes_table <- get(hes_tablename) 
    tabledata <- tabledata %>% union(hes_table)
  }
  if(exists(op_tablename)) {
    hes_table <- get(op_tablename) 
    tabledata <- tabledata %>% union(hes_table)
  }
  tabledata <- tabledata %>% group_by(patid) %>% summarise(date = min(date, na.rm =TRUE)) %>% ungroup() %>% analysis$cached(combined_tablename,indexes=c("patid","date"))
  assign(combined_tablename,tabledata)
  print(tabledata %>% count())
}

#
cm <- active_cohort_patids

for (i in conditions) {
  date_tablename <- paste0("first_",i, "_all")
  date_variablename <- paste0(i,"_diag_date")
  preindex_varname <- paste0("preindex_",i)
  diagnosis_date <- get(date_tablename)
  cm <- cm %>% left_join(diagnosis_date) %>% mutate(preindex = ifelse(!is.na(date) & date <= index.date, 1L, 0L)) %>% rename({{date_variablename}}:= date, {{preindex_varname}}:= preindex)
}

cm <- cm %>% analysis$cached(paste0(cohort.name,"_conditions_gp_hes"),unique_indexes="patid")


###CKD##########################################################################
##Identifying CKD5 by medcodes
#Cache raw ckd5 observations
ckd5_raw <- cprd$tables$observation %>% inner_join(codes$ckd5) %>% analysis$cached("all_ckd5_raw", indexes=c("patid","obsdate"))
ckd5_hes_raw <- cprd$tables$hesDiagnosisEpi %>% inner_join(codes$icd10_ckd5, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_ckd5_HES", indexes=c("patid","epistart"))
ckd5_op_raw <- cprd$tables$hesProceduresEpi %>% inner_join(codes$opcs4_ckd5,by=c("OPCS"="opcs4")) %>% analysis$cached("all_ckd5_op", indexes=c("patid","evdate"))

first_ckd5_gp <- active_cohort_patids %>% inner_join(ckd5_raw) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% filter(obsdate <= index.date) %>%
  group_by(patid) %>% summarise(date=min(obsdate,na.rm=TRUE))%>% ungroup()
first_ckd5_hes <- active_cohort_patids %>% inner_join(ckd5_hes_raw) %>% inner_join(cprd$tables$validDateLookup) %>% filter(epistart>=min_dob & epistart<=gp_ons_end_date) %>% filter(epistart <= index.date) %>%
  group_by(patid) %>% summarise(date=min(epistart,na.rm=TRUE))%>% ungroup()
first_ckd5_op <- active_cohort_patids %>% inner_join(ckd5_op_raw) %>% inner_join(cprd$tables$validDateLookup) %>% filter(evdate>=min_dob & evdate<=gp_ons_end_date) %>% filter(evdate <= index.date) %>%
  group_by(patid) %>% summarise(date=min(evdate,na.rm=TRUE))%>% ungroup()

#Find earliest date of ckd5 observations
earliest_ckd5 <- first_ckd5_gp %>% union(first_ckd5_hes) %>% union(first_ckd5_op) %>% group_by(patid) %>% summarise(ckd5_date=min(date,na.rm=TRUE)) %>% 
  ungroup() %>% select(patid,ckd5_date) %>% analysis$cached(paste0(cohort.name, "_earliest_ckd5"),indexes=c("patid","ckd5_date"))

#Working out CKD stage using eGFR (previously calculated using creatinine)
# Convert eGFR to CKD stage
clean_ckd <- all_egfr_cleaned %>% rename(egfr = testvalue, date = obsdate) %>%
  mutate(ckd_stage_coded= ifelse(egfr<15, 150L, ifelse(egfr<30, 140L, ifelse(egfr<45, 135L, ifelse(egfr<60, 130L, ifelse(egfr<90, 120L, ifelse(egfr>=90, 110L, NA))))))) %>%
  analysis$cached("all_clean_ckd_stages",indexes=c("patid","date","ckd_stage_coded"))

### A) Define period from current test until next test as having the ckd_stage of current test
#### Add in row labelling within each patient's values + max number of rows for each patient
ckd <- clean_ckd %>% group_by(patid) %>% dbplyr::window_order(date) %>% mutate(patid_row_id=row_number()) %>% mutate(patid_total_rows=max(patid_row_id,na.rm=TRUE)) %>% ungroup()

#### For rows where there is a next test, use this as end date; for last row, use start date as end date
ckd <- ckd %>% mutate(next_row=patid_row_id+1) %>% left_join(ckd,by=c("patid","next_row"="patid_row_id")) %>% 
  mutate(ckd_start=date.x, ckd_end=if_else(is.na(date.y),date.x,date.y), ckd_stage_coded=ckd_stage_coded.x) %>% select(patid,patid_row_id,ckd_stage_coded,ckd_start,ckd_end)

### B) Join together consecutive periods with the same ckd_stage
ckd_combined <- ckd %>% group_by(patid,ckd_stage_coded) %>% dbplyr::window_order(patid,ckd_stage_coded,patid_row_id) %>% mutate(lead_var=lead(ckd_start), cummax_var=cummax(ckd_end)) %>%
  mutate(compare=cumsum(lead_var>cummax_var)) %>% mutate(indx=ifelse(row_number()==1,0.0,lag(compare))) %>% ungroup() %>% group_by(patid,ckd_stage_coded,indx) %>%
  summarise(first_test_date=min(ckd_start,na.rm=TRUE), last_test_date=max(ckd_start,na.rm=TRUE), end_date=max(ckd_end,na.rm=TRUE), test_count=max(patid_row_id,na.rm=TRUE)-min(patid_row_id,na.rm=TRUE)+1) %>%
  ungroup() %>% analysis$cached("all_ckd_combined",indexes=c("test_count","first_test_date","last_test_date"))

ckd_combined %>% count()
ckd_combined %>% summarise(total=sum(test_count,na.rm=TRUE))

### C) Remove periods with 1 reading, or with multiple readings but <90 days between first and last test, and cache
ckd_combined_filtered <- ckd_combined %>% filter(test_count>1 & datediff(last_test_date,first_test_date)>=90) %>%
  analysis$cached("all_ckd_filtered",indexes=c("patid","ckd_stage_coded","first_test_date"))

### D) Just keep max stage per person and find date of onset
ckd_from_creatinine <- active_cohort_patids %>% inner_join(ckd_combined) %>% filter(first_test_date <= index.date) %>% group_by(patid) %>% summarise(ckd_stage=max(ckd_stage_coded,na.rm=TRUE)) %>%
  inner_join(ckd_combined_filtered,by=c("patid","ckd_stage"="ckd_stage_coded")) %>% group_by(patid,ckd_stage) %>%
  summarise(ckd_stage_onset=min(first_test_date,na.rm=TRUE)) %>% ungroup() %>%
  analysis$cached(paste0(cohort.name, "_ckd_creat"),indexes=c("patid","ckd_stage","ckd_stage_onset"))

##Combine ckd stage worked out from creatinine and ckd5 observations, if ckd5 medcode coded prior to index date then ckd5, otherwise take stage from creatinine
ckd_at_index <- ckd_from_creatinine %>% left_join(earliest_ckd5) %>% mutate(ckd_stage=ifelse(!is.na(ckd5_date),150L,ckd_stage), ckd_stage_onset=ifelse(!is.na(ckd5_date),ckd5_date,ckd_stage_onset)) %>%
  select(patid,ckd_stage,ckd_stage_onset)
extra_ckd5_only_people <- earliest_ckd5 %>% left_join(ckd_from_creatinine) %>% filter(is.na(ckd_stage)) %>% mutate(ckd_stage=150L) %>%
  select(patid,ckd_stage,ckd_stage_onset=ckd5_date)

ckd_at_index <- ckd_at_index %>% union(extra_ckd5_only_people) %>% select(patid, ckd_stage, ckd_stage_onset)

#For those with missing CKD stage from algorithm, use stage of most recent eGFR measurement in the last 2 years
ckd_recent_egfr <- biomarkers %>% select(patid, egfr_value) %>% mutate(recent_egfr_stage= ifelse(egfr_value<15, 150L, ifelse(egfr_value<30, 140L, ifelse(egfr_value<45, 135L, ifelse(egfr_value<60, 130L, ifelse(egfr_value<90, 120L, ifelse(egfr_value>=90, 110L, NA))))))) %>%
  select(patid, recent_egfr_stage)

ckd <- active_cohort_patids %>% left_join(ckd_at_index) %>% left_join(ckd_recent_egfr) %>% mutate(ckd_stage = ifelse(is.na(ckd_stage), recent_egfr_stage, ckd_stage)) %>% select(patid, ckd_stage, ckd_stage_onset) %>%
  mutate(ckd_stage = ifelse(ckd_stage == 110, "Stage 1", ifelse(ckd_stage == 120, "Stage 2", ifelse(ckd_stage == 130, "Stage 3a", ifelse(ckd_stage == 135, "Stage 3b", ifelse(ckd_stage == 140, "Stage 4", ifelse(ckd_stage == 150, "Stage 5", NA))))))) %>%
  analysis$cached(paste0(cohort.name, "_ckd_stages"),unique_indexes="patid",indexes=c("ckd_stage","ckd_stage_onset"))

###ACR##########################################################################

#Get raw urine albumin and creatinine and clean by removing missing test values, filtering to only most common/ correct units and cleaning dates
all_albumin_raw <- cprd$tables$observation %>% inner_join(codes$albumin_urine) %>% analysis$cached("all_albumin_raw", indexes=c("patid","obsdate", "testvalue", "numunitid"))
all_urinecreatinine_raw <- cprd$tables$observation %>% inner_join(codes$creatinine_urine) %>% analysis$cached("all_urine_creatinine_raw", indexes=c("patid","obsdate", "testvalue", "numunitid"))

albumin_cleaned <- all_albumin_raw %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% filter(year(obsdate)>=1990) %>% 
  filter(!is.na(testvalue)) %>% filter(numunitid == 183) %>% group_by(patid,obsdate) %>% summarise(testvalue=mean(testvalue)) %>% ungroup() %>% select(patid, obsdate, testvalue) %>% distinct()
urinecreatinine_cleaned <- all_urinecreatinine_raw %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% filter(year(obsdate)>=1990) %>% 
  filter(!is.na(testvalue)) %>% filter(numunitid == 218) %>% group_by(patid,obsdate) %>% summarise(testvalue=mean(testvalue)) %>% ungroup() %>% select(patid, obsdate, testvalue) %>% distinct()

#Combine albumin and creatinine and calculate acr, clean values
acr_calc <- albumin_cleaned %>% inner_join(urinecreatinine_cleaned, by = c("patid", "obsdate")) %>% mutate(testvalue = testvalue.x/testvalue.y) %>% select(patid, obsdate, testvalue) %>%
  clean_biomarker_values(testvalue, "acr") %>% rename(acr_calc_value = testvalue) %>% analysis$cached("all_acr_calc",indexes=c("patid","obsdate","acr_calc_value"))

#Combine coded and calculated acr values, take coded value unless missing, then take calculated value
acr_with_calc <- all_acr_cleaned %>% rename(acr_code_value = testvalue) %>% union(acr_calc) %>% group_by(patid,obsdate) %>% summarise(acr_code_value = max(acr_code_value), acr_calc_value = max(acr_calc_value)) %>%
  mutate(acr_value = ifelse(is.na(acr_code_value), acr_calc_value, acr_code_value)) %>% ungroup() %>% select(patid, acr_measure_date = obsdate, acr_value) %>% distinct() %>%
  analysis$cached("all_acr_combined", indexes=c("patid","acr_measure_date","acr_value"))

#Find most recent acr value
acr <- active_cohort_patids %>% inner_join(acr_with_calc) %>% filter(acr_measure_date <= index.date) %>% group_by(patid) %>% filter(acr_measure_date == max(acr_measure_date,na.rm=TRUE)) %>%
  analysis$cached(paste0(cohort.name, "_acr"), unique_indexes = "patid")

###RECENT PREVIOUS RESPIRATORY INFECTION HOSPITALISATION########################
two.years.earlier <- index.date - years(2)

resp_infection_hosps <- cprd$tables$hesDiagnosisHosp %>% inner_join(codes$icd10_respiratoryinfection, sql_on="LHS.ICD LIKE CONCAT(icd10,'%')") %>% analysis$cached("all_resp_infect_hosp", indexes = c("patid", "admidate"))

recent_resp_infect <- resp_infection_hosps %>% filter(admidate < index.date & admidate >= two.years.earlier) %>% filter(admidate != discharged) %>%
  left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>% group_by(patid) %>% 
  summarise(resp_infect_hosp_date = max(admidate)) %>% ungroup() %>% mutate(recent_hosp_resp_infect =1L) %>% analysis$cached(paste0(cohort.name,"_recent_resp_infect"), unique_indexes = "patid")

###RECENT PREVIOUS HOSPITALISATION FOR ANYTHING#################################
#Find hospitalisations for anything (excluding same day discharge and hospitalisations for respiratory infections) in last 2 years
recent_hosp_anything <- cprd$tables$hesDiagnosisHosp %>% anti_join(resp_infection_hosps %>% distinct(spno)) %>% filter(admidate < index.date & admidate >= two.years.earlier) %>% 
  filter(admidate != discharged) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>% 
  group_by(patid) %>% summarise(recent_hosp_date = max(admidate)) %>% ungroup() %>% mutate(recent_hosp_anything =1L) %>%
  analysis$cached(paste0(cohort.name,"_recent_hosp_anything"), unique_indexes = "patid")

################################################################################
###DIABETES MEDICATIONS#########################################################

six.months.earlier <- index.date - months(6)
##Insulin
#If ran setup then should just call this table
insulin_cleaned <- cprd$tables$drugIssue %>% inner_join(codes$insulin) %>% inner_join(cprd$tables$validDateLookup) %>% filter(!is.na(issuedate)) %>% 
  filter(issuedate>=min_dob & issuedate<=gp_ons_end_date) %>% select(patid,issuedate,dosageid,quantity,quantunitid,duration) %>% 
  analysis$cached("all_insulin_cleaned",indexes=c("patid","issuedate"))
#Find whether had any prescription for insulin in 6 months prior to and including index date
insulin_6months <- active_cohort_patids %>% left_join(insulin_cleaned) %>% filter(issuedate >= six.months.earlier & issuedate<=index.date) %>% select(patid) %>% 
  mutate(insulin_6months =1L) %>% distinct() %>% analysis$cached(paste0(cohort.name,"_insulin_6m"),unique_indexes="patid")

##OHAs
#Should call table from setup
ohas_cleaned <- cprd$tables$drugIssue %>% inner_join(cprd$tables$ohaLookup) %>% inner_join(cprd$tables$validDateLookup) %>% filter(!is.na(issuedate)) %>% 
  filter(issuedate>=min_dob & issuedate<=gp_ons_end_date) %>% select(patid,issuedate,dosageid,quantity,quantunitid,duration,INS,TZD,SU,DPP4,MFN,GLP1,Glinide,Acarbose,SGLT2) %>% 
  analysis$cached("all_ohas_cleaned",indexes=c("patid","issuedate","INS","TZD","SU","DPP4","MFN","GLP1","Glinide","Acarbose","SGLT2"))
ohas_6months <- active_cohort_patids %>% left_join(ohas_cleaned) %>% filter(issuedate >= six.months.earlier & issuedate<=index.date) %>% group_by(patid) %>% 
  summarise(INS = max(INS), TZD = max(TZD), SU = max(SU), DPP4 = max(DPP4), MFN = max(MFN), GLP1 = max(GLP1), Glinide = max(Glinide), Acarbose = max(Acarbose), SGLT2 = max(SGLT2)) %>% ungroup() %>% 
  mutate(OHA_6m =1L) %>% analysis$cached(paste0(cohort.name,"_ohas_6m"), unique_indexes= "patid")

#Combine diabetes medications
diabetes_medications <- active_cohort_patids %>% left_join(ohas_6months) %>% left_join(insulin_6months) %>% mutate(INS_6m = ifelse(insulin_6months ==1 | INS ==1, 1L, 0L), TZD_6m = TZD, SU_6m = SU, DPP4_6m = DPP4, MFN_6m = MFN, GLP1_6m = GLP1, Glinide_6m = Glinide, Acarbose_6m = Acarbose, SGLT2_6m = SGLT2) %>% select(patid, OHA_6m, INS_6m:SGLT2_6m) %>%
  mutate(INS_6m = ifelse(is.na(INS_6m), 0L, INS_6m), TZD_6m = ifelse(is.na(TZD_6m), 0L, TZD_6m), SU_6m = ifelse(is.na(SU_6m), 0L, SU_6m), DPP4_6m = ifelse(is.na(DPP4_6m), 0L, DPP4_6m), MFN_6m = ifelse(is.na(MFN_6m), 0L, MFN_6m), GLP1_6m = ifelse(is.na(GLP1_6m), 0L, GLP1_6m), Glinide_6m = ifelse(is.na(Glinide_6m), 0L, Glinide_6m), Acarbose_6m = ifelse(is.na(Acarbose_6m), 0L, Acarbose_6m), SGLT2_6m = ifelse(is.na(SGLT2_6m), 0L, SGLT2_6m), OHA_6m = ifelse(is.na(OHA_6m),0L,OHA_6m)) %>% 
  mutate(treatment_6m = ifelse(INS_6m == 1, "Insulin", ifelse(OHA_6m ==1 & INS_6m ==0, "OHA only", ifelse(INS_6m ==0 & OHA_6m ==0, "No treatment", NA)))) %>%
  analysis$cached(paste0(cohort.name,"_diabetesmeds_6m"), unique_indexes = "patid")

################################################################################
###OTHER MEDICATIONS############################################################
medications <- c("oralsteroids", "immunosuppressants", "labas", "ltras")

#Join codelist to drugs table and cache raw tables
for (i in medications) {
  print(i)
  tablename <- paste0("all_",i,"_raw")
  tabledata <- cprd$tables$drugIssue %>% inner_join(codes[[i]]) %>% analysis$cached(tablename,indexes=c("patid","issuedate"))
  assign(tablename,tabledata)
  print(tabledata %>% count())
}

#Clean issue dates and cache tables
for (i in medications) {
  print(i)
  raw_tablename <- paste0("all_",i,"_raw")
  clean_tablename <- paste0("all_",i,"_cleaned")
  tabledata <- get(raw_tablename) %>% inner_join(cprd$tables$validDateLookup) %>% filter(!is.na(issuedate)) %>% filter(issuedate>=min_dob & issuedate<=gp_ons_end_date) %>% 
    select(patid,issuedate,dosageid,quantity,quantunitid,duration) %>% analysis$cached(clean_tablename,indexes=c("patid", "issuedate"))
  assign(clean_tablename,tabledata)
  print(tabledata %>% count())
}

#Find where patients had medication prescribed in 6 months prior to index date
other_medications <- active_cohort_patids
six.months.earlier <- index.date - months(6)

for (i in medications) {
  cleaned_tablename <- paste0("all_",i,"_cleaned")
  new_variablename <- paste0(i,"_6m")
  last_6months <- get(cleaned_tablename) %>% filter(issuedate >= six.months.earlier & issuedate<=index.date) %>% distinct(patid) %>% mutate(med_6m =1)
  other_medications <- other_medications %>% left_join(last_6months) %>% mutate(med_6m = ifelse(is.na(med_6m), 0L, med_6m)) %>% rename({{new_variablename}}:= med_6m)
}
other_medications <- other_medications %>% analysis$cached(paste0(cohort.name,"_other_medications"),unique_indexes="patid")

################################################################################
###PULL TOGETHER################################################################
###Pull together all variables##################################################
active_cohort_all <- active_cohort %>% 
  left_join(deprivation) %>% mutate(imd_quintile= ifelse(is.na(imd_quintile), "Missing IMD", imd_quintile)) %>%
  left_join(ethnicity) %>% mutate(eth16 = ifelse(is.na(eth16), "17", eth16), eth5 = ifelse(is.na(eth5), "5", eth5)) %>% rename(eth5_gp = eth5, eth16_gp = eth16) %>%
  left_join(hes_ethnicity) %>% mutate(eth5 = ifelse(eth5_gp == "5", hes_eth5_cat, eth5_gp)) %>% select(-hes_eth5_cat) %>% mutate(eth5 = ifelse(is.na(eth5), "5", eth5)) %>%
  left_join(smoking) %>% mutate(smoking_status = ifelse(is.na(smoking_status), "Unknown smoking", smoking_status)) %>%
  left_join(alcohol) %>% mutate(alcohol_consumption = ifelse(is.na(alcohol_consumption), "Unknown alcohol", alcohol_consumption)) %>%
  left_join(region) %>% mutate(region = ifelse(is.na(region), "Missing region", region)) %>%
  left_join(cm) %>% mutate(number_complications = preindex_retinopathy + preindex_neuropathy + preindex_diabeticnephropathy) %>%
  left_join(ckd) %>% mutate(ckd_stage = ifelse(is.na(ckd_stage), "Missing CKD", ckd_stage)) %>% 
  left_join(acr) %>%  mutate(acr_cat = ifelse(acr_value <3, "A1", ifelse(acr_value >=3 & acr_value <= 30, "A2", ifelse(acr_value >30, "A3", NA)))) %>% mutate(acr_cat = ifelse(is.na(acr_cat), "Missing ACR", acr_cat)) %>%
  left_join(recent_resp_infect) %>% mutate(recent_hosp_resp_infect = ifelse(is.na(recent_hosp_resp_infect), 0L, recent_hosp_resp_infect)) %>%
  left_join(recent_hosp_anything) %>% mutate(recent_hosp_anything = ifelse(is.na(recent_hosp_anything), 0L, recent_hosp_anything)) %>%
  left_join(diabetes_medications) %>% 
  left_join(other_medications) %>%
  left_join(biomarkers) %>%
  mutate (hba1c_cat= ifelse(hba1c_value<48, "<48", ifelse(hba1c_value<53 & hba1c_value>=48, "48-53", ifelse(hba1c_value<64 & hba1c_value>=53, "53-64", ifelse(hba1c_value<75 & hba1c_value>=64, "64-75", ifelse(hba1c_value<86 & hba1c_value>=75, "75-86", ifelse(hba1c_value>=86, "86+", NA))))))) %>%
  mutate (bmi_cat= ifelse(bmi_value<18.5, "<18.5", ifelse(bmi_value<25 & bmi_value>=18.5, "18.5-24.9", ifelse(bmi_value<30 & bmi_value>=25, "25-29.9", ifelse(bmi_value<35 & bmi_value>=30, "30-34.9", ifelse(bmi_value<40 & bmi_value>=35, "35-39.9", ifelse(bmi_value>=40, "40+", NA))))))) %>%
  mutate (sbp_cat= ifelse(sbp_value <= 140, "<=140", ifelse(sbp_value >140, ">140", NA))) %>%
  mutate (totalcholesterol_cat= ifelse(totalcholesterol_value <= 5, "<=5", ifelse(totalcholesterol_value >5, ">5", NA))) %>%
  mutate(egfr_cat = ifelse(egfr_value<15, "<15", ifelse(egfr_value<30 & egfr_value>=15, "15-29", ifelse(egfr_value<45 & egfr_value>=30, "30-44", ifelse(egfr_value<60 & egfr_value>=45, "45-59", ifelse(egfr_value<90 & egfr_value>=60, "60-89", ifelse(egfr_value>=90, "90+", NA))))))) %>%
  mutate (hba1c_cat = ifelse(is.na(hba1c_cat), "Missing HbA1c", hba1c_cat), bmi_cat = ifelse(is.na(bmi_cat), "Missing BMI", bmi_cat), sbp_cat = ifelse(is.na(sbp_cat), "Missing SBP", sbp_cat), totalcholesterol_cat = ifelse(is.na(totalcholesterol_cat), "Missing TC", totalcholesterol_cat), egfr_cat = ifelse(is.na(egfr_cat), "Missing eGFR", egfr_cat)) %>%
  analysis$cached(paste0(cohort.name,"_cohort_all"), unique_indexes = "patid")


################################################################################
###END##########################################################################
rm(list=ls())
