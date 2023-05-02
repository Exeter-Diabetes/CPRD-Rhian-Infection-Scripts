################################################################################
###ANALYSIS#####################################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(lubridate)
library(survminer)
library(survival)
library(kableExtra)

#Load aurum package
library(aurum)

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Setting up/loading analysis test
analysis = cprd$analysis("Rhian_covid")


################################################################################
###COVID COHORT#################################################################
################################################################################

###Covid hospitalisations
cohort.name <- "feb2020"
infection <- "covid"
#Set index and end dates
index.date = as.Date("2020-02-01") #change these for different cohorts
end.date = as.Date("2020-10-31")

#Load cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name,"_cohort_all"), unique_indexes = "patid")

###HOSPITALISATION##############################################################
###Outcome: hospitalisation, any diagnoses in HES (excluding elective admissions: admimeth = 11/12/13)
hospitalisations <- cohort %>% select(patid) %>% left_join(cprd$tables$hesDiagnosisHosp) %>% filter(admidate >= index.date & admidate <= end.date) %>%
  inner_join(codes[[paste0("icd10_",infection)]], by=c("ICD"="icd10")) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>%
  mutate(same_day_discharge = ifelse(duration ==0, 1, 0)) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosps"), indexes = c("patid", "admidate"))

hospitalisations <- collect(hospitalisations)

################################################################################
outcome <- "hosp" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = hospitalisation_outcome, outcome_date = hosp_date)


###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)

#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age >=20 & age_at_index >=18)

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

##Table of code frequencies
cohort_outcome_codes <- cohort %>% filter(outcome ==1) %>% select(patid, outcome_date) %>% rename(admidate = outcome_date) %>% left_join(hospitalisations)
tab_covid_code_freq <- cohort_outcome_codes %>% select(patid, admidate, ICD) %>% distinct() %>% group_by(ICD) %>% summarise(freq = n()) %>% arrange(desc(freq))

covid_code_desc <- read_delim("exeter_icd10_covid.txt")
tab_covid_code_freq <- tab_covid_code_freq %>% left_join(covid_code_desc, c("ICD" = "ICD10")) %>% select(ICD10 = ICD, Description = Term_description, Frequency = freq)

kableExtra::kable(tab_covid_code_freq,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="12cm") %>%
  column_spec(3,width="3cm") %>%
  cat(.,file="STab3a_covid_outcome_code_frequencies.html")

################################################################################
###FLU COHORT###################################################################
################################################################################

###Flu hospitalisations
cohort.name <- "sep2016"
infection <- "influenza"
index.date = as.Date("2016-09-01") #change these for different cohorts
end.date = as.Date("2019-05-31")

#Load cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name,"_cohort_all"), unique_indexes = "patid")

###HOSPITALISATION##############################################################
###Outcome: hospitalisation, any diagnoses in HES (excluding elective admissions: admimeth = 11/12/13)
hospitalisations <- cohort %>% select(patid) %>% left_join(cprd$tables$hesDiagnosisHosp) %>% filter(admidate >= index.date & admidate <= end.date) %>%
  inner_join(codes[[paste0("icd10_",infection)]], by=c("ICD"="icd10")) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>%
  mutate(same_day_discharge = ifelse(duration ==0, 1, 0)) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosps"), indexes = c("patid", "admidate"))

hospitalisations <- collect(hospitalisations)

################################################################################
outcome <- "hosp" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = hospitalisation_outcome, outcome_date = hosp_date)

###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)

#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age >=20 & age_at_index >=18)

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

##Table of code frequencies
cohort_outcome_codes <- cohort %>% filter(outcome ==1) %>% select(patid, outcome_date) %>% rename(admidate = outcome_date) %>% left_join(hospitalisations)
tab_flu_code_freq <- cohort_outcome_codes %>% select(patid, admidate, ICD) %>% distinct() %>% group_by(ICD) %>% summarise(freq = n()) %>% arrange(desc(freq))

flu_code_desc <- read_delim("exeter_icd10_influenza.txt")
tab_flu_code_freq <- tab_flu_code_freq %>% left_join(flu_code_desc, c("ICD" = "ICD10")) %>% select(ICD10 = ICD, Description = Term_description, Frequency = freq)

kableExtra::kable(tab_flu_code_freq,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="12cm") %>%
  column_spec(3,width="3cm") %>%
  cat(.,file="STab3b_influenza_outcome_code_frequencies.html")

################################################################################
###PNEUMONIA COHORT#############################################################
################################################################################

###Pneumonia hospitalisations
cohort.name <- "sep2016"
infection <- "pneumonia"
index.date = as.Date("2016-09-01") #change these for different cohorts
end.date = as.Date("2019-05-31")

#Load cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name,"_cohort_all"), unique_indexes = "patid")

###HOSPITALISATION##############################################################
###Outcome: hospitalisation, any diagnoses in HES (excluding elective admissions: admimeth = 11/12/13)
hospitalisations <- cohort %>% select(patid) %>% left_join(cprd$tables$hesDiagnosisHosp) %>% filter(admidate >= index.date & admidate <= end.date) %>%
  inner_join(codes[[paste0("icd10_",infection)]], by=c("ICD"="icd10")) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>%
  mutate(same_day_discharge = ifelse(duration ==0, 1, 0)) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosps"), indexes = c("patid", "admidate"))

hospitalisations <- collect(hospitalisations)

################################################################################
outcome <- "hosp" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = hospitalisation_outcome, outcome_date = hosp_date)

###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)

#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age >=20 & age_at_index >=18)

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

##Table of code frequencies
cohort_outcome_codes <- cohort %>% filter(outcome ==1) %>% select(patid, outcome_date) %>% rename(admidate = outcome_date) %>% left_join(hospitalisations)
tab_pneumo_code_freq <- cohort_outcome_codes %>% select(patid, admidate, ICD) %>% distinct() %>% group_by(ICD) %>% summarise(freq = n()) %>% arrange(desc(freq))

pneumo_code_desc <- read_delim("exeter_icd10_pneumonia.txt")
tab_pneumo_code_freq <- tab_pneumo_code_freq %>% left_join(pneumo_code_desc, c("ICD" = "ICD10")) %>% select(ICD10 = ICD, Description = Term_description, Frequency = freq)

kableExtra::kable(tab_pneumo_code_freq,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  column_spec(1,width="3cm") %>%
  column_spec(2,width="12cm") %>%
  column_spec(3,width="3cm") %>%
  cat(.,file="STab3c_pneumonia_outcome_code_frequencies.html")

#END#
rm(list=ls())
