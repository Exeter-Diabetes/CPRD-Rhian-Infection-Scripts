################################################################################
###OUTCOMES#####################################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load tidyverse
library(tidyverse)

#Load aurum package
library(aurum)
library(EHRBiomarkr)

###Connecting to data analysis##################################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Connect to analysis
analysis = cprd$analysis("Rhian_covid")

#Set cohort/index dates
cohort.name <- "sep2016"
infection <- "influenza"
index.date <- as.Date("2016-09-01")
end.date <- as.Date("2019-05-31")

#Load cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name,"_cohort_all"), unique_indexes = "patid")

###HOSPITALISATION##############################################################
###Outcome: hospitalisation, any diagnoses in HES (excluding elective admissions: admimeth = 11/12/13)
hospitalisations <- cohort %>% select(patid) %>% left_join(cprd$tables$hesDiagnosisHosp) %>% filter(admidate >= index.date & admidate <= end.date) %>%
  inner_join(codes[[paste0("icd10_",infection)]], by=c("ICD"="icd10")) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>%
  mutate(same_day_discharge = ifelse(duration ==0, 1, 0)) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosps"), indexes = c("patid", "admidate"))

#Find first hospitalisation after index date
first_hospitalisation <- hospitalisations %>% select(patid, spno, admidate) %>% distinct() %>% group_by(patid) %>% summarise(hosp_date = min(admidate)) %>% mutate(hospitalisation_outcome =1) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosp_first"), indexes = c("patid", "hosp_date"))

##Hospitalisation with infection as primary diagnosis (for sensitivity analysis)
primary_hospitalisations <- cohort %>% select(patid) %>% left_join(cprd$tables$hesPrimaryDiagHosp) %>% filter(admidate >= index.date & admidate <= end.date) %>%
  inner_join(codes[[paste0("icd10_",infection)]], by=c("ICD_PRIMARY"="icd10")) %>% left_join(cprd$tables$hesHospital) %>% filter(admimeth != "11" & admimeth != "12" & admimeth != "13") %>%
  mutate(same_day_discharge = ifelse(duration ==0, 1, 0))

#First first primary diagnosis hospitalisation
primary_first_hospitalisation <- primary_hospitalisations %>% select(patid, spno, admidate) %>% distinct() %>% group_by(patid) %>% summarise(hosp_date_primary = min(admidate)) %>% mutate(primary_diag_hosp =1) %>% analysis$cached(paste0(cohort.name, "_", infection, "_prim_diag"))

###DEATH########################################################################
###Outcome: death, any cause in ONS, or a hospitalisation (as defined above) where the patient died (dismeth =4 or disdest =79)
#ONS deaths
deaths_ons <- cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause"="icd10")) %>% mutate(primary_death_cause =1) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause1"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause2"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause3"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause4"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause5"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause6"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause7"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause8"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause9"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause10"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause11"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause12"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause13"="icd10"))) %>%
  union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause14"="icd10"))) %>% union(cprd$tables$onsDeath %>% inner_join(codes[[paste0("icd10_",infection)]], by=c("cause15"="icd10"))) %>%
  analysis$cached(paste0(infection, "_deaths_any_ons"), indexes = "patid")

deaths_ons <- deaths_ons %>% group_by(patid) %>% mutate(primary_death_cause = max(primary_death_cause)) %>% ungroup() %>% select(patid, primary_death_cause) %>% distinct() %>% mutate(death_in_ONS =1)

#Deaths in HES (discharge method/ location = death)
HES_deaths <- hospitalisations %>% filter(dismeth==4 | disdest==79) %>% analysis$cached(paste0(cohort.name, "_", infection, "_hosp_deaths"), indexes = "patid")
HES_deaths <- HES_deaths %>% distinct(patid) %>% mutate(death_in_HES = 1)

###Add outcomes to cohort#######################################################
cohort <- cohort %>% left_join(first_hospitalisation) %>% mutate(hospitalisation_outcome = ifelse(is.na(hospitalisation_outcome), 0, hospitalisation_outcome)) %>%
  left_join(primary_first_hospitalisation) %>% mutate(primary_diag_hosp = ifelse(is.na(primary_diag_hosp), 0, primary_diag_hosp)) %>%
  left_join(deaths_ons) %>% mutate(death_in_ONS = ifelse(death_in_ONS ==1 & dod >= index.date & dod <= end.date & !is.na(dod), 1, 0)) %>% mutate(primary_death_cause = ifelse(primary_death_cause ==1 & death_in_ONS ==1, 1, 0)) %>%
  left_join(HES_deaths) %>% mutate(death_in_HES = ifelse(death_in_HES ==1 & dod >= index.date & dod <= end.date & !is.na(dod), 1, 0)) %>% mutate(death_outcome = ifelse(death_in_ONS==1 | death_in_HES ==1, 1, 0)) %>%
  mutate(primary_death_cause = ifelse(is.na(primary_death_cause), 0, primary_death_cause), death_in_ONS = ifelse(is.na(death_in_ONS), 0, death_in_ONS), death_in_HES = ifelse(is.na(death_in_HES), 0, death_in_HES), death_outcome = ifelse(is.na(death_outcome), 0, death_outcome)) %>%
  analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = c("hosp_date", "dod"))

###END##########################################################################
rm(list=ls())
