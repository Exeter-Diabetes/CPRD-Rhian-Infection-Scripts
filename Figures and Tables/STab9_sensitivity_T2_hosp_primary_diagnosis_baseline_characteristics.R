################################################################################
###Descriptive tables###########################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(lubridate)
library(tableone)
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
###COVID
###Set cohort###################################################################
cohort.name <- "feb2020"
infection <- "covid"
outcome <- "hosp_p" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = primary_diag_hosp, outcome_date = hosp_date_primary)
#Set index and end dates
index.date = as.Date("2020-02-01") #change these for different cohorts
end.date = as.Date("2020-10-31")

###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)
#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age_all >=20 & age_at_index >=18)
#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Setting variables to factors and setting reference category
#Gender
cohort$gender <- factor(cohort$gender)
levels(cohort$gender) = c("Male", "Female")
cohort$gender <- relevel(cohort$gender, ref="Female")
cohort$femalegender <- relevel(factor(cohort$gender), ref = "Male")
#Age
cohort$age_cat <- factor(cohort$age_cat)
#Ethnicity
cohort$eth5 <- factor(cohort$eth5)
levels(cohort$eth5) = c("White", "South Asian", "Black", "Other", "Mixed", "Unknown")
#Deprivation
cohort$imd_quintile <- factor(cohort$imd_quintile)
#Region
cohort$region <- factor(cohort$region)
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#Diabetes complications (as count of number)
cohort <- cohort %>% mutate(number_complications = ifelse(number_complications ==0, "0 complications", ifelse(number_complications ==1, "1 complication", ifelse(number_complications ==2, "2 complications", ifelse(number_complications ==3, "3 complications", NA)))))
cohort$number_complications <- factor(cohort$number_complications)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
#Smoking status
cohort$smoking_status <- factor(cohort$smoking_status)
#Recent respiratory infection
cohort$recent_hosp_resp_infect <- factor(cohort$recent_hosp_resp_infect)
levels(cohort$recent_hosp_resp_infect) <- c("No", "Yes")
#Recent hospitalisation for anything
cohort$recent_hosp_anything <- factor(cohort$recent_hosp_anything)
levels(cohort$recent_hosp_anything) <- c("No", "Yes")
#CKD Stage
cohort$ckd_stage <- factor(cohort$ckd_stage)
#ACR
cohort$acr_cat <- factor(cohort$acr_cat)
#Comorbidities
#Heart failure
cohort$preindex_heartfailure <- factor(cohort$preindex_heartfailure)
levels(cohort$preindex_heartfailure) <- c("No", "Yes")
#Myocardial infarction
cohort$preindex_myocardialinfarction <- factor(cohort$preindex_myocardialinfarction)
levels(cohort$preindex_myocardialinfarction) <- c("No", "Yes")
#Stroke
cohort$preindex_stroke <- factor(cohort$preindex_stroke)
levels(cohort$preindex_stroke) <- c("No", "Yes")
#Atrial fibrillation
cohort$preindex_af <- factor(cohort$preindex_af)
levels(cohort$preindex_af) <- c("No", "Yes")
#Angina
cohort$preindex_angina <- factor(cohort$preindex_angina)
levels(cohort$preindex_angina) <- c("No", "Yes")
#Asthma
cohort$preindex_asthma <- factor(cohort$preindex_asthma)
levels(cohort$preindex_asthma) <- c("No", "Yes")
#Chronic liver disease
cohort$preindex_cld <- factor(cohort$preindex_cld)
levels(cohort$preindex_cld) <- c("No", "Yes")
#COPD
cohort$preindex_copd <- factor(cohort$preindex_copd)
levels(cohort$preindex_copd) <- c("No", "Yes")
#Dementia
cohort$preindex_dementia <- factor(cohort$preindex_dementia)
levels(cohort$preindex_dementia) <- c("No", "Yes")
#Hypertension
cohort$preindex_hypertension <- factor(cohort$preindex_hypertension)
levels(cohort$preindex_hypertension) <- c("No", "Yes")
#Ischaemic heart disease
cohort$preindex_ihd <- factor(cohort$preindex_ihd)
levels(cohort$preindex_ihd) <- c("No", "Yes")
#Peripheral arterial disease
cohort$preindex_pad <- factor(cohort$preindex_pad)
levels(cohort$preindex_pad) <- c("No", "Yes")
#Revascularisation
cohort$preindex_revasc <- factor(cohort$preindex_revasc)
levels(cohort$preindex_revasc) <- c("No", "Yes")
#Solid organ transplant
cohort$preindex_solidorgantransplant <- factor(cohort$preindex_solidorgantransplant)
levels(cohort$preindex_solidorgantransplant) <- c("No", "Yes")
#Transient ischaemic attack
cohort$preindex_tia <- factor(cohort$preindex_tia)
levels(cohort$preindex_tia) <- c("No", "Yes")
#Other neurological condition
cohort$preindex_otherneuroconditions <- factor(cohort$preindex_otherneuroconditions)
levels(cohort$preindex_otherneuroconditions) <- c("No", "Yes")
#Haematological cancer
cohort$preindex_haem_cancer <- factor(cohort$preindex_haem_cancer)
levels(cohort$preindex_haem_cancer) <- c("No", "Yes")
#Solid cancer
cohort$preindex_solid_cancer <- factor(cohort$preindex_solid_cancer)
levels(cohort$preindex_solid_cancer) <- c("No", "Yes")
##Diabetes medications
cohort$treatment_6m <- factor(cohort$treatment_6m)

#Generate tableone
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
              "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
              "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
              "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

categorical_vars <-c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
                     "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
                     "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
                     "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1) #might also want to include deaths?

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_covid <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename ("2020 cohort" = All, "Covid-19 hospitalisations" = Hospitalisations)


################################################################################
###FLU
###Set cohort###################################################################
cohort.name <- "sep2016"
infection <- "influenza"
outcome <- "hosp_p" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = primary_diag_hosp, outcome_date = hosp_date_primary)
#Set index and end dates
index.date = as.Date("2016-09-01") #change these for different cohorts
end.date = as.Date("2019-05-31")

###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)
#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age_all >=20 & age_at_index >=18)
#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Setting variables to factors and setting reference category
#Gender
cohort$gender <- factor(cohort$gender)
levels(cohort$gender) = c("Male", "Female")
cohort$gender <- relevel(cohort$gender, ref="Female")
cohort$femalegender <- relevel(factor(cohort$gender), ref = "Male")
#Age
cohort$age_cat <- factor(cohort$age_cat)
#Ethnicity
cohort$eth5 <- factor(cohort$eth5)
levels(cohort$eth5) = c("White", "South Asian", "Black", "Other", "Mixed", "Unknown")
#Deprivation
cohort$imd_quintile <- factor(cohort$imd_quintile)
#Region
cohort$region <- factor(cohort$region)
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#Diabetes complications (as count of number)
cohort <- cohort %>% mutate(number_complications = ifelse(number_complications ==0, "0 complications", ifelse(number_complications ==1, "1 complication", ifelse(number_complications ==2, "2 complications", ifelse(number_complications ==3, "3 complications", NA)))))
cohort$number_complications <- factor(cohort$number_complications)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
#Smoking status
cohort$smoking_status <- factor(cohort$smoking_status)
#Recent respiratory infection
cohort$recent_hosp_resp_infect <- factor(cohort$recent_hosp_resp_infect)
levels(cohort$recent_hosp_resp_infect) <- c("No", "Yes")
#Recent hospitalisation for anything
cohort$recent_hosp_anything <- factor(cohort$recent_hosp_anything)
levels(cohort$recent_hosp_anything) <- c("No", "Yes")
#CKD Stage
cohort$ckd_stage <- factor(cohort$ckd_stage)
#ACR
cohort$acr_cat <- factor(cohort$acr_cat)
#Comorbidities
#Heart failure
cohort$preindex_heartfailure <- factor(cohort$preindex_heartfailure)
levels(cohort$preindex_heartfailure) <- c("No", "Yes")
#Myocardial infarction
cohort$preindex_myocardialinfarction <- factor(cohort$preindex_myocardialinfarction)
levels(cohort$preindex_myocardialinfarction) <- c("No", "Yes")
#Stroke
cohort$preindex_stroke <- factor(cohort$preindex_stroke)
levels(cohort$preindex_stroke) <- c("No", "Yes")
#Atrial fibrillation
cohort$preindex_af <- factor(cohort$preindex_af)
levels(cohort$preindex_af) <- c("No", "Yes")
#Angina
cohort$preindex_angina <- factor(cohort$preindex_angina)
levels(cohort$preindex_angina) <- c("No", "Yes")
#Asthma
cohort$preindex_asthma <- factor(cohort$preindex_asthma)
levels(cohort$preindex_asthma) <- c("No", "Yes")
#Chronic liver disease
cohort$preindex_cld <- factor(cohort$preindex_cld)
levels(cohort$preindex_cld) <- c("No", "Yes")
#COPD
cohort$preindex_copd <- factor(cohort$preindex_copd)
levels(cohort$preindex_copd) <- c("No", "Yes")
#Dementia
cohort$preindex_dementia <- factor(cohort$preindex_dementia)
levels(cohort$preindex_dementia) <- c("No", "Yes")
#Hypertension
cohort$preindex_hypertension <- factor(cohort$preindex_hypertension)
levels(cohort$preindex_hypertension) <- c("No", "Yes")
#Ischaemic heart disease
cohort$preindex_ihd <- factor(cohort$preindex_ihd)
levels(cohort$preindex_ihd) <- c("No", "Yes")
#Peripheral arterial disease
cohort$preindex_pad <- factor(cohort$preindex_pad)
levels(cohort$preindex_pad) <- c("No", "Yes")
#Revascularisation
cohort$preindex_revasc <- factor(cohort$preindex_revasc)
levels(cohort$preindex_revasc) <- c("No", "Yes")
#Solid organ transplant
cohort$preindex_solidorgantransplant <- factor(cohort$preindex_solidorgantransplant)
levels(cohort$preindex_solidorgantransplant) <- c("No", "Yes")
#Transient ischaemic attack
cohort$preindex_tia <- factor(cohort$preindex_tia)
levels(cohort$preindex_tia) <- c("No", "Yes")
#Other neurological condition
cohort$preindex_otherneuroconditions <- factor(cohort$preindex_otherneuroconditions)
levels(cohort$preindex_otherneuroconditions) <- c("No", "Yes")
#Haematological cancer
cohort$preindex_haem_cancer <- factor(cohort$preindex_haem_cancer)
levels(cohort$preindex_haem_cancer) <- c("No", "Yes")
#Solid cancer
cohort$preindex_solid_cancer <- factor(cohort$preindex_solid_cancer)
levels(cohort$preindex_solid_cancer) <- c("No", "Yes")
##Diabetes medications
cohort$treatment_6m <- factor(cohort$treatment_6m)

#Generate tableone
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
              "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
              "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
              "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

categorical_vars <-c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
                     "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
                     "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
                     "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1) #might also want to include deaths?

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_flu <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename ("2016 cohort" = All, "Influenza hospitalisations" = Hospitalisations)


################################################################################
###PNEUMONIA
###Set cohort###################################################################
cohort.name <- "sep2016"
infection <- "pneumonia"
outcome <- "hosp_p" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = primary_diag_hosp, outcome_date = hosp_date_primary)
#Set index and end dates
index.date = as.Date("2016-09-01") #change these for different cohorts
end.date = as.Date("2019-05-31")

###Calculate survival dates and times###########################################
cohort <- cohort %>% mutate(survival_date = end.date) %>% mutate(survival_date = ifelse(outcome ==1 & outcome_date < survival_date, outcome_date, survival_date)) %>%
  mutate(survival_date = ifelse(!is.na(regenddate) & regenddate < survival_date, regenddate, survival_date)) %>% mutate(survival_date = ifelse(!is.na(dod) & dod < survival_date, dod, survival_date)) %>% 
  mutate(outcome = ifelse(outcome ==1 & outcome_date == survival_date, 1, 0)) %>% mutate(survival_time = datediff(survival_date, index.date)) %>% 
  analysis$cached(paste0("surv_",infection, "_", outcome), unique_indexes = "patid", indexes = c("survival_date", "survival_time"))

#Collect
cohort <- collect(cohort)
#Filter cohort
cohort <- cohort %>% filter(diabetes_type == "type 2" & dm_diag_age_all >=20 & age_at_index >=18)
#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Setting variables to factors and setting reference category
#Gender
cohort$gender <- factor(cohort$gender)
levels(cohort$gender) = c("Male", "Female")
cohort$gender <- relevel(cohort$gender, ref="Female")
cohort$femalegender <- relevel(factor(cohort$gender), ref = "Male")
#Age
cohort$age_cat <- factor(cohort$age_cat)
#Ethnicity
cohort$eth5 <- factor(cohort$eth5)
levels(cohort$eth5) = c("White", "South Asian", "Black", "Other", "Mixed", "Unknown")
#Deprivation
cohort$imd_quintile <- factor(cohort$imd_quintile)
#Region
cohort$region <- factor(cohort$region)
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#Diabetes complications (as count of number)
cohort <- cohort %>% mutate(number_complications = ifelse(number_complications ==0, "0 complications", ifelse(number_complications ==1, "1 complication", ifelse(number_complications ==2, "2 complications", ifelse(number_complications ==3, "3 complications", NA)))))
cohort$number_complications <- factor(cohort$number_complications)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
#Smoking status
cohort$smoking_status <- factor(cohort$smoking_status)
#Recent respiratory infection
cohort$recent_hosp_resp_infect <- factor(cohort$recent_hosp_resp_infect)
levels(cohort$recent_hosp_resp_infect) <- c("No", "Yes")
#Recent hospitalisation for anything
cohort$recent_hosp_anything <- factor(cohort$recent_hosp_anything)
levels(cohort$recent_hosp_anything) <- c("No", "Yes")
#CKD Stage
cohort$ckd_stage <- factor(cohort$ckd_stage)
#ACR
cohort$acr_cat <- factor(cohort$acr_cat)
#Comorbidities
#Heart failure
cohort$preindex_heartfailure <- factor(cohort$preindex_heartfailure)
levels(cohort$preindex_heartfailure) <- c("No", "Yes")
#Myocardial infarction
cohort$preindex_myocardialinfarction <- factor(cohort$preindex_myocardialinfarction)
levels(cohort$preindex_myocardialinfarction) <- c("No", "Yes")
#Stroke
cohort$preindex_stroke <- factor(cohort$preindex_stroke)
levels(cohort$preindex_stroke) <- c("No", "Yes")
#Atrial fibrillation
cohort$preindex_af <- factor(cohort$preindex_af)
levels(cohort$preindex_af) <- c("No", "Yes")
#Angina
cohort$preindex_angina <- factor(cohort$preindex_angina)
levels(cohort$preindex_angina) <- c("No", "Yes")
#Asthma
cohort$preindex_asthma <- factor(cohort$preindex_asthma)
levels(cohort$preindex_asthma) <- c("No", "Yes")
#Chronic liver disease
cohort$preindex_cld <- factor(cohort$preindex_cld)
levels(cohort$preindex_cld) <- c("No", "Yes")
#COPD
cohort$preindex_copd <- factor(cohort$preindex_copd)
levels(cohort$preindex_copd) <- c("No", "Yes")
#Dementia
cohort$preindex_dementia <- factor(cohort$preindex_dementia)
levels(cohort$preindex_dementia) <- c("No", "Yes")
#Hypertension
cohort$preindex_hypertension <- factor(cohort$preindex_hypertension)
levels(cohort$preindex_hypertension) <- c("No", "Yes")
#Ischaemic heart disease
cohort$preindex_ihd <- factor(cohort$preindex_ihd)
levels(cohort$preindex_ihd) <- c("No", "Yes")
#Peripheral arterial disease
cohort$preindex_pad <- factor(cohort$preindex_pad)
levels(cohort$preindex_pad) <- c("No", "Yes")
#Revascularisation
cohort$preindex_revasc <- factor(cohort$preindex_revasc)
levels(cohort$preindex_revasc) <- c("No", "Yes")
#Solid organ transplant
cohort$preindex_solidorgantransplant <- factor(cohort$preindex_solidorgantransplant)
levels(cohort$preindex_solidorgantransplant) <- c("No", "Yes")
#Transient ischaemic attack
cohort$preindex_tia <- factor(cohort$preindex_tia)
levels(cohort$preindex_tia) <- c("No", "Yes")
#Other neurological condition
cohort$preindex_otherneuroconditions <- factor(cohort$preindex_otherneuroconditions)
levels(cohort$preindex_otherneuroconditions) <- c("No", "Yes")
#Haematological cancer
cohort$preindex_haem_cancer <- factor(cohort$preindex_haem_cancer)
levels(cohort$preindex_haem_cancer) <- c("No", "Yes")
#Solid cancer
cohort$preindex_solid_cancer <- factor(cohort$preindex_solid_cancer)
levels(cohort$preindex_solid_cancer) <- c("No", "Yes")
##Diabetes medications
cohort$treatment_6m <- factor(cohort$treatment_6m)

#Generate tableone
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
              "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
              "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
              "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

categorical_vars <-c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "number_complications", "bmi_cat", "smoking_status", 
                     "preindex_hypertension", "preindex_af", "preindex_angina", "preindex_myocardialinfarction", "preindex_revasc", "preindex_ihd", "preindex_heartfailure", "preindex_pad",
                     "recent_hosp_resp_infect", "recent_hosp_anything", "preindex_asthma", "preindex_copd", "preindex_tia", "preindex_stroke",  "preindex_dementia", "preindex_otherneuroconditions",
                     "preindex_haem_cancer", "preindex_solid_cancer", "preindex_solidorgantransplant", "preindex_cld", "ckd_stage", "acr_cat", "treatment_6m")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1) #might also want to include deaths?

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_pneumo <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename ("2016 cohort" = All, "Pneumonia hospitalisations" = Hospitalisations)

################################################################################
##Add together 
tabprint <- tabprint_covid %>% left_join(tabprint_flu) %>% left_join(tabprint_pneumo)

row_names <- data.frame(measure = c("n", "femalegender = Female (%)", "gender = Male (%)", "<40", "40-49","50-59", "60-69", "70-79", "80-89", "90+", "White", "South Asian", "Black", "Other", "Mixed", "Unknown",
               "1", "2", "3", "4", "5", "Missing IMD", "<1", "1-2", "3-5", "6-9", "10-14", "15-19", "20+", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c", "0 complications", "1 complication", "2 complications", "3 complications",
               "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI", "Active smoker", "Ex-smoker", "Non-smoker", "Unknown smoking", "preindex_hypertension = Yes (%)", "preindex_af = Yes (%)", "preindex_angina = Yes (%)",
               "preindex_myocardialinfarction = Yes (%)", "preindex_revasc = Yes (%)", "preindex_ihd = Yes (%)", "preindex_heartfailure = Yes (%)", "preindex_pad = Yes (%)", "recent_hosp_resp_infect = Yes (%)", "recent_hosp_anything = Yes (%)",
               "preindex_asthma = Yes (%)", "preindex_copd = Yes (%)", "preindex_tia = Yes (%)", "preindex_stroke = Yes (%)", "preindex_dementia = Yes (%)", "preindex_otherneuroconditions = Yes (%)", "preindex_haem_cancer = Yes (%)", "preindex_solid_cancer = Yes (%)",
               "preindex_solidorgantransplant = Yes (%)", "preindex_cld = Yes (%)", "Stage 1", "Stage 2", "Stage 3a", "Stage 3b", "Stage 4", "Stage 5", "Missing CKD", "A1", "A2", "A3", "Missing ACR",
               "Insulin", "OHA only", "No treatment"))

tabprint <- row_names %>% left_join(tabprint, by = "measure")

tab <- tabprint %>% 
  mutate(measure=ifelse(measure=="gender = Male (%)","Male", measure)) %>%
  mutate(measure=ifelse(measure=="femalegender = Female (%)","Female", measure)) %>%
  mutate(measure = ifelse(measure == "1", "1 (least deprived)", measure)) %>%
  mutate(measure = ifelse(measure == "5", "5 (most deprived)", measure)) %>%
  mutate(measure=ifelse(measure=="Missing IMD","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="Missing HbA1c","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="Missing BMI","Missing", measure)) %>%
  mutate(measure = ifelse(measure == "0 complications", "0", measure)) %>%
  mutate(measure = ifelse(measure == "1 complication", "1", measure)) %>%
  mutate(measure = ifelse(measure == "2 complications", "2", measure)) %>%
  mutate(measure = ifelse(measure == "3 complications", "3", measure)) %>%
  mutate(measure=ifelse(measure=="Unknown smoking","Unknown", measure)) %>%
  mutate(measure=ifelse(measure=="Missing ACR","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="Missing CKD","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="recent_hosp_resp_infect = Yes (%)","Respiratory infection", measure)) %>%
  mutate(measure=ifelse(measure=="recent_hosp_anything = Yes (%)","Anything else", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_heartfailure = Yes (%)","Heart failure", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_myocardialinfarction = Yes (%)","Previous myocardial infarction", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_stroke = Yes (%)","Previous stroke", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_af = Yes (%)","Atrial fibrillation", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_angina = Yes (%)","Angina", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_asthma = Yes (%)","Asthma", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_cld = Yes (%)","Chronic liver disease", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_copd = Yes (%)","Chronic obstructive pulmonary disease", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_dementia = Yes (%)","Dementia", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_hypertension = Yes (%)","Hypertension", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_ihd = Yes (%)","Other ischaemic heart disease", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_pad = Yes (%)","Peripheral arterial disease", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_revasc = Yes (%)","Previous cardiac revascularisation", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_solidorgantransplant = Yes (%)","Solid organ transplant", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_tia = Yes (%)","Previous transient ischaemic attack", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_otherneuroconditions = Yes (%)","Other neurological condition", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_haem_cancer = Yes (%)","Haematological cancer", measure)) %>%
  mutate(measure=ifelse(measure=="preindex_solid_cancer = Yes (%)","Solid cancer", measure)) %>%
  mutate(measure=ifelse(measure=="Insulin","Insulin prescription", measure)) %>%
  mutate(measure=ifelse(measure=="No treatment","No treatment", measure)) %>%
  mutate(measure=ifelse(measure=="OHA only","OHA prescription only", measure)) %>%
  rename(" " = measure)


#Outputting HTML table
kableExtra::kable(tab,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  pack_rows("Sex", 2,3, bold = TRUE) %>%
  pack_rows("Age group, years",4,10, bold=TRUE) %>% 
  pack_rows("Ethnicity",11,16, bold=TRUE)  %>% 
  pack_rows("Index of multiple deprivation quintile",17,22, bold=TRUE)  %>%  
  pack_rows("Duration of diagnosed diabetes, years",23,29, bold=TRUE)  %>%  
  pack_rows("HbA1c, mmol/mol",30,36, bold=TRUE)  %>% 
  pack_rows("Number of microvascular complications",37,40, bold=TRUE)  %>% 
  pack_rows("BMI, kg/m2",41,47, bold=TRUE)  %>%
  pack_rows("Smoking status",48,51, bold=TRUE)  %>% 
  pack_rows("Comorbidities",52,71, bold=TRUE)  %>% 
  pack_rows("Cardiovascular",52,59, bold=TRUE)  %>% 
  pack_rows("Recent hospitalisation",60,61, bold=TRUE)  %>% 
  pack_rows("Respiratory",62,63, bold=TRUE)  %>% 
  pack_rows("Neurological",64,67, bold=TRUE)  %>% 
  pack_rows("Oncological",68,69, bold=TRUE)  %>% 
  pack_rows("Other",70,71, bold=TRUE)  %>% 
  pack_rows("Chronic Kidney Disease (CKD) stage",72,78, bold=TRUE)  %>% 
  pack_rows("Albumin creatinine ratio (ACR) category",79,82, bold=TRUE)  %>% 
  pack_rows("Diabetes treatment (last 6 months)",83,85, bold=TRUE)  %>%
  column_spec(1,width="12cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="7cm") %>%
  column_spec(4,width="6cm") %>%
  column_spec(5,width="7cm") %>%
  column_spec(6,width="7cm") %>%
  cat(.,file="STab9_sensitivity_T2_hosp_primary_diagnosis_baseline_characteristics.html")

###END##########################################################################
rm(list=ls())


