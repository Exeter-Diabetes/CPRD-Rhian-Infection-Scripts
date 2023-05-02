################################################################################
###ANALYSIS#####################################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(lubridate)
library(survminer)
library(survival)
library(forestplot)

#Load aurum package
library(aurum)

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Setting up/loading analysis test
analysis = cprd$analysis("Rhian_covid")


###Set first cohort and outcome#################################################
cohort.name <- "sep2016"
infection <- "influenza"
outcome <- "hosp" #these will be pasted into output file names
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = hospitalisation_outcome, outcome_date = hosp_date)
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

#Check survival time distribution for those with outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

################################################################################
#Adding flu vacc
#Flu vaccination last two years
fluvacc_raw_1 <- cprd$tables$observation %>% inner_join(codes$fluvacc_stopflu_med) %>% analysis$cached("all_fluvacc_stopflu1", indexes= c("patid", "obsdate"))
fluvacc_raw_2 <- cprd$tables$drugIssue %>% inner_join(codes$fluvacc_stopflu_prod) %>% analysis$cached("all_fluvacc_stopflu2", indexes= c("patid", "issuedate"))

two.years.earlier <- index.date - years(2)

fluvacc_1_2years <- fluvacc_raw_1 %>% filter(!is.na(obsdate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(obsdate>=min_dob & obsdate<=gp_ons_end_date) %>% 
  filter(obsdate < index.date & obsdate >= two.years.earlier) %>% distinct(patid) %>% mutate(fluvacc_2years =1)
fluvacc_2_2years <- fluvacc_raw_2 %>% filter(!is.na(issuedate)) %>% inner_join(cprd$tables$validDateLookup) %>% filter(issuedate>=min_dob & issuedate<=gp_ons_end_date) %>% 
  filter(issuedate < index.date & issuedate >= two.years.earlier) %>% distinct(patid) %>% mutate(fluvacc_2years =1)

fluvacc <- collect(fluvacc_1_2years %>% union(fluvacc_2_2years))

cohort <- cohort %>% left_join(fluvacc, by = "patid") %>% mutate(fluvacc_2years = ifelse(is.na(fluvacc_2years), 0 , fluvacc_2years))
cohort %>% filter(fluvacc_2years ==1) %>% count()
cohort %>% filter(fluvacc_2years ==0) %>% count()

#Check survival time distribution for those with outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")
cohort %>% filter(outcome ==1 & fluvacc_2years ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")
cohort %>% filter(outcome ==1 & fluvacc_2years ==0) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")

################################################################################
#Setting variables to factors and setting reference category
#Gender
cohort$gender <- factor(cohort$gender)
levels(cohort$gender) = c("Male", "Female")
cohort$gender <- relevel(cohort$gender, ref="Female")
#Age
cohort$age_cat <- factor(cohort$age_cat)
cohort$age_cat <- relevel(cohort$age_cat, ref="60-69")
#Ethnicity
cohort$eth5 <- factor(cohort$eth5)
levels(cohort$eth5) = c("White", "South Asian", "Black", "Other", "Mixed", "Unknown")
cohort$eth5 <- relevel(cohort$eth5, ref="White")
#Deprivation
cohort$imd_quintile <- factor(cohort$imd_quintile)
cohort$imd_quintile <- relevel(cohort$imd_quintile, ref="1")
#Region
cohort$region <- factor(cohort$region)
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
cohort$duration_cat <- relevel(cohort$duration_cat, ref = "10-14")
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
cohort$hba1c_cat <- relevel(cohort$hba1c_cat, ref = "48-53")
#Diabetes complications (as count of number)
cohort$number_complications <- factor(cohort$number_complications)
cohort$number_complications <- relevel(cohort$number_complications, ref = "0")
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
cohort$bmi_cat <- relevel(cohort$bmi_cat, ref = "25-29.9")
#Smoking status
cohort$smoking_status <- factor(cohort$smoking_status)
cohort$smoking_status <- relevel(cohort$smoking_status, ref = "Non-smoker")
#Alcohol consumption
cohort$alcohol_consumption <- factor(cohort$alcohol_consumption)
cohort$alcohol_consumption <- relevel(cohort$alcohol_consumption, ref = "AlcoholConsumptionLevel0")
#Recent respiratory infection
cohort$recent_hosp_resp_infect <- factor(cohort$recent_hosp_resp_infect)
levels(cohort$recent_hosp_resp_infect) <- c("No", "Yes")
#Recent hospitalisation for anything
cohort$recent_hosp_anything <- factor(cohort$recent_hosp_anything)
levels(cohort$recent_hosp_anything) <- c("No", "Yes")
#CKD Stage
cohort$ckd_stage <- factor(cohort$ckd_stage)
cohort$ckd_stage <- relevel(cohort$ckd_stage, ref = "Stage 1")
#ACR
cohort$acr_cat <- factor(cohort$acr_cat)
cohort$acr_cat <- relevel(cohort$acr_cat, ref = "A1")
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
#Haematological cancer
cohort$preindex_haem_cancer <- factor(cohort$preindex_haem_cancer)
levels(cohort$preindex_haem_cancer) <- c("No", "Yes")
#Solid cancer
cohort$preindex_solid_cancer <- factor(cohort$preindex_solid_cancer)
levels(cohort$preindex_solid_cancer) <- c("No", "Yes")
#Other neurological conditions
cohort$preindex_otherneuroconditions <- factor(cohort$preindex_otherneuroconditions)
levels(cohort$preindex_otherneuroconditions) <- c("No", "Yes")
#Cystic fibrosis- will want to exclude these! Just a check
cohort$preindex_cysticfibrosis <- factor(cohort$preindex_cysticfibrosis)
levels(cohort$preindex_cysticfibrosis) <- c("No", "Yes")
#Pulmonary fibrosis
cohort$preindex_pulmonaryfibrosis <- factor(cohort$preindex_pulmonaryfibrosis)
levels(cohort$preindex_pulmonaryfibrosis) <- c("No", "Yes")
#Pulmonary hypertension
cohort$preindex_pulmonaryhypertension <- factor(cohort$preindex_pulmonaryhypertension)
levels(cohort$preindex_pulmonaryhypertension) <- c("No", "Yes")
#Bronchiectasis
cohort$preindex_bronchiectasis <- factor(cohort$preindex_bronchiectasis)
levels(cohort$preindex_bronchiectasis) <- c("No", "Yes")
##Diabetes medications
cohort$treatment_6m <- factor(cohort$treatment_6m)
cohort$treatment_6m <- relevel(cohort$treatment_6m, ref = "No treatment")
#Other medications
#Oral steroids
cohort$oralsteroids_6m <- factor(cohort$oralsteroids_6m)
levels(cohort$oralsteroids_6m) <- c("No", "Yes")
#Immunosuppressants
cohort$immunosuppressants_6m <- factor(cohort$immunosuppressants_6m)
levels(cohort$immunosuppressants_6m) <- c("No", "Yes")
#Long acting beta agonists
cohort$labas_6m <- factor(cohort$labas_6m)
levels(cohort$labas_6m) <- c("No", "Yes")
#Leukotriene receptor antagonists
cohort$ltras_6m <- factor(cohort$ltras_6m)
levels(cohort$ltras_6m) <- c("No", "Yes")


cohort_vaccinated <- cohort %>% filter(fluvacc_2years ==1)
cohort_unvaccinated <- cohort %>% filter(fluvacc_2years ==0)

################################################################################
###RUN MODELS AND BRING TOGETHER IN FOREST PLOT#################################
################################################################################

##Run models####################################################################
###Vaccinated
model_vaccinated <- coxph(Surv(survival_time, outcome) ~ gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat + number_complications + bmi_cat + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + 
                 preindex_heartfailure + preindex_pad + recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis +
                 preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions +
                 preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld + ckd_stage + acr_cat +  
                 treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort_vaccinated)
model_df_vaccinated <- data.frame(summary(model_vaccinated)$coefficients, summary(model_vaccinated)$conf.int)

model_df_vaccinated$row_name <- rownames(model_df_vaccinated)

model_df_vaccinated <- model_df_vaccinated %>% mutate(row_name = ifelse(row_name == "genderMale", "Male", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat<40", "<40", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat40-49", "40-49", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat50-59", "50-59", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat70-79", "70-79", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat80-89", "80-89", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat90+", "90+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5South Asian", "South Asian", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Black", "Black", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Other", "Other", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Mixed", "Mixed", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Unknown", "Unknown ethnicity", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile2", "2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile3", "3", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile4", "4", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile5", "5 (most deprived)", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintileMissing IMD", "Missing IMD", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionEast of England", "East of England", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionLondon", "London", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionMissing region", "Missing region", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionNorth East", "North East", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionNorth West", "North West", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth Central", "South Central", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth East Coast", "South East Coast", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth West", "South West", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionWest Midlands", "West Midlands", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionYorkshire And The Humber", "Yorkshire and the Humber", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat<1", "<1", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat1-2", "1-2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat15-19", "15-19", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat20+", "20+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat3-5", "3-5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat6-9", "6-9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat<48", "<48", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat53-64", "53-64", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat64-75", "64-75", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat75-86", "75-86", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat86+", "86+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_catMissing HbA1c", "Missing HbA1c", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications1", "1 complication", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications2", "2 complications", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications3", "3 complications", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat<18.5", "<18.5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat18.5-24.9", "18.5-24.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat30-34.9", "30-34.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat35-39.9", "35-39.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat40+", "40+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_catMissing BMI", "Missing BMI", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusActive smoker", "Active smoker", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusEx-smoker", "Ex-smoker", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusUnknown smoking", "Unknown smoking", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_hypertensionYes", "Hypertension", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_afYes", "Atrial fibrillation", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_anginaYes", "Angina", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_myocardialinfarctionYes", "Previous myocardial infarction", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_revascYes", "Previous cardiac revascularisation", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_ihdYes", "Other ischaemic heart disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_heartfailureYes", "Heart failure", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_padYes", "Peripheral arterial disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "recent_hosp_resp_infectYes", "Respiratory infection", row_name)) %>%
  mutate(row_name = ifelse(row_name == "recent_hosp_anythingYes", "Anything else", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_asthmaYes", "Asthma", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_copdYes", "Chronic obstructive pulmonary disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_bronchiectasisYes", "Bronchiectasis", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_pulmonaryfibrosisYes", "Pulmonary fibrosis", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_pulmonaryhypertensionYes", "Pulmonary hypertension", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_tiaYes", "Previous transient ischaemic attack", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_strokeYes", "Previous stroke", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_dementiaYes", "Dementia", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_otherneuroconditionsYes", "Other neurological condition", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_haem_cancerYes", "Haematological cancer", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_solid_cancerYes", "Solid cancer", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_solidorgantransplantYes", "Solid organ transplant", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_cldYes", "Chronic liver disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageMissing CKD", "Missing CKD", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 2", "Stage 2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 3a", "Stage 3a", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 3b", "Stage 3b", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 4", "Stage 4", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 5", "Stage 5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "acr_catA2", "A2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "acr_catA3", "A3", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "acr_catMissing ACR", "Missing ACR", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "treatment_6mInsulin", "Insulin prescription", row_name)) %>%
  mutate(row_name = ifelse(row_name == "treatment_6mOHA only", "OHA prescription only", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "oralsteroids_6mYes", "Oral steroids", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "immunosuppressants_6mYes", "Immunosuppressants", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "labas_6mYes", "Long acting beta agonists", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "ltras_6mYes", "Leukotriene receptor antagonists", row_name))


###Unvaccinated
model_unvaccinated <- coxph(Surv(survival_time, outcome) ~ gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat + number_complications + bmi_cat + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + 
                 preindex_heartfailure + preindex_pad + recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis +
                 preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions +
                 preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld + ckd_stage + acr_cat +  
                 treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort_unvaccinated)
model_df_unvaccinated <- data.frame(summary(model_unvaccinated)$coefficients, summary(model_unvaccinated)$conf.int)

model_df_unvaccinated$row_name <- rownames(model_df_unvaccinated)

model_df_unvaccinated <- model_df_unvaccinated %>% mutate(row_name = ifelse(row_name == "genderMale", "Male", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat<40", "<40", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat40-49", "40-49", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat50-59", "50-59", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat70-79", "70-79", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat80-89", "80-89", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat90+", "90+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5South Asian", "South Asian", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Black", "Black", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Other", "Other", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Mixed", "Mixed", row_name)) %>%
  mutate(row_name = ifelse(row_name == "eth5Unknown", "Unknown ethnicity", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile2", "2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile3", "3", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile4", "4", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintile5", "5 (most deprived)", row_name)) %>%
  mutate(row_name = ifelse(row_name == "imd_quintileMissing IMD", "Missing IMD", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionEast of England", "East of England", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionLondon", "London", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionMissing region", "Missing region", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionNorth East", "North East", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionNorth West", "North West", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth Central", "South Central", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth East Coast", "South East Coast", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionSouth West", "South West", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionWest Midlands", "West Midlands", row_name)) %>%
  mutate(row_name = ifelse(row_name == "regionYorkshire And The Humber", "Yorkshire and the Humber", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat<1", "<1", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat1-2", "1-2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat15-19", "15-19", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat20+", "20+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat3-5", "3-5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat6-9", "6-9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat<48", "<48", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat53-64", "53-64", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat64-75", "64-75", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat75-86", "75-86", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat86+", "86+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_catMissing HbA1c", "Missing HbA1c", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications1", "1 complication", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications2", "2 complications", row_name)) %>%
  mutate(row_name = ifelse(row_name == "number_complications3", "3 complications", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat<18.5", "<18.5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat18.5-24.9", "18.5-24.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat30-34.9", "30-34.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat35-39.9", "35-39.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat40+", "40+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_catMissing BMI", "Missing BMI", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusActive smoker", "Active smoker", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusEx-smoker", "Ex-smoker", row_name)) %>%
  mutate(row_name = ifelse(row_name == "smoking_statusUnknown smoking", "Unknown smoking", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_hypertensionYes", "Hypertension", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_afYes", "Atrial fibrillation", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_anginaYes", "Angina", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_myocardialinfarctionYes", "Previous myocardial infarction", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_revascYes", "Previous cardiac revascularisation", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_ihdYes", "Other ischaemic heart disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_heartfailureYes", "Heart failure", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_padYes", "Peripheral arterial disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "recent_hosp_resp_infectYes", "Respiratory infection", row_name)) %>%
  mutate(row_name = ifelse(row_name == "recent_hosp_anythingYes", "Anything else", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_asthmaYes", "Asthma", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_copdYes", "Chronic obstructive pulmonary disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_bronchiectasisYes", "Bronchiectasis", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_pulmonaryfibrosisYes", "Pulmonary fibrosis", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_pulmonaryhypertensionYes", "Pulmonary hypertension", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_tiaYes", "Previous transient ischaemic attack", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_strokeYes", "Previous stroke", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_dementiaYes", "Dementia", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_otherneuroconditionsYes", "Other neurological condition", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_haem_cancerYes", "Haematological cancer", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_solid_cancerYes", "Solid cancer", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_solidorgantransplantYes", "Solid organ transplant", row_name)) %>%
  mutate(row_name = ifelse(row_name == "preindex_cldYes", "Chronic liver disease", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageMissing CKD", "Missing CKD", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 2", "Stage 2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 3a", "Stage 3a", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 3b", "Stage 3b", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 4", "Stage 4", row_name)) %>%
  mutate(row_name = ifelse(row_name == "ckd_stageStage 5", "Stage 5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "acr_catA2", "A2", row_name)) %>%
  mutate(row_name = ifelse(row_name == "acr_catA3", "A3", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "acr_catMissing ACR", "Missing ACR", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "treatment_6mInsulin", "Insulin prescription", row_name)) %>%
  mutate(row_name = ifelse(row_name == "treatment_6mOHA only", "OHA prescription only", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "oralsteroids_6mYes", "Oral steroids", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "immunosuppressants_6mYes", "Immunosuppressants", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "labas_6mYes", "Long acting beta agonists", row_name)) %>% 
  mutate(row_name = ifelse(row_name == "ltras_6mYes", "Leukotriene receptor antagonists", row_name))


coefs1 <- model_df_vaccinated %>% select(row_name, mean = exp.coef., lower = lower..95, upper = upper..95, p_value = Pr...z..)
coefs2 <- model_df_unvaccinated %>% select(row_name, mean = exp.coef., lower = lower..95, upper = upper..95, p_value = Pr...z..)

#Adding in reference categories
row_names_with_ref <- data.frame(row_name = c("Female", "Male", "<40", "40-49","50-59", "60-69", "70-79", "80-89", "90+", "White", "South Asian", "Black", "Other", "Mixed", "Unknown ethnicity", "1 (least deprived)", "2", "3", "4", "5 (most deprived)", "Missing IMD",
                                              "East Midlands", "East of England","London","Missing region", "North East", "North West", "South Central", "South East Coast", "South West", "West Midlands", "Yorkshire and the Humber",
                                              "<1", "1-2", "3-5", "6-9", "10-14", "15-19", "20+", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c", "0 complications", "1 complication", "2 complications", "3 complications",
                                              "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI", "Active smoker", "Ex-smoker", "Non-smoker", "Unknown smoking", "Hypertension", "Atrial fibrillation", "Angina",
                                              "Previous myocardial infarction", "Previous cardiac revascularisation", "Other ischaemic heart disease", "Heart failure", "Peripheral arterial disease",
                                              "Respiratory infection", "Anything else", "Asthma", "Chronic obstructive pulmonary disease", "Bronchiectasis", "Pulmonary fibrosis", "Pulmonary hypertension", "Previous transient ischaemic attack", "Previous stroke", "Dementia", "Other neurological condition",
                                              "Haematological cancer", "Solid cancer", "Solid organ transplant", "Chronic liver disease", "Stage 1", "Stage 2", "Stage 3a", "Stage 3b", "Stage 4", "Stage 5", "Missing CKD",
                                              "A1", "A2", "A3", "Missing ACR", "Insulin prescription", "OHA prescription only", "No treatment", "Oral steroids", "Immunosuppressants", "Long acting beta agonists", "Leukotriene receptor antagonists"))

coefs1 <- row_names_with_ref %>% left_join(coefs1) %>% mutate(mean = ifelse(is.na(mean), 1, mean), lower = ifelse(is.na(lower), 1, lower), upper = ifelse(is.na(upper), 1, upper))
coefs2 <- row_names_with_ref %>% left_join(coefs2) %>% mutate(mean = ifelse(is.na(mean), 1, mean), lower = ifelse(is.na(lower), 1, lower), upper = ifelse(is.na(upper), 1, upper))

#Adding category names
row_names_with_groups <- data.frame(row_name = c("Sex", "Female", "Male", "Age group, years", "<40", "40-49","50-59", "60-69", "70-79", "80-89", "90+", "Ethnicity", "White", "South Asian", "Black", "Other", "Mixed", "Unknown ethnicity", "Index of multiple deprivation quintile*", "1 (least deprived)", "2", "3", "4", "5 (most deprived)",
                                                 "Duration of diabetes, years", "<1", "1-2", "3-5", "6-9", "10-14", "15-19", "20+", "HbA1c, mmol/mol", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c", "Number of microvascular complications", "0 complications", "1 complication", "2 complications", "3 complications",
                                                 "BMI, kg/m2", "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI", "Smoking status", "Active smoker", "Ex-smoker", "Non-smoker", "Unknown smoking", 
                                                 "Comorbidities", "Cardiovascular", "Hypertension", "Atrial fibrillation", "Angina", "Previous myocardial infarction", "Previous cardiac revascularisation", "Other ischaemic heart disease", "Heart failure", "Peripheral arterial disease", 
                                                 "Recent hospitalisation", "Respiratory infection", "Anything else", "Respiratory", "Asthma", "Chronic obstructive pulmonary disease", "Neurological", "Previous transient ischaemic attack", "Previous stroke", "Dementia", "Other neurological condition", "Oncological",
                                                 "Haematological cancer", "Solid cancer", "Other comorbidities", "Solid organ transplant", "Chronic liver disease", "Chronic Kidney Disease (CKD) stage", "Stage 1", "Stage 2", "Stage 3a", "Stage 3b", "Stage 4", "Stage 5", "Missing CKD", "Albumin creatinine ratio (ACR) category", "A1", "A2", "A3", "Missing ACR", 
                                                 "Diabetes treatment (last 6 months)", "Insulin prescription", "OHA prescription only", "No treatment"))

coefs1 <- row_names_with_groups %>% left_join(coefs1)
coefs2 <- row_names_with_groups %>% left_join(coefs2)

#Combine
coefs <- coefs1 %>% left_join(coefs2, by = "row_name")

#Tidy
coefs <- coefs %>% mutate(row_name = ifelse(row_name == "Unknown ethnicity", "Unknown", ifelse(row_name == "Missing IMD", "Missing", ifelse(row_name == "Missing region", "Missing", ifelse(row_name == "0 complications", "0", ifelse(row_name == "1 complication", "1", ifelse(row_name == "2 complications", "2", ifelse(row_name == "3 complications", "3", ifelse(row_name == "Missing HbA1c", "Missing", ifelse(row_name == "Missing BMI", "Missing", ifelse(row_name == "Unknown smoking", "Unknown", ifelse(row_name == "Missing CKD", "Missing", ifelse(row_name == "Other comorbidities", "Other", ifelse(row_name == "Missing ACR", "Missing", row_name))))))))))))))


################################################################################
###Forest plot

#Define forest plot x-axis intervals
tick <- c(0.25, 0.5, 1, 2, 5)

#Generate forest plot
fp <- coefs %>% forestplot(labeltext = row_name,
                           mean = c(mean.x, mean.y),
                           lower = c(lower.x, lower.y),
                           upper = c(upper.x, upper.y),
                           hrzl_lines = gpar(col="#444444"),
                           col=fpColors(box=c("#b2df8a","#33a02c"), lines=c("#b2df8a","#33a02c"), zero = "gray50"), # this changes the box/ line colour- use when have multiple infections plotted
                           xticks = tick,
                           boxsize = .2,
                           graphwidth = unit(10, 'cm'),
                           zero = 1,
                           xlog = TRUE,
                           ci.vertices = TRUE, #makes ends T
                           clip = c(0.25,5),
                           new_page = TRUE,
                           fn.ci_norm = c(fpDrawCircleCI, fpDrawDiamondCI), #can use this to change shape of box
                           lty.ci =1, #this is confidence interval line type
                           txt_gp = fpTxtGp(label = gpar(cex = 0.7), ticks  = gpar(cex = 0.7), xlab = gpar(cex = 0.7), legend = gpar(cex = 0.7), summary = gpar(cex = 0.7)),
                           xlab = "Hazard ratio",
                           is.summary = c(TRUE, rep(FALSE, 2), TRUE, rep(FALSE,7), TRUE, rep(FALSE,6), TRUE, rep(FALSE,5), TRUE, rep(FALSE,7), TRUE, rep(FALSE,7), TRUE, rep(FALSE,4), TRUE, rep(FALSE,7), TRUE, rep(FALSE,4), TRUE, TRUE, rep(FALSE,8), TRUE, rep(FALSE,2), TRUE, rep(FALSE,2), TRUE, rep(FALSE,4), TRUE, rep(FALSE,2), TRUE, rep(FALSE,2), TRUE, rep(FALSE,7), TRUE, rep(FALSE, 4), TRUE, rep(FALSE,3)),
                           legend = c("Vaccinated", "Unvaccinated"),
                           legend_args = fpLegend(pos = list(x=.79, y=1, align = "horizontal"))
)
fp

#Save as pdf
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig11_T2_flu_hosp_vaccination.pdf",width=10,height=14)
fp
dev.off()

#END#
rm(list=ls())
