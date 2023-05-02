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


################################################################################
###COVID COHORT#################################################################
################################################################################

###Set first cohort and outcome#################################################
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
cohort <- cohort %>% filter(diabetes_type == "type 1")

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Check survival time distribution for those with outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")

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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
cohort$duration_cat <- relevel(cohort$duration_cat, ref = "10-14")
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
cohort$hba1c_cat <- relevel(cohort$hba1c_cat, ref = "48-53")
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
cohort$bmi_cat <- relevel(cohort$bmi_cat, ref = "25-29.9")

#Assign new name
cohort_covid <- cohort

################################################################################
###FLU COHORT###################################################################
################################################################################

###Set next cohort and outcome##################################################
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
cohort <- cohort %>% filter(diabetes_type == "type 1")

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Check survival time distribution for those with outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")

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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
cohort$duration_cat <- relevel(cohort$duration_cat, ref = "10-14")
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
cohort$hba1c_cat <- relevel(cohort$hba1c_cat, ref = "48-53")
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
cohort$bmi_cat <- relevel(cohort$bmi_cat, ref = "25-29.9")


#Assign new name
cohort_flu <- cohort

################################################################################
###PNEUMONIA COHORT#############################################################
################################################################################

###Set next cohort and outcome##################################################
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
cohort <- cohort %>% filter(diabetes_type == "type 1")

#Exclude those with cystic fibrosis
cohort <- cohort %>% filter(is.na(cysticfibrosis_diag_date))

##Exclude people registered <1 year
index.date.minus1y <- index.date - years(1)
cohort <- cohort %>% filter(regstartdate <= index.date.minus1y)

#Exclude people with diabetes diagnosed during study
cohort <- cohort %>% filter(dm_diag_date <= index.date)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Check survival time distribution for those with outcome
cohort %>% filter(outcome ==1) %>% ggplot(aes(survival_time)) +  geom_histogram(binwidth = 5, colour= "white", fill="plum3")

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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
cohort$duration_cat <- relevel(cohort$duration_cat, ref = "10-14")
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
cohort$hba1c_cat <- relevel(cohort$hba1c_cat, ref = "48-53")
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)
cohort$bmi_cat <- relevel(cohort$bmi_cat, ref = "25-29.9")


#Assign new name
cohort_pneumo <- cohort

################################################################################
###RUN MODELS AND BRING TOGETHER IN FOREST PLOT#################################
################################################################################

##Run models####################################################################
#Covid
model_covid <- coxph(Surv(survival_time, outcome) ~ gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat + bmi_cat, data = cohort_covid)
model_df_covid <- data.frame(summary(model_covid)$coefficients, summary(model_covid)$conf.int)

model_df_covid$row_name <- rownames(model_df_covid)

model_df_covid <- model_df_covid %>% mutate(row_name = ifelse(row_name == "genderMale", "Male", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat<18", "<18", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat18-39", "18-39", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat40-49", "40-49", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat50-59", "50-59", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat70-79", "70-79", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat80+", "80+", row_name)) %>%
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
  mutate(row_name = ifelse(row_name == "duration_cat<5", "<5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat15-19", "15-19", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat20-24", "20-24", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat25-29", "25-29", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat30+", "30+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat5-9", "5-9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat<48", "<48", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat53-64", "53-64", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat64-75", "64-75", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat75-86", "75-86", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat86+", "86+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_catMissing HbA1c", "Missing HbA1c", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat<18.5", "<18.5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat18.5-24.9", "18.5-24.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat30-34.9", "30-34.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat35-39.9", "35-39.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat40+", "40+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_catMissing BMI", "Missing BMI", row_name))

#Flu
model_flu <- coxph(Surv(survival_time, outcome) ~ gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat + bmi_cat, data = cohort_flu)
model_df_flu <- data.frame(summary(model_flu)$coefficients, summary(model_flu)$conf.int)

model_df_flu$row_name <- rownames(model_df_flu)

model_df_flu <- model_df_flu %>% mutate(row_name = ifelse(row_name == "genderMale", "Male", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat<18", "<18", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat18-39", "18-39", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat40-49", "40-49", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat50-59", "50-59", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat70-79", "70-79", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat80+", "80+", row_name)) %>%
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
  mutate(row_name = ifelse(row_name == "duration_cat<5", "<5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat15-19", "15-19", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat20-24", "20-24", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat25-29", "25-29", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat30+", "30+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat5-9", "5-9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat<48", "<48", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat53-64", "53-64", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat64-75", "64-75", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat75-86", "75-86", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat86+", "86+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_catMissing HbA1c", "Missing HbA1c", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat<18.5", "<18.5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat18.5-24.9", "18.5-24.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat30-34.9", "30-34.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat35-39.9", "35-39.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat40+", "40+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_catMissing BMI", "Missing BMI", row_name))

#Pneumonia
model_pneumo <- coxph(Surv(survival_time, outcome) ~ gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat + bmi_cat, data = cohort_pneumo)
model_df_pneumo <- data.frame(summary(model_pneumo)$coefficients, summary(model_pneumo)$conf.int)

model_df_pneumo$row_name <- rownames(model_df_pneumo)

model_df_pneumo <- model_df_pneumo %>% mutate(row_name = ifelse(row_name == "genderMale", "Male", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat<18", "<18", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat18-39", "18-39", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat40-49", "40-49", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat50-59", "50-59", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat70-79", "70-79", row_name)) %>%
  mutate(row_name = ifelse(row_name == "age_cat80+", "80+", row_name)) %>%
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
  mutate(row_name = ifelse(row_name == "duration_cat<5", "<5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat15-19", "15-19", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat20-24", "20-24", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat25-29", "25-29", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat30+", "30+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "duration_cat5-9", "5-9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat<48", "<48", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat53-64", "53-64", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat64-75", "64-75", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat75-86", "75-86", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_cat86+", "86+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "hba1c_catMissing HbA1c", "Missing HbA1c", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat<18.5", "<18.5", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat18.5-24.9", "18.5-24.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat30-34.9", "30-34.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat35-39.9", "35-39.9", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_cat40+", "40+", row_name)) %>%
  mutate(row_name = ifelse(row_name == "bmi_catMissing BMI", "Missing BMI", row_name))

coefs1 <- model_df_covid %>% select(row_name, mean = exp.coef., lower = lower..95, upper = upper..95, p_value = Pr...z..)
coefs2 <- model_df_flu %>% select(row_name, mean = exp.coef., lower = lower..95, upper = upper..95, p_value = Pr...z..)
coefs3 <- model_df_pneumo %>% select(row_name, mean = exp.coef., lower = lower..95, upper = upper..95, p_value = Pr...z..)

#Adding in reference categories
row_names_with_ref <- data.frame(row_name = c("Female", "Male", "<18", "18-39", "40-49","50-59", "60-69", "70-79", "80+", "White", "South Asian", "Black", "Other", "Mixed", "Unknown ethnicity", "1 (least deprived)", "2", "3", "4", "5 (most deprived)", "Missing IMD",
                                              "East Midlands", "East of England","London","Missing region", "North East", "North West", "South Central", "South East Coast", "South West", "West Midlands", "Yorkshire and the Humber",
                                              "<5", "5-9", "10-14", "15-19", "20-24", "25-29", "30+", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c", "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI"))

coefs1 <- row_names_with_ref %>% left_join(coefs1) %>% mutate(mean = ifelse(is.na(mean), 1, mean), lower = ifelse(is.na(lower), 1, lower), upper = ifelse(is.na(upper), 1, upper))
coefs2 <- row_names_with_ref %>% left_join(coefs2) %>% mutate(mean = ifelse(is.na(mean), 1, mean), lower = ifelse(is.na(lower), 1, lower), upper = ifelse(is.na(upper), 1, upper))
coefs3 <- row_names_with_ref %>% left_join(coefs3) %>% mutate(mean = ifelse(is.na(mean), 1, mean), lower = ifelse(is.na(lower), 1, lower), upper = ifelse(is.na(upper), 1, upper))

#Adding category names
row_names_with_groups <- data.frame(row_name = c("Sex", "Female", "Male", "Age group, years*", "<18", "18-39", "40-49","50-59", "60-69", "70-79", "80+", "Ethnicity**", "White", "South Asian", "Black", "Other", "Mixed", "Index of multiple deprivation quintile***", "1 (least deprived)", "2", "3", "4", "5 (most deprived)",
                                                 "Duration of diabetes, years", "<5", "5-9", "10-14", "15-19", "20-24", "25-29", "30+", "HbA1c, mmol/mol", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c", "BMI, kg/m2", "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI"))

coefs1 <- row_names_with_groups %>% left_join(coefs1)
coefs2 <- row_names_with_groups %>% left_join(coefs2)
coefs3 <- row_names_with_groups %>% left_join(coefs3)

#Combine
coefs <- coefs1 %>% left_join(coefs2, by = "row_name") %>% left_join(coefs3, by = "row_name")

#Tidy
coefs <- coefs %>% mutate(row_name = ifelse(row_name == "Unknown ethnicity", "Unknown", ifelse(row_name == "Missing IMD", "Missing", ifelse(row_name == "Missing HbA1c", "Missing", ifelse(row_name == "Missing BMI", "Missing", row_name)))))

coefs <- coefs %>% mutate(lower.y = ifelse(upper.y == "Inf", NA, lower.y), mean.y = ifelse(upper.y == "Inf", NA, mean.y), p_value.y = ifelse(upper.y == "Inf", NA, p_value.y), upper.y = ifelse(upper.y == "Inf", NA, upper.y)) %>%
  mutate(lower.x = ifelse(upper.x == "Inf", NA, lower.x), mean.x = ifelse(upper.x == "Inf", NA, mean.x), p_value.x = ifelse(upper.x == "Inf", NA, p_value.x), upper.x = ifelse(upper.x == "Inf", NA, upper.x))

################################################################################
###All 3 infection plot

#Define forest plot x-axis intervals
tick <- c(0.01, 0.1, 0.25, 0.1, 0.5, 1, 2, 5, 10)

#Generate forest plot
fp <- coefs %>% forestplot(labeltext = row_name,
                           mean = c(mean.x, mean.y, mean),
                           lower = c(lower.x, lower.y, lower),
                           upper = c(upper.x, upper.y, upper),
                           hrzl_lines = gpar(col="#444444"),
                           col=fpColors(box=c("#fc8d62","#66c2a5", "#8da0cb"), lines=c("#fc8d62","#66c2a5", "#8da0cb"), zero = "gray50"), # this changes the box/ line colour- use when have multiple infections plotted
                           xticks = tick,
                           boxsize = .2,
                           graphwidth = unit(10, 'cm'),
                           clip = c(0.01,10),
                           zero = 1,
                           xlog = TRUE,
                           ci.vertices = TRUE, #makes ends T
                           new_page = TRUE,
                           fn.ci_norm = c(fpDrawNormalCI, fpDrawCircleCI, fpDrawDiamondCI), #can use this to change shape of box
                           lty.ci =1, #this is confidence interval line type
                           txt_gp = fpTxtGp(label = gpar(cex = 0.7), ticks  = gpar(cex = 0.7), xlab = gpar(cex = 0.7), legend = gpar(cex = 0.7), summary = gpar(cex = 0.7)),
                           xlab = "Hazard ratio",
                           is.summary = c(TRUE, rep(FALSE, 2), TRUE, rep(FALSE,7), TRUE, rep(FALSE,5), TRUE, rep(FALSE,5), TRUE, rep(FALSE,7), TRUE, rep(FALSE,7), TRUE, rep(FALSE,7)),
                           legend = c("Covid", "Influenza", "Pneumonia"),
                           legend_args = fpLegend(pos = list(x=.75, y=1, align = "horizontal"))
)
fp


#Save these to pdf
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig7_sensitivity_T1_hosp_primary_diagnosis_forest_plot.pdf",width=10,height=7)
fp
dev.off()


#END#
rm(list=ls())



