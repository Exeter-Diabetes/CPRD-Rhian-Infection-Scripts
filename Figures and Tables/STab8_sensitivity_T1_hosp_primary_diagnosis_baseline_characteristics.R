################################################################################
###Descriptive tables###########################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(tableone)
library(kableExtra)

#Load aurum package
library(aurum)

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")

#Setting up/loading analysis test
analysis = cprd$analysis("Rhian_covid")

###COVID########################################################################
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
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Setting variables to factors
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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)

#Generate table one for whole cohort
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

categorical_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1) 

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_covid <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename ("2020 cohort" = All, "Covid-19 hospitalisations" = Hospitalisations)

###INFLUENZA####################################################################
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
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Setting variables to factors
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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)

#Generate table one for whole cohort
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

categorical_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1)

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_flu <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename("2016 cohort" = All, "Influenza hospitalisations" = Hospitalisations)


###PNEUMONIA####################################################################
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
cohort <- cohort %>% filter(dm_diag_date_all <= index.date)

#Mean follow-up
mean(cohort$survival_time)
sd(cohort$survival_time)

#Edit age groups
cohort <- cohort %>% mutate(age_cat = ifelse(age_at_index<18, "<18", ifelse(age_at_index <40 & age_at_index>=18, "18-39", ifelse(age_at_index<50 & age_at_index>=40, "40-49", ifelse(age_at_index<60 & age_at_index>=50, "50-59", ifelse(age_at_index<70 & age_at_index>=60, "60-69", ifelse(age_at_index<80 & age_at_index>=70, "70-79", ifelse(age_at_index>=80, "80+", NA))))))))
#And duration
cohort <- cohort %>% mutate(duration_cat = ifelse(dm_duration_at_index <5, "<5", ifelse(dm_duration_at_index <10 & dm_duration_at_index >=5, "5-9", ifelse(dm_duration_at_index <15 & dm_duration_at_index>=10, "10-14", ifelse(dm_duration_at_index <20 & dm_duration_at_index>=15, "15-19", ifelse(dm_duration_at_index<25 & dm_duration_at_index>=20, "20-24", ifelse(dm_duration_at_index <30 & dm_duration_at_index>=25, "25-29", ifelse(dm_duration_at_index>=30, "30+", NA))))))))

#Setting variables to factors
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
#Diabetes duration
cohort$duration_cat <- factor(cohort$duration_cat)
#HbA1c
cohort$hba1c_cat <- factor(cohort$hba1c_cat)
#BMI
cohort$bmi_cat <- factor(cohort$bmi_cat)

#Generate table one for whole cohort
all_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

categorical_vars <- c("gender", "femalegender", "age_cat", "eth5", "imd_quintile", "duration_cat", "hba1c_cat", "bmi_cat")

tableone1 <- CreateTableOne(vars=all_vars,data=cohort,factorVars=categorical_vars, test=FALSE)

tabprint1 <-as_tibble(print(tableone1)) %>%
  add_column(measure=row.names(print(tableone1))) %>%
  mutate(measure = trimws(measure))


#Generate table one for just those with outcome
outcome_yes <- cohort %>% filter(outcome ==1)

tableone2 <- CreateTableOne(vars=all_vars,data=outcome_yes,factorVars=categorical_vars, test=FALSE)

tabprint2 <-as_tibble(print(tableone2)) %>%
  add_column(measure=row.names(print(tableone2))) %>%
  mutate(measure = trimws(measure))


#Combine
tabprint_pneumo <- tabprint1 %>% rename(All = Overall) %>% left_join(tabprint2, by = "measure") %>% select(measure, All, Hospitalisations = Overall) %>% filter(All != "") %>% mutate(Hospitalisations = ifelse(is.na(Hospitalisations), "0 ( 0.0)", Hospitalisations)) %>% rename("2016 cohort" = All, "Pneumonia hospitalisations" = Hospitalisations)

################################################################################
##Add together
tabprint <- tabprint_covid %>% left_join(tabprint_flu) %>% left_join(tabprint_pneumo)

row_names <- data.frame(measure = c("n", "femalegender = Female (%)", "gender = Male (%)", "<18", "18-39","40-49", "50-59", "60-69", "70-79", "80+", "White", "South Asian", "Black", "Other", "Mixed", "Unknown",
                                    "1", "2", "3", "4", "5", "Missing IMD", "<5", "5-9", "10-14", "15-19", "20-24", "25-29", "30+", "<48", "48-53", "53-64", "64-75", "75-86", "86+", "Missing HbA1c",
                                    "<18.5", "18.5-24.9", "25-29.9", "30-34.9", "35-39.9", "40+", "Missing BMI"))

tabprint <- row_names %>% left_join(tabprint, by = "measure")

tab <- tabprint %>% 
  mutate(measure=ifelse(measure=="gender = Male (%)","Male", measure)) %>%
  mutate(measure=ifelse(measure=="femalegender = Female (%)","Female", measure)) %>%
  mutate(measure = ifelse(measure == "1", "1 (least deprived)", measure)) %>%
  mutate(measure = ifelse(measure == "5", "5 (most deprived)", measure)) %>%
  mutate(measure=ifelse(measure=="Missing IMD","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="Missing HbA1c","Missing", measure)) %>%
  mutate(measure=ifelse(measure=="Missing BMI","Missing", measure)) %>%
  rename(" " = measure)

#Outputting HTML table
kableExtra::kable(tab,format="html", align="lll") %>% 
  kable_styling(bootstrap_options = c("striped","hover"), row_label_position = "lll") %>%
  pack_rows("Sex", 2,3, bold = TRUE) %>%
  pack_rows("Age, years",4,10, bold=TRUE) %>%  
  pack_rows("Ethnicity",11,16, bold=TRUE)  %>%  
  pack_rows("Index of multiple deprivation quintile",17,22, bold=TRUE)  %>% 
  pack_rows("Duration of diagnosed diabetes, years",23,29, bold=TRUE)  %>%  
  pack_rows("HbA1c, mmol/mol",30,36, bold=TRUE)  %>% 
  pack_rows("BMI, kg/m2",37,43, bold=TRUE)  %>%
  column_spec(1,width="12cm") %>%
  column_spec(2,width="6cm") %>%
  column_spec(3,width="7cm") %>%
  column_spec(4,width="6cm") %>%
  column_spec(5,width="7cm") %>%
  column_spec(5,width="7cm") %>%
  cat(.,file="STab8_sensitivity_T1_hosp_primary_diagnosis_baseline_characteristics.html")

###END##########################################################################
rm(list=ls())
