################################################################################
###ANALYSIS#####################################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(survival)
library(rms)
library(patchwork)
library(cowplot)

#Load aurum package
library(aurum)

###Connecting to data and setting up/connecting to analysis#####################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")

#Connect to analysis
analysis = cprd$analysis("Rhian_covid")


###COVID########################################################################
cohort.name <- "feb2020"
infection <- "covid"
outcome <- "hosp" #these will be pasted into output file names
infection.name <- "Covid-19"
#Set cohort
cohort <- cohort %>% analysis$cached(paste0(cohort.name, "_", infection, "_outcomes"), unique_indexes = "patid", indexes = "hosp_date")
#Only want patids with linkage
linkage_patids <- cprd$tables$patidsWithLinkage %>% filter(n_patid_hes <=20) %>% select(patid)
cohort <- cohort %>% inner_join(linkage_patids)
#Set outcome and outcome date variable
cohort <- cohort %>% mutate (outcome = hospitalisation_outcome, outcome_date = hosp_date)
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


################################################################################
#Set subgroups
subgroup1 <- "White"
subgroup2 <- "South Asian"
subgroup3 <- "Black"
cohort <- cohort %>% mutate(subgroup = eth5)
################################################################################

##Set reference to 30 for subgroup and then run model with interaction for each subgroup
dd <- datadist(cohort %>% filter(subgroup=="White") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_1 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_1)

predict_1 <- Predict(model_1, bmi_value = seq(18,50, by =1), subgroup ="White", ref.zero=TRUE, fun = exp)
predict_1_df <- as.data.frame(predict_1) %>% mutate(subgroup = paste0(subgroup1))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="South Asian") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_2 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_2)

predict_2 <- Predict(model_2, bmi_value = seq(18,50, by =1), subgroup ="South Asian", ref.zero=TRUE, fun = exp)
predict_2_df <- as.data.frame(predict_2) %>% mutate(subgroup = paste0(subgroup2))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="Black") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_3 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_3)

predict_3 <- Predict(model_3, bmi_value = seq(18,50, by =1), subgroup ="Black", ref.zero=TRUE, fun = exp)
predict_3_df <- as.data.frame(predict_3) %>% mutate(subgroup = paste0(subgroup3))

#Combine
subgroup_df <- predict_1_df %>% rbind(predict_2_df) %>% rbind(predict_3_df)

#Set values for 95%CI outside plot range to limit values
subgroup_df <- subgroup_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>3,3, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>3,3, upper)))

#Plot
plot<- ggplot(data=subgroup_df,aes(x=bmi_value, y=yhat, group = subgroup)) +
  geom_line(data=subgroup_df,aes(x=bmi_value, y=yhat, color = subgroup), size=0.75) +
  geom_ribbon(data=subgroup_df,aes(x=bmi_value,ymin=lower,ymax=upper, group = subgroup, fill = subgroup),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,3,0.5), limits = c(0.5,3)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = "Ethnicity", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size=10), legend.position = c(0.2,0.85), legend.background = element_blank())

#Paste subgroup names
cohort <- cohort %>% mutate(subgroup = ifelse(subgroup == "White", "White", ifelse(subgroup == "Black", "Black", ifelse(subgroup == "South Asian", "South Asian", NA))))


#Density line
density <- ggplot((cohort %>% filter(!is.na(subgroup))), aes(x = bmi_value, colour = subgroup)) +geom_density(position = "identity", size = 0.75) +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  theme(legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,1,1,1,"cm"),
        axis.line.x = element_line(color="grey50"))
density

plot <- plot_grid(plot, density, ncol = 1,align = "v", axis = "lr",
                  rel_heights = c(1,0.4), rel_widths = c(1,1))


covid_eth_plot <- plot
covid_eth_plot


###FLU##########################################################################
cohort.name <- "sep2016"
infection <- "influenza"
outcome <- "hosp" #these will be pasted into output file names
infection.name <- "Influenza"
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


################################################################################
#Set subgroups
subgroup1 <- "White"
subgroup2 <- "South Asian"
subgroup3 <- "Black"
cohort <- cohort %>% mutate(subgroup = eth5)
################################################################################

##Set reference to 30 for subgroup and then run model with interaction for each subgroup
dd <- datadist(cohort %>% filter(subgroup=="White") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_1 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,3) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_1)

predict_1 <- Predict(model_1, bmi_value = seq(18,50, by =1), subgroup ="White", ref.zero=TRUE, fun = exp)
predict_1_df <- as.data.frame(predict_1) %>% mutate(subgroup = paste0(subgroup1))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="South Asian") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_2 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,3) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_2)

predict_2 <- Predict(model_2, bmi_value = seq(18,50, by =1), subgroup ="South Asian", ref.zero=TRUE, fun = exp)
predict_2_df <- as.data.frame(predict_2) %>% mutate(subgroup = paste0(subgroup2))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="Black") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_3 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,3) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_3)

predict_3 <- Predict(model_3, bmi_value = seq(18,50, by =1), subgroup ="Black", ref.zero=TRUE, fun = exp)
predict_3_df <- as.data.frame(predict_3) %>% mutate(subgroup = paste0(subgroup3))

#Combine
subgroup_df <- predict_1_df %>% rbind(predict_2_df) %>% rbind(predict_3_df)

#Set values for 95%CI outside plot range to limit values
subgroup_df <- subgroup_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>3,3, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>3,3, upper)))

#Plot
plot<- ggplot(data=subgroup_df,aes(x=bmi_value, y=yhat, group = subgroup)) +
  geom_line(data=subgroup_df,aes(x=bmi_value, y=yhat, color = subgroup), size=0.75) +
  geom_ribbon(data=subgroup_df,aes(x=bmi_value,ymin=lower,ymax=upper, group = subgroup, fill = subgroup),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,3,0.5), limits = c(0.5,3)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = " ", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size=10), legend.position = c(0.2,0.85), legend.background = element_blank())

#Paste subgroup names
cohort <- cohort %>% mutate(subgroup = ifelse(subgroup == "White", "White", ifelse(subgroup == "Black", "Black", ifelse(subgroup == "South Asian", "South Asian", NA))))


#Density line
density <- ggplot((cohort %>% filter(!is.na(subgroup))), aes(x = bmi_value, colour = subgroup)) +geom_density(position = "identity", size = 0.75) +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  theme(legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,1,1,1,"cm"),
        axis.line.x = element_line(color="grey50"))
density

plot <- plot_grid(plot, density, ncol = 1,align = "v", axis = "lr",
                  rel_heights = c(1,0.4), rel_widths = c(1,1))

flu_eth_plot <- plot
flu_eth_plot

###PNEUMONIA####################################################################
cohort.name <- "sep2016"
infection <- "pneumonia"
outcome <- "hosp" #these will be pasted into output file names
infection.name <- "Pneumonia"
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


################################################################################
#Set subgroups
subgroup1 <- "White"
subgroup2 <- "South Asian"
subgroup3 <- "Black"
cohort <- cohort %>% mutate(subgroup = eth5)
################################################################################

##Set reference to 30 for subgroup and then run model with interaction for each subgroup
dd <- datadist(cohort %>% filter(subgroup=="White") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_1 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_1)

predict_1 <- Predict(model_1, bmi_value = seq(18,50, by =1), subgroup ="White", ref.zero=TRUE, fun = exp)
predict_1_df <- as.data.frame(predict_1) %>% mutate(subgroup = paste0(subgroup1))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="South Asian") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_2 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_2)

predict_2 <- Predict(model_2, bmi_value = seq(18,50, by =1), subgroup ="South Asian", ref.zero=TRUE, fun = exp)
predict_2_df <- as.data.frame(predict_2) %>% mutate(subgroup = paste0(subgroup2))

#Repeat for other subgroups
dd <- datadist(cohort %>% filter(subgroup=="Black") %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 30
model_3 <- cph(Surv(survival_time,outcome) ~ subgroup * rcs(bmi_value,5) + gender + age_cat + imd_quintile + region + duration_cat + hba1c_cat + number_complications + smoking_status +
                 preindex_hypertension + preindex_af + preindex_angina + preindex_myocardialinfarction + preindex_revasc + preindex_ihd + preindex_heartfailure + preindex_pad +
                 recent_hosp_resp_infect + recent_hosp_anything + preindex_asthma + preindex_copd + preindex_bronchiectasis + preindex_pulmonaryfibrosis + preindex_pulmonaryhypertension + 
                 preindex_tia + preindex_stroke + preindex_dementia + preindex_otherneuroconditions + preindex_haem_cancer + preindex_solid_cancer + preindex_solidorgantransplant + preindex_cld +
                 ckd_stage + acr_cat + treatment_6m + oralsteroids_6m + immunosuppressants_6m + labas_6m + ltras_6m, data = cohort, x= T, y=T)

anova(model_3)

predict_3 <- Predict(model_3, bmi_value = seq(18,50, by =1), subgroup ="Black", ref.zero=TRUE, fun = exp)
predict_3_df <- as.data.frame(predict_3) %>% mutate(subgroup = paste0(subgroup3))

#Combine
subgroup_df <- predict_1_df %>% rbind(predict_2_df) %>% rbind(predict_3_df)

#Set values for 95%CI outside plot range to limit values
subgroup_df <- subgroup_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>3,3, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>3,3, upper)))

#Plot
plot<- ggplot(data=subgroup_df,aes(x=bmi_value, y=yhat, group = subgroup)) +
  geom_line(data=subgroup_df,aes(x=bmi_value, y=yhat, color = subgroup), size=0.75) +
  geom_ribbon(data=subgroup_df,aes(x=bmi_value,ymin=lower,ymax=upper, group = subgroup, fill = subgroup),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,3,0.5), limits = c(0.5,3)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = " ", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) +
  theme(legend.title = element_blank(), legend.text = element_text(size=10), legend.position = c(0.2,0.85), legend.background = element_blank())

#Paste subgroup names
cohort <- cohort %>% mutate(subgroup = ifelse(subgroup == "White", "White", ifelse(subgroup == "Black", "Black", ifelse(subgroup == "South Asian", "South Asian", NA))))


#Density line
density <- ggplot((cohort %>% filter(!is.na(subgroup))), aes(x = bmi_value, colour = subgroup)) +geom_density(position = "identity", size = 0.75) +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  theme(legend.title = element_blank(), panel.background = element_blank(), legend.position = "none") +
  theme(axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin = margin(0,1,1,1,"cm"),
        axis.line.x = element_line(color="grey50"))
density

plot <- plot_grid(plot, density, ncol = 1,align = "v", axis = "lr",
                  rel_heights = c(1,0.4), rel_widths = c(1,1))

pneumo_eth_plot <- plot
pneumo_eth_plot

################################################################################
#Combine all plots
ethnicity_plot <- covid_eth_plot + flu_eth_plot + pneumo_eth_plot

#Save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig3b_T2_hosp_BMI_splines_ethnicity.pdf",width=16,height=6)
ethnicity_plot
dev.off()


#END#
rm(list=ls())


