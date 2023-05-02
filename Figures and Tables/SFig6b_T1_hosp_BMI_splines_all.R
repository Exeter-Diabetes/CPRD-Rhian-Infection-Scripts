################################################################################
###ANALYSIS#####################################################################
################################################################################

.libPaths("C:/Users/rh530/OneDrive - University of Exeter/R/win-library/4.1")

#Load packages
library(tidyverse)
library(lubridate)
library(survminer)
library(survival)
library(patchwork)
library(rms)
library(cowplot)

#Load aurum package
library(aurum)

###Connecting to data and setting up analysis###################################
#Initialise connection
cprd = CPRDData$new(cprdEnv = "test-remote",cprdConf = "C:/Users/rh530/.aurum.yaml")
codesets = cprd$codesets()
codes = codesets$getAllCodeSetVersion(v = "31/10/2021")

#Setting up/loading analysis test
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


################################################################################

##Set reference to 26 (around median)
dd <- datadist(cohort %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 26
model <- cph(Surv(survival_time,outcome) ~ rcs(bmi_value,5) + gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat, data = cohort, x= T, y=T)
predict <- Predict(model, bmi_value = seq(18,50, by =1), ref.zero=TRUE, fun = exp)
predict_df <- as.data.frame(predict)

#Set values for 95%CI outside range to limit values
predict_df <- predict_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>4,4, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>4,4, upper)))

#Plot
plot<- ggplot(data=predict_df,aes(x=bmi_value, y=yhat)) +
  geom_line(data=predict_df,aes(x=bmi_value, y=yhat), size=0.75) +
  geom_ribbon(data=predict_df,aes(x=bmi_value,ymin=lower,ymax=upper),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,4,0.5), limits = c(0.5,4)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = "All", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) 

#Density line
density <- ggplot(cohort, aes(x = bmi_value)) +geom_density(position = "identity", size = 0.75) +
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

covid_plot <- plot

covid_plot

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


################################################################################

##Set reference to 26
dd <- datadist(cohort %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 26
model <- cph(Surv(survival_time,outcome) ~ rcs(bmi_value,3) + gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat, data = cohort, x= T, y=T)
predict <- Predict(model, bmi_value = seq(18,50, by =1), ref.zero=TRUE, fun = exp)
predict_df <- as.data.frame(predict)

#Set values for 95%CI outside range to limit values
predict_df <- predict_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>4,4, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>4,4, upper)))

#Plot
plot<- ggplot(data=predict_df,aes(x=bmi_value, y=yhat)) +
  geom_line(data=predict_df,aes(x=bmi_value, y=yhat), size=0.75) +
  geom_ribbon(data=predict_df,aes(x=bmi_value,ymin=lower,ymax=upper),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,4,0.5), limits = c(0.5,4)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = " ", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) 

#Desnity line
density <- ggplot(cohort, aes(x = bmi_value)) +geom_density(position = "identity", size = 0.75) +
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

flu_plot <- plot


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


################################################################################

##Set reference to 26
dd <- datadist(cohort %>% select(!cysticfibrosis_diag_date))
options(datadist="dd")
dd$limits["Adjust to","bmi_value"] <- 26
model <- cph(Surv(survival_time,outcome) ~ rcs(bmi_value,5) + gender + age_cat + eth5 + imd_quintile + region + duration_cat + hba1c_cat, data = cohort, x= T, y=T)
predict <- Predict(model, bmi_value = seq(18,50, by =1), ref.zero=TRUE, fun = exp)
predict_df <- as.data.frame(predict)

#Set values for 95%CI outside range to limit values
predict_df <- predict_df %>% mutate(lower = ifelse(lower<0.5, 0.5, ifelse(lower>4,4, lower)), upper = ifelse(upper<0.5, 0.5, ifelse(upper>4,4, upper)))

#Plot
plot<- ggplot(data=predict_df,aes(x=bmi_value, y=yhat)) +
  geom_line(data=predict_df,aes(x=bmi_value, y=yhat), size=0.75) +
  geom_ribbon(data=predict_df,aes(x=bmi_value,ymin=lower,ymax=upper),alpha=0.2) + 
  geom_hline(yintercept =1, linetype = "dashed") +
  scale_x_continuous(breaks = seq(20,50,10), limits = c(18,50)) +
  scale_y_continuous(breaks = seq(0.5,4,0.5), limits = c(0.5,4)) +
  labs(x = "BMI (kg/m2)", y = "Hazard ratio", title = " ", subtitle = paste0(infection.name)) +
  theme(axis.line = element_line(colour =  "grey50" ),
        plot.title = element_text(size=10, face = "bold"), plot.subtitle = element_text(size=10), plot.margin = margin(1,1,0,1,"cm"),
        axis.text = element_text(size = 10), axis.title = element_text(size = 10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank()) 

#Density line
density <- ggplot(cohort, aes(x = bmi_value)) +geom_density(position = "identity", size = 0.75) +
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

pneumo_plot <- plot

################################################################################
all_plot <- covid_plot + flu_plot + pneumo_plot
################################################################################


#Plot next to each other and save
pdf.options(reset = TRUE, onefile = TRUE)
pdf("SFig6b_T1_hosp_BMI_splines_all.pdf",width=16,height=6)
all_plot
dev.off()


#END#
rm(list=ls())
