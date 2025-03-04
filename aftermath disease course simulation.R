#### Cohort ####

source("aftermath disease course data analysis.R")

cohort_params <- list(
  N = 5000,
  symptom_duration_mean_reported = 17, 
  symptom_duration_sd_reported = 14, 
  symptom_duration_timescale = 1/4,
  symptom_underestimation_factor = 0.5, # making this up, will vary from 0 to ??
  recurrence_time_mean = 245,
  recurrence_time_sd = 175,
  incidence_18mo = 0.082,
  proportion_micro_pos = 0.04275/0.08055, 
  proportion_ever_subclinical = 3/4, # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = 1,
  duration_subclinical_cv = 1,
  auc = 0.7
)

attach(cohort_params)

cohort <- data.frame(ID = 1:N)


#### Ability to predict individual risk ####
# simulate an ROC curve with desired AUC (will be able to get data from Samyra instead, for prediction model with desired covariates)
# https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc, from Zelen and Severo Handbook of Mathematical Functions 1964
t <- sqrt(log(1/(1-auc)^2))
z <- t-((2.515517 + 0.802853*t + 0.0103328*t^2) / 
          (1 + 1.432788*t + 0.189269*t^2 + 0.001308*t^3))
d <- z*sqrt(2)

n <- 10000
x <- c(rnorm(n/2, mean = 0), rnorm(n/2, mean = d))
y <- c(rep(0, n/2), rep(1, n/2))

roc.curve <- roc(y, x)



#### Symptomatic TB onset ####

symptom_onset_mean = recurrence_time_mean - symptom_duration_mean_reported/symptom_underestimation_factor

# timing
scale <- recurrence_time_sd^2/symptom_onset_mean
shape <- symptom_onset_mean^2/recurrence_time_sd^2

# total events
incidence_18mo_micropos <- incidence_18mo * proportion_micro_pos
# 8.2% within 18 months; extend out indefinitely
incidence_total = 0.082 / pgamma(18*30, shape = shape, scale = scale)

cohort$TB <- rbinom(prob = incidence_total, size = 1, n = N)
cohort$symptom_onset [cohort$TB == 1] <- rgamma(n = sum(cohort$TB == 1), shape = shape, scale = scale)

# sputum+ pulmonary?
cohort <- cohort %>% mutate(pulmonary_with_micro = case_when(
  TB == 1 ~ rbinom(prob = incidence_18mo_micropos/incidence_18mo, size = 1, n = n()),
  TRUE ~ 0))

#### Timing of routine diagnosis ####
# Use negative binomial to increase variance (trial data mean duration 17 days, sd 14 days)
# And scale by (t/120)^1/4 for mean and var 
# Size = (mu^2)/(sd^2 - mu), varies by same fator as mean and var.

cohort <- cohort %>% mutate(diagnosis_routine = case_when(
  TB == 1 ~ symptom_onset + rnbinom(n = n(), 
                                    size = symptom_duration_mean_reported^2/(symptom_duration_sd_reported^2) * 
                                      (symptom_onset/120)^(symptom_duration_timescale),
                                    mu = symptom_duration_mean_reported/symptom_underestimation_factor) * 
                                        (symptom_onset/120)^(symptom_duration_timescale),
  TRUE ~ NA_real_))

cohort %>% filter(TB==1) %>% ggplot(aes(x=symptom_onset, y=diagnosis_routine)) + geom_point()
cohort %>% filter(TB==1) %>% ggplot(aes(x=symptom_onset, y=diagnosis_routine - symptom_onset)) + geom_point()
cohort %>% filter(TB==1, diagnosis_routine <= 30*18) %>% summarise(mean(diagnosis_routine - symptom_onset), sd(diagnosis_routine - symptom_onset))

#### Subclinical TB ####

mean_duration_subclinical <- symptom_duration_mean_reported/symptom_underestimation_factor*
  duration_ratio_subclinical_symptomatic/proportion_ever_subclinical

cohort <- cohort %>% mutate(
  ever_subclinical_sputumpos = case_when(pulmonary_with_micro == 1 ~ rbinom(prob = proportion_ever_subclinical, size = 1, n = n()),
                                      TRUE ~ 0),
  # In these people, sputum onset precedes symptom onset by a gamma-distributed time
  sputum_onset = case_when(ever_subclinical_sputumpos == 1 ~ 
                             pmax(0, symptom_onset - mean_duration_subclinical* rgamma(scale = duration_subclinical_cv^2, 
                                                                                       shape = 1/(duration_subclinical_cv^2), n = n())),
                           ever_subclinical_sputumpos == 0 & pulmonary_with_micro == 1 ~ 
                             runif(n = n(), min = symptom_onset, max = diagnosis_routine),
                           TRUE ~ NA_real_))

# Check the number who are subclinical sputum+ from the start:

cohort %>% filter(TB==1) %>% summarise(mean(sputum_onset == 0, na.rm = TRUE), n()) # around 5% of all recurrences
cohort %>% filter(TB==1) %>% summarise(sum(sputum_onset == 0, na.rm = TRUE)/N, n()) # around 0.3% of the overall cohort -- both c/w crtriumph (or a little low).

#### Prediction ####

# assign risk rankings based on sensitivity and specificity of roc.curve for TB:
# First, assign scores to the TB cases such that a given prediction score P captures P% of the TB cases
# And then, for the non-cases, assign with probability (1-spec) at a 
cohort$TB_risk_rank[cohort$TB == 1] <- (1:sum(cohort$TB == 1))/sum(cohort$TB == 1)
cohort$nonTB_risk_rank[cohort$TB == 0] <- (1:sum(cohort$TB == 0))/sum(cohort$TB == 0)
for (i in 1:N)
{TB_status <- cohort$TB[i]
  if (TB_status == 1)
  {cohort$risk_score[i] <- cohort$TB_risk_rank[i]}
  else
  {cohort$risk_score[i] <- 1 - roc.curve$specificities[which.min(roc.curve$sensitivities > cohort$nonTB_risk_rank[i])]}
}
cohort$risk_score

# checking:
cohort %>% group_by(TB) %>% summarise(mean(risk_score), mean(risk_score > 0.9), mean(risk_score > 0.5), mean(risk_score > 0.1))



#### Examine outputs #### 
mean(cohort$TB)

colors <- c("Symptomatic" = "blue", "Sputum+" = "red")
fillcolors <- c("Micro+ Pulmonary" = "red", "Clinical or extrapulmonary" = "gray")
ggplot(cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
         mutate(newID = row_number(),
                TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                                   TRUE ~ "Clinical or extrapulmonary"))) + 
  geom_segment(aes(x = sputum_onset, xend =  diagnosis_routine, y = newID - 0.2, col='Sputum+'), size=0.8) + 
  geom_segment(aes(x = symptom_onset, xend = diagnosis_routine, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
  scale_color_manual("Time with TB", values = colors) + 
  scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
  geom_point(aes(x = diagnosis_routine, y = newID, fill = TBtype), pch=21) + 
  ylab("Cases arranged by timing of routine diagnosis") + 
  xlab("Days since prior treatment completion") +
  ggtitle("Timing of TB recurrence detection (routine)\nand preceding symptomatic and/or infectious time") + 
  xlim(0,30*24)


# Plot sputum+ only:
ggplot(cohort %>% filter(TB==1, pulmonary_with_micro == 1) %>% arrange(diagnosis_routine) %>% 
         mutate(newID = row_number(),
                TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                                   TRUE ~ "Clinical or extrapulmonary"))) + 
  geom_segment(aes(x = sputum_onset, xend =  diagnosis_routine, y = newID - 0.2, col='Sputum+'), size=0.8) + 
  geom_segment(aes(x = symptom_onset, xend = diagnosis_routine, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
  scale_color_manual("Time with TB", values = colors) + 
  scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
  geom_point(aes(x = diagnosis_routine, y = newID, fill = TBtype), pch=21) + 
  ylab("Cases arranged by timing of routine diagnosis") + 
  xlab("Days since prior treatment completion") +
  ggtitle("Timing of TB recurrence detection (routine)\nand preceding symptomatic and/or infectious time") + 
  xlim(0,30*24)

table(cohort$TB, cohort$risk_score > 0.5); prop.table(table(cohort$TB, cohort$risk_score > 0.5), margin=1)
table(cohort$TB, cohort$risk_score > 0.9); prop.table(table(cohort$TB, cohort$risk_score > 0.9), margin=1)

