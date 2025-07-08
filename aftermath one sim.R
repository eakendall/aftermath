#### Cohort ####

source("aftermath cohort sim functions.R")

cohort_params <- list(
  N = 1000,
  
  recurrence_time_mean = 192, # diagnosis minus estimated symptom duration, among those diagnosed by 18mo
  recurrence_time_shapefactor = 1,
  
  incidence_18mo = 0.092,
  proportion_micro_pos = 0.04392/0.0822,
  auc = 0.69,

  symptom_duration_mean_reported = 17, 
  symptom_duration_sd_reported = 14*1.2, # Increased >1x; long tails perhaps more common in trial nonparticipants/losses, or perhaps underreported; will explore sensitivity
  symptom_duration_timescale = 0,
  symptom_underestimation_factor = 0.33, # 1 means no underestimation, 0.33 means 3x longer than reported
  
  proportion_ever_subclinical = 3/4, # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = 1,
  duration_subclinical_cv = 1,
    subclinical_baseline_amongTB_max = 0.20, 
    subclinical_6m_amongcohort_min = 0.003, 
    subclinical_6m_amongcohort_max = 0.020, 
  
  coverage_phone = 0.7,
  coverage_home_reduction = 0.85,
  sensitivity_symptoms_home = 0.8,
  sensitivity_symptoms_phone_reduction = 0.65,
  success_sputum_home = 0.9,
  success_sputum_phone_reduction = 0.8,
  
  home_visit_passive_detection_impact = 1, # time by which symptomatic time is shortened after initial home visit, assumed.
  aftermath_counseling_passive_detection_impact = 1, #0.9 was assumed, less than the 20% reduction we'll model for intentional counseling
  intentional_counseling_passive_detection_impact = 1 #assumed 
  
)

cohort <- create_cohort(cohort_params)

#### Examine outputs #### 
mean(cohort$TB)
check_subclinical(cohort, cohort_params) # will reject sampled paramters if this is false

mean(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(diagnosis_routine) %>% unlist(.))
median(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(diagnosis_routine) %>% unlist(.))
sd(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(diagnosis_routine) %>% unlist(.))

mean(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(symptom_onset) %>% unlist(.))
median(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(symptom_onset) %>% unlist(.))
sd(cohort %>% filter(TB==1, diagnosis_routine <= 540) %>% select(symptom_onset) %>% unlist(.))

colors <- c("Symptomatic" = "blue", "Sputum+" = "red")
fillcolors <- c("Micro+ Pulmonary" = "red", "Clinical or extrapulmonary" = "gray")
ggplot(cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
         mutate(newID = row_number(),
                TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                                   TRUE ~ "Clinical or extrapulmonary"))) + 
  geom_segment(aes(x = sputum_onset/30, xend =  diagnosis_routine/30, y = newID - 0.2, col='Sputum+'), size=0.8) + 
  geom_segment(aes(x = symptom_onset/30, xend = diagnosis_routine/30, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
  scale_color_manual("Time with TB", values = colors) + 
  scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
  geom_point(aes(x = diagnosis_routine/30, y = newID, fill = TBtype), pch=21) + 
  ylab("Cases arranged by timing of routine diagnosis") + 
  xlab("Months since prior treatment completion") +
  # ggtitle("Timing of TB recurrence detection (routine)\nand preceding symptomatic and/or infectious time") + 
  xlim(0,30) + 
  theme_bw()

# # with black dots at "NAAT+"
# colors <- c("Symptomatic" = "blue", "NAAT+" = "red")
# ggplot(cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
#          mutate(newID = row_number(),
#                 TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
#                                    TRUE ~ "Clinical or extrapulmonary"))) + 
#   geom_segment(aes(x = sputum_onset/30, xend =  diagnosis_routine/30, y = newID - 0.2, col='NAAT+'), size=0.8) + 
#   geom_segment(aes(x = symptom_onset/30, xend = diagnosis_routine/30, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
#   scale_color_manual("Time with TB", values = colors) + 
#   geom_point(aes(x = diagnosis_routine/30, y = newID, fill = TBtype), pch=21) + 
#   ylab("Cases arranged by timing of routine diagnosis") + 
#   xlab("Months since prior treatment completion") +
#   # ggtitle("Timing of TB recurrence detection (routine)\nand preceding symptomatic and/or infectious time") + 
#   # x axis breaks at 6, 12, 18, 24 mo
#   scale_x_continuous(breaks = seq(0, 30, 6), limits = c(0,30)) +
#   theme_bw()


# # Plot sputum+ only:
# ggplot(cohort %>% filter(TB==1, pulmonary_with_micro == 1) %>% arrange(diagnosis_routine) %>% 
#          mutate(newID = row_number(),
#                 TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
#                                    TRUE ~ "Clinical or extrapulmonary"))) + 
#   geom_segment(aes(x = sputum_onset, xend =  diagnosis_routine, y = newID - 0.2, col='NAAT+'), size=0.8) + 
#   geom_segment(aes(x = symptom_onset, xend = diagnosis_routine, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
#   scale_color_manual("Time with TB", values = colors) + 
#   scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
#   geom_point(aes(x = diagnosis_routine, y = newID), pch=21) + 
#   ylab("Cases arranged by timing of routine diagnosis") + 
#   xlab("Days since prior treatment completion") +
#   ggtitle("Timing of TB recurrence detection (routine)\nand preceding symptomatic and/or infectious time") + 
#   xlim(0,30*30)

table(cohort$TB, cohort$risk_score > 0.5); prop.table(table(cohort$TB, cohort$risk_score > 0.5), margin=1)
table(cohort$TB, cohort$risk_score > 0.9); prop.table(table(cohort$TB, cohort$risk_score > 0.9), margin=1)

table(cohort$TB, cohort$risk_quantile > 0.5); prop.table(table(cohort$TB, cohort$risk_quantile > 0.5), margin=1)
table(cohort$TB, cohort$risk_quantile > 0.9); prop.table(table(cohort$TB, cohort$risk_quantile > 0.9), margin=1)




# plot an illustartive cohort with a screening timepoint

design_illustrative = data.frame(
  "counseling_coverage" = 0, 
  "timing_months" = c(5),
  "target_coverage" = c(1),
  "screening_method" = c("symptoms"),
  "screening_location" = c("home")
) %>% 
  arrange(timing_months)

intervention_parameters <- list(
  coverage = list("phone" = cohort_params$coverage_phone, 
                  "home" = cohort_params$coverage_phone * cohort_params$coverage_home_reduction),
  sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home, 
                              "phone" = cohort_params$sensitivity_symptoms_home * cohort_params$sensitivity_symptoms_phone_reduction),
  success_sputum = list("home" = cohort_params$success_sputum_home, 
                        "phone" = cohort_params$success_sputum_home * cohort_params$success_sputum_phone_reduction))


cohort_screened <- apply_intervention(cohort, design_illustrative, intervention_parameters, cohort_params)

pdf("Aftermath model Fig2.pdf", width = 6, height = 9)
plot_screening(apply_intervention(cohort_screened,
                                  cohort_params = cohort_params,
                                  intervention_parameters = intervention_parameters,
                                  design = design_illustrative),
               screening_design = design_illustrative, colorfill = F)
dev.off()
cohort %>% filter(diagnosis_routine <= 30*30) %>% count(n())
