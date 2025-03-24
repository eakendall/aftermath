#### Cohort ####

source("aftermath cohort simulation functions.R")

cohort_params <- list(
  N = 10000,
  symptom_duration_mean_reported = 17, 
  symptom_duration_sd_reported = 14*1.5, # Increased 1.5x; long tails perhaps more common in trial nonparticipants/losses, or perhaps underreported; will explore sensitivity
  symptom_duration_timescale = 1/4,
  symptom_underestimation_factor = 0.5, # 1 means no underestimation, 0.33 means 3x longer than reported
  recurrence_time_mean = 216, # diagnosis minus reported symptom duration
  recurrence_time_sd = 151,
  incidence_18mo = 0.082,
  proportion_micro_pos = 0.04275/0.08055, 
  proportion_ever_subclinical = 3/4, # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = 1,
  duration_subclinical_cv = 1,
  auc = 0.7,
  home_visit_passive_detection_impact = 0.8 # time by which symptomatic time is shortened after initial home visit, assumed.
    # all aftermath participants are assumed to have experienced this due to initial counseling (in addition to underestimation factor above).
)

cohort <- create_cohort(cohort_params)

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
  xlim(0,30*24) + 
  ggthemes::theme_fivethirtyeight()


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

table(cohort$TB, cohort$risk_quantile > 0.5); prop.table(table(cohort$TB, cohort$risk_quantile > 0.5), margin=1)
table(cohort$TB, cohort$risk_quantile > 0.9); prop.table(table(cohort$TB, cohort$risk_quantile > 0.9), margin=1)



# Check the number who are subclinical sputum+ from the start:

cohort %>% filter(TB==1) %>% summarise(mean(sputum_onset == 0, na.rm = TRUE), n()) # around 5% of all recurrences
cohort %>% filter(TB==1) %>% summarise(sum(sputum_onset == 0, na.rm = TRUE)/N, n()) # around 0.3% of the overall cohort -- both c/w crtriumph (or a little low).

# Check the number who are subclinical at 6mo postrx:
cohort %>% filter(TB==1) %>% count(u = diagnosis_routine >= 30*6, s = sputum_onset <=30*6) %>%
  filter(s, u) %>% reframe(n/N) 
# target is 11/861 = 0.013 cx+, so slightly fewer (~0.09 = 0.13*0.7?) detectable by Xpert. 

