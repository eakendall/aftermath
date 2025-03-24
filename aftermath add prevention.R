#### Prevention ####

# Model a single preventive intervention (whichever you think will be most cost effective), with a:
# target % population, 
# estimated coverage of intended populatoin
# estimated effectiveness (relative risk of recurrence; timing assumed not to change)

coverage_prev <- 0.5 # of those targeted
efficacy_prev <- 0.7 # relative risk of TB

apply_prevention <- function(
    cohort, 
    target_coverage = 1 # (targeting highest risk first)
)
{
  cohort <- cohort %>% mutate(prevented = case_when(
    TB == 1 & risk_score > target_coverage ~ 
      rbinom(n = n(), size = 1, prob = coverage_prev * efficacy_prev),
    TB == 1 ~ 0, 
    TRUE ~ NA))
  return(cohort)  
}

cohort <- apply_prevention(cohort, 
                           target_coverage = 0.5)


# Plot who was detected or prevented

ggplot(cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
         mutate(newID = row_number(),
                TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                                   TRUE ~ "Clinical or extrapulmonary"),
                detected = !is.na(detection_timing))) + 
  geom_vline(data = screening_design, aes(xintercept = timing_months*30, linetype = screening_method)) + 
  geom_segment(aes(x = sputum_onset, xend =  diagnosis_routine, y = newID - 0.2, col='Sputum+', linetype = as.factor(prevented), alpha = as.factor(!prevented)), size=0.8) + 
  geom_segment(aes(x = symptom_onset, xend = diagnosis_routine, y = newID + 0.2, col='Symptomatic', linetype = as.factor(prevented), alpha = as.factor(!prevented)), size=0.8) + 
  scale_color_manual("Time with TB", values = colors) + 
  scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
  geom_point(aes(x = diagnosis_routine, y = newID, fill = TBtype), pch=21) + 
  geom_point(data = cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% mutate(newID = row_number()) %>% filter(prevented == 1),
             aes(x = diagnosis_routine, y = newID), col = 'black', pch=4, size=4) + 
  ylab("Cases arranged by timing of routine diagnosis") + 
  xlab("Days since prior treatment completion") +
  xlim(0,30*24) + 
  geom_point(aes(x = detection_timing, y = newID), pch=8, col = "green", size=3) + 
  theme_minimal()




time_with_tb <- function(cohort)
{
  # add time 
  outcomes <- cohort %>% mutate(
    symptom_days_soc = diagnosis_routine - symptom_onset,
    sputum_days_soc = diagnosis_routine - sputum_onset,
    symptom_days_screening = case_when(detection_timing > symptom_onset & 
                                         detection_timing < diagnosis_routine ~ 
                                         diagnosis_routine - detection_timing,
                                       TRUE ~ sputum_days_soc),
    sputum_days_screening = case_when(detection_timing > sputum_onset &
                                        detection_timing < diagnosis_routine ~
                                        diagnosis_routine - detection_timing,
                                      TRUE ~ sputum_days_soc),
    symptom_days_prevention = case_when(prevented == 1 ~ 0,
                                        TRUE ~ symptom_days_soc),
    sputum_days_prevention = case_when(prevented == 1 ~ 0,
                                       TRUE ~ sputum_days_soc),
    symptom_days_both = case_when(prevented == 1 ~ 0,
                                  TRUE ~ symptom_days_screening),
    sputum_days_both = case_when(prevented == 1 ~ 0,
                                 TRUE ~ sputum_days_screening))
  results <- outcomes %>% summarise(
    symptom_days_soc = sum(symptom_days_soc, na.rm = TRUE),
    symptom_days_screening = sum(symptom_days_screening, na.rm = TRUE),
    symptom_days_prevention = sum(symptom_days_prevention, na.rm = TRUE),
    symptom_days_both = sum(symptom_days_both, na.rm = TRUE),
    sputum_days_soc = sum(sputum_days_soc, na.rm = TRUE),
    sputum_days_screening = sum(sputum_days_screening, na.rm = TRUE),
    sputum_days_prevention = sum(sputum_days_prevention, na.rm = TRUE),
    sputum_days_both = sum(sputum_days_both, na.rm = TRUE)
  ) %>% 
    pivot_longer(names_sep = "_", cols = everything(), names_to = c("outcome", "unit", "scenario"))
  return(results)
}
