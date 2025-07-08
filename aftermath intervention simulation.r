source("aftermath cohort simulation functions.R")

#### Screening #### 

screening_design <- data.frame(
  "timing_months" = c(3, 3, 6, 12),
  "target_coverage" = c(0.5, 1, 1, 1),
  "screening_method" = c("micro", "symptoms", "symptoms", "symptoms"),
  "screening_location" = c("home", "phone", "phone", "phone")
) %>% 
  arrange(timing_months)

intervention_parameters <- list(
  coverage = list("home" = 0.81, "phone" = 0.95),
  sensitivity_symptoms = list("home" = 0.8, "phone" = 0.5), #* Aye to explore data, but we'll have to guess in the end
  success_sputum = list("home" = 0.9, "phone" = 0.7))

cohort <- create_cohort(cohort_params)

temp <- get_screening_coverage(cohort, screening_design, cohort_params)
covered_screening <- temp$covered
cohort <- temp$cohort

# apply df of screening rounds
for (r in 1:nrow(screening_design)) {
  cohort <- apply_screening_round(cohort, 
                            covered_screening[,r],
                            timing_months = screening_design$timing_months[r], 
                            screening_method = screening_design$screening_method[r], 
                            screening_location = screening_design$screening_location[r],
                            intervention_parameters)
}


# rseults
cohort %>% filter(TB==1) %>% mutate(detected = !is.na(detection_timing)) %>% 
  summarise(mean(detected))

plot_screening(cohort, )

time_with_tb(cohort)    


time_with_tb(cohort)  %>% filter(outcome == "symptom") %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
  summarise(screening/soc)
time_with_tb(cohort)  %>% filter(outcome == "sputum") %>% pivot_wider(names_from = "scenario", values_from = "value") %>%
  summarise(screening/soc)


costs(screening_design, cohort, covered_screening)


# cost effectiveness: cost per month of symptomatic or infectious TB prevented
sum(costs(screening_design, cohort, covered_screening))/
  (time_with_tb(cohort)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
     mutate(months_averted = (soc - screening)/30) %>% select(months_averted))
# per "symptomatic" then "sputum+" month averted

# Combine with an estimate of DALYs per case detection from my recent analysis (tho that's CXR and not recurrent TB) to estimate ICER?

# Mortality? (Do we have the data to estimate this? What would we need to assume?)

plot_screening(cohort, screening_design)
