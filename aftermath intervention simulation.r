# As interventions, consider:
# Category A: intervention to prevent recurrence. E.g. extended treatment, nutrition, ?? 
# (For modeling, can just pick the one we think will be most cost effective.)
# 
# Category B: early detection. Symptom or micro based. Optimize timing based on observed diagnosis timing, 
# estimated symptom duration (explore variance and underreporting in sensitivity analysis), and estimated presymptomatic duration (again explore duration and variance). Other unknowns include infectivity, symptom, mortality correlations. 
# 
# Model risk prediction with error, for curve of PCT of recurrences vs PCT of total pop you could target. 
# For each intervention, consider different levels of targeting.


#### Screening #### 

# Model each roune with a:
  # timing (months post treatment), 
  # target % population.
  # screening method (symptom, micro, both), and
  # screening location (home, telephonic), with
    ## estimated sensitivity for symptoms and
    ## ability to get sputum for micro, and
    ## estimated coverage as % of target population.

screening_design <- data.frame(
  "timing_months" = c(2, 6, 10),
  "target_coverage" = c(1, 1, 1),
  "screening_method" = c("both", "symptoms", "symptoms"),
  "screening_location" = c("home", "phone", "phone")
)

intervention_parameters <- list(
  coverage = list("home" = 0.95, "phone" = 0.81),
  sensitivity_symptoms = list("home" = 0.8, "phone" = 0.5), #* Aye to explore data, but we'll have to guess in the end
  success_sputum = list("home" = 0.9, "phone" = 0.7))

# Identify the dates when any cases are detected by screening
apply_screening <- function(
    cohort, 
    timing_months = 6,
    target_coverage = 0.5, # (targeting highest risk first)
    screening_method = "symptoms", # (symptoms, micro, or both)
    screening_location = "home", # (home, phone)
    intervention_parameters
)
{
  if(!("detection_timing" %in% colnames(cohort))) cohort$detection_timing <- NA
  
  # symptom screening
    cohort <- cohort %>% mutate(
      detected_current_round = case_when(
        TB != 1 ~ 0,
        risk_score > target_coverage ~ 0, 
        timing_months*30 >= diagnosis_routine ~ 0,
        screening_method %in% c("symptoms", "both") &
          timing_months*30 > symptom_onset ~ 
            rbinom(n = n(), size = 1, prob = intervention_parameters$coverage[[screening_location]] * 
                                            intervention_parameters$sensitivity_symptoms[[screening_location]]),
        
        screening_method %in% c("micro", "both") & 
          timing_months*30 > sputum_onset ~ 
            rbinom(n = n(), size = 1, prob = intervention_parameters$coverage[[screening_location]] * 
                                            intervention_parameters$success_sputum[[screening_location]]),
        TRUE ~ 0),
      detection_timing = case_when(
            detected_current_round == 1 & 
              (is.na(detection_timing) | timing_months*30 < detection_timing) ~ timing_months*30,
            TRUE ~ detection_timing)
        )

return(cohort %>% select(-detected_current_round))  
}


# apply df of screening rounds
for (r in 1:nrow(screening_design)) {
  cohort <- apply_screening(cohort, 
                            timing_months = screening_design$timing_months[r], 
                            target_coverage = screening_design$target_coverage[r], 
                            screening_method = screening_design$screening_method[r], 
                            screening_location = screening_design$screening_location[r],
                            intervention_parameters)
}


# Plot who was detected
cohort %>% filter(TB==1) %>% mutate(detected = !is.na(detection_timing)) %>% 
  summarise(mean(detected))

ggplot(cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
         mutate(newID = row_number(),
                TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                                   TRUE ~ "Clinical or extrapulmonary"),
                detected = !is.na(detection_timing))) + 
  geom_vline(data = screening_design, aes(xintercept = timing_months*30, linetype = screening_method)) + 
  geom_segment(aes(x = sputum_onset, xend =  diagnosis_routine, y = newID - 0.2, col='Sputum+'), size=0.8) + 
  geom_segment(aes(x = symptom_onset, xend = diagnosis_routine, y = newID + 0.2, col='Symptomatic'), size=0.8) + 
  scale_color_manual("Time with TB", values = colors) + 
  scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
  geom_point(aes(x = diagnosis_routine, y = newID, fill = TBtype), pch=21) + 
  ylab("Cases arranged by timing of routine diagnosis") + 
  xlab("Days since prior treatment completion") +
  xlim(0,30*24) + 
  geom_point(aes(x = detection_timing, y = newID), pch=8, col = "green", size=3) + 
  theme_minimal()


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




#### Outcomes of interest #### 

# Cases detected (is this valuable in itself?)
 ## Cases detected "early" by intervention
# Incident recurrences 
  # Incident recurrences "prevented" by intervention
# Time with recurrent TB
 ## Time sputum+ (associated with transmission)
 ## Time symptomatic (associated with lung damage and with morality risk??)


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
  )
  return(results)
}

time_with_tb(cohort)    
# To get uncertainty, could create a large chort and bootstrap sample N, or could run medium-sized cohort repeatedly, 
# or could capture parameter uncertainty by probabilistic sampling. Probably want to do both, so sample parameters for each of many largish cohort runs.

time_with_tb(cohort)  %>% summarise(symptom_days_screening/symptom_days_soc)
time_with_tb(cohort)  %>% summarise(sputum_days_screening/sputum_days_soc)
time_with_tb(cohort)  %>% summarise(sputum_days_prevention/sputum_days_soc)


# Mortality? (Do we have the data to estimate this? What would we need to assume)

# Costs
  ## Cost per case detected early or prevented (combine these? Though not equally valuable?)
  ## Cost per month with TB prevented 
      ### Overall
      ### Symptomatic
      ### Sputum+

costs <- function(screening_design, cohort)
{
  initial_contact_cost = list(home = 3, phone = 1)
  sputum_test_cost = 15
  treatment_cost = 1000
  prevention_cost = 500
  
  symptom_prevalence_nonTB = 0.05 #** Will get from Aye's analysis of symptoms after retraining
  
  # count units of:
  # total people contacted
  # number of people with symptomatic TB, sputum+ asymptomatic TB, and non-TB symptoms
  # number of people with TB detected by screening
    
  # total people with initial contact:
  contacted <- screening_design$target_coverage * N * # targeted
    unlist(intervention_parameters$coverage[screening_design$screening_location])  # included
  contact_cost <- contacted * unlist(initial_contact_cost[screening_design$screening_location]) # cost
  
  # people with symptomatic TB, who are targeted and reached:
  sxtb <- colSums(cohort$TB * sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset < s) * 
                    sapply(screening_design$target_coverage, function(x) cohort$risk_score < x) *
                    unlist(intervention_parameters$coverage[screening_design$screening_location]) *
                 sapply(screening_design$timing_months*30, function(s) s <= cohort$diagnosis_routine), na.rm=T)

  # people with sputum+ asymptomatic TB, who are targeted and reached:
  asxtb <- colSums(cohort$TB * sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset > s & cohort$sputum_onset < s) * 
                     sapply(screening_design$target_coverage, function(x) cohort$risk_score < x) *
                     unlist(intervention_parameters$coverage[screening_design$screening_location]) *
            sapply(screening_design$timing_months*30, function(s) s <= cohort$diagnosis_routine), na.rm=T)
  
  # people with non-TB but with symptoms (never TB, or after diagnosis, or not yet sputum+ or symptom+:
  nontbsx <- colSums((cohort$TB ==0 |
            sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset > s & (cohort$sputum_onset > s | is.na(cohort$sputum_onset))) |
            sapply(screening_design$timing_months*30, function(s) s > cohort$diagnosis_routine)) & 
          sapply(screening_design$target_coverage, function(x) cohort$risk_score < x), na.rm=T) *
    unlist(intervention_parameters$coverage[screening_design$screening_location]) *
    symptom_prevalence_nonTB
  
  # sputum tests, 
  # if universal testing: 
  sputa <-contacted * unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
    screening_design$screening_method %in% c("micro", "both") + 
  # if symptom-based testing: 
    (sxtb + nontbsx) * unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
    (ifelse(screening_design$screening_method == "symptoms",1,0))
   

  
  # home visits (including a visit if symptom+ on phone call?)
  # sputum tests
  # treatments (incremental?)
  # patient costs?
  return(costs)
}



# cost effectiveneess
  costs <- costs %>% mutate(
    cost_detected = cost_per_case_detected * nrow(cohort) * mean(detected),
    cost_prevented = cost_per_month_prevented * nrow(cohort) * mean(symptom_days_prevention),
    cost_prevented_symptomatic = cost_per_month_prevented_symptomatic * nrow(cohort) * mean(symptom_days_prevention),
    cost_prevented_sputum = cost_per_month_prevented_sputum * nrow(cohort) * mean(sputum_days_prevention)
