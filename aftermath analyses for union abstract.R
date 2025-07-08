library(binom)
library(clipr)

source("aftermath cohort simulation functions.R")

#### incorporate uncertainty ####
# Sample cohort_params from posterior distributions
cohort_param_ranges <- list(
  N = 100000, 
  symptom_duration_mean_reported = 17,
  symptom_duration_sd_reported = 14*c(1, 2), 
  symptom_duration_timescale = 0, # increase this to 1/4 in a sensitivity analysis
  symptom_underestimation_factor = c(0.33, 0.67), 
  recurrence_time_mean = c((233-17*2) - 17*2, (233-17*2) + 17*2), # diagnosis minus symptom duration
  # recurrence_time_sd = 151,
  recurrence_time_shapefactor = c(0.5, 1.5), # variation in (mean/sd)^2, compared to aftermath ratio where it was (199/151)^2
  incidence_18mo = c(0.065, 0.100),
  proportion_micro_pos = c(0.42, 0.64), #binom.exact(n = 90, x = 0.04275/0.08055*90)
  proportion_ever_subclinical = c(0.6, 0.9), # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2), 
  duration_subclinical_cv = c(0.5, 1.5),
  auc = c(0.6, 0.8),

  coverage_phone = c(0.9, 1),
  coverage_home_reduction = c(0.75, 0.95), #0.81/0.95 = 0.85
  sensitivity_symptoms_home = c(0.7, 0.9),
  sensitivity_symptoms_phone_reduction = c(0.5, 0.75),
  success_sputum_home = c(0.85, 0.95),
  success_sputum_phone_reduction = c(0.6, 1.0),
  
  home_visit_passive_detection_impact = c(1), # except in sensitivity analysis
  
  # also sample costs??
  initial_contact_cost_home = c(5, 15), 
  initial_contact_cost_phone = c(1, 2), 
  sputum_test_cost = c(14,18)
  
)
  
N_samples <- 500

cohort_params <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 
  
# For each sample, simulate a cohort and describe key features: 
 # month with highest incidence, month with highest undiagnosed prevalence (symptomatic and sputum+, and overall), and 
# remaining expected duration of prevalent disease at each of 3, 6, 9, 12, 15, and 18mo
cohort_features <- 
  data.frame(month_with_highest_incidence = numeric(N_samples),
             month_with_highest_prevalence_sx = numeric(N_samples),
             month_with_highest_prevalence_inf = numeric(N_samples),
             month_with_highest_prevalence_total = numeric(N_samples),
             remaining_duration_3mo = numeric(N_samples),
             remaining_duration_6mo = numeric(N_samples),
             remaining_duration_9mo = numeric(N_samples),
             remaining_duration_12mo = numeric(N_samples),
             remaining_duration_15mo = numeric(N_samples),
             remaining_duration_18mo = numeric(N_samples))

for (n in 1:N_samples)
{
  # setup
  cohort <- create_cohort(cohort_params[n,])
  
  cohort_features[n, "month_with_highest_incidence"] <- 
    (cohort %>% filter(TB==1) %>% 
       mutate(month = floor(symptom_onset/30)) %>% group_by(month) %>% summarise(n = n()) %>% 
       arrange(desc(n)) %>% select(month)) %>%  .[1,1]
  
  cohort_features[n, "month_with_highest_prevalence_sx"] <-
    ((1:20))[which.max(sapply(30 * (1:20), function(x) sum(cohort$symptom_onset < x & cohort$diagnosis_routine > x, na.rm=T)))]

  cohort_features[n, "month_with_highest_prevalence_inf"] <-
    ((1:20))[which.max(sapply(30 * (1:20), function(x) sum(cohort$sputum_onset < x & cohort$diagnosis_routine > x, na.rm=T)))]

  cohort_features[n, "month_with_highest_prevalence_total"] <-
    ((1:20))[which.max(sapply(30 * (1:20), function(x) 
      sum((cohort$sputum_onset < x | cohort$symptom_onset < x) & cohort$diagnosis_routine > x, na.rm=T)))]
  
  cohort_features[n, "remaining_duration_3mo"] <-
    cohort %>% filter(diagnosis_routine > 3*30 & (sputum_onset < 3*30 | symptom_onset < 3*30)) %>% 
      summarise(mean(diagnosis_routine - 3*30, na.rm = TRUE))
  cohort_features[n, "remaining_duration_6mo"] <-
    cohort %>% filter(diagnosis_routine > 6*30 & (sputum_onset < 6*30 | symptom_onset < 6*30)) %>% 
      summarise(mean(diagnosis_routine - 6*30, na.rm = TRUE))
  cohort_features[n, "remaining_duration_9mo"] <-
    cohort %>% filter(diagnosis_routine > 9*30 & (sputum_onset < 9*30 | symptom_onset < 9*30)) %>% 
      summarise(mean(diagnosis_routine - 9*30, na.rm = TRUE))
  cohort_features[n, "remaining_duration_12mo"] <-
    cohort %>% filter(diagnosis_routine > 12*30 & (sputum_onset < 12*30 | symptom_onset < 12*30)) %>% 
      summarise(mean(diagnosis_routine - 12*30, na.rm = TRUE))
  cohort_features[n, "remaining_duration_15mo"] <-
    cohort %>% filter(diagnosis_routine > 15*30 & (sputum_onset < 15*30 | symptom_onset < 15*30)) %>% 
      summarise(mean(diagnosis_routine - 15*30, na.rm = TRUE))
  cohort_features[n, "remaining_duration_18mo"] <-
    cohort %>% filter(diagnosis_routine > 18*30 & (sputum_onset < 18*30 | symptom_onset < 18*30)) %>% 
      summarise(mean(diagnosis_routine - 18*30, na.rm = TRUE))
}

cohort_features[,] %>% summarise_all(list(median))
cohort_features[,] %>% summarise_all(function(x) quantile(x, 0.25))
cohort_features[,] %>% summarise_all(function(x) quantile(x, 0.75)) # for the abstract
    
# For each sample, simulate a cohort and estimate the impact of the aftermath home intervention:

screening_design_aftermath <- 
  screening_design <- data.frame(
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# identify same-cost options with greater predicted impact:

## earlier screening
screening_design_earlier <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 6, 9),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## more frequent but targeted. Could include the highest ~0.62 risk, if only the first is home and home costs 3x as much?
#i.e. 3+1+1 vs (5/8)*(3+5*1)
screening_design_targeted <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 6, 9, 12, 15, 18),
    "target_coverage" = rep(0.625,6),
    "screening_method" = rep("symptoms", 6),
    "screening_location" = c("home", rep("phone",5))
  ) %>% 
  arrange(timing_months)

## sputum at initial visit for highest risk
# (should modify to include confirmatory testing or false positives?)
# if cost of a sputum test + a home visit is ~6x a phone call, 
#  then can do sputum+home visit for highest x risk where x*6 + (1-x)*1 = 3 --> x = 40%
screening_design_sputum <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 6, 12, 18),
    "target_coverage" = c(0.4, 1, 1, 1),
    "screening_method" = c("both", "symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# aftermath for a fraction, plus sputum at first visit (which increases total cost by (2*3+1+1)/(3+1+1)=1.6x, so target top 62%)
screening_design_moretargeted <- 
  screening_design <- data.frame(
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1,1,1) * 0.625,
    "screening_method" = c("both", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# will also explore aftermath but with impact on subsequent passive diagnosis.
# Need to use "SOC" for aftermath, since that's the SOC with no shortening. 

results_aftermath <- 
  results_earlier <- 
  results_targeted <- 
  results_sputum <- 
  results_moretargeted <- 
  results_home_improves_passive <- 
            data.frame(detections = numeric(N_samples),
                      mean_symptom_days = numeric(N_samples),
                      symptomatic_months_soc = numeric(N_samples),
                      infectious_months_soc = numeric(N_samples),
                      symptomatic_months_averted = numeric(N_samples), 
                      infectious_months_averted = numeric(N_samples), 
                      cost = numeric(N_samples),
                      cost_per_symptomatic_month_averted = numeric(N_samples),
                      cost_per_infectious_months_averted = numeric(N_samples))

for (n in 1:N_samples)
{

  # setup
  cohort <- create_cohort(cohort_params[n,])

  intervention_parameters <- list(
    coverage = list("phone" = cohort_params$coverage_phone[n], 
                    "home" = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]),
    sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home[n], 
                                "phone" = cohort_params$sensitivity_symptoms_home[n] * cohort_params$sensitivity_symptoms_phone_reduction[n]),
    success_sputum = list("home" = cohort_params$success_sputum_home[n], 
                          "phone" = cohort_params$success_sputum_home[n] * cohort_params$success_sputum_phone_reduction[n]))
  
  # intervention
  
  screened <- get_screening_coverage(cohort, screening_design_aftermath, cohort_params[n,])
  covered_screening <- screened$covered
  cohort_screened <- screened$cohort
  
  for (r in 1:nrow(screening_design_aftermath)) {
    cohort_screened <- apply_screening_round(cohort_screened, 
                              covered_screening[,r],
                              timing_months = screening_design_aftermath$timing_months[r], 
                              screening_method = screening_design_aftermath$screening_method[r], 
                              screening_location = screening_design_aftermath$screening_location[r],
                              intervention_parameters)
  }
  
  # alternative interventions
  screened_earlier <- get_screening_coverage(cohort, screening_design_earlier, cohort_params[n,])
  covered_screening_earlier <- screened_earlier$covered
  cohort_screened_earlier <- screened_earlier$cohort
  
  for (r in 1:nrow(screening_design_earlier)) {
    cohort_screened_earlier <- apply_screening_round(cohort_screened_earlier, 
                                             covered_screening_earlier[,r],
                                             timing_months = screening_design_earlier$timing_months[r], 
                                             screening_method = screening_design_earlier$screening_method[r], 
                                             screening_location = screening_design_earlier$screening_location[r],
                                             intervention_parameters)
  }
  
  screened_targeted <- get_screening_coverage(cohort, screening_design_targeted, cohort_params[n,])
  covered_screening_targeted <- screened_targeted$covered
  cohort_screened_targeted <- screened_targeted$cohort
  
  for (r in 1:nrow(screening_design_targeted)) {
    cohort_screened_targeted <- apply_screening_round(cohort_screened_targeted, 
                                                     covered_screening_targeted[,r],
                                                     timing_months = screening_design_targeted$timing_months[r], 
                                                     screening_method = screening_design_targeted$screening_method[r], 
                                                     screening_location = screening_design_targeted$screening_location[r],
                                                     intervention_parameters)
  }
  
  screened_sputum <- get_screening_coverage(cohort, screening_design_sputum, cohort_params[n,])
  covered_screening_sputum <- screened_sputum$covered
  cohort_screened_sputum <- screened_sputum$cohort
  
  for (r in 1:nrow(screening_design_sputum)) {
    cohort_screened_sputum <- apply_screening_round(cohort_screened_sputum, 
                                                     covered_screening_sputum[,r],
                                                     timing_months = screening_design_sputum$timing_months[r], 
                                                     screening_method = screening_design_sputum$screening_method[r], 
                                                     screening_location = screening_design_sputum$screening_location[r],
                                                     intervention_parameters)
  }
  
  screened_moretargeted <- get_screening_coverage(cohort, screening_design_moretargeted, cohort_params[n,])
  covered_screening_moretargeted <- screened_moretargeted$covered
  cohort_screened_moretargeted <- screened_moretargeted$cohort
  
  for (r in 1:nrow(screening_design_moretargeted)) {
    cohort_screened_moretargeted <- apply_screening_round(cohort_screened_moretargeted, 
                                                      covered_screening_moretargeted[,r],
                                                      timing_months = screening_design_moretargeted$timing_months[r], 
                                                      screening_method = screening_design_moretargeted$screening_method[r], 
                                                      screening_location = screening_design_moretargeted$screening_location[r],
                                                      intervention_parameters)
  }
  
  
  # home visit improves passive
  cohort_params_alt <- cohort_params[n,]; cohort_params_alt["home_visit_passive_detection_impact"] <- 0.8
  
  screened_alt <- get_screening_coverage(cohort, screening_design_aftermath, cohort_params_alt)
  covered_screening_alt <- screened_alt$covered
  cohort_screened_alt <- screened_alt$cohort
  
  for (r in 1:nrow(screening_design_aftermath)) {
    cohort_screened_alt <- apply_screening_round(cohort_screened_alt, 
                                             covered_screening_alt[,r],
                                             timing_months = screening_design_aftermath$timing_months[r], 
                                             screening_method = screening_design_aftermath$screening_method[r], 
                                             screening_location = screening_design_aftermath$screening_location[r],
                                             intervention_parameters)
  }
  
  
  
  # impact
  (detected <- sum(!is.na(cohort_screened$detection_timing)))
  (symptomdays <- ((time_with_tb(cohort_screened) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened$TB))[1,])
  (time <- time_with_tb(cohort_screened) %>% filter(scenario == "soc") %>% 
    mutate(months = value/30) %>%
    tibble::column_to_rownames('outcome') %>% select(months))
  (impact <- (time_with_tb(cohort_screened)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
      mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_earlier <- sum(!is.na(cohort_screened_earlier$detection_timing)))
  (symptomdays_earlier <- ((time_with_tb(cohort_screened_earlier) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened_earlier$TB))[1,])
  (time_earlier <- time_with_tb(cohort_screened_earlier) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_earlier <- (time_with_tb(cohort_screened_earlier)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_targeted <- sum(!is.na(cohort_screened_targeted$detection_timing)))
  (symptomdays_targeted <- ((time_with_tb(cohort_screened_targeted) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened_targeted$TB))[1,])
  (time_targeted <- time_with_tb(cohort_screened_targeted) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_targeted <- (time_with_tb(cohort_screened_targeted)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                        mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_sputum <- sum(!is.na(cohort_screened_sputum$detection_timing)))
  (symptomdays_sputum <- ((time_with_tb(cohort_screened_sputum) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened_sputum$TB))[1,])
  (time_sputum <- time_with_tb(cohort_screened_sputum) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_sputum <- (time_with_tb(cohort_screened_sputum)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                        mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_moretargeted <- sum(!is.na(cohort_screened_moretargeted$detection_timing)))
  (symptomdays_moretargeted <- ((time_with_tb(cohort_screened_moretargeted) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened_moretargeted$TB))[1,])
  (time_moretargeted <- time_with_tb(cohort_screened_moretargeted) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_moretargeted <- (time_with_tb(cohort_screened_moretargeted)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                         mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_alt <- sum(!is.na(cohort_screened_alt$detection_timing)))
  (symptomdays_alt <- ((time_with_tb(cohort_screened_alt) %>% filter(scenario == "screening") %>% select(value))/sum(cohort_screened_alt$TB))[1,])
  (time_alt <- time_with_tb(cohort_screened) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_alt <- (time_with_tb(cohort_screened_alt)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                       mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  # costs
  (cost <- sum(costs(screening_design_aftermath, cohort_screened, covered_screening, 
                     cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  (cost_earlier <- sum(costs(screening_design_earlier, cohort_screened_earlier, covered_screening_earlier,
                             cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  (cost_targeted <- sum(costs(screening_design_targeted, cohort_screened_targeted, covered_screening_targeted,
                              cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  (cost_sputum <- sum(costs(screening_design_sputum, cohort_screened_sputum, covered_screening_sputum,
                            cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  (cost_moretargeted <- sum(costs(screening_design_moretargeted, cohort_screened_moretargeted, covered_screening_moretargeted,
                              cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  (cost_alt <- sum(costs(screening_design_aftermath, cohort_screened_alt, covered_screening_alt,
                         cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
  
  # cost-benefit ratios
  cost/impact
  

  results_aftermath[n,] <- c(detected, symptomdays, time$months[1], time$months[2], impact$months_averted[1], impact$months_averted[2], cost,
                    cost/impact$months_averted[1], cost/impact$months_averted[2])
  results_earlier[n,] <- c(detected_earlier, symptomdays_earlier, time_earlier$months[1], time_earlier$months[2], 
                           impact_earlier$months_averted[1], impact_earlier$months_averted[2], cost_earlier,
                             cost_earlier/impact_earlier$months_averted[1], cost_earlier/impact_earlier$months_averted[2])
  results_targeted[n,] <- c(detected_targeted, symptomdays_targeted, time_targeted$months[1], time_targeted$months[2], 
                           impact_targeted$months_averted[1], impact_targeted$months_averted[2], cost_targeted,
                           cost_targeted/impact_targeted$months_averted[1], cost_targeted/impact_targeted$months_averted[2])
  results_sputum[n,] <- c(detected_sputum, symptomdays_sputum, time_sputum$months[1], time_sputum$months[2], 
                           impact_sputum$months_averted[1], impact_sputum$months_averted[2], cost_sputum,
                           cost_sputum/impact_sputum$months_averted[1], cost_sputum/impact_sputum$months_averted[2])
  results_moretargeted[n,] <- c(detected_moretargeted, symptomdays_moretargeted, time_moretargeted$months[1], time_moretargeted$months[2], 
                            impact_moretargeted$months_averted[1], impact_moretargeted$months_averted[2], cost_moretargeted,
                            cost_moretargeted/impact_moretargeted$months_averted[1], cost_moretargeted/impact_moretargeted$months_averted[2])
  results_home_improves_passive[n,] <- c(detected_alt, symptomdays_alt, time_alt$months[1], time_alt$months[2], 
                          impact_alt$months_averted[1], impact_alt$months_averted[2], cost_alt,
                          cost_alt/impact_alt$months_averted[1], cost_alt/impact_alt$months_averted[2])
}


# collate results for a table:

results_aftermath$intervention <- "aftermath"
results_earlier$intervention <- "earlier"
results_targeted$intervention <- "targeted"
results_sputum$intervention <- "sputum"
results_moretargeted$intervention <- "moretargeted"
results_home_improves_passive$intervention <- "home_visit_improves_passive"

allresults <- rbind(results_aftermath, results_earlier, results_targeted, results_sputum, results_moretargeted, results_home_improves_passive)
allresults %>% group_by(intervention) %>% mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
                                                 inf_reduction = infectious_months_averted/infectious_months_soc*100, 
                                                 cost_per_sx_month_averted = cost/symptomatic_months_averted,
                                                 cost_per_inf_month_averted = cost/infectious_months_averted) %>%
  select(detections, mean_symptom_days, sx_reduction, inf_reduction, cost_per_sx_month_averted, cost_per_inf_month_averted) %>%
  summarise_all(list(mean, 
                text = function(x) paste0(round(median(x),0), " (", round(quantile(x, 0.25),0), ", ", round(quantile(x, 0.75),0), ")"))) %>% 
  write_clip()







# summarize results for aftermath intervention
summary(results_aftermath$cost_per_symptomatic_month_averted)
summary(results_aftermath$cost_per_infectious_months_averted)

summary(results_earlier$cost_per_symptomatic_month_averted)
summary(results_earlier$cost_per_infectious_months_averted)

summary(results_targeted$cost_per_symptomatic_month_averted)
summary(results_targeted$cost_per_infectious_months_averted)

summary(results_sputum$cost_per_symptomatic_month_averted)
summary(results_sputum$cost_per_infectious_months_averted)

summary(results_moretargeted$cost_per_symptomatic_month_averted)
summary(results_moretargeted$cost_per_infectious_months_averted)

summary(results_home_improves_passive$cost_per_symptomatic_month_averted)
summary(results_home_improves_passive$cost_per_infectious_months_averted)


# sensitivty analysis
# for each parameter in cohort_param_ranges, compare results$cost_per_symptomatic_case_averted 
# between the top 10% and bottom 10% of values for that paramter
parameters_varied <- names(cohort_param_ranges)[lapply(cohort_param_ranges, length)>1]
oneways <- array(NA, dim=c(2, length(parameters_varied)))
dimnames(oneways) = list(c("top", "bottom"), parameters_varied)
for (p in parameters_varied) {
  oneways["top", p] <- mean(results_aftermath[cohort_params[p] >= quantile(unlist(cohort_params[p]), 0.8),]$cost_per_symptomatic_month_averted)
  oneways["bottom", p] <- mean(results_aftermath[cohort_params[p] <= quantile(unlist(cohort_params[p]), 0.2),]$cost_per_symptomatic_month_averted)
}
# plot tornado diagram, centered at mean of results$cost_per_symptomatic_month_averted
oneways %>%
  as.data.frame() %>%
  rownames_to_column("top_or_bottom") %>%
  pivot_longer(cols = -top_or_bottom, names_to = "parameter", values_to = "cost_per_symptomatic_month_averted") %>%
  ggplot(aes(x = parameter, 
             col = top_or_bottom)) +
  geom_errorbar(aes(ymin = mean(results$cost_per_symptomatic_month_averted), ymax=cost_per_symptomatic_month_averted)) +
  coord_flip() +
  geom_hline(yintercept = mean(results$cost_per_symptomatic_month_averted), linetype = "dashed") + 
  # specify text labels in color legend
  labs(color = "Quantile of parameter values") + 
  scale_color_manual(labels = c("Lowest 20%", "Highest 20%"), values=c("red","blue"))
  
  
#############





