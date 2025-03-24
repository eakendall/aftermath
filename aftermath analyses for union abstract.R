library(binom)
library(clipr)

source("aftermath cohort simulation functions.R")

#### incorporate uncertainty ####
# Sample cohort_params from posterior distributions
cohort_param_ranges <- list(
  N = 20000, 
  symptom_duration_mean_reported = 17,
  symptom_duration_sd_reported = 14*c(1, 2), 
  symptom_duration_timescale = 0, # increase this to 1/4 in a sensitivity analysis
  symptom_underestimation_factor = c(0.33, 0.67), 
  recurrence_time_mean = c(233 - 17*3, 233 + 17*3), # diagnosis minus reported symptom duration
  recurrence_time_sd = 151,
  incidence_18mo = c(0.065, 0.100),
  proportion_micro_pos = c(0.42, 0.53), #binom.exact(n = 90, x = 0.04275/0.08055*90)
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
  success_sputum_phone_reduction = c(0.7, 0.85),
  
  home_visit_passive_detection_impact = 1 # except in sensitivity analysis
  
)
  
N_samples <- 100

cohort_params <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 
  
# For each sample, simulate a cohort and estimate the impact of the aftermath home intervention:

screening_design_aftermath <- 
  screening_design <- data.frame(
    "timing_months" = c(6, 12),
    "target_coverage" = c(1, 1),
    "screening_method" = c("symptoms", "symptoms"),
    "screening_location" = c("home", "phone")
  ) %>% 
  arrange(timing_months)

# identify same-cost options with greater predicted impact:

## earlier screening
screening_design_earlier <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 9),
    "target_coverage" = c(1, 11),
    "screening_method" = c("symptoms", "symptoms"),
    "screening_location" = c("home", "phone")
  ) %>% 
  arrange(timing_months)

## more frequent but targeted
screening_design_targeted <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 6, 9, 12),
    "target_coverage" = rep(0.5,4),
    "screening_method" = rep("symptoms", 4),
    "screening_location" = c("home", rep("phone",3))
  ) %>% 
  arrange(timing_months)

## sputum at initial visit for highest risk
# (should modify to include confirmatory testing or false positives?)
screening_design_sputum <- 
  screening_design <- data.frame(
    "timing_months" = c(3, 3, 9),
    "target_coverage" = c(0.5, 1, 1),
    "screening_method" = c("sputum", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# will also explore aftermath but with impact on subsequent passive diagnosis

results_aftermath <- 
  results_earlier <- 
  results_targeted <- 
  results_sputum <- 
  results_home_improves_passive <- 
            data.frame(detections = numeric(N_samples),
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
  (time <- time_with_tb(cohort_screened) %>% filter(scenario == "soc") %>% 
    mutate(months = value/30) %>%
    tibble::column_to_rownames('outcome') %>% select(months))
  (impact <- (time_with_tb(cohort_screened)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
      mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_earlier <- sum(!is.na(cohort_screened_earlier$detection_timing)))
  (time_earlier <- time_with_tb(cohort_screened_earlier) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_earlier <- (time_with_tb(cohort_screened_earlier)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_targeted <- sum(!is.na(cohort_screened_targeted$detection_timing)))
  (time_targeted <- time_with_tb(cohort_screened_targeted) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_targeted <- (time_with_tb(cohort_screened_targeted)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                        mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_sputum <- sum(!is.na(cohort_screened_sputum$detection_timing)))
  (time_sputum <- time_with_tb(cohort_screened_sputum) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_sputum <- (time_with_tb(cohort_screened_sputum)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                        mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  (detected_alt <- sum(!is.na(cohort_screened_alt$detection_timing)))
  (time_alt <- time_with_tb(cohort_screened_alt) %>% filter(scenario == "soc") %>% 
      mutate(months = value/30) %>%
      tibble::column_to_rownames('outcome') %>% select(months))
  (impact_alt <- (time_with_tb(cohort_screened_alt)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                       mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted))
  
  # costs
  (cost <- sum(costs(screening_design_aftermath, cohort_screened, covered_screening)))
  (cost_earlier <- sum(costs(screening_design_earlier, cohort_screened_earlier, covered_screening_earlier)))
  (cost_targeted <- sum(costs(screening_design_targeted, cohort_screened_targeted, covered_screening_targeted)))
  (cost_sputum <- sum(costs(screening_design_sputum, cohort_screened_sputum, covered_screening_sputum)))
  (cost_alt <- sum(costs(screening_design_aftermath, cohort_screened_alt, covered_screening_alt)))
  
  # cost-benefit ratios
  cost/impact
  

  results_aftermath[n,] <- c(detected, time$months[1], time$months[2], impact$months_averted[1], impact$months_averted[2], cost,
                    cost/impact$months_averted[1], cost/impact$months_averted[2])
  results_earlier[n,] <- c(detected_earlier, time_earlier$months[1], time_earlier$months[2], 
                           impact_earlier$months_averted[1], impact_earlier$months_averted[2], cost_earlier,
                             cost_earlier/impact_earlier$months_averted[1], cost_earlier/impact_earlier$months_averted[2])
  results_targeted[n,] <- c(detected_targeted, time_targeted$months[1], time_targeted$months[2], 
                           impact_targeted$months_averted[1], impact_targeted$months_averted[2], cost_targeted,
                           cost_targeted/impact_targeted$months_averted[1], cost_targeted/impact_targeted$months_averted[2])
  results_sputum[n,] <- c(detected_sputum, time_sputum$months[1], time_sputum$months[2], 
                           impact_sputum$months_averted[1], impact_sputum$months_averted[2], cost_sputum,
                           cost_sputum/impact_sputum$months_averted[1], cost_sputum/impact_sputum$months_averted[2])
  results_home_improves_passive[n,] <- c(detected_alt, time_alt$months[1], time_alt$months[2], 
                          impact_alt$months_averted[1], impact_alt$months_averted[2], cost_alt,
                          cost_alt/impact_alt$months_averted[1], cost_alt/impact_alt$months_averted[2])
}


# collate results for a table:

results_aftermath$intervention <- "aftermath"
results_earlier$intervention <- "earlier"
results_targeted$intervention <- "targeted"
results_sputum$intervention <- "sputum"
results_home_improves_passive$intervention <- "home_visit_improves_passive"

allresults <- rbind(results_aftermath, results_earlier, results_targeted, results_sputum, results_home_improves_passive)
allresults %>% group_by(intervention) %>% mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc, 
                                                 inf_reduction = infectious_months_averted/infectious_months_soc, 
                                                 cost_per_sx_month_averted = cost/symptomatic_months_averted,
                                                 cost_per_inf_month_averted = cost/infectious_months_averted) %>%
  select(detections, sx_reduction, inf_reduction, cost_per_sx_month_averted, cost_per_inf_month_averted) %>%
  summarise_all(list(mean, 
                text = function(x) paste0(round(median(x),2), " (", round(quantile(x, 0.25),2), ", ", round(quantile(x, 0.75),2), ")"))) %>% 
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

summary(results_alt$cost_per_symptomatic_month_averted)
summary(results_alt$cost_per_infectious_months_averted)


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





