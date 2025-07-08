library(binom)
library(clipr)

source("aftermath cohort sim functions.R")
N_cohort = 1000
N_samples <- 20

#### Incorporate uncertainty: Sample cohort_params from posterior distributions ####
cohort_param_ranges <- list(
  N = N_cohort, # number of participants in each cohort sample
  
  recurrence_time_mean = c(158, 225), # diagnosis minus symptom duration, of those diagnosed by day 540
  # recurrence_time_sd = 179,
  recurrence_time_shapefactor = c(0.5, 1.5), # variation in (mean/sd)^2, compared to aftermath ratio where it was (213/179)^2
  incidence_18mo = c(0.074, 0.110),
  proportion_micro_pos = c(0.43, 0.63), #binom.agresti.coull(n = 88, x = 0.04392/0.0822*88)
  auc = c(0.56, 0.77),
  
  symptom_duration_mean_reported = 17,
  symptom_duration_sd_reported = 14*c(1, 1.5), 
  symptom_duration_timescale = 0, # increase this to 1/4 in a sensitivity analysis
  symptom_underestimation_factor = c(0.25, 0.5), 
  
  proportion_ever_subclinical = c(0.6, 0.9), # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2), 
  duration_subclinical_cv = c(0.5, 1.5),
  subclinical_baseline_amongTB_max = 0.20, 
  subclinical_6m_amongcohort_min = 0.005, 
  subclinical_6m_amongcohort_max = 0.016, 
  
  coverage_phone = c(0.9, 1.0), 
  coverage_home_reduction = c(0.75, 0.95), #0.81/0.95 = 0.85
  sensitivity_symptoms_home = c(0.85, 0.95),
  sensitivity_symptoms_phone_reduction = c(0.6, 0.82),
  success_sputum_home = c(0.85, 0.95),
  success_sputum_phone_reduction = c(0.6, 1.0),
  
  home_visit_passive_detection_impact = c(1), # except in sensitivity analysis
  aftermath_counseling_passive_detection_impact = c(0.9,1), 
  intentional_counseling_passive_detection_impact = c(0.8), 
  intentional_counseling_passive_detection_duration = 180, # days before effect of counseling on care seeking "wears off"
  prevention_efficacy = 0.6,
  
  # also sample costs (for allocating intervention coverage in cost-equivalent way)
  initial_contact_cost_home = c(1,5), 
  initial_contact_cost_phone_factor = c(1/14, 1/8), 
  sputum_test_cost = c(2,4),
  prevention_cost = c(7.24, 9.14)
)
  


cohort_params <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 
  
##### For each sample, simulate a cohort and describe key features: ####
 # month with highest incidence, month with highest undiagnosed prevalence (symptomatic and sputum+, and overall), and 
# remaining expected duration of prevalent disease at each of 3, 6, 9, 12, 15, and 18mo
run_cohort_features <- FALSE
if (run_cohort_features)
{  cohort_features <- 
    data.frame(
            accepted_subclinical = numeric(N_samples),
            month_with_highest_incidence = numeric(N_samples),
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
    
    cohort_features[n, "accepted_subclinical"] <- check_subclinical(cohort, cohort_params[n,])
    
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
  
  head(cohort_features)
  cohort_features[,] %>% summarise_all(median)
  cohort_features[,] %>% summarise_all(mean)
  cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) median(x, na.rm=T))
  cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.25, na.rm=T))
  cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.75, na.rm=T)) # for the abstract
}


#### More intervnetion setup ####

intervention_names = c("guidelines", "earlier", "frequent", "sputum", "sputummoretargeted", "counseling", "prevention")

results <- list()
for (name in intervention_names) {
  results[[name]] <- 
    data.frame(detections = numeric(N_samples),
               mean_symptom_days = numeric(N_samples),
               symptomatic_months_soc = numeric(N_samples),
               infectious_months_soc = numeric(N_samples),
               symptomatic_months_averted = numeric(N_samples), 
               infectious_months_averted = numeric(N_samples))

}

screening_design_guidelines <- 
  data.frame(
    "counseling_coverage" = 0, 
    "prevention_coverage" = 0,
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## earlier screening
screening_design_earlier <- 
  data.frame(
    "counseling_coverage" = 0, 
    "prevention_coverage" = 0,
    "timing_months" = c(3, 6, 12),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# counseling, with passive detection impact (assumed same cost as 12 +18 month visits)
screening_design_counseling <- 
  data.frame(
    "counseling_coverage" = 1,
    "prevention_coverage" = 0,
    "timing_months" = c(6),
    "target_coverage" = c(1),
    "screening_method" = c("symptoms"),
    "screening_location" = c("home")
  ) %>% 
  arrange(timing_months)


# For the rest, we'll define within the loop because target coverage depends on cost

#### For each sample, simulate a cohort and estimate the impact of the aftermath intervention: ####
# Running as loop, could parallelize for efficiency ***

for (n in 1:N_samples)
{
  #### More parameter-specific setup ####
  cohort <- create_cohort(cohort_params[n,])
  # if doesn't validate on subclinical TB, abort and fill NAs
  if (!check_subclinical(cohort, cohort_params[n,])) {
    # change the nth element of each item in results to NA
    for (name in intervention_names) 
      results[[name]][n,] <- rep(NA, ncol(results[[name]]))
    next
  }
  
  intervention_parameters <- list(
    coverage = list("phone" = cohort_params$coverage_phone[n], 
                    "home" = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]),
    sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home[n], 
                                "phone" = cohort_params$sensitivity_symptoms_home[n] * cohort_params$sensitivity_symptoms_phone_reduction[n]),
    success_sputum = list("home" = cohort_params$success_sputum_home[n], 
                          "phone" = cohort_params$success_sputum_home[n] * cohort_params$success_sputum_phone_reduction[n]))
  
  #### Define additional cost-equivalent intervention options ####

  ## more frequent but targeted. 
  targeting_frequent = (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) / # cost of aftermath version
                  (1 + 5*cohort_params$initial_contact_cost_phone_factor[n])
  screening_design_frequent <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(3, 6, 9, 12, 15, 18),
      "target_coverage" = rep(targeting_frequent, 6),
      "screening_method" = rep("symptoms", 6),
      "screening_location" = c("home", rep("phone",5))
    ) %>% 
    arrange(timing_months)
  
  ## sputum at 6mo for highest risk, in place of 12 and 18mo visits for all (the untargeted still get a 6mo visit)
  # (should modify to include confirmatory testing or false positives?)
  # (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) = #cost of aftermath version
  #                 (cohort_params$initial_contact_cost_home[n] * (1 + targeting*2*cohort_params$initial_contact_cost_phone_factor[n]) +
  #                    targeting*cohort_params$sputum_test_cost[n])
  targeting_sputum =  (cohort_params$initial_contact_cost_home[n] * (2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
                  (cohort_params$initial_contact_cost_home[n] * 2*cohort_params$initial_contact_cost_phone_factor[n] +
                     cohort_params$sputum_test_cost[n])
    screening_design_sputum <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(6, 6, 12, 18),
      "target_coverage" = c(targeting_sputum, 1, targeting_sputum, targeting_sputum),
      "screening_method" = c("sputum", "symptoms", "symptoms", "symptoms"),
      "screening_location" = c("home", "home", "phone", "phone")
    ) %>% 
    arrange(timing_months)
  
  # aftermath for a fraction, plus sputum at first visit (and the others get nothing)
    targeting_sputummoretargeted =  (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
      (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) +
         cohort_params$sputum_test_cost[n])
  screening_design_sputummoretargeted <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(6, 12, 18),
      "target_coverage" = c(1,1,1) * targeting_sputummoretargeted,
      "screening_method" = c("both", "symptoms", "symptoms"),
      "screening_location" = c("home", "phone", "phone")
    ) %>% 
    arrange(timing_months)


  # prevention of recurrence
  targeting_prevention = 
    (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
    (cohort_params$prevention_cost[n])
  screening_design_prevention <- data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = targeting_prevention,
      "timing_months" = c(6),
      "target_coverage" = c(0),
      "screening_method" = c("symptoms"),
      "screening_location" = c("home")
    ) %>% 
    arrange(timing_months)
  
  

  ##### Run the interventions ####
  outputs <- list()
  for (intervention in intervention_names)
    outputs[[intervention]] <- apply_intervention(cohort, get(paste0("screening_design_", intervention)), 
                                                  intervention_parameters, cohort_params[n,])  
  
  # impact
  (detected <- lapply(outputs, function(x) sum(!is.na(x$detection_timing))))
  (symptomdays <- lapply(outputs, function(x) ((time_with_tb(x) %>% filter(scenario == "screening") %>% select(value))/sum(x$TB))[1,]))
  (time <- lapply(outputs, function(x) time_with_tb(x) %>% filter(scenario == "soc") %>% 
                    mutate(months = value/30) %>%
                    tibble::column_to_rownames('outcome') %>% select(months)))
  (impact <- lapply(outputs, function(x) (time_with_tb(x)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                                            mutate(months_averted = (soc - screening)/30)) %>% 
                      tibble::column_to_rownames('outcome')  %>% select(months_averted)))
  
  
  
  for (name in intervention_names)
    results[[name]] [n,] <- c(detected[[name]], symptomdays[[name]], time[[name]]$months[1], time[[name]]$months[2], 
                              impact[[name]]$months_averted[1], impact[[name]]$months_averted[2])
  # , 
  # cost,
  # cost/impact$months_averted[1], cost/impact$months_averted[2])
}


#### Look at results ####
# collate results for a table:
# for all dataframe elements of list "results", report the mean and 25th and 75th percentiles of each column

lapply(results, function(x) 
  # mean of each columns of table
  x %>% 
    mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
           inf_reduction = infectious_months_averted/infectious_months_soc*100) %>%
    select(detections, mean_symptom_days, sx_reduction, inf_reduction) %>%
    
    summarise_all(function(y) paste0(round(median(y, na.rm = T),0), " (", 
                                     round(quantile(y, 0.25, na.rm = T),0), ", ", 
                                     round(quantile(y, 0.75, na.rm = T),0), ")"))) %>%
  # convert list elements to rows of a single dataframe, with intervention_names as rownames
  bind_rows() %>%
  mutate(intervention = intervention_names) %>%
  # make intervention the first column, and mean_symptom_days the second column
  select(intervention, mean_symptom_days, everything()) %>%
  write_clip() 

  







# summarize results for aftermath intervention
lapply(results, function(x) summary(x$cost_per_symptomatic_month_averted))
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





