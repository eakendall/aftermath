library(tidyverse)
library(binom)
library(sensitivity)
library(conflicted)

conflicts_prefer(dplyr::summarize)

source("aftermath cohort sim functions.R")

N_cohort <- 10000
N_samples <- 5000

set.seed(12345)

#### Load fixed empirical inputs ####

fixed_empirical_inputs <- readRDS("outputs/fixed_empirical_inputs.rds")

#### Parameter ranges ####
# Parameters that depend on true symptom duration / symptom onset timing are NOT
# estimated here. Instead, draw uncertainty multipliers to apply downstream.

cohort_param_ranges <- list(
  N = N_cohort,
  
  #### To be applied downstream after Weibull fitting ####
  recurrence_time_mean_multiplier = c(0.8, 1.2),
  recurrence_time_cv_multiplier = c(0.8, 1.2),
  
  #### Empirical / semi-empirical parameters ####
  incidence_18mo_multiplier = c(0.8, 1.2), # now refers to incidence of diagnosis by 540d not onset. Wider uncertainty than estimated from aftermath directly, but that's appropriate because of simplifications in weibull fitting.
  proportion_micro_pos = c(0.42, 0.64),
  # binom.agresti.coull(fixed_empirical_inputs$micropos_n, fixed_empirical_inputs$recurrence_n)
  auc = 0.69, #c(0.56, 0.77), # not using for manuscript
  
  #### Reported symptom duration distribution ####
  symptom_duration_meanlog_reported = 2.76,
  # fixed_empirical_inputs$reported_symptom_duration_recurrence$meanlog
  symptom_duration_sdlog_reported = c(0.6, 0.8),
  # fixed_empirical_inputs$reported_symptom_duration_recurrence$sdlog
  
  ####Symptom_underestimation_factor ####
  # Reported duration is this fraction of true duration.
  # true_duration = reported_duration / reported_fraction_of_true_symptom_duration
  reported_fraction_of_true_symptom_duration = c(1/6, 1/2),
  ### Aftermath faster care-seeking factor ### 
  programmatic_symptom_duration_factor = c(1, 2),
  
  #### Subclinical natural history ####
  proportion_ever_subclinical = c(0.6, 0.9),
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2),
  duration_subclinical_cv = c(0.5, 1.5),
  subclinical_baseline_amongTB_max = 0.20,
  subclinical_6m_amongcohort_min = 0.004,
  subclinical_6m_amongcohort_max = 0.017,
  
  #### Intervention parameters ####
  coverage_phone = c(0.9, 1.0),
  coverage_home_reduction = c(0.75, 0.95),
  sensitivity_symptoms_home = c(0.85, 0.95),
  sensitivity_symptoms_phone_reduction = c(0.6, 0.82),
  # fixed_empirical_inputs$relative_symptom_reporting_phone
  success_sputum_home = 1, #c(0.85, 0.95), # assuming can swab tongue if participating
  success_sputum_phone_reduction = c(0.6, 1.0),
  
  home_visit_passive_detection_impact = 1,
  intentional_counseling_passive_detection_impact = 0.8,
  intentional_counseling_passive_detection_duration = 180,
  prevention_efficacy = 0.6,
  
  #### Costs and non-TB symptoms ####
  initial_contact_cost_home = c(2, 4),
  initial_contact_cost_home_vs_phone_factor = c(4, 12),
  sputum_test_cost = c(5, 9),
  symptom_prevalence_nontb = c(0.1, 0.3),
  prevention_cost = c(7.24, 9.14),
  
  #### Mortality ####
  case_fatality = c(0.05, 0.15)
)

#### Draw parameters ####

draw_param <- function(x, n) {
  if (length(x) == 2) {
    runif(n, min = x[1], max = x[2])
  } else {
    rep(x, n)
  }
}

cohort_params <- lapply(
  cohort_param_ranges,
  draw_param,
  n = N_samples
) %>%
  as.data.frame() %>%
  mutate(draw = row_number()) %>%
  relocate(draw)

#### Derived parameters that do NOT require Weibull fitting ####

cohort_params <- cohort_params %>%
  mutate(
    initial_contact_cost_phone =
      initial_contact_cost_home / initial_contact_cost_home_vs_phone_factor,
    
    coverage_home =
      coverage_phone * coverage_home_reduction,
    
    sensitivity_symptoms_phone =
      sensitivity_symptoms_home * sensitivity_symptoms_phone_reduction,
    
    success_sputum_phone =
      success_sputum_home * success_sputum_phone_reduction
  )

#### Save ####

dir.create("outputs", showWarnings = FALSE)

write_csv(
  cohort_params,
  "outputs/probabilistic_parameter_draws_pre_weibull.csv"
)

saveRDS(
  cohort_params,
  "outputs/probabilistic_parameter_draws_pre_weibull.rds"
)