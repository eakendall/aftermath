library(tidyverse)
library(binom)
library(sensitivity)
library(conflicted)

conflicts_prefer(dplyr::summarize)

N_cohort <- 50000
N_samples <- 5000 # changed temporarily from 5000

date <- "20260529"

set.seed(12345)

#### Load fixed empirical inputs ####

fixed_empirical_inputs <- readRDS("outputs/fixed_empirical_inputs.rds")

#### Parameter ranges ####

cohort_param_ranges <- list(
  N = N_cohort,
  
  #### Empirical / semi-empirical parameters ####
  incidence_18mo_multiplier = c(0.8, 1.2),
  proportion_micro_pos = c(0.42, 0.64),
  auc = 0.69,
  
  #### Reported symptom duration distribution ####
  symptom_duration_meanlog_reported = 2.76,
  symptom_duration_sdlog_reported = c(0.6, 0.8),
  
  #### Symptom duration scaling ####
  # Reported duration is this fraction of true duration.
  # true_duration = reported_duration / reported_fraction_of_true_symptom_duration
  reported_fraction_of_true_symptom_duration = c(1/6, 1/2),
  programmatic_symptom_duration_factor = c(1, 2),
  
  #### Natural history: first detectable recurrent-TB state ####
  # A small proportion of micropositive recurrent TB is already NAAT+/sputum+ at treatment completion.
  # This will be assigned explicitly in create_cohort(), rather than arising from pmax().
  proportion_micropos_subclinical_at_eot = c(0.00, 0.10),
  
  # Among micropositive recurrent TB episodes, probability that the first detectable state
  # is sputum+/NAAT+ before symptom-screen positivity.
  proportion_micropos_sputum_first = c(0.60, 0.80),
  
  # Among sputum-first episodes, duration from sputum+/NAAT+ onset to symptom onset.
  # Mean is set relative to mean true symptomatic duration downstream.
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2),
  duration_subclinical_cv = c(0.5, 1.5),
  
  #### Subclinical prevalence calibration criteria ####
  # Six-month: symptom-negative, NAAT+/sputum+ recurrent TB prevalence
  # among the full post-TB cohort.
  subclinical_6m_amongcohort_min = 0.004,
  subclinical_6m_amongcohort_max = 0.017,
  
  #### Intervention parameters ####
  coverage_phone = c(0.9, 1.0),
  coverage_home_reduction = c(0.75, 0.95),
  sensitivity_symptoms_home = c(0.7, 0.9),
  sensitivity_symptoms_phone_reduction = c(0.6, 0.82),
  success_sputum_home = c(0.9, 1.0),
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
  paste0(
    "outputs/probabilistic_parameter_draws_pre_weibull_",
    date,
    ".csv"
  )
)

saveRDS(
  cohort_params,
  paste0(
    "outputs/probabilistic_parameter_draws_pre_weibull_",
    date,
    ".rds"
  )
)