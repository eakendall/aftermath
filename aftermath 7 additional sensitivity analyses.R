# plot consistency of simulated cohorts with aftermath data

# plot consistency of simulated cohorts with aftermath data

# time to diagnosis

# plot survival_dataset_truncated with confidence ribbon for time to event
library(ggsurvfit)
library(tidycmprsk)
library(survival)
library(stats)

# compare with cumulative incidence curves -- will need to derive each separately and overlay:
data <- read.csv("../Data June 2025 from Aye/variables_trialdata1806.csv")
data <- data %>% 
  filter(!(end_reason=="death" & tb_related_death == "Unknown")) %>%
  mutate(
    recurrence = case_when(
      end_reason == "TB recurrence" ~ 1,
      end_reason == "death" & tb_related_death == "Yes" ~ 1,
      end_reason == "death" & tb_related_death == "Probable" ~ 1,
      end_reason == "death" & tb_related_death == "No" ~ 0,
      end_reason == "completion" ~ 0,
      end_reason == "LFU" ~ 0,
      TRUE ~ 0),
    end_days = case_when(
      end_reason == "death" & recurrence == 1 ~ interval(txcompl_date, death_date)/days(1),
      TRUE ~ txcompl_endreason_days
    )
  )

data <- data %>% mutate(micropos = case_when(ev_micro_test == "Microbiological confirmation" & ev_TBtype == "Pulmonary TB (PTB)" ~ 1, TRUE ~ 0),
                        clindx = case_when(ev_micro_test == "Clinical confirmation" | ev_TBtype == "Extra pulmonary TB (EPTB)" ~ 1 , TRUE ~ 0))

survival_dataset <- data %>% 
  select(record_id, term_reason, end_reason, txcompl_endreason_days, micropos, clindx) %>% 
  mutate(
    event = case_when(
      end_reason == "TB recurrence" ~ 1,
      TRUE ~ 0
    ),
    time = txcompl_endreason_days
  )
survival_dataset_truncated <- survival_dataset; 
survival_dataset_truncated$time[survival_dataset$time > 540] <- 540
survival_dataset_truncated$event[survival_dataset$time > 540] <- 0

survival_dataset_truncated$eventtype <- case_when(
  survival_dataset_truncated$micropos == 1 & survival_dataset_truncated$event == 1 ~ 1,
  survival_dataset_truncated$clindx == 1 & survival_dataset_truncated$event == 1 ~ 2,
  TRUE ~ 0
)

cuminc_fit <- cuminc(Surv(time, as.factor(event)) ~ 1, data = survival_dataset_truncated)

(dataplot <- dataplotbase <- ggcuminc(cuminc_fit, outcome = 1) +
    add_confidence_interval() +
    scale_ggsurvfit())


# overlay simulated timing
# For each set of cohort_params, simulate a cohort and plot the distribution of diagnosis_routing
set.seed(12345)

for (n in 1: min(500, nrow(cohort_params))) {
  params <- cohort_params[n, ]
  
  cohort <- create_cohort(params)
  
  if(check_subclinical(cohort, params)) 
  {
    cohort$diagnosis_routine[is.na(cohort$diagnosis_routine)] <- 10000 # set to beyond plot's x axis limit if no recurrence, so this will be cumulative % across full cohort not just recurrences
    
    # plot cumulative distribution of diagnosis_routine 
    dataplot <- dataplot + 
      stat_ecdf(data= cohort, aes(x = diagnosis_routine), col = "blue", alpha = 0.1) + coord_cartesian(xlim = c(0, 18*30), ylim= c(0, 0.12))
  }
}

dataplot + coord_cartesian(xlim = c(0,30*18), ylim=c(0,0.12)) + 
  scale_x_continuous(breaks = seq(0, 18*30, by = 3*30), labels = seq(0, 18, by = 3)) + 
  xlab("Months since treatment completion") + 
  ylab("Cumulative notifications of recurrent TB") 


#### Compare timing of clinical and micro diagnoses: ####

ggplot(data %>% filter(end_reason == "TB recurrence"), 
       aes(x = txcompl_endreason_days/30, color = as.factor(micropos))) +
  stat_ecdf() +
  xlab("Months to TB recurrence diagnosis") +
  ylab("Cumulative proportion of identified TB recurrences during Aftermath study") +
  # ggtitle("Timing of TB recurrence diagnosis by diagnostic method") + 
  # set color legend to have title of "TB type" and levels of "Micro+ pulmonary" for 1 and "other" for 0
  scale_color_discrete(name = "TB type", labels = c("Micro+ pulmonary", "Other")) + 
  theme_minimal()
# Looks pretty similar.




#### Supplement: retained vs rejected simulations (Fig Sx) ####
# Compare sampled parameters and key natural-history outcomes
# between simulations retained vs rejected by check_subclinical()

sampled_parameters <- cohort_params %>%
  summarise(across(everything(), ~ n_distinct(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "n_distinct") %>%
  filter(n_distinct > 1) %>%
  pull(parameter)

exclude_parameters <- c(
  "draw", "N", "case_fatality", 
  "recurrence_time_mean_multiplier", "recurrence_time_cv_multiplier", "incidence_18mo_multiplier",
  "coverage_phone", "coverage_home_reduction",
  "sensitivity_symptoms_home", "sensitivity_symptoms_phone_reduction",
  "success_sputum_home", "success_sputum_phone_reduction",
  "initial_contact_cost_home", "initial_contact_cost_home_vs_phone_factor",
  "sputum_test_cost", "symptom_prevalence_nontb",
  "initial_contact_cost_phone", "coverage_home",
  "sensitivity_symptoms_phone", "success_sputum_phone",
  "recurrence_shape", "recurrence_scale", "recurrence_time_mean", "recurrence_time_cv",
  "probability_dx540_given_recur", "probability_ever_recur", "prevention_cost",
  "objective_value",
  "mean_symptom_duration_by_540_sim", "median_symptom_duration_by_540_sim", 
  "prop_onset_le_7_among_dx540_sim",  "prop_onset_le_30_among_dx540_sim"
)

sampled_parameters <- setdiff(sampled_parameters, exclude_parameters)

nice_names <- c(
  incidence_18mo = "18-month recurrence incidence",
  probability_ever_recur = "Probability of eventual recurrence",
  proportion_micro_pos = "Proportion NAAT-positive at diagnosis",
  symptom_duration_sdlog_reported = "SD log reported symptom duration",
  reported_fraction_of_true_symptom_duration = "Reported fraction of true symptom duration",
  programmatic_symptom_duration_factor = "Increase in symptom duration under programmatic conditions",
  proportion_ever_subclinical = "Proportion with asymptomatic  NAAT+ period",
  duration_ratio_subclinical_symptomatic = "Asymptomatic:symptomatic NAAT+ time ratio",
  duration_subclinical_cv = "Coefficient of variation, asymptomatic TB duration",
  subclinical_baseline_amongTB_max = "Maximum subclinical prevalence at treatment completion",
  subclinical_6m_amongcohort_min = "Minimum 6-month asymptomatic TB prevalence",
  subclinical_6m_amongcohort_max = "Maximum 6-month asymptomatic TB prevalence",
  median_dx_by_540_sim= "Median days to diagnosis (if diagnosed by 18mo)",
  p25_dx_by_540_sim = "25th quantile, days to diagnosis (if diagnosed by 18mo)",
  p75_dx_by_540_sim = "75th quantile, days to diagnosis (if diagnosed by 18mo)",
  prop_dx_le_90_among_dx540_sim = "Proportion diagnosed by day 90 (if diagnosed by 18mo)",
  prop_dx_le_360_among_dx540_sim = "Proportion diagnosed by day 360 (if diagnosed by 18mo)"
  # subclinical_baseline_among_micropos =   "Baseline subclinical prevalence among micropositive recurrences",
  # subclinical_6mo_amongcohort = "6-month subclinical prevalence"
)

key_outcomes <- c(
  "cum_notifs_18mo",
  "prev_sx_6mo",
  "prev_inf_6mo",
  "prev_subclinical_6mo",
  "cum_months_sx_24mo",
  "cum_months_inf_24mo",
  "cum_months_subclinical_24mo"
)

accepted_status <- tibble(
  sim_id = seq_len(nrow(cohort_params)),
  accepted_subclinical = cohort_features$accepted_subclinical
)

param_compare_data <- cohort_params %>%
  mutate(sim_id = row_number()) %>%
  select(sim_id, all_of(sampled_parameters)) %>%
  pivot_longer(
    cols = -sim_id,
    names_to = "quantity",
    values_to = "value"
  ) %>%
  mutate(type = "Sampled parameter") %>%
  bind_rows(
    cohort_features %>%
      mutate(sim_id = row_number()) %>%
      select(sim_id, all_of(key_outcomes)) %>%
      pivot_longer(
        cols = -sim_id,
        names_to = "quantity",
        values_to = "value"
      ) %>%
      mutate(type = "Model outcome")
  ) %>%
  left_join(accepted_status, by = "sim_id")

param_compare_summary <- param_compare_data %>%
  group_by(type, quantity, accepted_subclinical) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q25 = quantile(value, 0.25, na.rm = TRUE),
    q75 = quantile(value, 0.75, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  param_compare_summary,
  "outputs/subclinical_acceptance_parameter_outcome_summary.csv"
)

nice_order <- nice_names[names(nice_names) %in% sampled_parameters] %>%
  unname() %>%
  stringr::str_wrap(width = 22)

param_compare_data <- param_compare_data %>%
  mutate(
    nice_quantity = recode(quantity, !!!nice_names),
    nice_quantity = stringr::str_wrap(nice_quantity, width = 22),
    nice_quantity = factor(nice_quantity, levels = nice_order),
    accepted_label = ifelse(accepted_subclinical, "Retained", "Rejected")
  )
#maintain ordering of nice_names in plot
param_compare_plot <- param_compare_data %>%
  
  filter(type == "Sampled parameter") %>%
  ggplot(aes(x = accepted_label, y = value, color = accepted_label)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(. ~ nice_quantity, scales = "free_y") +
  xlab("Retained by subclinical calibration") +
  ylab("Value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

ggsave(
  "outputs/subclinical_acceptance_parameter_outcome_distributions.pdf",
  param_compare_plot,
  width = 14,
  height = 10
)



#### Quintile violin sensitivity plots ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

prcc_outcomes <- tibble(
  sim_id = seq_len(nrow(cohort_params)),
  
  guidelines_sx_reduction =
    results$guidelines$symptomatic_months_averted /
    results$guidelines$symptomatic_months_soc * 100,
  
  guidelines_inf_reduction =
    results$guidelines$infectious_months_averted /
    results$guidelines$infectious_months_soc * 100,
  
  earlier_three_vs_guidelines_sx =
    results$earlier_three$symptomatic_months_averted /
    results$guidelines$symptomatic_months_averted,
  
  earlier_three_vs_guidelines_inf =
    results$earlier_three$infectious_months_averted /
    results$guidelines$infectious_months_averted,
  
  earlier_three_sputum_vs_guidelines_sx =
    results$earlier_three_sputum$symptomatic_months_averted /
    results$guidelines$symptomatic_months_averted,
  
  earlier_three_sputum_vs_guidelines_inf =
    results$earlier_three_sputum$infectious_months_averted /
    results$guidelines$infectious_months_averted
)

outcome_nice_names <- c(
  guidelines_sx_reduction = "Guidelines vs SOC:\n% symptomatic time averted",
  guidelines_inf_reduction = "Guidelines vs SOC:\n% infectious time averted",
  earlier_three_vs_guidelines_sx = "3/6/9m symptoms vs guidelines:\nrelative symptomatic time averted",
  earlier_three_vs_guidelines_inf = "3/6/9m symptoms vs guidelines:\nrelative infectious time averted",
  earlier_three_sputum_vs_guidelines_sx = "3/6/9m + 3m micro vs guidelines:\nrelative symptomatic time averted",
  earlier_three_sputum_vs_guidelines_inf = "3/6/9m + 3m micro vs guidelines:\nrelative infectious time averted"
)

parameters_to_plot <- prcc_results %>%
  pull(variable) %>%
  unique()
parameters_to_plot2 <- prcc2_results %>%
  pull(variable) %>%
  unique()
parameters_to_plot_combined <- union(parameters_to_plot, parameters_to_plot2)


plot_data <- cohort_params %>%
  mutate(sim_id = row_number()) %>%
  select(sim_id, all_of(parameters_to_plot_combined)) %>%
  pivot_longer(
    cols = -sim_id,
    names_to = "parameter",
    values_to = "parameter_value"
  ) %>%
  group_by(parameter) %>%
  mutate(
    quintile = ntile(parameter_value, 5),
    quintile_group = case_when(
      quintile == 1 ~ "Lowest quintile",
      quintile == 5 ~ "Highest quintile",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(quintile_group)) %>%
  left_join(prcc_outcomes, by = "sim_id") %>%
  pivot_longer(
    cols = starts_with("guidelines") | starts_with("earlier_three"),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  mutate(
    parameter_label = recode(parameter, !!!param_names, .default = parameter),
    parameter_label = str_wrap(parameter_label, width = 22),
    outcome_label = recode(outcome, !!!outcome_nice_names, .default = outcome),
    outcome_label = str_wrap(outcome_label, width = 18),
    quintile_group = factor(
      quintile_group,
      levels = c("Highest quintile", "Lowest quintile")
    )
  )

density_data <- plot_data %>%
  group_by(parameter_label, outcome_label, quintile_group) %>%
  group_modify(~{
    if (nrow(.x) < 5 || length(unique(.x$value)) < 2) return(tibble())
    d <- density(.x$value, na.rm = TRUE)
    tibble(y = d$x, dens = d$y)
  }) %>%
  ungroup() %>%
  group_by(parameter_label, outcome_label, quintile_group) %>%
  mutate(
    dens_scaled = dens / max(dens, na.rm = TRUE) * 0.42,
    x = case_when(
      quintile_group == "Highest quintile" ~ 1 - dens_scaled,
      quintile_group == "Lowest quintile" ~ 1 + dens_scaled
    )
  ) %>%
  ungroup()

box_data <- plot_data %>%
  mutate(
    x_box = case_when(
      quintile_group == "Highest quintile" ~ 0.92,
      quintile_group == "Lowest quintile" ~ 1.08
    )
  )

quintile_sensitivity_plot <- ggplot() +
  geom_polygon(
    data = density_data,
    aes(
      x = x,
      y = y,
      group = interaction(parameter_label, outcome_label, quintile_group),
      fill = quintile_group
    ),
    alpha = 0.75,
    color = "grey40"
  ) +
  geom_boxplot(
    data = box_data,
    aes(
      x = x_box,
      y = value,
      group = interaction(parameter_label, outcome_label, quintile_group),
      fill = quintile_group
    ),
    width = 0.05,
    outlier.shape = NA,
    alpha = 0.55
  ) +
  facet_grid(outcome_label ~ parameter_label, scales = "free_y") +
  scale_fill_manual(
    values = c(
      "Highest quintile" = "steelblue4",
      "Lowest quintile" = "skyblue3"
    ),
    name = "Parameter value"
  ) +
  scale_x_continuous(limits = c(0.5, 1.5), breaks = NULL) +
  xlab("") +
  ylab("Outcome value") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 7),
    strip.text.y = element_text(size = 8)
  )

ggsave(
  "outputs/quintile_sensitivity_prcc_outcomes_half_violin.pdf",
  quintile_sensitivity_plot,
  width = 18,
  height = 10
)

#### Compare intervention outputs with vs without subclinical calibration filter ####

key_interventions <- c(
  "guidelines",
  "earlier_three",
  "earlier_three_sputum",
  "four_visits_36912",
  "frequent"
)

output_nice_names <- c(
  proportion_detected = "Proportion of recurrences detected",
  sx_reduction = "Symptomatic months averted",
  inf_reduction = "Infectious months averted",
  cost_per_patient = "Cost per TB survivor",
  cost_per_detection = "Cost per case detected"
)

intervention_nice_names <- c(
  guidelines = "Guidelines\n6, 12, 18m",
  earlier_three = "3, 6, 9m\nsymptoms",
  earlier_three_sputum = "3, 6, 9m\n+ 3m micro",
  four_visits_36912 = "3, 6, 9, 12m\nsymptoms",
  frequent = "3, 6, 9, 12, 15, 18m\nsymptoms"
)

prep_results_for_filter_comparison <- function(results_obj, filter_label) {
  bind_rows(
    lapply(names(results_obj), function(intervention) {
      results_obj[[intervention]] %>%
        mutate(
          intervention = intervention,
          filter_label = filter_label,
          sim_id = row_number(),
          proportion_detected = detections / recurrences,
          sx_reduction = symptomatic_months_averted / symptomatic_months_soc,
          inf_reduction = infectious_months_averted / infectious_months_soc,
          cost_per_patient = cost / N_cohort,
          cost_per_detection = cost / detections
        )
    })
  )
}

filter_compare_data <- bind_rows(
  prep_results_for_filter_comparison(results_unfiltered, "All simulations"),
  prep_results_for_filter_comparison(results_filtered, "Retained simulations")
) %>%
  filter(intervention %in% key_interventions) %>%
  select(
    filter_label, sim_id, intervention,
    proportion_detected, sx_reduction, inf_reduction,
    cost_per_patient, cost_per_detection
  ) %>%
  pivot_longer(
    cols = c(proportion_detected, sx_reduction, inf_reduction,
             cost_per_patient, cost_per_detection),
    names_to = "outcome",
    values_to = "value"
  ) %>%
  mutate(
    intervention_label = recode(intervention, !!!intervention_nice_names),
    outcome_label = recode(outcome, !!!output_nice_names),
    outcome_label = stringr::str_wrap(outcome_label, width = 28),
    filter_label = factor(filter_label, levels = c("All simulations", "Retained simulations"))
  )

filter_compare_summary <- filter_compare_data %>%
  group_by(filter_label, intervention, outcome) %>%
  summarise(
    median = median(value, na.rm = TRUE),
    q025 = quantile(value, 0.025, na.rm = TRUE),
    q975 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  filter_compare_summary,
  "outputs/filter_vs_no_filter_key_outputs_summary.csv"
)

filter_compare_plot <- filter_compare_data %>%
  ggplot(aes(x = filter_label, y = value, fill = filter_label)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  facet_grid(outcome_label ~ intervention_label, scales = "free_y") +
  xlab("") +
  ylab("Value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8)
  )

ggsave(
  "outputs/filter_vs_no_filter_key_outputs_violin.pdf",
  filter_compare_plot,
  width = 13,
  height = 9
)