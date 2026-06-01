# plot consistency of simulated cohorts with aftermath data

library(tidyverse)
library(lubridate)
library(survival)

#### Empirical data: cumulative hazard to recurrent TB diagnosis/death by 540d ####

data <- read.csv("../Data June 2025 from Aye/variables_trialdata1806.csv")

data <- data %>%
  filter(!(end_reason == "death" & tb_related_death == "Unknown")) %>%
  mutate(
    recurrence = case_when(
      end_reason == "TB recurrence" ~ 1,
      end_reason == "death" & tb_related_death %in% c("Yes", "Probable") ~ 1,
      TRUE ~ 0
    ),
    end_days = case_when(
      end_reason == "death" & recurrence == 1 ~
        as.numeric(interval(ymd(txcompl_date), ymd(death_date)) / days(1)),
      TRUE ~ txcompl_endreason_days
    ),
    micropos = case_when(
      ev_micro_test == "Microbiological confirmation" &
        ev_TBtype == "Pulmonary TB (PTB)" ~ 1,
      TRUE ~ 0
    ),
    clindx = case_when(
      ev_micro_test == "Clinical confirmation" |
        ev_TBtype == "Extra pulmonary TB (EPTB)" ~ 1,
      TRUE ~ 0
    )
  )

survival_dataset_truncated <- data %>%
  transmute(
    record_id,
    micropos,
    clindx,
    time = pmin(end_days, 540),
    event = if_else(recurrence == 1 & end_days <= 540, 1, 0)
  ) %>%
  filter(is.finite(time), !is.na(event))

cumhaz_fit <- survfit(
  Surv(time, event == 1) ~ 1,
  data = survival_dataset_truncated
)

cumhaz_summary <- summary(cumhaz_fit)

cumhaz_data <- tibble(
  time = cumhaz_summary$time,
  surv = cumhaz_summary$surv,
  surv_lower = cumhaz_summary$lower,
  surv_upper = cumhaz_summary$upper,
  cumhaz = -log(surv),
  cumhaz_lower = -log(surv_upper),
  cumhaz_upper = -log(surv_lower)
)

emp_540 <- summary(cumhaz_fit, times = 540, extend = TRUE) %>%
  with(tibble(
    time = time,
    cumhaz = -log(surv),
    cumhaz_lower = -log(upper),
    cumhaz_upper = -log(lower),
    cuminc = 1 - surv,
    cuminc_lower = 1 - upper,
    cuminc_upper = 1 - lower
  ))

print(emp_540)

#### Simulated traces: cumulative hazard to routine diagnosis ####

make_sim_cumhaz <- function(params, sim_id, limit_days = 540) {
  cohort <- create_cohort(params)
  N_full <- nrow(cohort)
  
  cohort %>%
    filter(TB == 1, is.finite(diagnosis_routine), diagnosis_routine <= limit_days) %>%
    arrange(diagnosis_routine) %>%
    mutate(
      sim_id = sim_id,
      time = diagnosis_routine,
      cuminc = row_number() / N_full,
      cumhaz = -log(1 - cuminc)
    ) %>%
    select(sim_id, time, cumhaz)
}

set.seed(12345)

n_to_plot <- min(100, nrow(cohort_params_main))

sim_cumhaz_data <- map_dfr(
  seq_len(n_to_plot),
  function(n) {
    params <- cohort_params_main[n, ]
    cohort_trace <- tryCatch(
      make_sim_cumhaz(params, sim_id = n, limit_days = 540),
      error = function(e) NULL
    )
    cohort_trace
  }
)

#### Plot ####

dataplot <- ggplot() +
  geom_ribbon(
    data = cumhaz_data,
    aes(x = time, ymin = cumhaz_lower, ymax = cumhaz_upper),
    alpha = 0.25
  ) +
  geom_step(
    data = cumhaz_data,
    aes(x = time, y = cumhaz),
    linewidth = 0.8
  ) +
  geom_step(
    data = sim_cumhaz_data,
    aes(x = time, y = cumhaz, group = sim_id),
    color = "blue",
    alpha = 0.2,
    linewidth = 0.5
  ) +
  coord_cartesian(xlim = c(0, 540), ylim = c(0, 0.13)) +
  scale_x_continuous(
    breaks = seq(0, 540, by = 90),
    labels = seq(0, 18, by = 3)
  ) +
  theme_minimal() +
  xlab("Months since treatment completion") +
  ylab("Cumulative hazard of recurrent TB notification or TB-attributed death")

dataplot


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
# between simulations retained vs rejected

sampled_parameters <- cohort_params %>%
  summarise(across(everything(), ~ n_distinct(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "parameter", values_to = "n_distinct") %>%
  filter(n_distinct > 1) %>%
  pull(parameter)

exclude_parameters <- c(
  "draw", "N", "case_fatality", "incidence_18mo_multiplier",
  "coverage_phone", "coverage_home_reduction",
  "sensitivity_symptoms_home", "sensitivity_symptoms_phone_reduction",
  "success_sputum_home", "success_sputum_phone_reduction",
  "initial_contact_cost_home", "initial_contact_cost_home_vs_phone_factor",
  "sputum_test_cost", "symptom_prevalence_nontb",
  "initial_contact_cost_phone", "coverage_home",
  "sensitivity_symptoms_phone", "success_sputum_phone","prevention_cost",
  "recurrence_shape", "recurrence_scale", "recurrence_time_mean", "recurrence_time_cv",
  "probability_dx540_given_recur", "probability_ever_recur", 
  # "median_dx_by_540_sim"  ,                    
  "p25_dx_by_540_sim"      ,                   
  "p75_dx_by_540_sim"       ,                  
  "prop_dx_le_90_among_dx540_sim",             
  "prop_dx_le_360_among_dx540_sim",
  "objective_value",
  "mean_symptom_duration_by_540_sim",
  "median_symptom_duration_by_540_sim",
  "prop_first_event_le_7_among_dx540_sim",
  "prop_first_event_le_30_among_dx540_sim",
  "prop_symptom_onset_le_7_among_dx540_sim",
  # "prop_symptom_onset_le_30_among_dx540_sim",
  "prop_sputum_first_by_540_sim",
  "mean_subclinical_duration_by_540_sim",
  "optim_convergence_code",                    
  "valid_weibull_numeric",                     
  "valid_weibull_fit"
)

sampled_parameters <- setdiff(sampled_parameters, exclude_parameters)

nice_names <- c(
  incidence_18mo = "18-month recurrence incidence",
  # probability_ever_recur = "Probability of eventual recurrence",
  proportion_micro_pos = "Proportion NAAT-positive at diagnosis",
  symptom_duration_sdlog_reported = "SD log of reported symptom duration",
  reported_fraction_of_true_symptom_duration = "Reported fraction of true symptom duration",
  programmatic_symptom_duration_factor = "Increase in symptom duration under programmatic conditions",
  proportion_micropos_sputum_first = "Proportion of micro+ recurrences with asymptomatic NAAT+ period",
  proportion_micropos_subclinical_at_eot = "Proportion of micro+ recurrences NAAT+ at treatment completion",
  duration_ratio_subclinical_symptomatic = "Asymptomatic:symptomatic NAAT+ time ratio",
  duration_subclinical_cv = "Coefficient of variation, asymptomatic TB duration",
  median_dx_by_540_sim= "Median days to diagnosis (if diagnosed by 18mo)",
  p25_dx_by_540_sim = "25th quantile, days to diagnosis (if diagnosed by 18mo)",
  p75_dx_by_540_sim = "75th quantile, days to diagnosis (if diagnosed by 18mo)",
  prop_dx_le_90_among_dx540_sim = "Proportion diagnosed by day 90 (if diagnosed by 18mo)",
  prop_dx_le_360_among_dx540_sim = "Proportion diagnosed by day 360 (if diagnosed by 18mo)",
  prop_symptom_onset_le_30_among_dx540_sim = "Proportion of symptom onset times occurring before day 30"
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
  "cum_months_subclinical_24mo",
  "proportion_sputum_first",
  "proportion_subclinical_at_eot",
  "median_time_to_first_event"
)

accepted_status <- tibble(
  sim_id = seq_len(nrow(cohort_params)),
  accepted = keep_index
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
  group_by(type, quantity, accepted) %>%
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
  "outputs/acceptance_parameter_outcome_summary.csv"
)

nice_order <- nice_names[names(nice_names) %in% sampled_parameters] %>%
  unname() %>%
  stringr::str_wrap(width = 22)

param_compare_data <- param_compare_data %>%
  mutate(
    nice_quantity = recode(quantity, !!!nice_names),
    nice_quantity = stringr::str_wrap(nice_quantity, width = 22),
    nice_quantity = factor(nice_quantity, levels = nice_order),
    accepted_label = ifelse(accepted, "Retained", "Rejected")
  )
#maintain ordering of nice_names in plot
param_compare_plot <- param_compare_data %>%
  
  filter(type == "Sampled parameter") %>%
  ggplot(aes(x = accepted_label, y = value, color = accepted_label)) +
  geom_violin(trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(. ~ nice_quantity, scales = "free_y") +
  xlab("Retained after calibration") +
  ylab("Value") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 8)
  )

ggsave(
  "outputs/acceptance_parameter_outcome_distributions.pdf",
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
  sim_id = seq_len(nrow(cohort_params_main)),
  
  guidelines_sx_reduction =
    results_main$guidelines$symptomatic_months_averted /
    results_main$guidelines$symptomatic_months_soc * 100,
  
  guidelines_inf_reduction =
    results_main$guidelines$infectious_months_averted /
    results_main$guidelines$infectious_months_soc * 100,
  
  earlier_three_vs_guidelines_sx =
    results_main$earlier_three$symptomatic_months_averted /
    results_main$guidelines$symptomatic_months_averted,
  
  earlier_three_vs_guidelines_inf =
    results_main$earlier_three$infectious_months_averted /
    results_main$guidelines$infectious_months_averted,
  
  earlier_three_sputum_vs_earlier_three_sx =
    results_main$earlier_three_sputum$symptomatic_months_averted /
    results_main$earlier_three$symptomatic_months_averted,
  
  earlier_three_sputum_vs_earlier_three_inf =
    results_main$earlier_three_sputum$infectious_months_averted /
    results_main$earlier_three$infectious_months_averted
)

outcome_nice_names <- c(
  guidelines_sx_reduction = "Guidelines vs SOC:\n% symptomatic time averted",
  guidelines_inf_reduction = "Guidelines vs SOC:\n% infectious time averted",
  earlier_three_vs_guidelines_sx = "3/6/9m symptoms vs guidelines:\nrelative symptomatic time averted",
  earlier_three_vs_guidelines_inf = "3/6/9m symptoms vs guidelines:\nrelative infectious time averted",
  earlier_three_sputum_vs_earlier_three_sx = "3/6/9m + 3m micro vs 3/6/9m symptoms:\nrelative symptomatic time averted",
  earlier_three_sputum_vs_earlier_three_inf = "3/6/9m + 3m micro vs 3/6/9m symptoms:\nrelative infectious time averted"
)

parameters_to_plot <- prcc_results %>%
  pull(variable) %>%
  unique()
parameters_to_plot2 <- prcc2_results %>%
  pull(variable) %>%
  unique()
parameters_to_plot_combined <- union(parameters_to_plot, parameters_to_plot2)


plot_data <- cohort_params_main %>%
  mutate(sim_id = row_number()) %>%
  select(sim_id, all_of(parameters_to_plot_combined)) %>%
  pivot_longer(
    cols = -sim_id,
    names_to = "parameter",
    values_to = "parameter_value"
  ) %>%
  group_by(parameter) %>%
  mutate(
    quintile = ntile(parameter_value, 10),
    quintile_group = case_when(
      quintile == 1 ~ "Lowest decile",
      quintile == 5 ~ "Highest decile",
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
      levels = c("Highest decile", "Lowest decile")
    )
  )

density_data <- plot_data %>%
  group_by(parameter_label, outcome_label, quintile_group) %>%
  group_modify(~{
    if (nrow(.x) < 10 || length(unique(.x$value)) < 2) return(tibble())
    d <- density(.x$value, na.rm = TRUE)
    tibble(y = d$x, dens = d$y)
  }) %>%
  ungroup() %>%
  group_by(parameter_label, outcome_label, quintile_group) %>%
  mutate(
    dens_scaled = dens / max(dens, na.rm = TRUE) * 0.42,
    x = case_when(
      quintile_group == "Highest decile" ~ 1 - dens_scaled,
      quintile_group == "Lowest decile" ~ 1 + dens_scaled
    )
  ) %>%
  ungroup()

box_data <- plot_data %>%
  mutate(
    x_box = case_when(
      quintile_group == "Highest decile" ~ 0.92,
      quintile_group == "Lowest decile" ~ 1.08
    )
  )

decile_sensitivity_plot <- ggplot() +
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
      "Highest decile" = "steelblue4",
      "Lowest decile" = "skyblue3"
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
  decile_sensitivity_plot,
  width = 18,
  height = 10
)

# compare deciles of a given parameter.
# parameter: 1/symptom_underestimation_factor
# outcome: symptomatic_months_averted/symptomatic_months_soc
# scenario: earlier_three_sputum
# sim numbers considered:
results_main$earlier_three_sputum %>%
  reframe(symptomatic_reduction = symptomatic_months_averted/symptomatic_months_soc*100) %>%
  bind_cols(cohort_params_main) %>%
  group_by(decile = ntile(1/reported_fraction_of_true_symptom_duration, 10)) %>%
  summarize(median_reduction = median(symptomatic_reduction, na.rm=T),
            lci_reduction = quantile(symptomatic_reduction, 0.025, na.rm=T),
            uci_reduction = quantile(symptomatic_reduction, 0.975, na.rm=T))



#### Compare outcomes for retained vs rejected simulations ####

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
  earlier_three = "Symptoms at\n3, 6, & 9m",
  earlier_three_sputum = "Symptoms at 3, 6, & 9m\n+ 3m micro",
  four_visits_36912 = "Symptoms at\n3, 6, 9, & 12m",
  frequent = "Symptoms at 3, 6, 9, 12, 15, & 18m"
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

results_rejected <- lapply(
  results_unfiltered,
  function(x) x[!keep_index, ]
)

filter_compare_data <- bind_rows(
  prep_results_for_filter_comparison(results_rejected, "Rejected simulations"),
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
    intervention = factor(intervention, levels = key_interventions),
    outcome = factor(outcome, levels = names(output_nice_names)),
    intervention_label = factor(
      recode(as.character(intervention), !!!intervention_nice_names),
      levels = unname(intervention_nice_names[key_interventions])
    ),
    outcome_label = factor(
      stringr::str_wrap(recode(as.character(outcome), !!!output_nice_names), width = 22),
      levels = stringr::str_wrap(unname(output_nice_names), width = 22)
    ),
    filter_label = factor(
      filter_label,
      levels = c("Rejected simulations", "Retained simulations")
    )
  )

filter_compare_plot <- filter_compare_data %>%
  ggplot(aes(x = filter_label, y = value, fill = filter_label)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.12, outlier.shape = NA) +
  facet_grid(
    outcome_label ~ intervention_label,
    scales = "free_y",
    switch = "both"
  ) +
  xlab("Intervention") +
  ylab("Outcome") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text.x = element_text(size = 8),
    strip.text.y = element_text(size = 8),
    strip.placement = "outside",
    strip.background = element_blank(),
    axis.title.x = element_text(face = "bold", margin = margin(t = 12)),
    axis.title.y = element_text(face = "bold", margin = margin(r = 12))
  )

ggsave(
  "outputs/filter_vs_no_filter_key_outputs_violin.pdf",
  filter_compare_plot,
  width = 13,
  height = 10
)

