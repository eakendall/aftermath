library(tidyverse)
library(binom)
library(clipr)
library(sensitivity)
library(conflicted)
library(future)
library(future.apply)
library(progressr)

handlers(global = TRUE)
handlers("txtprogressbar")

conflicts_prefer(dplyr::select, dplyr::filter, dplyr::summarize)

#### Settings ####

run_cohort_features <- TRUE
run_interventions <- TRUE
apply_subclinical_filter <- FALSE

# Start conservative for memory; increase if stable
plan(multisession, workers = 3)

script_start_time <- Sys.time()

cohort_params <- readRDS(
  paste0("outputs/cohort_params_final_", date, ".rds")
)

N_samples <- nrow(cohort_params)
N_cohort <- cohort_params$N[1]

#### Screening designs ####

intervention_names <- c(
  "guidelines",
  "3m_sx", "6m_sx", "12m_sx", "18m_sx",
  "3m_swab", "6m_swab", "12m_swab", "18m_swab",
  "earlier_two", "earlier_three", "earlier_three_sputum",
  "frequent", "frequent_sputum", "risk_targeted",
  "four_visits_36912", "five_visits_3691215",
  "four_visits_sputum", "five_visits_sputum"
)

screening_design_guidelines <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(6, 12, 18),
  target_coverage = c(1, 1, 1),
  screening_method = c("symptoms", "symptoms", "symptoms"),
  screening_location = c("home", "phone", "phone")
) %>% arrange(timing_months)

screening_design_3m_sx <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = 3,
  target_coverage = 1,
  screening_method = "symptoms",
  screening_location = "home"
) %>% arrange(timing_months)

screening_design_6m_sx <- screening_design_3m_sx; screening_design_6m_sx$timing_months <- 6
screening_design_12m_sx <- screening_design_3m_sx; screening_design_12m_sx$timing_months <- 12
screening_design_18m_sx <- screening_design_3m_sx; screening_design_18m_sx$timing_months <- 18

screening_design_3m_swab <- screening_design_3m_sx; screening_design_3m_swab$screening_method <- "both"
screening_design_6m_swab <- screening_design_6m_sx; screening_design_6m_swab$screening_method <- "both"
screening_design_12m_swab <- screening_design_12m_sx; screening_design_12m_swab$screening_method <- "both"
screening_design_18m_swab <- screening_design_18m_sx; screening_design_18m_swab$screening_method <- "both"

screening_design_earlier_two <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6),
  target_coverage = c(1, 1),
  screening_method = c("symptoms", "symptoms"),
  screening_location = c("home", "phone")
) %>% arrange(timing_months)

screening_design_earlier_three <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6, 9),
  target_coverage = c(1, 1, 1),
  screening_method = c("symptoms", "symptoms", "symptoms"),
  screening_location = c("home", "phone", "phone")
) %>% arrange(timing_months)

screening_design_earlier_three_sputum <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6, 9),
  target_coverage = c(1, 1, 1),
  screening_method = c("both", "symptoms", "symptoms"),
  screening_location = c("home", "phone", "phone")
) %>% arrange(timing_months)

screening_design_frequent <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6, 9, 12, 15, 18),
  target_coverage = rep(1, 6),
  screening_method = rep("symptoms", 6),
  screening_location = c("home", rep("phone", 5))
) %>% arrange(timing_months)

screening_design_frequent_sputum <- screening_design_frequent
screening_design_frequent_sputum$screening_method[1] <- "both"

screening_design_risk_targeted <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 3, 6, 9, 12, 15, 18),
  target_coverage = c(0.5, 1, 1, 1, 0.5, 0.5, 0.5),
  screening_method = c("sputum", "symptoms", "symptoms", "symptoms", "symptoms", "symptoms", "symptoms"),
  screening_location = c("home", "home", "phone", "phone", "phone", "phone", "phone")
) %>% arrange(timing_months)

screening_design_four_visits_36912 <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6, 9, 12),
  target_coverage = rep(1, 4),
  screening_method = rep("symptoms", 4),
  screening_location = c("home", rep("phone", 3))
) %>% arrange(timing_months)

screening_design_five_visits_3691215 <- data.frame(
  counseling_coverage = 0,
  prevention_coverage = 0,
  timing_months = c(3, 6, 9, 12, 15),
  target_coverage = rep(1, 5),
  screening_method = rep("symptoms", 5),
  screening_location = c("home", rep("phone", 4))
) %>% arrange(timing_months)

screening_design_four_visits_sputum <- screening_design_four_visits_36912
screening_design_four_visits_sputum$screening_method[1] <- "both"

screening_design_five_visits_sputum <- screening_design_five_visits_3691215
screening_design_five_visits_sputum$screening_method[1] <- "both"

screening_designs <- mget(paste0("screening_design_", intervention_names))

#### Helper: one cohort feature run ####

run_one_cohort_feature <- function(n) {
  cohort <- create_cohort(cohort_params[n, ])
  
  subclinical_baseline_among_micropos <-
    sum(cohort$subclinical_at_eot == 1, na.rm = TRUE) /
    sum(cohort$pulmonary_with_micro == 1, na.rm = TRUE)
  
  subclinical_6mo_amongcohort <-
    sum(cohort$sputum_onset < 180 & cohort$symptom_onset >= 180,
        na.rm = TRUE) / nrow(cohort)
  
  accepted_subclinical <- check_subclinical(cohort, cohort_params[n, ])
  
  intervention_parameters <- list(
    coverage = list(
      phone = cohort_params$coverage_phone[n],
      home = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]
    ),
    sensitivity_symptoms = list(
      home = cohort_params$sensitivity_symptoms_home[n],
      phone = cohort_params$sensitivity_symptoms_home[n] *
        cohort_params$sensitivity_symptoms_phone_reduction[n]
    ),
    success_sputum = list(
      home = cohort_params$success_sputum_home[n],
      phone = cohort_params$success_sputum_home[n] *
        cohort_params$success_sputum_phone_reduction[n]
    )
  )
  
  outputs <- apply_intervention(
    cohort,
    screening_design_guidelines,
    intervention_parameters,
    cohort_params[n, ]
  )
  
  cohort_row <- data.frame(
    accepted_subclinical = accepted_subclinical,
    subclinical_baseline_among_micropos = subclinical_baseline_among_micropos,
    subclinical_6mo_amongcohort = subclinical_6mo_amongcohort,
    fails_6mo_low =
      subclinical_6mo_amongcohort <= cohort_params$subclinical_6m_amongcohort_min[n],
    fails_6mo_high =
      subclinical_6mo_amongcohort >= cohort_params$subclinical_6m_amongcohort_max[n],
    
    month_with_highest_incidence =
      cohort %>% filter(TB == 1) %>%
      mutate(month = floor(symptom_onset / 30)) %>%
      count(month, sort = TRUE) %>%
      slice(1) %>%
      pull(month),
    
    month_with_highest_prevalence_sx =
      (1:20)[which.max(sapply(30 * (1:20), function(x)
        sum(cohort$symptom_onset < x & cohort$diagnosis_routine > x, na.rm = TRUE)))],
    
    month_with_highest_prevalence_inf =
      (1:20)[which.max(sapply(30 * (1:20), function(x)
        sum(cohort$sputum_onset < x & cohort$diagnosis_routine > x, na.rm = TRUE)))],
    
    month_with_highest_prevalence_total =
      (1:20)[which.max(sapply(30 * (1:20), function(x)
        sum((cohort$sputum_onset < x | cohort$symptom_onset < x) &
              cohort$diagnosis_routine > x, na.rm = TRUE)))],
    
    remaining_duration_3mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 90 &
          (cohort$sputum_onset < 90 | cohort$symptom_onset < 90)
      ] - 90, na.rm = TRUE),
    
    remaining_duration_6mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 180 &
          (cohort$sputum_onset < 180 | cohort$symptom_onset < 180)
      ] - 180, na.rm = TRUE),
    
    remaining_duration_9mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 270 &
          (cohort$sputum_onset < 270 | cohort$symptom_onset < 270)
      ] - 270, na.rm = TRUE),
    
    remaining_duration_12mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 360 &
          (cohort$sputum_onset < 360 | cohort$symptom_onset < 360)
      ] - 360, na.rm = TRUE),
    
    remaining_duration_15mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 450 &
          (cohort$sputum_onset < 450 | cohort$symptom_onset < 450)
      ] - 450, na.rm = TRUE),
    
    remaining_duration_18mo =
      mean(cohort$diagnosis_routine[
        cohort$diagnosis_routine > 540 &
          (cohort$sputum_onset < 540 | cohort$symptom_onset < 540)
      ] - 540, na.rm = TRUE),
    
    cum_notifs_6mo =
      sum(cohort$TB == 1 & cohort$diagnosis_routine <= 180, na.rm = TRUE) / N_cohort,
    cum_notifs_18mo =
      sum(cohort$TB == 1 & cohort$diagnosis_routine <= 540, na.rm = TRUE) / N_cohort,
    cum_notifs_24mo =
      sum(cohort$TB == 1 & cohort$diagnosis_routine <= 720, na.rm = TRUE) / N_cohort,
    cum_notifs_overall =
      sum(cohort$TB == 1, na.rm = TRUE) / N_cohort,
    
    median_time_to_onset = median(cohort$symptom_onset, na.rm = TRUE),
    median_time_to_first_event = median(cohort$first_event_time, na.rm = TRUE),
    median_time_to_diagnosis = median(cohort$diagnosis_routine, na.rm = TRUE),
    
    proportion_sputum_first =
      sum(cohort$sputum_first == 1, na.rm = TRUE) /
      sum(cohort$pulmonary_with_micro == 1, na.rm = TRUE),
    
    proportion_subclinical_at_eot =
      sum(cohort$subclinical_at_eot == 1, na.rm = TRUE) /
      sum(cohort$pulmonary_with_micro == 1, na.rm = TRUE),
    
    mean_symptom_duration = mean(cohort$diagnosis_routine - cohort$symptom_onset, na.rm = TRUE),
    sd_symptom_duration = sd(cohort$diagnosis_routine - cohort$symptom_onset, na.rm = TRUE),
    proportion_micro_pos =
      sum(cohort$pulmonary_with_micro == 1, na.rm = TRUE) /
      sum(cohort$TB == 1, na.rm = TRUE),
    proportion_with_subclinical =
      sum(cohort$sputum_onset < cohort$symptom_onset, na.rm = TRUE) /
      sum(cohort$pulmonary_with_micro == 1, na.rm = TRUE),
    
    mean_subclinical_duration =
      mean(cohort$symptom_onset[cohort$sputum_onset < cohort$symptom_onset] -
             cohort$sputum_onset[cohort$sputum_onset < cohort$symptom_onset],
           na.rm = TRUE),
    
    sd_subclinical_duration =
      sd(cohort$symptom_onset[cohort$sputum_onset < cohort$symptom_onset] -
           cohort$sputum_onset[cohort$sputum_onset < cohort$symptom_onset],
         na.rm = TRUE),
    
    prev_sx_6mo =
      sum(cohort$symptom_onset < 180 & cohort$diagnosis_routine >= 180,
          na.rm = TRUE) / N_cohort,
    prev_inf_6mo =
      sum(cohort$sputum_onset < 180 & cohort$diagnosis_routine >= 180,
          na.rm = TRUE) / N_cohort,
    prev_subclinical_6mo =
      sum(cohort$sputum_onset < 180 & cohort$symptom_onset >= 180,
          na.rm = TRUE) / N_cohort,
    
    prev_sx_12mo =
      sum(cohort$symptom_onset < 360 & cohort$diagnosis_routine >= 360,
          na.rm = TRUE) / N_cohort,
    prev_inf_12mo =
      sum(cohort$sputum_onset < 360 & cohort$diagnosis_routine >= 360,
          na.rm = TRUE) / N_cohort,
    prev_subclinical_12mo =
      sum(cohort$sputum_onset < 360 & cohort$symptom_onset >= 360,
          na.rm = TRUE) / N_cohort,
    
    cum_months_sx =
      sum(cohort$diagnosis_routine - cohort$symptom_onset, na.rm = TRUE) / 30,
    cum_months_inf =
      sum(cohort$diagnosis_routine - cohort$sputum_onset, na.rm = TRUE) / 30,
    cum_months_subclinical =
      sum(pmax(cohort$symptom_onset - cohort$sputum_onset, 0), na.rm = TRUE) / 30,
    
    cum_months_sx_24mo =
      sum(pmax(pmin(cohort$diagnosis_routine, 720) -
                 pmin(cohort$symptom_onset, 720), 0),
          na.rm = TRUE) / 30,
    cum_months_inf_24mo =
      sum(pmax(pmin(cohort$diagnosis_routine, 720) -
                 pmin(cohort$sputum_onset, 720), 0),
          na.rm = TRUE) / 30,
    cum_months_subclinical_24mo =
      sum(pmax(pmin(cohort$symptom_onset, 720) -
                 pmin(cohort$sputum_onset, 720), 0),
          na.rm = TRUE) / 30,
    
    duration_diagnosed_before_6 =
      mean(cohort$diagnosis_routine[cohort$TB == 1 & cohort$diagnosis_routine <= 180] -
             cohort$symptom_onset[cohort$TB == 1 & cohort$diagnosis_routine <= 180],
           na.rm = TRUE) / 30,
    
    duration_diagnosed_before_9 =
      mean(cohort$diagnosis_routine[cohort$TB == 1 & cohort$diagnosis_routine <= 270] -
             cohort$symptom_onset[cohort$TB == 1 & cohort$diagnosis_routine <= 270],
           na.rm = TRUE) / 30,
    
    duration_diagnosed_after_6 =
      mean(cohort$diagnosis_routine[cohort$TB == 1 & cohort$diagnosis_routine > 180] -
             cohort$symptom_onset[cohort$TB == 1 & cohort$diagnosis_routine > 180],
           na.rm = TRUE) / 30,
    
    duration_diagnosed_after_9 =
      mean(cohort$diagnosis_routine[cohort$TB == 1 & cohort$diagnosis_routine > 270] -
             cohort$symptom_onset[cohort$TB == 1 & cohort$diagnosis_routine > 270],
           na.rm = TRUE) / 30
  )
  
  cascade_row <- data.frame(
    accepted_subclinical = accepted_subclinical,
    cumulative_incidence = sum(cohort$TB == 1),
    TB_beyond_6mo = sum(cohort$diagnosis_routine > 180 & cohort$TB == 1),
    TB_beyond6_before18 =
      sum(cohort$diagnosis_routine > 180 &
            (cohort$symptom_onset < 540 | cohort$sputum_onset < 540) &
            cohort$TB == 1, na.rm = TRUE),
    TB_at_visit =
      sum(cohort$TB == 1 &
            (((cohort$symptom_onset < 180 | cohort$sputum_onset < 180) &
                cohort$diagnosis_routine > 180) |
               ((cohort$symptom_onset < 360 | cohort$sputum_onset < 360) &
                  cohort$diagnosis_routine > 360) |
               ((cohort$symptom_onset < 540 | cohort$sputum_onset < 540) &
                  cohort$diagnosis_routine > 540)), na.rm = TRUE),
    symptomatic_at_visit =
      sum(cohort$TB == 1 &
            ((cohort$symptom_onset < 180 & cohort$diagnosis_routine > 180) |
               (cohort$symptom_onset < 360 & cohort$diagnosis_routine > 360) |
               (cohort$symptom_onset < 540 & cohort$diagnosis_routine > 540)),
          na.rm = TRUE),
    detected_and_linked =
      sum(outputs$TB == 1 & !is.na(outputs$detection_timing)),
    total_time_of_linked =
      sum((outputs$diagnosis_routine - outputs$symptom_onset)[
        outputs$TB == 1 & !is.na(outputs$detection_timing)
      ]),
    remaining_time_of_linked =
      sum((outputs$detection_timing - outputs$symptom_onset)[
        outputs$TB == 1 & !is.na(outputs$detection_timing)
      ]),
    duration_soc_detected_61218 =
      mean(outputs$diagnosis_routine[outputs$TB == 1 & !is.na(outputs$detection_timing)] -
             outputs$symptom_onset[outputs$TB == 1 & !is.na(outputs$detection_timing)],
           na.rm = TRUE) / 30,
    contribution_soc_detected_61218 =
      sum(outputs$diagnosis_routine[outputs$TB == 1 & !is.na(outputs$detection_timing)] -
            outputs$symptom_onset[outputs$TB == 1 & !is.na(outputs$detection_timing)],
          na.rm = TRUE) /
      sum(outputs$diagnosis_routine[outputs$TB == 1] -
            outputs$symptom_onset[outputs$TB == 1],
          na.rm = TRUE),
    duration_notdetected_61218 =
      mean(outputs$diagnosis_routine[outputs$TB == 1 & is.na(outputs$detection_timing)] -
             outputs$symptom_onset[outputs$TB == 1 & is.na(outputs$detection_timing)],
           na.rm = TRUE) / 30
  )
  
  rm(cohort, outputs, intervention_parameters)
  invisible(gc())
  
  list(cohort = cohort_row, cascade = cascade_row)
}

#### Run cohort features ####

if (run_cohort_features) {
  with_progress({
    p <- progressor(along = seq_len(N_samples))
    
    feature_list <- future_lapply(
      seq_len(N_samples),
      function(n) {
        p(sprintf("cohort features: %s/%s", n, N_samples))
        run_one_cohort_feature(n)
      },
      future.seed = TRUE
    )
  })
  
  cohort_features <- bind_rows(lapply(feature_list, `[[`, "cohort"))
  cascade_features <- bind_rows(lapply(feature_list, `[[`, "cascade"))
  
  rm(feature_list)
  invisible(gc(full = TRUE))
  
  saveRDS(cohort_features, file = paste0("outputs/cohort_features_", date, ".rds"))
  saveRDS(cascade_features, file = paste0("outputs/cascade_features_", date, ".rds"))
}

# cohort features diagnostic: 

cohort_features %>%
  summarise(
    accepted = mean(accepted_subclinical),
    median_sputum_first = median(proportion_sputum_first, na.rm = TRUE),
    median_subclinical_at_eot = median(proportion_subclinical_at_eot, na.rm = TRUE),
    median_first_event = median(median_time_to_first_event, na.rm = TRUE),
    median_symptom_onset = median(median_time_to_onset, na.rm = TRUE),
    median_6mo_subclinical = median(subclinical_6mo_amongcohort, na.rm = TRUE)
  )

#### Helper: one intervention set ####

run_one_intervention_set <- function(n) {
  cohort <- create_cohort(cohort_params[n, ])
  
  if (apply_subclinical_filter && !check_subclinical(cohort, cohort_params[n, ])) {
    empty <- data.frame(
      recurrences = NA_real_,
      detections = NA_real_,
      mean_symptom_days = NA_real_,
      symptomatic_months_soc = NA_real_,
      infectious_months_soc = NA_real_,
      symptomatic_months_averted = NA_real_,
      infectious_months_averted = NA_real_,
      cost = NA_real_
    )
    return(setNames(rep(list(empty), length(intervention_names)), intervention_names))
  }
  
  intervention_parameters <- list(
    coverage = list(
      phone = cohort_params$coverage_phone[n],
      home = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]
    ),
    sensitivity_symptoms = list(
      home = cohort_params$sensitivity_symptoms_home[n],
      phone = cohort_params$sensitivity_symptoms_home[n] *
        cohort_params$sensitivity_symptoms_phone_reduction[n]
    ),
    success_sputum = list(
      home = cohort_params$success_sputum_home[n],
      phone = cohort_params$success_sputum_home[n] *
        cohort_params$success_sputum_phone_reduction[n]
    )
  )
  
  out <- vector("list", length(intervention_names))
  names(out) <- intervention_names
  
  for (intervention in intervention_names) {
    design <- screening_designs[[paste0("screening_design_", intervention)]]
    
    output <- apply_intervention(
      cohort,
      design,
      intervention_parameters,
      cohort_params[n, ]
    )
    
    this_cost <- costs(
      design,
      cohort,
      cohort_params[n, ],
      intervention_parameters
    )
    
    tb_time <- time_with_tb(output, limit_days = 720)
    
    time_soc <- tb_time %>%
      filter(scenario == "soc") %>%
      mutate(months = value / 30) %>%
      select(outcome, months)
    
    impact <- tb_time %>%
      pivot_wider(names_from = "scenario", values_from = "value") %>%
      mutate(months_averted = (soc - screening) / 30) %>%
      select(outcome, months_averted)
    
    out[[intervention]] <- data.frame(
      recurrences = sum(output$TB),
      detections = sum(!is.na(output$detection_timing)),
      mean_symptom_days =
        (tb_time %>% filter(scenario == "screening") %>% pull(value))[1] /
        sum(output$TB),
      symptomatic_months_soc =
        time_soc$months[time_soc$outcome == "symptom"],
      infectious_months_soc =
        time_soc$months[time_soc$outcome == "sputum"],
      symptomatic_months_averted =
        impact$months_averted[impact$outcome == "symptom"],
      infectious_months_averted =
        impact$months_averted[impact$outcome == "sputum"],
      cost = this_cost
    )
    
    rm(output, this_cost, tb_time, time_soc, impact)
  }
  
  rm(cohort, intervention_parameters)
  invisible(gc())
  
  out
}

#### Run interventions ####

if (run_interventions) {
  with_progress({
    p <- progressor(along = seq_len(N_samples))
    
    results_by_sim <- future_lapply(
      seq_len(N_samples),
      function(n) {
        p(sprintf("interventions: %s/%s", n, N_samples))
        run_one_intervention_set(n)
      },
      future.seed = TRUE
    )
  })
  
  results <- vector("list", length(intervention_names))
  names(results) <- intervention_names
  
  for (intervention in intervention_names) {
    results[[intervention]] <- bind_rows(
      lapply(results_by_sim, `[[`, intervention)
    )
  }
  
  rm(results_by_sim)
  invisible(gc(full = TRUE))
  
  if (apply_subclinical_filter) {
    saveRDS(results, paste0("outputs/results_aftermath_", date, ".RDS"))
  } else {
    saveRDS(results, paste0("outputs/results_no_subclinical_filter_", date, ".RDS"))
  }
}

