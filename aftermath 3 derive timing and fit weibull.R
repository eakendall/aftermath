#### 03_fit_recurrence_timing_fast.R ####

library(tidyverse)
library(conflicted)
library(future)
library(future.apply)
library(progressr)

handlers(global = TRUE)
handlers("txtprogressbar")

plan(multisession, workers = 2)

conflicts_prefer(dplyr::select, dplyr::filter, dplyr::summarize)

set.seed(12345)

#### Run controls ####

test_mode <- FALSE
n_test_draws <- 10

n_sim_fit <- 3000
n_sim_final <- 20000

diagnosis_horizon_days <- 540

#### Load inputs ####

fixed_empirical_inputs <- readRDS("outputs/fixed_empirical_inputs.rds")
data <- fixed_empirical_inputs$data

cohort_params <- readRDS(
  paste0(
    "outputs/probabilistic_parameter_draws_pre_weibull_",
    date,
    ".rds"
  )
)

if (test_mode) {
  cohort_params <- cohort_params %>% slice(1:n_test_draws)
  message("TEST MODE: running ", nrow(cohort_params), " parameter draws.")
} else {
  message("FULL MODE: running ", nrow(cohort_params), " parameter draws.")
}

#### Empirical targets ####

empirical_diagnosis_hazard_540 <-
  fixed_empirical_inputs$diagnosis_cumhaz_540$cumhaz

observed_diagnosis_times <- data %>%
  filter(
    recurrence == 1,
    txcompl_endreason_days > 0,
    txcompl_endreason_days <= diagnosis_horizon_days
  ) %>%
  pull(txcompl_endreason_days)

observed_targets <- tibble(
  target_median_dx = median(observed_diagnosis_times, na.rm = TRUE),
  target_p25_dx = quantile(observed_diagnosis_times, 0.25, na.rm = TRUE),
  target_p75_dx = quantile(observed_diagnosis_times, 0.75, na.rm = TRUE),
  target_prop_dx_le_90 = mean(observed_diagnosis_times <= 90, na.rm = TRUE),
  target_prop_dx_le_360 = mean(observed_diagnosis_times <= 360, na.rm = TRUE),
  n_observed_dx_by_540 = length(observed_diagnosis_times)
)

print(observed_targets)

#### Helper functions ####

simulate_recurrence_course <- function(
    n,
    recurrence_shape,
    recurrence_scale,
    symptom_duration_meanlog_reported,
    symptom_duration_sdlog_reported,
    reported_fraction_of_true_symptom_duration,
    programmatic_symptom_duration_factor,
    proportion_micropos_sputum_first,
    duration_ratio_subclinical_symptomatic,
    duration_subclinical_cv
) {
  
  first_event_time <- rweibull(
    n,
    shape = recurrence_shape,
    scale = recurrence_scale
  )
  
  reported_symptom_duration <- rlnorm(
    n,
    meanlog = symptom_duration_meanlog_reported,
    sdlog = symptom_duration_sdlog_reported
  )
  
  true_symptom_duration <-
    reported_symptom_duration /
    reported_fraction_of_true_symptom_duration *
    programmatic_symptom_duration_factor
  
  mean_true_symptom_duration <- mean(true_symptom_duration, na.rm = TRUE)
  
  mean_subclinical_duration <-
    mean_true_symptom_duration *
    duration_ratio_subclinical_symptomatic
  
  sputum_first <- rbinom(
    n = n,
    size = 1,
    prob = proportion_micropos_sputum_first
  )
  
  subclinical_duration <- ifelse(
    sputum_first == 1,
    mean_subclinical_duration *
      rgamma(
        n = n,
        shape = 1 / duration_subclinical_cv^2,
        scale = duration_subclinical_cv^2
      ),
    0
  )
  
  sputum_onset <- ifelse(
    sputum_first == 1,
    first_event_time,
    NA_real_
  )
  
  symptom_onset <- first_event_time + subclinical_duration
  
  diagnosis_time <- symptom_onset + true_symptom_duration
  
  tibble(
    first_event_time = first_event_time,
    sputum_first = sputum_first,
    sputum_onset = sputum_onset,
    subclinical_duration = subclinical_duration,
    symptom_onset = symptom_onset,
    true_symptom_duration = true_symptom_duration,
    diagnosis_time = diagnosis_time
  )
}

summarize_dx <- function(sim_data, horizon = 540) {
  diagnosis_time <- sim_data$diagnosis_time
  dx540 <- diagnosis_time[diagnosis_time <= horizon]
  
  if (length(dx540) < 50) {
    return(tibble(
      median_dx_by_540 = NA_real_,
      p25_dx_by_540 = NA_real_,
      p75_dx_by_540 = NA_real_,
      prop_dx_le_90_among_dx540 = NA_real_,
      prop_dx_le_360_among_dx540 = NA_real_,
      probability_dx540_given_recur = mean(diagnosis_time <= horizon),
      mean_symptom_duration_by_540 = NA_real_,
      median_symptom_duration_by_540 = NA_real_,
      mean_subclinical_duration_by_540 = NA_real_,
      prop_first_event_le_7_among_dx540 = NA_real_,
      prop_first_event_le_30_among_dx540 = NA_real_,
      prop_symptom_onset_le_7_among_dx540 = NA_real_,
      prop_symptom_onset_le_30_among_dx540 = NA_real_,
      prop_sputum_first_by_540 = NA_real_
    ))
  }
  
  dx540_index <- diagnosis_time <= horizon
  
  tibble(
    median_dx_by_540 = median(dx540),
    p25_dx_by_540 = quantile(dx540, 0.25),
    p75_dx_by_540 = quantile(dx540, 0.75),
    prop_dx_le_90_among_dx540 = mean(dx540 <= 90),
    prop_dx_le_360_among_dx540 = mean(dx540 <= 360),
    probability_dx540_given_recur = mean(diagnosis_time <= horizon),
    
    mean_symptom_duration_by_540 =
      mean(sim_data$true_symptom_duration[dx540_index], na.rm = TRUE),
    median_symptom_duration_by_540 =
      median(sim_data$true_symptom_duration[dx540_index], na.rm = TRUE),
    
    mean_subclinical_duration_by_540 =
      mean(sim_data$subclinical_duration[dx540_index], na.rm = TRUE),
    
    prop_first_event_le_7_among_dx540 =
      mean(sim_data$first_event_time[dx540_index] <= 7, na.rm = TRUE),
    prop_first_event_le_30_among_dx540 =
      mean(sim_data$first_event_time[dx540_index] <= 30, na.rm = TRUE),
    
    prop_symptom_onset_le_7_among_dx540 =
      mean(sim_data$symptom_onset[dx540_index] <= 7, na.rm = TRUE),
    prop_symptom_onset_le_30_among_dx540 =
      mean(sim_data$symptom_onset[dx540_index] <= 30, na.rm = TRUE),
    
    prop_sputum_first_by_540 =
      mean(sim_data$sputum_first[dx540_index] == 1, na.rm = TRUE)
  )
}

weibull_mean <- function(shape, scale) {
  scale * gamma(1 + 1 / shape)
}

weibull_cv <- function(shape, scale) {
  m <- weibull_mean(shape, scale)
  s <- scale * sqrt(
    gamma(1 + 2 / shape) -
      gamma(1 + 1 / shape)^2
  )
  s / m
}

#### Fit one draw ####

fit_one_draw <- function(draw_row, i, n_total) {
  
  message(
    "[", i, "/", n_total, "] draw=", draw_row$draw,
    " | fraction=", round(draw_row$reported_fraction_of_true_symptom_duration, 3),
    " | sdlog=", round(draw_row$symptom_duration_sdlog_reported, 3),
    " | prog factor=", round(draw_row$programmatic_symptom_duration_factor, 2),
    " | sputum first=", round(draw_row$proportion_micropos_sputum_first, 2)
  )
  
  incidence_18mo <-
    empirical_diagnosis_hazard_540 *
    draw_row$incidence_18mo_multiplier
  
  objective <- function(log_params) {
    
    shape <- exp(log_params[1])
    scale <- exp(log_params[2])
    
    sim_data <- simulate_recurrence_course(
      n = n_sim_fit,
      recurrence_shape = shape,
      recurrence_scale = scale,
      symptom_duration_meanlog_reported =
        draw_row$symptom_duration_meanlog_reported,
      symptom_duration_sdlog_reported =
        draw_row$symptom_duration_sdlog_reported,
      reported_fraction_of_true_symptom_duration =
        draw_row$reported_fraction_of_true_symptom_duration,
      programmatic_symptom_duration_factor =
        draw_row$programmatic_symptom_duration_factor,
      proportion_micropos_sputum_first =
        draw_row$proportion_micropos_sputum_first,
      duration_ratio_subclinical_symptomatic =
        draw_row$duration_ratio_subclinical_symptomatic,
      duration_subclinical_cv =
        draw_row$duration_subclinical_cv
    )
    
    sim_targets <- summarize_dx(
      sim_data,
      horizon = diagnosis_horizon_days
    )
    
    if (
      any(!is.finite(unlist(sim_targets))) ||
      sim_targets$probability_dx540_given_recur <= 0
    ) {
      return(1e20)
    }
    
    sim_iqr <- sim_targets$p75_dx_by_540 -
      sim_targets$p25_dx_by_540
    
    obs_iqr <- observed_targets$target_p75_dx -
      observed_targets$target_p25_dx
    
    error_median <-
      (sim_targets$median_dx_by_540 -
         observed_targets$target_median_dx)^2
    
    error_iqr <-
      0.5 * (sim_iqr - obs_iqr)^2
    
    error_p90 <-
      5000 *
      (sim_targets$prop_dx_le_90_among_dx540 -
         observed_targets$target_prop_dx_le_90)^2
    
    error_p360 <-
      5000 *
      (sim_targets$prop_dx_le_360_among_dx540 -
         observed_targets$target_prop_dx_le_360)^2
    
    error_median + error_iqr + error_p90 + error_p360
  }
  
  fit <- optim(
    par = log(c(1.2, 200)),
    fn = objective,
    method = "Nelder-Mead",
    control = list(maxit = 300, reltol = 1e-6)
  )
  
  shape <- exp(fit$par[1])
  scale <- exp(fit$par[2])
  
  sim_data_final <- simulate_recurrence_course(
    n = n_sim_final,
    recurrence_shape = shape,
    recurrence_scale = scale,
    symptom_duration_meanlog_reported =
      draw_row$symptom_duration_meanlog_reported,
    symptom_duration_sdlog_reported =
      draw_row$symptom_duration_sdlog_reported,
    reported_fraction_of_true_symptom_duration =
      draw_row$reported_fraction_of_true_symptom_duration,
    programmatic_symptom_duration_factor =
      draw_row$programmatic_symptom_duration_factor,
    proportion_micropos_sputum_first =
      draw_row$proportion_micropos_sputum_first,
    duration_ratio_subclinical_symptomatic =
      draw_row$duration_ratio_subclinical_symptomatic,
    duration_subclinical_cv =
      draw_row$duration_subclinical_cv
  )
  
  final_targets <- summarize_dx(
    sim_data_final,
    horizon = diagnosis_horizon_days
  )
  
  probability_dx540_given_recur <-
    final_targets$probability_dx540_given_recur
  
  probability_ever_recur <-
    incidence_18mo / probability_dx540_given_recur
  
  message(
    "   fitted shape=", round(shape, 3),
    " scale=", round(scale, 1),
    " P(dx<=540|recur)=", round(probability_dx540_given_recur, 3),
    " P(ever recur)=", round(probability_ever_recur, 3),
    " mean sx dur(dx540)=", round(final_targets$mean_symptom_duration_by_540, 1),
    " mean subclin dur(dx540)=", round(final_targets$mean_subclinical_duration_by_540, 1),
    " obj=", round(fit$value, 1)
  )
  
  tibble(
    incidence_18mo = incidence_18mo,
    
    recurrence_shape = shape,
    recurrence_scale = scale,
    recurrence_time_mean = weibull_mean(shape, scale),
    recurrence_time_cv = weibull_cv(shape, scale),
    
    probability_dx540_given_recur = probability_dx540_given_recur,
    probability_ever_recur = probability_ever_recur,
    
    median_dx_by_540_sim =
      final_targets$median_dx_by_540,
    p25_dx_by_540_sim =
      final_targets$p25_dx_by_540,
    p75_dx_by_540_sim =
      final_targets$p75_dx_by_540,
    prop_dx_le_90_among_dx540_sim =
      final_targets$prop_dx_le_90_among_dx540,
    prop_dx_le_360_among_dx540_sim =
      final_targets$prop_dx_le_360_among_dx540,
    
    mean_symptom_duration_by_540_sim =
      final_targets$mean_symptom_duration_by_540,
    median_symptom_duration_by_540_sim =
      final_targets$median_symptom_duration_by_540,
    mean_subclinical_duration_by_540_sim =
      final_targets$mean_subclinical_duration_by_540,
    
    prop_first_event_le_7_among_dx540_sim =
      final_targets$prop_first_event_le_7_among_dx540,
    prop_first_event_le_30_among_dx540_sim =
      final_targets$prop_first_event_le_30_among_dx540,
    prop_symptom_onset_le_7_among_dx540_sim =
      final_targets$prop_symptom_onset_le_7_among_dx540,
    prop_symptom_onset_le_30_among_dx540_sim =
      final_targets$prop_symptom_onset_le_30_among_dx540,
    prop_sputum_first_by_540_sim =
      final_targets$prop_sputum_first_by_540,
    
    target_median_dx =
      observed_targets$target_median_dx,
    target_p25_dx =
      observed_targets$target_p25_dx,
    target_p75_dx =
      observed_targets$target_p75_dx,
    target_prop_dx_le_90 =
      observed_targets$target_prop_dx_le_90,
    target_prop_dx_le_360 =
      observed_targets$target_prop_dx_le_360,
    
    objective_value = fit$value,
    optim_convergence_code = fit$convergence,
    optim_message = ifelse(is.null(fit$message), NA_character_, fit$message),
    converged = fit$convergence == 0,
    
    valid_calibration =
      is.finite(probability_ever_recur) &&
      probability_ever_recur > 0 &&
      probability_ever_recur <= 1
  )
}

#### Timing test: one draw ####

test_i <- 1
test_draw_row <- cohort_params %>% slice(test_i)

system.time({
  test_fit <- fit_one_draw(
    test_draw_row,
    i = test_i,
    n_total = nrow(cohort_params)
  )
})

test_fit

#### Run all draws ####

n_total <- nrow(cohort_params)

with_progress({
  p <- progressor(along = seq_len(n_total))
  
  fit_results <- future_lapply(
    seq_len(n_total),
    function(i) {
      p(sprintf("fitting draw %s/%s", i, n_total))
      
      draw_row <- cohort_params %>% slice(i)
      
      bind_cols(
        draw_row,
        fit_one_draw(
          draw_row,
          i = i,
          n_total = n_total
        )
      )
    },
    future.seed = TRUE
  )
})


cohort_params_final <- bind_rows(fit_results)

cohort_params_final <- cohort_params_final %>%
  mutate(
    valid_weibull_fit =
      valid_calibration &
      is.finite(objective_value) &
      probability_ever_recur > 0 &
      probability_ever_recur <= 1 &
      objective_value <= quantile(objective_value, 0.99, na.rm = TRUE)
  )

rm(fit_results)
invisible(gc(full = TRUE))

#### Save outputs ####

write_csv(
  cohort_params_final,
  paste0("outputs/cohort_params_final_", date, ".csv")
)

saveRDS(
  cohort_params_final,
  paste0("outputs/cohort_params_final_", date, ".rds")
)

write_csv(
  observed_targets,
  paste0("outputs/observed_diagnosis_timing_targets_", date, ".csv")
)

message(
  "Done. Saved outputs/cohort_params_final_",
  date,
  ".rds and .csv."
)
