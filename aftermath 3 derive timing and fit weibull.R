#### 03_fit_recurrence_timing_fast.R ####

library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::select, dplyr::filter, dplyr::summarize)

set.seed(12345)

#### Run controls ####

test_mode <- FALSE
n_test_draws <- 10

n_sim_fit <- 10000
n_sim_final <- 50000

diagnosis_horizon_days <- 540

#### Load inputs ####

fixed_empirical_inputs <- readRDS("outputs/fixed_empirical_inputs.rds")
data <- fixed_empirical_inputs$data

cohort_params <- readRDS("outputs/probabilistic_parameter_draws_pre_weibull.rds")

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
    reported_fraction_of_true_symptom_duration
) {
  onset <- rweibull(
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
    reported_symptom_duration / reported_fraction_of_true_symptom_duration
  
  diagnosis_time <- onset + true_symptom_duration
  
  diagnosis_time
}

summarize_dx <- function(diagnosis_time, horizon = 540) {
  dx540 <- diagnosis_time[diagnosis_time <= horizon]
  
  if (length(dx540) < 50) {
    return(tibble(
      median_dx_by_540 = NA_real_,
      p25_dx_by_540 = NA_real_,
      p75_dx_by_540 = NA_real_,
      prop_dx_le_90_among_dx540 = NA_real_,
      prop_dx_le_360_among_dx540 = NA_real_,
      probability_dx540_given_recur = mean(diagnosis_time <= horizon)
    ))
  }
  
  tibble(
    median_dx_by_540 = median(dx540),
    p25_dx_by_540 = quantile(dx540, 0.25),
    p75_dx_by_540 = quantile(dx540, 0.75),
    prop_dx_le_90_among_dx540 = mean(dx540 <= 90),
    prop_dx_le_360_among_dx540 = mean(dx540 <= 360),
    probability_dx540_given_recur = mean(diagnosis_time <= horizon)
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
    " | sdlog=", round(draw_row$symptom_duration_sdlog_reported, 3)
  )
  
  incidence_18mo <-
    empirical_diagnosis_hazard_540 *
    draw_row$incidence_18mo_multiplier
  
  objective <- function(log_params) {
    
    shape <- exp(log_params[1])
    scale <- exp(log_params[2])
    
    diagnosis_time <- simulate_recurrence_course(
      n = n_sim_fit,
      recurrence_shape = shape,
      recurrence_scale = scale,
      symptom_duration_meanlog_reported =
        draw_row$symptom_duration_meanlog_reported,
      symptom_duration_sdlog_reported =
        draw_row$symptom_duration_sdlog_reported,
      reported_fraction_of_true_symptom_duration =
        draw_row$reported_fraction_of_true_symptom_duration
    )
    
    sim_targets <- summarize_dx(
      diagnosis_time,
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
    control = list(maxit = 80)
  )
  
  shape <- exp(fit$par[1])
  scale <- exp(fit$par[2])
  
  diagnosis_time_final <- simulate_recurrence_course(
    n = n_sim_final,
    recurrence_shape = shape,
    recurrence_scale = scale,
    symptom_duration_meanlog_reported =
      draw_row$symptom_duration_meanlog_reported,
    symptom_duration_sdlog_reported =
      draw_row$symptom_duration_sdlog_reported,
    reported_fraction_of_true_symptom_duration =
      draw_row$reported_fraction_of_true_symptom_duration
  )
  
  final_targets <- summarize_dx(
    diagnosis_time_final,
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
    converged = fit$convergence == 0,
    
    valid_calibration =
      is.finite(probability_ever_recur) &&
      probability_ever_recur > 0 &&
      probability_ever_recur <= 1
  )
}

#### Run all draws ####

n_total <- nrow(cohort_params)

cohort_params_final <- map_dfr(
  seq_len(n_total),
  function(i) {
    draw_row <- cohort_params %>% slice(i)
    
    bind_cols(
      draw_row,
      fit_one_draw(
        draw_row,
        i = i,
        n_total = n_total
      )
    )
  }
)

#### Save outputs ####

write_csv(
  cohort_params_final,
  "outputs/cohort_params_final.csv"
)

saveRDS(
  cohort_params_final,
  "outputs/cohort_params_final.rds"
)

write_csv(
  observed_targets,
  "outputs/observed_diagnosis_timing_targets.csv"
)

message("Done. Saved outputs/cohort_params_final.rds and .csv.")