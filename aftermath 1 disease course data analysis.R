library(tidyverse)
library(survival)
library(readxl)
library(conflicted)

conflicts_prefer(dplyr::select, dplyr::filter)

#### Load trial data ####

data <- read.csv("../Data June 2025 from Aye/variables_trialdata1806.csv")
dictionary <- readxl::read_excel("../Data Jan 2025 from Aye/variables_trialdata1512_codebook.xlsx")

#### Define recurrence outcome ####
# This does not depend on assumptions about true symptom duration.

data <- data %>%
  filter(!(end_reason == "death" & tb_related_death == "Unknown")) %>%
  mutate(
    recurrence = case_when(
      end_reason == "TB recurrence" ~ 1,
      end_reason == "death" & tb_related_death == "Yes" ~ 1,
      end_reason == "death" & tb_related_death == "Probable" ~ 1,
      end_reason == "death" & tb_related_death == "No" ~ 0,
      end_reason == "completion" ~ 0,
      end_reason == "LFU" ~ 0,
      TRUE ~ 0
    ),
    end_days = case_when(
      end_reason == "death" & recurrence == 1 ~ interval(txcompl_date, death_date) / days(1),
      TRUE ~ txcompl_endreason_days
    )
  )

#### Classify micropositive pulmonary versus other recurrent TB ####

data <- data %>%
  mutate(
    micropos = case_when(
      ev_micro_test == "Microbiological confirmation" &
        ev_TBtype == "Pulmonary TB (PTB)" ~ 1,
      TRUE ~ 0
    ),
    clindx = case_when(
      ev_micro_test == "Clinical confirmation" |
        ev_TBtype == "Extra pulmonary TB (EPTB)" ~ 1, # unconfirmed TB deaths aren't classified as either
      TRUE ~ 0
    )
  )

#### Survival-estimated cumulative hazard of diagnosed recurrent TB by day 540 ####

diagnosis_survival_dataset <- data %>%
  select(record_id, end_reason, recurrence, end_days, txcompl_endreason_days) %>%
  mutate(
    event = recurrence == 1,
    time = end_days
  )

diagnosis_survival_dataset_truncated <- diagnosis_survival_dataset %>%
  mutate(
    event = if_else(time > 540, FALSE, event),
    time = pmin(time, 540)
  )

diagnosis_survfit_540 <- survival::survfit(
  Surv(time, event) ~ 1,
  data = diagnosis_survival_dataset_truncated
)

diagnosis_cumhaz_540_summary <- summary(
  diagnosis_survfit_540,
  times = 540,
  extend = TRUE
)

diagnosis_cumhaz_540 <- tibble(
  time = diagnosis_cumhaz_540_summary$time,
  cumhaz = diagnosis_cumhaz_540_summary$cumhaz,
  cumhaz_lower = diagnosis_cumhaz_540_summary$lower,
  cumhaz_upper = diagnosis_cumhaz_540_summary$upper
)

#### Fixed empirical summaries ####

fixed_empirical_inputs <- list(
  data = data,
  dictionary = dictionary,
  
  n_total = nrow(data),
  cumulative_recurrence = mean(data$recurrence),
  cumulative_micropos = mean(data$micropos),
  cumulative_clindx = mean(data$clindx),
  cumulative_any_diagnosed = mean(data$clindx | data$micropos),
  
  recurrence_n = sum(data$recurrence == 1, na.rm = TRUE),
  micropos_n = sum(data$micropos == 1, na.rm = TRUE),
  
  diagnosis_survival_dataset_truncated = diagnosis_survival_dataset_truncated,
  diagnosis_cumhaz_540 = diagnosis_cumhaz_540,
  
  reported_symptom_duration_recurrence = data %>%
    filter(recurrence == 1) %>%
    summarise(
      mean = mean(ev_sym_durmax, na.rm = TRUE),
      meanlog = mean(log(ev_sym_durmax), na.rm = TRUE),
      median = median(ev_sym_durmax, na.rm = TRUE),
      sd = sd(ev_sym_durmax, na.rm = TRUE),
      sdlog = sd(log(ev_sym_durmax), na.rm = TRUE),
      n_nonmissing = sum(!is.na(ev_sym_durmax))
    ),
  
  reported_symptom_duration_micropos = data %>%
    filter(micropos == 1) %>%
    summarise(
      mean = mean(ev_sym_durmax, na.rm = TRUE),
      median = median(ev_sym_durmax, na.rm = TRUE),
      sd = sd(ev_sym_durmax, na.rm = TRUE),
      n_nonmissing = sum(!is.na(ev_sym_durmax))
    ),

  diagnosis_timing_recurrence_by_540 = data %>%
    filter(recurrence == 1, txcompl_endreason_days <= 540) %>%
    summarise(
      mean_diagnosis_day = mean(txcompl_endreason_days, na.rm = TRUE),
      median_diagnosis_day = median(txcompl_endreason_days, na.rm = TRUE),
      sd_diagnosis_day = sd(txcompl_endreason_days, na.rm = TRUE),
      n = n()
    ),
  
  diagnosis_timing_micropos_by_540 = data %>%
    filter(micropos == 1, txcompl_endreason_days <= 540) %>%
    summarise(
      mean_diagnosis_day = mean(txcompl_endreason_days, na.rm = TRUE),
      median_diagnosis_day = median(txcompl_endreason_days, na.rm = TRUE),
      sd_diagnosis_day = sd(txcompl_endreason_days, na.rm = TRUE),
      n = n()
    ),
  
  ### Symptom reporting in home vs phone screening
  relative_symptom_reporting_phone = 
    data %>% group_by(arm) %>% 
    pivot_longer(cols = ends_with("sympwith_maxdur"), 
                 names_to = "month",
                 values_to = "symdur",
                 names_pattern = "m(\\d+)_sympwith_maxdur") %>%
    mutate(symptoms = !is.na(symdur)) %>%
    select(symdur, symptoms, arm, month) %>%
    summarise(mean(symptoms)) %>% 
    pivot_wider(names_from = arm, values_from = `mean(symptoms)`) %>%
    reframe(relative_reporting = `Telephonic` / `Home Visit`)
)
  
dir.create("outputs", showWarnings = FALSE)

saveRDS(
  fixed_empirical_inputs,
  file = "outputs/fixed_empirical_inputs.rds"
)

write_csv(
  fixed_empirical_inputs$reported_symptom_duration_recurrence,
  "outputs/reported_symptom_duration_recurrence.csv"
)

write_csv(
  fixed_empirical_inputs$reported_symptom_duration_micropos,
  "outputs/reported_symptom_duration_micropos.csv"
)