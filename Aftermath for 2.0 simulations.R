library(tidyverse)
library(survival)
library(survminer)
library(Temporal)
library(binom)
library(clipr)
library(sensitivity)
library(pROC)
library(conflicted)
conflicts_prefer(dplyr::select, dplyr::filter)


# Define parameters 

N_cohort = 10000 # balancing run time vs stochastic effects
N_samples <- 5000 # for capturing paramter uncertainty and stochasticity

set.seed(12345)

#### Incorporate uncertainty: Sample cohort_params from posterior distributions ####
cohort_param_ranges <- list(
  N = N_cohort, # number of participants in each cohort sample
  
  recurrence_shape = 0.75,
  recurrence_scale = 406, 
  recurrence_time_mean = c(150, 218), # mean of diagnosis time minus symptom duration, of those diagnosed by day 540
  # Based on empirical mean of 184, and will multiply times fitted weibull parameters of shape 0.84, scale 406 
  recurrence_time_shapefactor = c(0.9, 1.1), # variation in (mean/sd)^2, compared to aftermath calculation
  incidence_18mo = c(0.074, 0.110), # cumulative incidence of *symptoms* by 18 months
  proportion_micro_pos = c(0.43, 0.63), #binom.confint(0.489*96, 96)
  auc = c(0.56, 0.77),
  
  symptom_duration_meanlog_reported = 2.68,
  symptom_duration_sdlog_reported = c(0.6, 0.8), 
  symptom_duration_timescale = 0, #c(0.125, 0.375), # exponential increase in mean and sd, 240d vs 120d
  symptom_underestimation_factor = c(0.25, 0.5), 
  
  proportion_ever_subclinical = c(0.6, 0.9), # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2), 
  duration_subclinical_cv = c(0.5, 1.5),
  subclinical_baseline_amongTB_max = 0.20, 
  subclinical_6m_amongcohort_min = 0.004, 
  subclinical_6m_amongcohort_max = 0.017, 
  
  coverage_phone = c(0.9, 1.0), 
  coverage_home_reduction = c(0.75, 0.95), #0.81/0.95 = 0.85
  sensitivity_symptoms_home = c(0.85, 0.95),
  sensitivity_symptoms_phone_reduction = c(0.6, 0.82),
  success_sputum_home = c(0.85, 0.95),
  success_sputum_phone_reduction = c(0.6, 1.0),
  
 
  # new addition for aftermath 2.0 simulations: a proportion of TB ends in death rather than notification; both are at the time of "diagnosis_routine"
  case_fatality = c(0.05, 0.15)
)

cohort_param_samples <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 

# Simulate cohort without intervention

create_cohort <- function(cohort_params)
{
  with(cohort_params, {
    cohort <- data.frame(ID = 1:N)
    
    #### Symptomatic TB onset ####
    
    # timing
    
    # # if I change the mean while keeping the sd proportional, the scale should change in proportion to the mean, and the shape shouldn't change. So:
    scale <- recurrence_scale * recurrence_time_mean/(184)
    shape <- recurrence_shape *recurrence_time_shapefactor
    
    # total events (micro+ symptom onset)
    incidence_18mo_micropos <- incidence_18mo * proportion_micro_pos
    # within 18 months; extend out indefinitely
    incidence_total = incidence_18mo / pweibull(18*30, shape = shape, scale = scale)
    
    cohort$TB <- rbinom(prob = incidence_total, size = 1, n = N)
    cohort$TBdeath <- ifelse(cohort$TB == 1, rbinom(prob = case_fatality, size = 1, n = sum(cohort$TB == 1)),
                             NA)
    cohort$symptom_onset [cohort$TB == 1] <- rweibull(n = sum(cohort$TB == 1), shape = shape, scale = scale)
    
    # sputum+ pulmonary?
    cohort <- cohort %>% mutate(pulmonary_with_micro = case_when(
      TB == 1 ~ rbinom(prob = incidence_18mo_micropos/incidence_18mo, size = 1, n = n()),
      TRUE ~ 0))
    
    #### Timing of routine diagnosis ####
    # Use lognormal to capture trial distribution of reported symptom durations (trial data mean duration 17 days, sd 14 days)
    # And scale by (t/mean onset time)^1/4 for mean and var (so by 1/2 for sd)
    
    cohort <- cohort %>% mutate(diagnosis_routine = case_when(
      TB == 1 ~ symptom_onset + 
        (1/symptom_underestimation_factor * 
           (symptom_onset/recurrence_time_mean)^(symptom_duration_timescale) * 
           rlnorm(n = n(), meanlog = symptom_duration_meanlog_reported,
                  sdlog = symptom_duration_sdlog_reported)),
      TRUE ~ NA_real_))
    
    
    #### Subclinical TB ####
    
    # need to calculate the duration of symptoms, so that we can back-calculate the mean duration of subclinical TB 
    symptom_duration_mean_reported <- exp(symptom_duration_meanlog_reported + symptom_duration_sdlog_reported^2/2)
    
    mean_duration_subclinical <- symptom_duration_mean_reported/symptom_underestimation_factor*
      duration_ratio_subclinical_symptomatic/proportion_ever_subclinical
    
    cohort <- cohort %>% mutate(
      ever_subclinical_sputumpos = case_when(pulmonary_with_micro == 1 ~ rbinom(prob = proportion_ever_subclinical, size = 1, n = n()),
                                             TRUE ~ 0),
      # In these people, sputum onset precedes symptom onset by a gamma-distributed time. For the rest, it happens at a random time during the symptomatic period. 
      sputum_onset = case_when(ever_subclinical_sputumpos == 1 ~ 
                                 pmax(0, symptom_onset - mean_duration_subclinical* rgamma(scale = duration_subclinical_cv^2, 
                                                                                           shape = 1/(duration_subclinical_cv^2), n = n())),
                               ever_subclinical_sputumpos == 0 & pulmonary_with_micro == 1 ~ 
                                 runif(n = n(), min = symptom_onset, max = diagnosis_routine),
                               TRUE ~ NA_real_))
    
    #### Risk prediction ####
    
    roc.curve <- simulate_auc(auc = auc)
    
    # assign risk rankings based on sensitivity and specificity of roc.curve for TB:
    # First, assign scores to the TB cases such that a given prediction score cutoff of P proportion 1-P of the TB cases
    # And then, for the non-cases, assign with probability (1-spec) at a 
    cohort$TB_risk_rank[cohort$TB == 1] <- (1:sum(cohort$TB == 1))/sum(cohort$TB == 1)
    cohort$nonTB_risk_rank[cohort$TB == 0] <- (1:sum(cohort$TB == 0))/sum(cohort$TB == 0)
    for (i in 1:N)
    {TB_status <- cohort$TB[i]
    if (TB_status == 1)
    {cohort$risk_score[i] <- cohort$TB_risk_rank[i]}
    else
    {cohort$risk_score[i] <- 1 - roc.curve$specificities[which.min(roc.curve$sensitivities > cohort$nonTB_risk_rank[i])]}
    }
    
    # but actually, we may want to target a certain % of patients, using risk quantiles: 
    cohort <- cohort %>% arrange(risk_score) %>% mutate(risk_quantile = row_number()/n())
    
    return(cohort)
  })
}

# AUC helper function
# simulate an ROC curve with desired AUC,using approach from
# https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc
simulate_auc <- function(auc)
{
  t <- sqrt(log(1/(1-auc)^2))
  z <- t-((2.515517 + 0.802853*t + 0.0103328*t^2) / 
            (1 + 1.432788*t + 0.189269*t^2 + 0.001308*t^3))
  d <- z*sqrt(2)
  
  n <- 10000
  x <- c(rnorm(n/2, mean = 0), rnorm(n/2, mean = d))
  y <- c(rep(0, n/2), rep(1, n/2))
  
  roc.curve <- roc(y, x)
  return(roc.curve)
}

# validate the proportion of subclinical cases at baseline and at 6 months, and reject if inconsistent with Report data:
# sputum-positive TB from baseline (from the start of post-treatment follow-up) in up to 9/(9+67) = 12% [6-20%] of recurrences 
# cross-sectional prevalence of asymptomatic TB at 6 months of between 0.5% and 1.6%. 
check_subclinical <- function(cohort, cohort_params)
{
  # check the proportion of subclinical cases at baseline:
  subclinical_baseline_amongTB <- cohort %>% filter(TB == 1) %>% summarize(sum(sputum_onset==0, na.rm=T))/
    cohort %>% filter(TB == 1) %>% summarize(n())
  subclinical_6mo_amongcohort <- cohort %>% summarize(sum(sputum_onset<180 & symptom_onset >= 180, na.rm=T))/
    cohort %>% summarize(n())
  if (subclinical_baseline_amongTB < cohort_params$subclinical_baseline_amongTB_max & 
      subclinical_6mo_amongcohort > cohort_params$subclinical_6m_amongcohort_min & 
      subclinical_6mo_amongcohort < cohort_params$subclinical_6m_amongcohort_max)
    return(TRUE) else
      return(FALSE)
}

### Apply screening interventions ###
# Model each round with:
# timing (months post treatment), 
# target % population.
# screening method (symptom, micro, both), and
# screening location (home, telephonic), with
## estimated sensitivity for symptoms and
## ability to get sputum for micro, and
## estimated coverage as % of target population.

# This function just says who's contacted by screening, and makes baseline adjustments to later TB based on counseling and (in commented out section) prevention.
get_screening_coverage <- function(cohort, screening_design, cohort_params) {
  
  # Track screening contacts made
  covered_screening <- 
    mapply(function(loc, cov) sapply(cohort$risk_quantile, function(risk) 
      (risk >= (1-cov)) * (rbinom(n = 1, size = 1, prob = intervention_parameters$coverage[[loc]]))), 
      screening_design$screening_location, 
      screening_design$target_coverage)
  
  
  return(list("covered" = covered_screening, 
              "cohort" = cohort))
}



# Identify the dates when any cases are detected by screening

apply_screening_round <- function(
    cohort, 
    covered_screening_column,
    timing_months = 6,
    screening_method = "symptoms", # (symptoms, micro, or both)
    screening_location = "home", # (home, phone)
    intervention_parameters
)
{
  if(!("detection_timing" %in% colnames(cohort))) cohort$detection_timing <- NA
  
  # symptom screening
  cohort$screened_current_round <- covered_screening_column
  
  cohort <- cohort %>% mutate(
    detected_current_round = case_when(
      TB != 1 ~ 0, # no TB
      screened_current_round==0 ~ 0, # not screened
      timing_months*30 >= diagnosis_routine ~ 0, # already diagnosed
      screening_method  == "both" & 
        timing_months*30 >= symptom_onset & timing_months*30 >= sputum_onset  ~ # both sx+ and sput+
        rbinom(n = n(), size = 1, prob = 1 - (1 - intervention_parameters$success_sputum[[screening_location]])*
                 (1-intervention_parameters$sensitivity_symptoms[[screening_location]])), # assuming independent probabilities
      screening_method  == "both" &
        timing_months*30 >= symptom_onset ~ # sx+ only
        rbinom(n = n(), size = 1, prob = intervention_parameters$sensitivity_symptoms[[screening_location]]),
      screening_method  == "both" &
        timing_months*30 >= sputum_onset ~ # sput+ only
        rbinom(n = n(), size = 1, prob = intervention_parameters$success_sputum[[screening_location]]),
      # symptom screening only
      screening_method =="symptoms" &
        timing_months*30 >= symptom_onset ~ 
        rbinom(n = n(), size = 1, prob = intervention_parameters$sensitivity_symptoms[[screening_location]]),
      # micro screening only
      screening_method == "micro" & 
        timing_months*30 >= sputum_onset ~ 
        rbinom(n = n(), size = 1, prob = intervention_parameters$success_sputum[[screening_location]]),
      TRUE ~ 0),
    
    detection_timing = case_when(
      detected_current_round == 1 & 
        (is.na(detection_timing) | timing_months*30 < detection_timing) ~ timing_months*30,
      TRUE ~ detection_timing)
  )
  return(cohort %>% select(-c(screened_current_round, detected_current_round)))  
}



# function we'll use to apply all the intervention screening encounters
apply_intervention <- function(cohort, design, intervention_parameters, cohort_params) {
  # apply screening rounds
  screened <- get_screening_coverage(cohort, design, cohort_params)
  covered_screening <- screened$covered
  cohort_screened <- screened$cohort
  
  for (r in 1:nrow(design)) {
    cohort_screened <- apply_screening_round(cohort_screened, 
                                             covered_screening[,r],
                                             timing_months = design$timing_months[r], 
                                             screening_method = design$screening_method[r], 
                                             screening_location = design$screening_location[r],
                                             intervention_parameters)
  }
  
  return(cohort_screened)
}


# Visualize the cohort and screening detections
plot_screening <- function(cohort, screening_design, colorfill = TRUE)
{
  colors <- c("Symptomatic" = "royalblue", "NAAT+" = "red")
  fillcolors <- c("Micro+ Pulmonary" = "red", "Clinical or extrapulmonary" = "gray")
  # linetypes <- c("symptoms" = "dotted", "micro" = "dashed", "both" = "dotdash")
  linetypes <- c("Not averted" = "solid", "Averted by screening" = "11")
  plotdata <- cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
    mutate(newID = row_number(),
           TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                              TRUE ~ "Clinical or extrapulmonary"),
           detected = !is.na(detection_timing))
  plot <- ggplot(plotdata) + 
    geom_vline(data = screening_design, aes(xintercept = timing_months, linetype = screening_method),
               linetype = "dotdash") + 
    # add small horizontal text "symptoms screening" along each vertical line
    geom_text(data = screening_design, aes(x = timing_months, y = 0), label = "symptom screening", 
              angle = 90, size = 3, hjust= 0, vjust=-0.5) + 
    guides(colour = guide_legend(reverse=T)) +
    scale_color_manual("Time with TB", values = colors) + 
    # scale_linetype_manual("Screening Method", values = linetypes) + 
    ylab("Cases arranged by timing of routine diagnosis") + 
    xlab("Months since prior treatment completion") +
    scale_linetype_manual("Outcome", values = linetypes) + 
    # if detection_timing is NA, plot full segment
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  diagnosis_routine/30, y = newID + 0.15, col='NAAT+',
                     lty = "Not averted"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  detection_timing/30, y = newID + 0.15, col='NAAT+',
                     lty = "Not averted"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & sputum_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID + 0.15, col='NAAT+',
                     lty = "Averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  diagnosis_routine/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Not averted"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  detection_timing/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Not averted"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & symptom_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Averted by screening"), size=0.8) + 
    geom_point(aes(x = diagnosis_routine/30, y = newID), fill = 'black', pch=21) + 
    geom_point(aes(x = detection_timing/30, y = newID), pch=10, col = "green2", size=2) +
    scale_x_continuous(breaks = seq(0, 30, 6), limits = c(0,30)) +
    theme_minimal() +
    # overlay legend on bottom right corner of plot
    theme(legend.position = c(0.8, 0.3), 
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "white", color = NA))
  
  
  if(colorfill) plot <- plot + 
    scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
    geom_point(aes(x = diagnosis_routine/30, y = newID, fill = TBtype), pch=21)  
  
  plot <- plot + guides(color = guide_legend(order = 1))
  
  
  return(plot)
  
}


# Count time with symptomatic and tiem with infectious TB, with and without the intervention
time_with_tb <- function(cohort, limit_days = NULL)
{
  if (!is.null(limit_days)) cohort <- cohort %>% mutate(
    diagnosis_routine = pmin(diagnosis_routine, limit_days),
    symptom_onset = pmin(symptom_onset, limit_days),
    sputum_onset = pmin(sputum_onset, limit_days),
    detection_timing = pmin(detection_timing, limit_days)
  )
  # add time 
  outcomes <- cohort %>% mutate(
    symptom_days_soc = diagnosis_routine - symptom_onset,
    sputum_days_soc = diagnosis_routine - sputum_onset,
    symptom_days_screening = case_when(detection_timing >= symptom_onset & 
                                         detection_timing <= diagnosis_routine ~ 
                                         detection_timing - symptom_onset,
                                       detection_timing < symptom_onset ~ 0, # don't need this for sxs but do for sputum
                                       TRUE ~ diagnosis_routine - symptom_onset),
    sputum_days_screening = case_when(detection_timing >= sputum_onset &
                                        detection_timing <= diagnosis_routine ~
                                        detection_timing - sputum_onset,
                                      detection_timing < sputum_onset ~ 0,
                                      TRUE ~ diagnosis_routine - sputum_onset))
  results <- outcomes %>% summarise(
    symptom_days_soc = sum(symptom_days_soc, na.rm = TRUE),
    symptom_days_screening = sum(symptom_days_screening, na.rm = TRUE),
    sputum_days_soc = sum(sputum_days_soc, na.rm = TRUE),
    sputum_days_screening = sum(sputum_days_screening, na.rm = TRUE)
  ) %>% 
    pivot_longer(names_sep = "_", cols = everything(), names_to = c("outcome", "unit", "scenario"))
  return(results)
}



#### Outcomes of interest #### 

# Cases detected (will need to add. original model assumes all with be detected eventually, but now we'll model some as deaths (or resolutions) rather than notificatoins)
## Cases detected "early" by intervention
# Time with recurrent TB
## Time sputum+ (associated with transmission)
## Time symptomatic (associated with lung damage and with morality risk??)


##### Run the above model

n <- 1

cohort <- create_cohort(cohort_param_samples[n,])

screening_design_aftermath2 <- 
  data.frame(
    "counseling_coverage" = 0, 
    "prevention_coverage" = 0,
    "timing_months" = c(3, 6, 12),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("both", "both", "both"),
    "screening_location" = c("home", "home", "home")
  ) %>% 
  arrange(timing_months)

# If running within a loop, uncomment lines below to skip and generate NAs for any simulation that doesn't match TBDM/CTRiumph data on subclinical TB
if (!check_subclinical(cohort, cohort_param_samples[n,])) {
    # # change the nth element of each item in results to NA
  # results[n,] <- rep(NA, ncol(results))
  # next
  print("Error: Fails subclinical validation checks")
}

intervention_parameters <- list(
  coverage = list("phone" = cohort_param_samples$coverage_phone[n], 
                  "home" = cohort_param_samples$coverage_phone[n] * cohort_param_samples$coverage_home_reduction[n]),
  sensitivity_symptoms = list("home" = cohort_param_samples$sensitivity_symptoms_home[n], 
                              "phone" = cohort_param_samples$sensitivity_symptoms_home[n] * cohort_param_samples$sensitivity_symptoms_phone_reduction[n]),
  success_sputum = list("home" = cohort_param_samples$success_sputum_home[n], 
                        "phone" = cohort_param_samples$success_sputum_home[n] * cohort_param_samples$success_sputum_phone_reduction[n]))


##### Run the interventions ####
outputs <- apply_intervention(cohort, screening_design_aftermath2, 
                                  intervention_parameters, cohort_param_samples[n,])  

# impact
mean(outputs$TB)
(recurred_ever = sum(outputs$TB))
(recurred_within_18m = sum(outputs$symptom_onset <= 30*18 & outputs$TB == 1)) # symptomatic
(diagnosed_within_18m_soc = sum(outputs$diagnosis_routine <= 30*18 & outputs$TB == 1 & outputs$TBdeath != 1)) # the rest are either still undiagnosed or died before diagnosis

(detected_by_screening <- sum(!is.na(outputs$detection_timing)))
(tb_time_first_18m_soc <- time_with_tb(outputs, limit_days = 30*18) %>% filter(scenario == "soc") %>% 
                  mutate(months = value/30) %>%
                  tibble::column_to_rownames('outcome') %>% select(months))
(tb_time_first_18m_averted <- time_with_tb(outputs, limit_days = 30*18)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                                          mutate(months_averted = (soc - screening)/30) %>% 
                    tibble::column_to_rownames('outcome')  %>% select(months_averted))

# time to event
outputs$event_soc <- outputs$diagnosis_routine <= 18*30 & outputs$TB == 1 & outputs$TBdeath != 1 # routine TB notification
outputs$event_screening <- (!is.na(outputs$detection_timing) & outputs$detection_timing <= 18*30) | outputs$event_soc # routine TB notification or detection through screening

outputs$timing_soc <- pmin(outputs$diagnosis_routine, 18*30) # notification or TB death, or end of follow up
outputs$timing_screening <- pmin(ifelse(is.na(outputs$detection_timing), 18*30, outputs$detection_timing), outputs$diagnosis_routine, 18*30)

mean(outputs$event_soc)
mean(outputs$event_screening)
summary(outputs$timing_soc[outputs$event_soc == 1])
summary(outputs$timing_screening[outputs$event_screening == 1])

km_soc <- survfit(Surv(timing_soc, event_soc) ~ 1, data = outputs)
km_screening <- survfit(Surv(timing_screening, event_screening) ~ 1, data = outputs)
plot(km_soc, col = "black")
lines(km_screening, col = "blue")
