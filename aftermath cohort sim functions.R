library(tidyverse)
library(survival)
library(survminer)
library(Temporal)
library(pROC)
library(conflicted)
conflicts_prefer(dplyr::select, dplyr::filter)

create_cohort <- function(cohort_params)
{
  with(cohort_params, {
    cohort <- data.frame(ID = 1:N)
  
    #### Symptomatic TB onset ####
    
    # timing
    # used python script to estimate gamma parameters
  
    # # if I change the mean while keeping the sd proportional, the scale should change in proportion to the mean the shape shouldn't change. So:
    scale <- 113.88529634851781 * recurrence_time_mean/(192)
    shape <- 1.6859068628739582 *recurrence_time_shapefactor
    # And I can vary the shape (mean/sd)^2 as another parameter.
    
    # plot_densities <- data.frame(symptom_onset = rgamma(shape = shape, scale = scale, n = 10000))
    # ggplot(plot_densities, aes(x=symptom_onset/30)) + 
    #   geom_density() + xlim(0,24) + xlab("Months")
    
    # total events
    incidence_18mo_micropos <- incidence_18mo * proportion_micro_pos
    # within 18 months; extend out indefinitely
    incidence_total = incidence_18mo / pgamma(18*30, shape = shape, scale = scale)
    
    cohort$TB <- rbinom(prob = incidence_total, size = 1, n = N)
    cohort$symptom_onset [cohort$TB == 1] <- rgamma(n = sum(cohort$TB == 1), shape = shape, scale = scale)
    
    # sputum+ pulmonary?
    cohort <- cohort %>% mutate(pulmonary_with_micro = case_when(
      TB == 1 ~ rbinom(prob = incidence_18mo_micropos/incidence_18mo, size = 1, n = n()),
      TRUE ~ 0))
    
    #### Timing of routine diagnosis ####
    # Use negative binomial to increase variance (trial data mean duration 17 days, sd 14 days)
    # And scale by (t/120)^1/4 for mean and var 
    # An alternative parametrization (often used in ecology) is by the mean mu (see above), 
    # and size, the dispersion parameter, where prob = size/(size+mu). The variance is mu + mu^2/size in this parametrization.
    # Var = mu + mu^2/size --> size = (mu^2)/(sd^2 - mu), varies by same factor as mean and var.
    
    
    cohort <- cohort %>% mutate(diagnosis_routine = case_when(
      TB == 1 ~ symptom_onset + rnbinom(n = n(), 
                                        size = symptom_duration_mean_reported^2/(symptom_duration_sd_reported^2) * 
                                          (symptom_onset/120)^(symptom_duration_timescale),
                                        mu = symptom_duration_mean_reported/symptom_underestimation_factor/home_visit_passive_detection_impact * 
                                          (symptom_onset/120)^(symptom_duration_timescale)),
            TRUE ~ NA_real_)) %>%
      mutate(diagnosis_routine_original = diagnosis_routine,
             TB_original = TB,
             symptom_onset_original = symptom_onset,
             sputum_onset_original = sputum_onset)
    
    
    # plot_densities <- plot_densities %>% 
    #   mutate(diagnosis_routine = rnbinom(n = 10000, 
    #                size = symptom_duration_mean_reported^2/(symptom_duration_sd_reported^2) * 
    #                  (symptom_onset/120)^(symptom_duration_timescale),
    #                mu = symptom_duration_mean_reported/symptom_underestimation_factor) * 
    #          (symptom_onset/120)^(symptom_duration_timescale))
    # ggplot(plot_densities) + geom_density(aes(x = symptom_onset/30, col = "symptom onset")) + 
    #   # geom_density(aes(x = diagnosis_routine), col='blue')+ 
    #   geom_density(aes(x = (diagnosis_routine + symptom_onset)/30, col = "routine diagnosis")) +
    #   xlim(0,24) + xlab("Months") +
    #   guides(colour = guide_legend(reverse=T))  
    # 
    # cohort %>% filter(TB==1) %>% ggplot(aes(x=symptom_onset, y=diagnosis_routine)) + geom_point()
    # cohort %>% filter(TB==1) %>% ggplot(aes(x=symptom_onset, y=diagnosis_routine - symptom_onset)) + geom_point()
    # cohort %>% filter(TB==1, diagnosis_routine <= 30*18) %>% summarise(mean(diagnosis_routine - symptom_onset), sd(diagnosis_routine - symptom_onset))
    
    #### Subclinical TB ####
    
    mean_duration_subclinical <- symptom_duration_mean_reported/symptom_underestimation_factor*
      duration_ratio_subclinical_symptomatic/proportion_ever_subclinical
    
    cohort <- cohort %>% mutate(
      ever_subclinical_sputumpos = case_when(pulmonary_with_micro == 1 ~ rbinom(prob = proportion_ever_subclinical, size = 1, n = n()),
                                             TRUE ~ 0),
      # In these people, sputum onset precedes symptom onset by a gamma-distributed time
      sputum_onset = case_when(ever_subclinical_sputumpos == 1 ~ 
                                 pmax(0, symptom_onset - mean_duration_subclinical* rgamma(scale = duration_subclinical_cv^2, 
                                                                                           shape = 1/(duration_subclinical_cv^2), n = n())),
                               ever_subclinical_sputumpos == 0 & pulmonary_with_micro == 1 ~ 
                                 runif(n = n(), min = symptom_onset, max = diagnosis_routine),
                               TRUE ~ NA_real_))
    
    
    #### Prediction ####
    
    # simulate an ROC curve with desired AUC,using approach from
    # https://stats.stackexchange.com/questions/422926/generate-synthetic-data-given-auc
    t <- sqrt(log(1/(1-auc)^2))
    z <- t-((2.515517 + 0.802853*t + 0.0103328*t^2) / 
              (1 + 1.432788*t + 0.189269*t^2 + 0.001308*t^3))
    d <- z*sqrt(2)
    
    n <- 10000
    x <- c(rnorm(n/2, mean = 0), rnorm(n/2, mean = d))
    y <- c(rep(0, n/2), rep(1, n/2))
    
    roc.curve <- roc(y, x)
    
    
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

# validate the proportion of subclinical cases at baseline and at 6 months, and reject if inconsistent with Report data:
# sputum-positive TB from baseline (from the start of post-treatment follow-up) in up to 9/(9+67) = 12% [6-20%] of recurrences 
# cross-sectional prevalence of asymptomatic TB at 6 months of between 0.5% and 1.6%. 

check_subclinical <- function(cohort, cohort_params)
{
  # check the proportion of subclinical cases at baseline:
  subclinical_baseline_amongTB <- cohort %>% filter(TB == 1) %>% summarize(sum(sputum_onset==0, na.rm=T))/
    cohort %>% filter(TB == 1) %>% summarize(n())
  subclinical_6mo_amongcohort <- cohort %>% summarize(sum(sputum_onset<180 & symptom_onset > 180, na.rm=T))/
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

# This function just says who's contacted by screening, and makes baseline adjustments to later TB based on counseling and prevention.
get_screening_coverage <- function(cohort, screening_design, cohort_params) {
  
  # Track screening contacts made
  covered_screening <- 
    mapply(function(loc, cov) sapply(cohort$risk_quantile, function(risk) 
      (risk >= (1-cov)) * (rbinom(n = 1, size = 1, prob = intervention_parameters$coverage[[loc]]))), 
      screening_design$screening_location, 
      screening_design$target_coverage)
  
  # Throughout, if get counseling about recurrence up front, adjust timing by intentional_counseling_passive_detection_impact
  # for all with symptom onset before intentional_counseling_passive_detection_duration. 
  
  if (screening_design$counseling_coverage[1] > 0)
  {
    cohort <- cohort %>% mutate(
      diagnosis_routine = case_when(
        (TB == 1) & (risk_quantile > (1 - screening_design$counseling_coverage[1])) & (diagnosis_routine > symptom_onset) & 
          (symptom_onset < cohort_params$intentional_counseling_passive_detection_duration) ~ 
          diagnosis_routine - (diagnosis_routine - symptom_onset)*(1-cohort_params$intentional_counseling_passive_detection_impact),
        TRUE ~ diagnosis_routine))
  }
  
  # And apply prevention here too, as changing whether TB will develop. 
  if (screening_design$prevention_coverage[1] > 0)
  {
    cohort <- cohort %>% mutate(
      TB = case_when(
        (TB == 1) & (risk_quantile > (1 - screening_design$prevention_coverage[1])) ~ 
          rbinom(n = n(), size = 1, prob = cohort_params$prevention_efficacy),
        TRUE ~ TB),
      symptom_onset = case_when(TB == 0 & !is.na(symptom_onset) ~ NA, TRUE ~ symptom_onset),
      sputum_onset = case_when(TB == 0 & !is.na(sputum_onset) ~ NA, TRUE ~ sputum_onset),
      diagnosis_routine = case_when(TB == 0 & !is.na(diagnosis_routine) ~ NA, TRUE ~ diagnosis_routine),
      pulmonary_with_micro = case_when(TB == 0 & pulmonary_with_micro == 1 ~ 0, TRUE ~ pulmonary_with_micro),
      ever_subclinical_sputumpos = case_when(TB == 0 & ever_subclinical_sputumpos == 1 ~ 0, TRUE ~ ever_subclinical_sputumpos))
  }
  
  # # And Once screened at home for the first time, adjust future routine diagnosis timing by home_visit_passive_detection_impact
  # ## identify first home contact
  # if (sum(screening_design$screening_location == "home") > 0) {
  #   if (sum(screening_design$screening_location == "home") == 1) {
  #     first_home_contact <- ifelse(covered_screening[,screening_design$screening_location == "home"] == 1, 
  #                                  screening_design$timing_months[screening_design$screening_location == "home"] * 30, 
  #                                  NA)
  #   } else {
  #     first_home_contact <- 
  #       (screening_design$timing_months[screening_design$screening_location == "home"] * 30)[
  #         apply(covered_screening[,screening_design$screening_location == "home"], 1, function(x) which(x == 1)[1])
  #       ]
  #     
  #   }
  # } else first_home_contact <- rep(NA, nrow(cohort))
  # ## and adjust later care-seeking accordingly 
  # cohort <- cohort %>% mutate(
  #   diagnosis_routine = case_when(
  #     is.na(first_home_contact) ~ diagnosis_routine,
  #     symptom_onset > first_home_contact ~ diagnosis_routine - (diagnosis_routine - symptom_onset)*(1-cohort_params$home_visit_passive_detection_impact),
  #     TRUE ~ diagnosis_routine),
  #   sputum_onset = case_when(
  #     sputum_onset > diagnosis_routine ~ NA,
  #     TRUE ~ sputum_onset))
  
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
      TB != 1 ~ 0,
      screened_current_round==0 ~ 0, 
      timing_months*30 >= diagnosis_routine ~ 0,
      screening_method %in% c("symptoms", "both") &
        timing_months*30 > symptom_onset ~ 
        rbinom(n = n(), size = 1, prob = intervention_parameters$sensitivity_symptoms[[screening_location]]),
      # micro screening
      screening_method %in% c("micro", "both") & 
        timing_months*30 > sputum_onset ~ 
        rbinom(n = n(), size = 1, prob = intervention_parameters$success_sputum[[screening_location]]),
      TRUE ~ 0),
    
    detection_timing = case_when(
      detected_current_round == 1 & 
        (is.na(detection_timing) | timing_months*30 < detection_timing) ~ timing_months*30,
      TRUE ~ detection_timing)
  )
  return(cohort %>% select(-c(screened_current_round, detected_current_round)))  
}

 

# function we'll use to apply the same steps to all 
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


plot_screening <- function(cohort, screening_design, colorfill = TRUE)
{
  colors <- c("Symptomatic" = "blue", "NAAT+" = "red")
  fillcolors <- c("Micro+ Pulmonary" = "red", "Clinical or extrapulmonary" = "gray")
  linetypes <- c("symptoms" = "dotted", "micro" = "dashed", "both" = "dotdash")
  plotdata <- cohort %>% filter(TB==1) %>% arrange(diagnosis_routine) %>% 
    mutate(newID = row_number(),
           TBtype = case_when(pulmonary_with_micro == 1 ~ "Micro+ Pulmonary",
                              TRUE ~ "Clinical or extrapulmonary"),
           detected = !is.na(detection_timing))
  plot <- ggplot(plotdata) + 
    geom_vline(data = screening_design, aes(xintercept = timing_months, linetype = screening_method)) + 
    scale_color_manual("Time with TB", values = colors) + 
    scale_linetype_manual("Screening Method", values = linetypes) + 
    ylab("Cases arranged by timing of routine diagnosis") + 
    xlab("Months since prior treatment completion") +
    # if detection_timing is NA, plot full segment
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  diagnosis_routine/30, y = newID + 0.15, col='NAAT+'), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  detection_timing/30, y = newID + 0.15, col='NAAT+'), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & sputum_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID + 0.15, col='NAAT+'), size=0.8, lty="11") + 
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  diagnosis_routine/30, y = newID - 0.15, col='Symptomatic'), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  detection_timing/30, y = newID - 0.15, col='Symptomatic'), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & symptom_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID - 0.15, col='Symptomatic'), size=0.8, lty="11") + 
    geom_point(aes(x = diagnosis_routine/30, y = newID), fill = 'black', pch=21) + 
    geom_point(aes(x = detection_timing/30, y = newID), pch=10, col = "green2", size=2) +
    scale_x_continuous(breaks = seq(0, 30, 6), limits = c(0,30)) +
    theme_minimal() +
    # reverse order of color legend
    guides(colour = guide_legend(reverse=T),
           # remove lty legend
           linetype = "none") +
    # overlay legend on bottom right corner of plot
    theme(legend.position = c(0.8, 0.2), 
          legend.background = element_rect(fill = "transparent"),
          legend.box.background = element_rect(fill = "white"))
  
  if(colorfill) plot <- plot + 
    scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
    geom_point(aes(x = diagnosis_routine/30, y = newID, fill = TBtype), pch=21)  
    
  return(plot)
  
}

  
#### Outcomes of interest #### 

# Cases detected (is this valuable in itself? this simulation assumes all with be detected eventually, or deaths aren't observed)
## Cases detected "early" by intervention
# Time with recurrent TB
## Time sputum+ (associated with transmission)
## Time symptomatic (associated with lung damage and with morality risk??)


time_with_tb <- function(cohort)
{
  # add time 
  outcomes <- cohort %>% mutate(
    symptom_days_soc = diagnosis_routine_original - symptom_onset_original,
    sputum_days_soc = diagnosis_routine_original - sputum_onset_original,
    symptom_days_screening = case_when(detection_timing > symptom_onset & 
                                         detection_timing < diagnosis_routine ~ 
                                         diagnosis_routine - detection_timing,
                                       TRUE ~ diagnosis_routine - symptom_onset),
    sputum_days_screening = case_when(detection_timing > sputum_onset &
                                        detection_timing < diagnosis_routine ~
                                        diagnosis_routine - detection_timing,
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


# Costs
## Cost per case detected early or prevented (combine these? Though not equally valuable?)
## Cost per month with TB prevented 
### Overall
### Symptomatic
### Sputum+

costs <- function(screening_design, cohort, covered_screening,
                  initial_contact_cost_home = 15, 
                  initial_contact_cost_phone = 5, #** from hojoon's team, need to talk with them about uncertainty
                  sputum_test_cost = 16 # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0214675 
)
{
  initial_contact_cost = list(home = initial_contact_cost_home, phone = initial_contact_cost_phone)
  # treatment_cost = 1000
  
  symptom_prevalence_nonTB = 0.05 #** Will get from Aye's analysis of symptoms after retraining. 
  # Or take from prev survey: 27% of all people with past TB had a positive symptom screen. !!
  
  # count units of:
  # total people contacted
  # number of people with symptomatic TB, sputum+ asymptomatic TB, and non-TB symptoms
  # number of people with TB detected by screening
  
  # total people with initial contact:
  contacted <- colSums(covered_screening)
  contact_cost <- contacted * unlist(initial_contact_cost[screening_design$screening_location]) # cost
  
  # people with symptomatic TB, who are targeted and reached:
  sxtb <- colSums(cohort$TB * sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset < s) * 
                    covered_screening *
                    sapply(screening_design$timing_months*30, function(s) s <= cohort$diagnosis_routine), na.rm=T)
  
  # people with sputum+ asymptomatic TB, who are targeted and reached:
  asxtb <- colSums(cohort$TB * sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset > s & cohort$sputum_onset < s) * 
                     covered_screening *
                     sapply(screening_design$timing_months*30, function(s) s <= cohort$diagnosis_routine), na.rm=T)
  
  # people with non-TB but with symptoms (never TB, or after diagnosis, or not yet sputum+ or symptom+:
  nontbsx <- colSums((cohort$TB ==0 |
                        sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset > s & (cohort$sputum_onset > s | is.na(cohort$sputum_onset))) |
                        sapply(screening_design$timing_months*30, function(s) s > cohort$diagnosis_routine)) & 
                       covered_screening) *
    symptom_prevalence_nonTB
  
  # sputum tests, 
  # if universal testing: 
  sputa <- contacted * unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
    screening_design$screening_method %in% c("micro", "both") + 
    # if symptom-based testing: 
    (sxtb + nontbsx) * unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
    (ifelse(screening_design$screening_method == "symptoms",1,0))
  
  costs <- 
    # home visits (including a visit if symptom+ on phone call?)
    # sputum tests
    contact_cost + sputa * sputum_test_cost
  # treatments (incremental?)
  # patient costs?
  return(costs)
}

