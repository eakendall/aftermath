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
    
    scale <- recurrence_scale
    shape <- recurrence_shape
    
    # total events (micro+ symptom onset)
    cohort$TB <- rbinom(
      n = N,
      size = 1,
      prob = probability_ever_recur
    )
    
    cohort$symptom_onset [cohort$TB == 1] <- rweibull(n = sum(cohort$TB == 1), shape = shape, scale = scale)
    
    # sputum+ pulmonary?
    cohort <- cohort %>%
      mutate(
        pulmonary_with_micro = case_when(
          TB == 1 ~ rbinom(n = n(), size = 1, prob = proportion_micro_pos),
          TRUE ~ 0
        )
      )    
    #### Timing of routine diagnosis ####
    # Use lognormal to capture trial distribution of reported symptom durations (trial data mean duration 17 days, sd 14 days)
    
    # Under routine conditions, the time to routine diagnosis increases by factor programmatic_symptom_duration_factor
    
    cohort <- cohort %>%
      mutate(
        raw_symptom_duration = case_when(
          TB == 1 ~
            1 / reported_fraction_of_true_symptom_duration /
            home_visit_passive_detection_impact *
            programmatic_symptom_duration_factor *
            rlnorm(
              n = n(),
              meanlog = symptom_duration_meanlog_reported,
              sdlog = symptom_duration_sdlog_reported
            ),
          TRUE ~ NA_real_
        ),
        
        symptom_duration = case_when(
          TB == 1 ~
            pmin(
              raw_symptom_duration,
              max_symptom_duration_fraction_of_onset_time * symptom_onset
            ),
          TRUE ~ NA_real_
        ),
        
        diagnosis_routine = case_when(
          TB == 1 ~ symptom_onset + symptom_duration,
          TRUE ~ NA_real_
        )
      )
     
    #### Subclinical TB ####
    
    symptom_duration_mean_reported <- exp(
      symptom_duration_meanlog_reported +
        symptom_duration_sdlog_reported^2 / 2
    )
    
    true_symptom_duration_mean <-
      symptom_duration_mean_reported /
      reported_fraction_of_true_symptom_duration
    
    mean_sputumpos_symptomatic_duration_among_micropos <-
      true_symptom_duration_mean *
      (
        proportion_ever_subclinical +
          0.5 * (1 - proportion_ever_subclinical)
      )
    
    base_mean_duration_subclinical <-
      mean_sputumpos_symptomatic_duration_among_micropos *
      duration_ratio_subclinical_symptomatic /
      proportion_ever_subclinical
    
    cohort <- cohort %>%
      mutate(
        ever_subclinical_sputumpos = case_when(
          pulmonary_with_micro == 1 ~
            rbinom(
              prob = proportion_ever_subclinical,
              size = 1,
              n = n()
            ),
          TRUE ~ 0
        ))
    
    cohort <- cohort %>% mutate(sputumpos_at_eot = 0)
    
    eligible_eot <- which(cohort$ever_subclinical_sputumpos == 1)
    
    n_eot <- round(length(eligible_eot) * proportion_subclinical_sputumpos_at_eot)
    
    if (length(eligible_eot) > 0 && n_eot > 0) {
      eot_weights <- 1 / pmax(cohort$symptom_onset[eligible_eot], 1)
      
      eot_selected <- sample(
        eligible_eot,
        size = min(n_eot, length(eligible_eot)),
        prob = eot_weights,
        replace = FALSE
      )
      
      cohort$sputumpos_at_eot[eot_selected] <- 1
    }
    
    
     cohort <- cohort %>% mutate(
        raw_subclinical_duration = case_when(
          ever_subclinical_sputumpos == 1 &
            sputumpos_at_eot == 0 ~
            base_mean_duration_subclinical *
            rgamma(
              n = n(),
              shape = 1 / duration_subclinical_cv^2,
              scale = duration_subclinical_cv^2
            ),
          TRUE ~ NA_real_
        ),
        
        capped_subclinical_duration = case_when(
          ever_subclinical_sputumpos == 1 &
            sputumpos_at_eot == 0 ~
            pmin(
              raw_subclinical_duration,
              max_subclinical_fraction_of_presymptom_time * symptom_onset
            ),
          TRUE ~ NA_real_
        ),
        
        sputum_onset = case_when(
          ever_subclinical_sputumpos == 1 &
            sputumpos_at_eot == 1 ~
            0,
          
          ever_subclinical_sputumpos == 1 &
            sputumpos_at_eot == 0 ~
            symptom_onset - capped_subclinical_duration,
          
          ever_subclinical_sputumpos == 0 &
            pulmonary_with_micro == 1 ~
            runif(
              n = n(),
              min = symptom_onset,
              max = diagnosis_routine
            ),
          
          TRUE ~ NA_real_
        )
      )
     
     cohort <- cohort %>%
       mutate(
         diagnosis_routine_original = diagnosis_routine,
         TB_original = TB,
         symptom_onset_original = symptom_onset,
         sputum_onset_original = sputum_onset
       )
    
    #### Prediction ####
    
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
  if (subclinical_6mo_amongcohort > cohort_params$subclinical_6m_amongcohort_min & 
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
  colors <- c("Symptomatic" = "royalblue", "NAAT+" = "red")
  fillcolors <- c("Micro+ Pulmonary" = "red", "Clinical or extrapulmonary" = "gray")
  pointcolors <- c("Successful" = "green2")
  # linetypes <- c("symptoms" = "dotted", "micro" = "dashed", "both" = "dotdash")
  linetypes <- c("Not averted by screening" = "solid", "Averted by screening" = "11")
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
    ylab("Recurrences arranged by timing of routine diagnosis") + 
    xlab("Months since prior treatment completion") +
    scale_linetype_manual("Time with TB", values = linetypes) + 
    # if detection_timing is NA, plot full segment
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  diagnosis_routine/30, y = newID + 0.15, col='NAAT+',
                     lty = "Not averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = sputum_onset/30, xend =  detection_timing/30, y = newID + 0.15, col='NAAT+',
                     lty = "Not averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & sputum_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID + 0.15, col='NAAT+',
                     lty = "Averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  diagnosis_routine/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Not averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing)), 
                 aes(x = symptom_onset/30, xend =  detection_timing/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Not averted by screening"), size=0.8) + 
    geom_segment(data = plotdata %>% filter(!is.na(detection_timing) & symptom_onset < detection_timing), 
                 aes(x = detection_timing/30, xend = diagnosis_routine/30, y = newID - 0.15, col='Symptomatic',
                     lty = "Averted by screening"), size=0.8) + 
    geom_point(aes(x = detection_timing/30, y = newID), pch=10, col = "green2", size=2) +
    geom_point(aes(x = diagnosis_routine/30, y = newID), fill = 'black', pch=21) + 
    # geom_point(aes(x = detection_timing/30, y = newID), pch=10, col = "green2", size=2) +
    scale_x_continuous(breaks = seq(0, 30, 6), limits = c(0,30)) +
    theme_minimal() +
    # overlay legend on bottom right corner of plot
    theme(legend.position = c(0.8, 0.3), 
          legend.background = element_rect(fill = "transparent", color = NA),
          legend.box.background = element_rect(fill = "white", color = NA)) 
  
  if(colorfill) plot <- plot + 
    scale_fill_manual("Routine TB diagnosis", values = fillcolors) + 
    geom_point(aes(x = diagnosis_routine/30, y = newID, fill = TBtype), pch=21) 

  library(ggnewscale)
 plot <-  plot + new_scale_color() + 
    scale_color_manual("Detection through screening", values = pointcolors) + 
    geom_point(data = plotdata %>% filter(detected == 1), aes(x = detection_timing/30, y = newID, col = "Successful"), pch=10, size=2) 
  
  plot <- plot + guides(color = guide_legend(order = 1))
  
  
  return(plot)
  
}


#### Outcomes of interest #### 

# Cases detected (is this valuable in itself? this simulation assumes all with be detected eventually, or deaths aren't observed)
## Cases detected "early" by intervention
# Time with recurrent TB
## Time sputum+ (associated with transmission)
## Time symptomatic (associated with lung damage and with morality risk??)


time_with_tb <- function(cohort, limit_days = NULL)
{
  if (!is.null(limit_days)) cohort <- cohort %>% mutate(
    diagnosis_routine_original = pmin(diagnosis_routine_original, limit_days),
    diagnosis_routine = pmin(diagnosis_routine, limit_days),
    symptom_onset_original = pmin(symptom_onset_original, limit_days),
    symptom_onset = pmin(symptom_onset, limit_days),
    sputum_onset_original = pmin(sputum_onset_original, limit_days),
    sputum_onset = pmin(sputum_onset, limit_days),
    detection_timing = pmin(detection_timing, limit_days)
  )
  # add time 
  outcomes <- cohort %>% mutate(
    symptom_days_soc = diagnosis_routine_original - symptom_onset_original,
    sputum_days_soc = diagnosis_routine_original - sputum_onset_original,
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


# Costs
## Cost per case detected early or prevented (combine these? Though not equally valuable?)
## Cost per month with TB prevented 
### Overall
### Symptomatic
### Sputum+

costs <- function(screening_design, cohort, params, intervention_parameters)
{
  
  initial_contact_cost = list(home = params$initial_contact_cost_home, 
                              phone = params$initial_contact_cost_home/params$initial_contact_cost_home_vs_phone_factor)
  # treatment_cost = 1000
  naat_cost = params$sputum_test_cost
  symptom_prevalence_nontb = params$symptom_prevalence_nontb # in prevalence survey 27% of all people with past TB had a positive symptom screen. much lower in aftermath. will assume this is independent of TB risk (if currently TB negative).
  
  # count units of:
  # total people contacted
  # number of people with symptomatic TB, sputum+ asymptomatic TB, and non-TB symptoms
  # number of people with TB detected by screening
  
  # total people with initial contact:
  # contacted <- screening_design$target_coverage * 
  #                        params$coverage_phone * 
  #                        ifelse(screening_design$screening_location=="home", params$coverage_home_reduction, 1)
  contacted_individual <- # will create matrix of individuals (rows) x screening rounds (columns)
    # apply over screening rounds, 1:nrow(screening_design):
    mapply(function(loc, cov) sapply(cohort$risk_quantile, function(risk) 
      (risk >= (1-cov)) * (rbinom(n = 1, size = 1, prob = intervention_parameters$coverage[[loc]]))), 
      screening_design$screening_location, 
      screening_design$target_coverage)
    
  contact_cost <- rowSums(unlist(initial_contact_cost[screening_design$screening_location]) * t(contacted_individual))
  # contact_cost/N_cohort/screening_design$target_coverage # check
  
  sxtb <- sapply(screening_design$timing_months*30, function(s) cohort$symptom_onset < s & s <= cohort$diagnosis_routine)
  asxtb <- sapply(screening_design$timing_months*30, function(s) !is.na(cohort$sputum_onset) & cohort$sputum_onset < s & cohort$symptom_onset >= s) 
  nontbsx <- (!sxtb | is.na(sxtb)) & (!asxtb | is.na(asxtb)) * rbinom(n = nrow(cohort), size = 1, prob = symptom_prevalence_nontb)
  # colSums(sxtb[which(cohort$TB==1),])/sum(cohort$TB==1) # check, prev by time point
  # colSums(asxtb[which(cohort$TB==1),])/sum(cohort$TB==1) # check, prev by time point, ~half of sxtb bc sxtb included micro neg
  # mean(nontbsx)

  # sputum tests, 
  # symptom-driven sputum testing
  sputumprob_sx <- t(unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
                      unlist(screening_design$screening_method %in% c("symptoms", "both")) * t((sxtb|nontbsx) * contacted_individual )) 
  # colSums(sputumprob_sx[which(cohort$TB==1),])/sum(cohort$TB==1) # check
  # universal sputum testing (if people get multiple will reduce to 1)
    sputumprob_universal <- t(unlist(intervention_parameters$success_sputum[screening_design$screening_location]) * 
                      unlist(screening_design$screening_method %in% c("sputum", "both")) * t(contacted_individual))
  # sample each as binomial:
    sputum_sx <- matrix(
      rbinom(n = length(sputumprob_sx), size = 1, prob = sputumprob_sx),
      nrow = nrow(sputumprob_sx),
      ncol = ncol(sputumprob_sx)
    )
    sputum_universal <- matrix(
      rbinom(n = length(sputumprob_universal), size = 1, prob = sputumprob_universal),
      nrow = nrow(sputumprob_universal),
      ncol = ncol(sputumprob_universal)
    )
    sputum <- pmax(sputum_sx, sputum_universal, na.rm = T) # if doing both in same column, get max 1
    # if an element of screening_design$timing_months is the same as the element before, only count sputum tests in that column of sputum if the column before is 0 or NA
    if (nrow(screening_design) > 1) 
      for (i in 2:nrow(screening_design))
      {
        if (screening_design$timing_months[i] == screening_design$timing_months[i-1])
        {
          sputum[,i] <- sputum[,i] * (1 - sputum[,i-1])
        }
      }
    
  
  costs <- 
    sum(
    # home visits (including a visit if symptom+ on phone call?)
    # sputum tests
    contact_cost + colSums(sputum * params$sputum_test_cost, na.rm = T)
  # treatments (incremental?)
  # patient costs?
    )
  return(costs)
}

