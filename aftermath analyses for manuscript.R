library(binom)
library(clipr)

source("aftermath cohort sim functions.R")
N_cohort = 10000
N_samples <- 5000

#### Incorporate uncertainty: Sample cohort_params from posterior distributions ####
cohort_param_ranges <- list(
  N = N_cohort, # number of participants in each cohort sample
  
  recurrence_time_mean = c(158, 225), # diagnosis minus symptom duration, of those diagnosed by day 540
  # recurrence_time_sd = 179,
  recurrence_time_shapefactor = c(0.5, 1.5), # variation in (mean/sd)^2, compared to aftermath ratio where it was (213/179)^2
  incidence_18mo = c(0.074, 0.110),
  proportion_micro_pos = c(0.43, 0.63), #binom.agresti.coull(n = 88, x = 0.04392/0.0822*88)
  auc = c(0.56, 0.77),
  
  symptom_duration_mean_reported = 17,
  symptom_duration_sd_reported = 14*c(1, 1.5), 
  symptom_duration_timescale = 0, # increase this to 1/4 in a sensitivity analysis
  symptom_underestimation_factor = c(0.25, 0.5), 
  
  proportion_ever_subclinical = c(0.6, 0.9), # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2), 
  duration_subclinical_cv = c(0.5, 1.5),
  subclinical_baseline_amongTB_max = 0.20, 
  subclinical_6m_amongcohort_min = 0.005, 
  subclinical_6m_amongcohort_max = 0.016, 
  
  coverage_phone = c(0.9, 1.0), 
  coverage_home_reduction = c(0.75, 0.95), #0.81/0.95 = 0.85
  sensitivity_symptoms_home = c(0.85, 0.95),
  sensitivity_symptoms_phone_reduction = c(0.6, 0.82),
  success_sputum_home = c(0.85, 0.95),
  success_sputum_phone_reduction = c(0.6, 1.0),
  
  home_visit_passive_detection_impact = c(1), # except in sensitivity analysis
  aftermath_counseling_passive_detection_impact = c(0.9,1), 
  intentional_counseling_passive_detection_impact = c(0.8), 
  intentional_counseling_passive_detection_duration = 180, # days before effect of counseling on care seeking "wears off"
  prevention_efficacy = 0.6,
  
  # also sample costs (for allocating intervention coverage in cost-equivalent way)
  initial_contact_cost_home = c(1,5), 
  initial_contact_cost_phone_factor = c(1/14, 1/8), 
  sputum_test_cost = c(2,4),
  prevention_cost = c(7.24, 9.14)
)
  


cohort_params <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 
  
##### For each sample, simulate a cohort and describe key features: ####
 # month with highest incidence, month with highest undiagnosed prevalence (symptomatic and sputum+, and overall), and 
# remaining expected duration of prevalent disease at each of 3, 6, 9, 12, 15, and 18mo
run_cohort_features <- TRUE
if (run_cohort_features)
{  cohort_features <- 
    data.frame(
            accepted_subclinical = numeric(N_samples),
            month_with_highest_incidence = numeric(N_samples),
             month_with_highest_prevalence_sx = numeric(N_samples),
             month_with_highest_prevalence_inf = numeric(N_samples),
             month_with_highest_prevalence_total = numeric(N_samples),
             remaining_duration_3mo = numeric(N_samples),
             remaining_duration_6mo = numeric(N_samples),
             remaining_duration_9mo = numeric(N_samples),
             remaining_duration_12mo = numeric(N_samples),
             remaining_duration_15mo = numeric(N_samples),
             remaining_duration_18mo = numeric(N_samples),
            cum_notifs_6mo = numeric(N_samples),
            cum_notifs_18mo = numeric(N_samples),
            cum_notifs_overall = numeric(N_samples),
            median_time_to_onset = numeric(N_samples),
            median_time_to_diagnosis = numeric(N_samples),
            mean_symptom_duration = numeric(N_samples),
            sd_symptom_duration = numeric(N_samples),
            proportion_micro_pos = numeric(N_samples),
            proportion_with_subclinical = numeric(N_samples),
            mean_subclinical_duration = numeric(N_samples),
            sd_subclinical_duration = numeric(N_samples),
            prev_sx_6mo = numeric(N_samples),
            prev_inf_6mo = numeric(N_samples),
            prev_sx_12mo = numeric(N_samples),
            prev_inf_12mo = numeric(N_samples),
            cum_months_sx = numeric(N_samples),
            cum_months_inf = numeric(N_samples))
  
    # For visual, for guidelines intervention, quantify for a cascade:
    cascade_features <- 
      data.frame(
        accepted_subclinical = numeric(N_samples),
        cumulative_incidence = numeric(N_samples),
        TB_beyond_6mo = numeric(N_samples),
        TB_beyond6_before18 = numeric(N_samples),
        TB_at_visit = numeric(N_samples),
        symptomatic_at_visit = numeric(N_samples),
        detected_and_linked = numeric(N_samples),
        total_time_of_linked = numeric(N_samples),
        remaining_time_of_linked = numeric(N_samples)
      )
  

  for (n in 1:N_samples)
  {
    # setup
    cohort <- create_cohort(cohort_params[n,])
    
    cohort_features[n, "accepted_subclinical"] <- 
      cascade_features[n, "accepted_subclinical"] <-
        check_subclinical(cohort, cohort_params[n,])
    
    cohort_features[n, "month_with_highest_incidence"] <- 
      (cohort %>% filter(TB==1) %>% 
         mutate(month = floor(symptom_onset/30)) %>% group_by(month) %>% summarise(n = n()) %>% 
         arrange(desc(n)) %>% select(month)) %>%  .[1,1]
    
    cohort_features[n, "month_with_highest_prevalence_sx"] <-
      ((1:20))[which.max(sapply(30 * (1:20), function(x) sum(cohort$symptom_onset < x & cohort$diagnosis_routine > x, na.rm=T)))]
  
    cohort_features[n, "month_with_highest_prevalence_inf"] <-
      ((1:20))[which.max(sapply(30 * (1:20), function(x) sum(cohort$sputum_onset < x & cohort$diagnosis_routine > x, na.rm=T)))]
  
    cohort_features[n, "month_with_highest_prevalence_total"] <-
      ((1:20))[which.max(sapply(30 * (1:20), function(x) 
        sum((cohort$sputum_onset < x | cohort$symptom_onset < x) & cohort$diagnosis_routine > x, na.rm=T)))]
    
    cohort_features[n, "remaining_duration_3mo"] <-
      cohort %>% filter(diagnosis_routine > 3*30 & (sputum_onset < 3*30 | symptom_onset < 3*30)) %>% 
        summarise(mean(diagnosis_routine - 3*30, na.rm = TRUE))
    cohort_features[n, "remaining_duration_6mo"] <-
      cohort %>% filter(diagnosis_routine > 6*30 & (sputum_onset < 6*30 | symptom_onset < 6*30)) %>% 
        summarise(mean(diagnosis_routine - 6*30, na.rm = TRUE))
    cohort_features[n, "remaining_duration_9mo"] <-
      cohort %>% filter(diagnosis_routine > 9*30 & (sputum_onset < 9*30 | symptom_onset < 9*30)) %>% 
        summarise(mean(diagnosis_routine - 9*30, na.rm = TRUE))
    cohort_features[n, "remaining_duration_12mo"] <-
      cohort %>% filter(diagnosis_routine > 12*30 & (sputum_onset < 12*30 | symptom_onset < 12*30)) %>% 
        summarise(mean(diagnosis_routine - 12*30, na.rm = TRUE))
    cohort_features[n, "remaining_duration_15mo"] <-
      cohort %>% filter(diagnosis_routine > 15*30 & (sputum_onset < 15*30 | symptom_onset < 15*30)) %>% 
        summarise(mean(diagnosis_routine - 15*30, na.rm = TRUE))
    cohort_features[n, "remaining_duration_18mo"] <-
      cohort %>% filter(diagnosis_routine > 18*30 & (sputum_onset < 18*30 | symptom_onset < 18*30)) %>% 
        summarise(mean(diagnosis_routine - 18*30, na.rm = TRUE))
    
    cohort_features[n, "cum_notifs_6mo"] <- 
      sum(cohort$TB==1 & cohort$diagnosis_routine <= 6*30, na.rm=T)/N_cohort
    cohort_features[n, "cum_notifs_18mo"] <-
      sum(cohort$TB==1 & cohort$diagnosis_routine <= 18*30, na.rm=T)/N_cohort
    cohort_features[n, "cum_notifs_overall"] <-
      sum(cohort$TB==1, na.rm=T)/N_cohort
    cohort_features[n, "median_time_to_onset"] <-
      median(cohort$symptom_onset, na.rm=T)
    cohort_features[n, "median_time_to_diagnosis"] <-
      median(cohort$diagnosis_routine, na.rm=T)
    cohort_features[n, "mean_symptom_duration"] <-
      mean(cohort$diagnosis_routine - cohort$sputum_onset, na.rm=T)
    cohort_features[n, "sd_symptom_duration"] <-
      sd(cohort$diagnosis_routine - cohort$sputum_onset, na.rm=T)
    cohort_features[n, "proportion_micro_pos"] <-
      sum(cohort$pulmonary_with_micro==1, na.rm=T)/sum(cohort$TB==1, na.rm=T)
    cohort_features[n, "proportion_with_subclinical"] <-
      sum(cohort$sputum_onset < cohort$symptom_onset, na.rm=T)/sum(cohort$pulmonary_with_micro==1, na.rm=T)
    cohort_features[n, "mean_subclinical_duration"] <-
      mean(cohort$symptom_onset[cohort$sputum_onset < cohort$symptom_onset] - 
             cohort$sputum_onset[cohort$sputum_onset < cohort$symptom_onset], na.rm=T)
    cohort_features[n, "sd_subclinical_duration"] <-
      sd(cohort$symptom_onset[cohort$sputum_onset < cohort$symptom_onset] - 
           cohort$sputum_onset[cohort$sputum_onset < cohort$symptom_onset], na.rm=T)
    cohort_features[n, "prev_sx_6mo"] <-
      sum(cohort$symptom_onset < 6*30 & cohort$diagnosis_routine >= 6*30, na.rm=T)/N_cohort
    cohort_features[n, "prev_inf_6mo"] <-
      sum(cohort$sputum_onset < 6*30 & cohort$diagnosis_routine >= 6*30, na.rm=T)/N_cohort
    cohort_features[n, "prev_sx_12mo"] <-
      sum(cohort$symptom_onset < 12*30 & cohort$diagnosis_routine >= 12*30, na.rm=T)/N_cohort
    cohort_features[n, "prev_inf_12mo"] <-
      sum(cohort$sputum_onset < 12*30 & cohort$diagnosis_routine >= 12*30, na.rm=T)/N_cohort
    cohort_features[n, "cum_months_sx"] <-
      sum(cohort$diagnosis_routine - cohort$symptom_onset, na.rm=T)/30
    cohort_features[n, "cum_months_inf"] <-
      sum(cohort$diagnosis_routine - cohort$sputum_onset, na.rm=T)/30
      
    
    cascade_features[n, "cumulative_incidence"] <- 
      sum(cohort$TB==1)
    cascade_features[n, "TB_beyond_6mo"] <-
      sum(cohort$diagnosis_routine > 6*30 & cohort$TB==1)
    cascade_features[n, "TB_beyond6_before18"] <-
      sum(cohort$diagnosis_routine > 6*30 & (cohort$symptom_onset < 18*30 | cohort$sputum_onset < 18*30) & cohort$TB==1, na.rm=T)
    cascade_features[n, "TB_at_visit"] <-
      sum(cohort$TB==1 & (((cohort$symptom_onset < 6*30 | cohort$sputum_onset < 6*30) & cohort$diagnosis_routine > 6*30)  | 
                          ((cohort$symptom_onset < 12*30 | cohort$sputum_onset < 12*30) & cohort$diagnosis_routine > 12*30) |
                          ((cohort$symptom_onset < 18*30 | cohort$sputum_onset < 18*30) & cohort$diagnosis_routine > 18*30)), na.rm=T)
    cascade_features[n, "symptomatic_at_visit"] <-
      sum(cohort$TB==1 & ((cohort$symptom_onset < 6*30 & cohort$diagnosis_routine > 6*30)  | 
                          (cohort$symptom_onset < 12*30 & cohort$diagnosis_routine > 12*30) |
                          (cohort$symptom_onset < 18*30 & cohort$diagnosis_routine > 18*30)), na.rm=T)
    
    
    screening_design_guidelines <- 
      data.frame(
        "counseling_coverage" = 0, 
        "prevention_coverage" = 0,
        "timing_months" = c(6, 12, 18),
        "target_coverage" = c(1, 1, 1),
        "screening_method" = c("symptoms", "symptoms", "symptoms"),
        "screening_location" = c("home", "phone", "phone")
      ) %>% 
      arrange(timing_months)
    
    intervention_parameters <- list(
      coverage = list("phone" = cohort_params$coverage_phone[n], 
                      "home" = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]),
      sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home[n], 
                                  "phone" = cohort_params$sensitivity_symptoms_home[n] * cohort_params$sensitivity_symptoms_phone_reduction[n]),
      success_sputum = list("home" = cohort_params$success_sputum_home[n], 
                            "phone" = cohort_params$success_sputum_home[n] * cohort_params$success_sputum_phone_reduction[n]))
    
    outputs <- apply_intervention(cohort, screening_design_guidelines,
                                                  intervention_parameters, cohort_params[n,]) 
  
    cascade_features[n, "detected_and_linked"] <-
      sum(outputs$TB==1 & !is.na(outputs$detection_timing))
    cascade_features[n, "total_time_of_linked"] <-
      sum((outputs$diagnosis_routine - outputs$symptom_onset)[outputs$TB==1 & !is.na(outputs$detection_timing)])
    cascade_features[n, "remaining_time_of_linked"] <-
      sum((outputs$detection_timing - outputs$symptom_onset)[outputs$TB==1 & !is.na(outputs$detection_timing)])
    
  }        
}
  
#### Evaluate cohort characteristics and natural history results ####
head(cohort_features)
cohort_features[,] %>% summarise_all(median)
cohort_features[,] %>% summarise_all(mean)

sum(cohort_features$accepted_subclinical) # proportion of cohort that accepted subclinical TB screening
mean(cohort_features$accepted_subclinical) # proportion of cohort that accepted subclinical TB screening

cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) median(x, na.rm=T))
cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.025, na.rm=T))
cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.975, na.rm=T)) 
cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.25, na.rm=T))
cohort_features[cohort_features$accepted_subclinical==1,] %>% summarise_all(function(x) quantile(x, 0.75, na.rm=T)) # for the abstract

  
  # For visual, for guidelines intervention, quantify for a cascade:
  #   * cumulative incidence, 
  #   * proportion with onset after 6 months, 
  #   * proportion of those with onset before 18 months, 
  #   * proportion of those with TB at time of avisit (6, 12, or 18), 
  #   * proportion of those symptomatic at visit, 
  #   * proportion of those detected and linked, 
  #  * proportion of those people's' disease courses (mean) averted through detection
  cascade_features[,] %>% filter(accepted_subclinical==1) %>% summarise_all(median)  
  # organize data for a cascade plot, with each step as a % of total and as a % of previous step
  # * all TB recurrence 
  # * proportion with TB before 6mo, within 6-18mo, and after 18M (3 bars, add up to first bar, spaced accordingly)
  # within 6-18mo bar, proportions with asymptomatic TB at a study visit, sx at a study visit, and all bwteen visits
  # of those with TB during a study visit, proportion detected and linked vs not
  # of those detected and linked, proportion of their time with TB alredy occurred vs averted.
  
  cascade_features$n <- 1:N_samples
  # Make a cascade bar graph, where there are bars from:
  plotdata <- 
    cascade_features %>%
    filter(accepted_subclinical==1) %>%
      mutate(
         recur_start = 0,
         recur_end = cumulative_incidence,
         recur_label = "Develops recurrent TB",
         recur_proportion = cumulative_incidence/N_cohort,
         recur_group = "All recurrences",
         
         before_6mo_start = 0,
         before_6mo_end = cumulative_incidence - TB_beyond_6mo,
         before_6mo_label = "Occurs too early (<6mo)",
         before_6mo_proportion = (cumulative_incidence - TB_beyond_6mo)/cumulative_incidence, 
         before_6mo_group = "Timing of recurrence",
         
         between618mo_start = cumulative_incidence - TB_beyond_6mo,
         between618mo_end = cumulative_incidence - (TB_beyond_6mo - TB_beyond6_before18),
         between618mo_label = "Active between 6-18 months",
         between618mo_proportion = (TB_beyond6_before18)/cumulative_incidence,
         between618mo_group = "Timing of recurrence",

         after18mo_start = cumulative_incidence - (TB_beyond_6mo - TB_beyond6_before18),
         after18mo_end = cumulative_incidence,
         after18mo_label = "Occurs too late (>18mo)",
         after18mo_proportion = (TB_beyond_6mo - TB_beyond6_before18)/cumulative_incidence,
         after18mo_group = "Timing of recurrence",

         asxatvisit_start = between618mo_start,
         asxatvisit_end = between618mo_start + (TB_at_visit - symptomatic_at_visit),
         asxatvisit_label = "Asymptomatic at visit",
         asxatvisit_proportion = (TB_at_visit - symptomatic_at_visit)/(TB_beyond6_before18),
         asxatvisit_group = "Status at time of visit",
         
         sxatvisit_start = between618mo_start + (TB_at_visit - symptomatic_at_visit),
         sxatvisit_end = between618mo_start + TB_at_visit,
         sxatvisit_label = "Symptomatic at visit",
         sxatvisit_proportion = symptomatic_at_visit/(TB_beyond6_before18),
         sxatvisit_group = "Status at time of visit",

         notbatvisit_start = between618mo_start + TB_at_visit,
         notbatvisit_end = between618mo_end,
         notbatvisit_label = "No TB at visit",
         notbatvisit_proportion = (TB_beyond6_before18 - TB_at_visit)/(TB_beyond6_before18),
         notbatvisit_group = "Status at time of visit",
         
         undetected_start = sxatvisit_start,
         undetected_end = sxatvisit_start + (symptomatic_at_visit - detected_and_linked),
         undetected_label = "Missed",
         undetected_proportion = (symptomatic_at_visit - detected_and_linked)/(symptomatic_at_visit),
         undetected_group = "Outcome of visit",
         
         detectedlinked_start = sxatvisit_start + (symptomatic_at_visit - detected_and_linked),
         detectedlinked_end = sxatvisit_start + symptomatic_at_visit,
         detectedlinked_label = "Found early",
         detectedlinked_proportion = detected_and_linked/(symptomatic_at_visit),
         detectedlinked_group = "Outcome of visit",
         
         unaverted_start = detectedlinked_start,
         unaverted_end = detectedlinked_start + detected_and_linked*(total_time_of_linked - remaining_time_of_linked)/total_time_of_linked,
         unaverted_label = "Already elapsed",
         unaverted_proportion = (total_time_of_linked - remaining_time_of_linked)/(total_time_of_linked),
         unaverted_group = "Time with TB",
         
         averted_start = detectedlinked_start + detected_and_linked*(total_time_of_linked - remaining_time_of_linked)/total_time_of_linked,
         averted_end = detectedlinked_start + detected_and_linked,
         averted_label = "Averted",
         averted_proportion = remaining_time_of_linked/(total_time_of_linked),
         averted_group = "Time with TB")
  
  # lengthen data, using portion before last "_" as one variable name and portion after as another
    plotdata <- plotdata %>% pivot_longer(
      cols = c(recur_start:averted_group), 
      names_to = c("outcome", ".value"), 
      names_pattern = "(.*)_(.*)"
    ) %>% select(n, outcome, start, end, label, proportion, group,
                 cumulative_incidence, TB_beyond6_before18, TB_at_visit, symptomatic_at_visit, detected_and_linked)
         
  # bar chart, with x axis label and separate color for each group
    # We'll plot averages, with error bars at the top for uncertainty. 
    # This means we'll need to summarize the data first, then plot.
    # And within each group, we'll have the uncertainty for each bar reflect uncertainty in the proportion with a given outcome, 
    # as a proportion of the relevant grup denominator (N_cohort, cumulative_incidence, TB_beyond6_before18, TB_at_visit, symptomatic_at_visit, and detected_linked)
    
    # keep outcomes in order 
  plotdata$outcome <- factor(plotdata$outcome, 
                             levels = c("recur", "before_6mo", "between618mo", "after18mo", 
                                        "asxatvisit", "sxatvisit", "notbatvisit",
                                        "undetected", "detectedlinked", 
                                        "unaverted", "averted"))
  
  plotdata$group <- factor(plotdata$group, 
                             levels = c("All recurrences", 
                                        "Timing of recurrence",
                                        "Status at time of visit",
                                        "Outcome of visit",
                                        "Time with TB"))
  
  plotdata$outcomex <- as.numeric(plotdata$group) + as.numeric(plotdata$outcome)
  labelmapping <- plotdata %>% count(outcomex, label) %>% select(outcomex, label)
  pathmapping <- cbind(c("recur", "before_6mo", "between618mo", "after18mo", 
                         "asxatvisit", "sxatvisit", "notbatvisit",
                         "undetected", "detectedlinked", 
                         "unaverted", "averted"), 
                       c(T, F, T, F, 
                           F, T, F,
                           F, T,
                           F, T))
  groups <- c(1, 2, 2, 2, 
              3, 3, 3, 4, 4, 5, 5)
  items <- (1:length(groups))-1
  # in plotdata, assign "path" by looking up outcome in pathmapping
  plotdata$path <- as.logical(pathmapping[plotdata$outcome, 2])
  
  axistext <- data.frame(text = c("Occurrence", "Timing", "Status",  "Detection", "Time with TB"),
                      position = c(2, 5, 9, 12.5, 15.5))
    
  
  plotdata <- plotdata %>% mutate(denominator = case_when(
    group == "All recurrences" ~ N_cohort,
    group == "Timing of recurrence" ~ cumulative_incidence,
    group == "Status at time of visit" ~ TB_beyond6_before18,
    group == "Outcome of visit" ~ symptomatic_at_visit,
    group == "Time with TB" ~ detected_and_linked 
  ))
  
  
  
  
  ####### cascade figure, for guidelines intervention ####
  ggplot(
      # summarize means for each outcome
      plotdata %>% group_by(outcomex, group) %>% 
        summarise(start = mean(start), end = mean(end), label = first(label), 
                  proportion = mean(proportion), cumulative_incidence = first(cumulative_incidence),
                  TB_beyond6_before18 = mean(TB_beyond6_before18),
                  TB_at_visit = mean(TB_at_visit),
                  symptomatic_at_visit = mean(symptomatic_at_visit),
                  detected_and_linked = mean(detected_and_linked),
                  path = ifelse(first(path), "Yes", "No"),
                  text = paste0(round(100*proportion,0),"%")),
      aes(x = outcomex, ymin = start/N_cohort, ymax = end/N_cohort, col = group)) + 
    geom_linerange(size = 10, aes(alpha = factor(path))) + 
    scale_alpha_discrete(range = c(0.3, 0.9), guide = guide_legend(override.aes = list(fill = "black"))) +
    geom_label(aes(x = outcomex, y = (start + (end-start)/2)/N_cohort, label=text)) +
    # error on top of bars, within a group, is based on quantiles of that bar size as a proportion of denoninator,
    # so calculate the quantiles of the proportion and then multiply by mean size and add to mean start.
    geom_errorbar(data = plotdata %>% group_by(outcomex, group) %>%
                    mutate(endmin = mean(start)/N_cohort + mean(denominator)/N_cohort*quantile((end-start)/denominator, 0.025),
                          endmax = mean(start)/N_cohort + mean(denominator)/N_cohort*quantile((end-start)/denominator, 0.975)), 
                  aes(ymin = endmin, ymax = endmax), col = "black", alpha=0.04, width = 0.3) +
    
    theme_minimal() +
  # y axis as percentages, and labeled "Proportion of cohort"
    ylab("Proportion of treatment-completing cohort") + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    xlab("Step in recurrent-TB case-finding cascade") +
    # rotate x axis lablels, and replace them with variable "label"
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.title = element_text(face="bold")) + 
    scale_x_continuous(breaks=labelmapping$outcomex,
                         labels=labelmapping$label) +
    guides(col=guide_legend(title="Portion of cascade"),
           alpha = guide_legend(title="Leads to detection\nthrough screening?", reverse = T)) +
    # add extra text in bold along x axis:
    annotate("text", x = axistext$position, y = -0.003, label = axistext$text, fontface=2)
    
    
    



#### More intervnetion setup ####

intervention_names = c("guidelines", "earlier", "frequent", "sputum", "sputummoretargeted", "counseling", "prevention")

results <- list()
for (name in intervention_names) {
  results[[name]] <- 
    data.frame(recurrences = numeric(N_samples),
               coverage = numeric(N_samples), # of targeted portion of intervention
               detections = numeric(N_samples), # by intervention
               mean_symptom_days = numeric(N_samples),
               symptomatic_months_soc = numeric(N_samples),
               infectious_months_soc = numeric(N_samples),
               symptomatic_months_averted = numeric(N_samples), 
               infectious_months_averted = numeric(N_samples))

}

screening_design_guidelines <- 
  data.frame(
    "counseling_coverage" = 0, 
    "prevention_coverage" = 0,
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## earlier screening
screening_design_earlier <- 
  data.frame(
    "counseling_coverage" = 0, 
    "prevention_coverage" = 0,
    "timing_months" = c(3, 6, 12),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

# counseling, with passive detection impact (assumed same cost as 12 +18 month visits)
screening_design_counseling <- 
  data.frame(
    "counseling_coverage" = 1,
    "prevention_coverage" = 0,
    "timing_months" = c(6),
    "target_coverage" = c(1),
    "screening_method" = c("symptoms"),
    "screening_location" = c("home")
  ) %>% 
  arrange(timing_months)


# For the rest, we'll define within the loop because target coverage depends on cost

#### For each sample, simulate a cohort and estimate the impact of the aftermath intervention and alternatives: ####
# Running as loop, could parallelize for efficiency ***

for (n in 1:N_samples)
{
  print(n)
  #### More parameter-specific setup ####
  cohort <- create_cohort(cohort_params[n,])
  # if doesn't validate on subclinical TB, abort and fill NAs
  if (!check_subclinical(cohort, cohort_params[n,])) {
    # change the nth element of each item in results to NA
    for (name in intervention_names) 
      results[[name]][n,] <- rep(NA, ncol(results[[name]]))
    next
  }
  
  intervention_parameters <- list(
    coverage = list("phone" = cohort_params$coverage_phone[n], 
                    "home" = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]),
    sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home[n], 
                                "phone" = cohort_params$sensitivity_symptoms_home[n] * cohort_params$sensitivity_symptoms_phone_reduction[n]),
    success_sputum = list("home" = cohort_params$success_sputum_home[n], 
                          "phone" = cohort_params$success_sputum_home[n] * cohort_params$success_sputum_phone_reduction[n]))
  
  #### Define additional cost-equivalent intervention options ####

  ## more frequent but targeted. 
  targeting_frequent = (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) / # cost of aftermath version
                  (1 + 5*cohort_params$initial_contact_cost_phone_factor[n])
  screening_design_frequent <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(3, 6, 9, 12, 15, 18),
      "target_coverage" = rep(targeting_frequent, 6),
      "screening_method" = rep("symptoms", 6),
      "screening_location" = c("home", rep("phone",5))
    ) %>% 
    arrange(timing_months)
  
  ## sputum at 6mo for highest risk, in place of 12 and 18mo visits for all (the untargeted still get a 6mo visit)
  # (should modify to include confirmatory testing or false positives?)
  # (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) = #cost of aftermath version
  #                 (cohort_params$initial_contact_cost_home[n] * (1 + targeting*2*cohort_params$initial_contact_cost_phone_factor[n]) +
  #                    targeting*cohort_params$sputum_test_cost[n])
  targeting_sputum =  (cohort_params$initial_contact_cost_home[n] * (2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
                  (cohort_params$initial_contact_cost_home[n] * 2*cohort_params$initial_contact_cost_phone_factor[n] +
                     cohort_params$sputum_test_cost[n])
    screening_design_sputum <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(6, 6, 12, 18),
      "target_coverage" = c(targeting_sputum, 1, targeting_sputum, targeting_sputum),
      "screening_method" = c("sputum", "symptoms", "symptoms", "symptoms"),
      "screening_location" = c("home", "home", "phone", "phone")
    ) %>% 
    arrange(timing_months)
  
  # aftermath for a fraction, plus sputum at first visit (and the others get nothing)
    targeting_sputummoretargeted =  (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
      (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) +
         cohort_params$sputum_test_cost[n])
  screening_design_sputummoretargeted <- 
    data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = 0,
      "timing_months" = c(6, 12, 18),
      "target_coverage" = c(1,1,1) * targeting_sputummoretargeted,
      "screening_method" = c("both", "symptoms", "symptoms"),
      "screening_location" = c("home", "phone", "phone")
    ) %>% 
    arrange(timing_months)


  # prevention of recurrence
  targeting_prevention = 
    (cohort_params$initial_contact_cost_home[n] * (1 + 2*cohort_params$initial_contact_cost_phone_factor[n]) ) / 
    (cohort_params$prevention_cost[n])
  screening_design_prevention <- data.frame(
      "counseling_coverage" = 0,
      "prevention_coverage" = targeting_prevention,
      "timing_months" = c(6),
      "target_coverage" = c(0),
      "screening_method" = c("symptoms"),
      "screening_location" = c("home")
    ) %>% 
    arrange(timing_months)
  
  

  ##### Run the interventions ####
  outputs <- list()
  for (intervention in intervention_names)
  {
    coverage[[intervention]] <- ifelse(intervention %in% c("frequent", "sputum", "sputummoretargeted", "prevention"), get(paste0("targeting_", intervention)), NA) 
    outputs[[intervention]] <- apply_intervention(cohort, get(paste0("screening_design_", intervention)), 
                                                  intervention_parameters, cohort_params[n,])  
  }
  
  # impact
  (recurred = lapply(outputs, function(x) sum(x$TB)))
  
  (detected <- lapply(outputs, function(x) sum(!is.na(x$detection_timing))))
  (symptomdays <- lapply(outputs, function(x) ((time_with_tb(x) %>% filter(scenario == "screening") %>% select(value))/sum(x$TB))[1,]))
  (time <- lapply(outputs, function(x) time_with_tb(x) %>% filter(scenario == "soc") %>% 
                    mutate(months = value/30) %>%
                    tibble::column_to_rownames('outcome') %>% select(months)))
  (impact <- lapply(outputs, function(x) (time_with_tb(x)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
                                            mutate(months_averted = (soc - screening)/30)) %>% 
                      tibble::column_to_rownames('outcome')  %>% select(months_averted)))
  
  
  
  for (name in intervention_names)
    results[[name]] [n,] <- c(recurred[[name]], coverage[[name]], 
                              detected[[name]], symptomdays[[name]], time[[name]]$months[1], time[[name]]$months[2], 
                              impact[[name]]$months_averted[1], impact[[name]]$months_averted[2])
    
  # cost,
  # cost/impact$months_averted[1], cost/impact$months_averted[2])
}

saveRDS(results, "results_aftermath__20250716.RDS")
results <- readRDS("results_aftermath__20250716.RDS")

#### Look at results ####
# collate results for a table:
# for all dataframe elements of list "results", report the mean and 25th and 75th percentiles of each column

cbind(
  lapply(results, function(x) 
  # mean of each columns of table
  x %>% 
    mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
           inf_reduction = infectious_months_averted/infectious_months_soc*100,
           coverage_percent = coverage * 100,
           proportion_detected_percent = detections/recurrences*100) %>%
    select(coverage_percent, proportion_detected_percent, sx_reduction, inf_reduction) %>%
  summarise_all(  function(y) paste0(round(median(y, na.rm = T),0), "% (", 
                     round(quantile(y, 0.025, na.rm = T),0), "-", 
                     round(quantile(y, 0.975, na.rm = T),0), "%)")))  %>% bind_rows() %>% mutate(intervention = intervention_names),
  lapply(results, function(x)
  # mean of each columns of table
  x %>% 
    mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
           inf_reduction = infectious_months_averted/infectious_months_soc*100,
           coverage_percent = coverage * 100,
           proportion_detected_percent = detections/recurrences*100) %>%
    select(recurrences, detections,  mean_symptom_days) %>%
    summarise_all(  
      function(y) paste0(round(median(y, na.rm = T),0), " (", 
                         round(quantile(y, 0.025, na.rm = T),0), ", ", 
                         round(quantile(y, 0.975, na.rm = T),0), ")"))) %>% bind_rows() %>% mutate(intervention = intervention_names) %>%
    select(-intervention)
  ) %>%
  # make intervention the first column, and mean_symptom_days the second column
  select(intervention, coverage_percent, recurrences, detections, proportion_detected_percent, sx_reduction, inf_reduction) %>%
  filter(intervention != "sputum") %>% 
  write_clip() 





# old code::
# summarize results for aftermath intervention
lapply(results, function(x) summary(x$cost_per_symptomatic_month_averted))
summary(results_aftermath$cost_per_infectious_months_averted)

summary(results_earlier$cost_per_symptomatic_month_averted)
summary(results_earlier$cost_per_infectious_months_averted)

summary(results_targeted$cost_per_symptomatic_month_averted)
summary(results_targeted$cost_per_infectious_months_averted)

summary(results_sputum$cost_per_symptomatic_month_averted)
summary(results_sputum$cost_per_infectious_months_averted)

summary(results_moretargeted$cost_per_symptomatic_month_averted)
summary(results_moretargeted$cost_per_infectious_months_averted)

summary(results_home_improves_passive$cost_per_symptomatic_month_averted)
summary(results_home_improves_passive$cost_per_infectious_months_averted)


# sensitivty analysis
# for each parameter in cohort_param_ranges, compare results$cost_per_symptomatic_case_averted 
# between the top 10% and bottom 10% of values for that paramter
parameters_varied <- names(cohort_param_ranges)[lapply(cohort_param_ranges, length)>1]
oneways <- array(NA, dim=c(2, length(parameters_varied)))
dimnames(oneways) = list(c("top", "bottom"), parameters_varied)
for (p in parameters_varied) {
  oneways["top", p] <- mean(results_aftermath[cohort_params[p] >= quantile(unlist(cohort_params[p]), 0.8),]$cost_per_symptomatic_month_averted)
  oneways["bottom", p] <- mean(results_aftermath[cohort_params[p] <= quantile(unlist(cohort_params[p]), 0.2),]$cost_per_symptomatic_month_averted)
}
# plot tornado diagram, centered at mean of results$cost_per_symptomatic_month_averted
oneways %>%
  as.data.frame() %>%
  rownames_to_column("top_or_bottom") %>%
  pivot_longer(cols = -top_or_bottom, names_to = "parameter", values_to = "cost_per_symptomatic_month_averted") %>%
  ggplot(aes(x = parameter, 
             col = top_or_bottom)) +
  geom_errorbar(aes(ymin = mean(results$cost_per_symptomatic_month_averted), ymax=cost_per_symptomatic_month_averted)) +
  coord_flip() +
  geom_hline(yintercept = mean(results$cost_per_symptomatic_month_averted), linetype = "dashed") + 
  # specify text labels in color legend
  labs(color = "Quantile of parameter values") + 
  scale_color_manual(labels = c("Lowest 