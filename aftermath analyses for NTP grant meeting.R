library(binom)
library(clipr)
library(cowplot)
library(patchwork)

source("aftermath cohort simulation functions.R")

#### incorporate uncertainty ####
# Sample cohort_params from posterior distributions
cohort_param_ranges <- list(
  N = 100000, 
  symptom_duration_mean_reported = 17,
  symptom_duration_sd_reported = 14*c(1, 2), 
  symptom_duration_timescale = 0, # increase this to 1/4 in a sensitivity analysis
  symptom_underestimation_factor = c(0.33, 0.67), 
  recurrence_time_mean = c((233-17*2) - 17*2, (233-17*2) + 17*2), # diagnosis minus symptom duration
  # recurrence_time_sd = 151,
  recurrence_time_shapefactor = c(0.5, 1.5), # variation in (mean/sd)^2, compared to aftermath ratio where it was (199/151)^2
  incidence_18mo = c(0.065, 0.100),
  proportion_micro_pos = c(0.42, 0.64), #binom.exact(n = 90, x = 0.04275/0.08055*90)
  proportion_ever_subclinical = c(0.6, 0.9), # proportion of sputum+ TB that is sputum+ before symptom screen+
  # (of those that will be sputum+ when routinely diagnosed based on symptoms; the rest become sputum+ somewhere between symptom onset and routine diagnosis)
  duration_ratio_subclinical_symptomatic = c(0.8, 1.2), 
  duration_subclinical_cv = c(0.5, 1.5),
  auc = c(0.6, 0.8),

  coverage_phone = c(0.6,0.8),
  coverage_home_reduction = c(0.75, 0.95), #0.81/0.95 = 0.85
  sensitivity_symptoms_home = c(0.7, 0.9),
  sensitivity_symptoms_phone_reduction = c(0.5, 0.75),
  success_sputum_home = c(0.85, 0.95),
  success_sputum_phone_reduction = c(0.7, 0.9),
  
  home_visit_passive_detection_impact = c(1), 
  aftermath_counseling_passive_detection_impact = c(0.9,1), 
  intentional_counseling_passive_detection_impact = c(0.8), 
  
  
  # also sample costs??
  initial_contact_cost_home = c(10, 20), 
  initial_contact_cost_phone = c(3,7), 
  sputum_test_cost = c(14,18)
  
)


## Instead of home visit improving passive, we want thebaseline counseling to improve passive. Make counseling variable a part of screening intervention parameters. And then incorporate into simluation

  
N_samples <- 1000

cohort_params <- lapply(cohort_param_ranges, function(x) if(length(x)==2) runif(N_samples, min = x[1], max = x[2]) else rep(x, N_samples)) %>%
  as.data.frame() 
  
# For each sample, simulate a cohort and describe key features: 
 # month with highest incidence, month with highest undiagnosed prevalence (symptomatic and sputum+, and overall), and 
# remaining expected duration of prevalent disease at each of 3, 6, 9, 12, 15, and 18mo
cohort_features <- 
  data.frame(month_with_highest_incidence = numeric(N_samples),
             month_with_highest_prevalence_sx = numeric(N_samples),
             month_with_highest_prevalence_inf = numeric(N_samples),
             month_with_highest_prevalence_total = numeric(N_samples),
             remaining_duration_3mo = numeric(N_samples),
             remaining_duration_6mo = numeric(N_samples),
             remaining_duration_9mo = numeric(N_samples),
             remaining_duration_12mo = numeric(N_samples),
             remaining_duration_15mo = numeric(N_samples),
             remaining_duration_18mo = numeric(N_samples))

for (n in 1:N_samples)
{
  # setup
  cohort <- create_cohort(cohort_params[n,])
  
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
}

cohort_features[,] %>% summarise_all(function(x) median(x, na.rm = T))
cohort_features[,] %>% summarise_all(function(x) quantile(x, 0.25, na.rm = T))
cohort_features[,] %>% summarise_all(function(x) quantile(x, 0.75, na.rm = T)) # for the abstract
    
# For each sample, simulate a cohort and estimate the impact of the aftermath home intervention:

# for all, will assume impact on subsequent passive diagnosis.
# but will also run aftermath with no such impact (use "SOC" for aftermath, since that's the SOC with no shortening)

#### Define interventions. ####

screening_design_guidelines <- 
  screening_design <- data.frame(
    "counseling_coverage" = 0, 
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

screening_design_counseling <- 
  screening_design <- data.frame(
    "counseling_coverage" = 1, 
    "timing_months" = c(6, 12, 18),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## add 3m symptom home visit
screening_design_addearly <- 
  screening_design <- data.frame(
    counseling_coverage = 1,
    "timing_months" = c(3, 6, 12, 18),
    "target_coverage" = c(1, 1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## sputum at initial visit 
screening_design_addsputum <- 
  screening_design <- data.frame(
    counseling_coverage = 1,
    "timing_months" = c(3, 6, 12, 18),
    "target_coverage" = c(1, 1, 1, 1),
    "screening_method" = c("both", "symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## earlier screening
screening_design_earlier <- 
  screening_design <- data.frame(
    "counseling_coverage" = 1, 
    "timing_months" = c(3, 9, 15),
    "target_coverage" = c(1, 1, 1),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)


## limit to high risk
screening_design_limited <- 
  screening_design <- data.frame(
    "counseling_coverage" = 0, 
    "timing_months" = c(6, 12, 18),
    "target_coverage" = rep(0.5, 3),
    "screening_method" = c("symptoms", "symptoms", "symptoms"),
    "screening_location" = c("home", "phone", "phone")
  ) %>% 
  arrange(timing_months)

## incrase frequency and add sputum, each in targeted way 
screening_design_intensified <- 
  screening_design <- data.frame(
    "counseling_coverage" = 0.5, 
    "timing_months" = c(3, 6, 9, 12, 15, 18),
    "target_coverage" = c(0.5,1,0.5,1,0.5,1),
    "screening_method" = c("both", rep("symptoms", 5)),
    "screening_location" = c("home", "home", rep("phone",4))
  ) %>% 
  arrange(timing_months)

intervention_names = c("guidelines", "counseling", "addearly", "addsputum", "earlier", "limited", "intensified")

##### Run interventions ####
results <- list()
for (name in intervention_names) {
  results[[name]] <- 
            data.frame(detections = numeric(N_samples),
                      mean_symptom_days = numeric(N_samples),
                      symptomatic_months_soc = numeric(N_samples),
                      infectious_months_soc = numeric(N_samples),
                      symptomatic_months_averted = numeric(N_samples), 
                      infectious_months_averted = numeric(N_samples))
  # , 
  #                     cost = numeric(N_samples),
  #                     cost_per_symptomatic_month_averted = numeric(N_samples),
  #                     cost_per_infectious_months_averted = numeric(N_samples))
}




########## apply interventions to all sampled parameter sets #######

for (n in 1:N_samples)
{

  # setup
  cohort <- create_cohort(cohort_params[n,])

  intervention_parameters <- list(
    coverage = list("phone" = cohort_params$coverage_phone[n], 
                    "home" = cohort_params$coverage_phone[n] * cohort_params$coverage_home_reduction[n]),
    sensitivity_symptoms = list("home" = cohort_params$sensitivity_symptoms_home[n], 
                                "phone" = cohort_params$sensitivity_symptoms_home[n] * cohort_params$sensitivity_symptoms_phone_reduction[n]),
    success_sputum = list("home" = cohort_params$success_sputum_home[n], 
                          "phone" = cohort_params$success_sputum_home[n] * cohort_params$success_sputum_phone_reduction[n]))
  
  # interventions
  outputs <- list()
  for (intervention in intervention_names)
    outputs[[intervention]] <- apply_intervention(cohort, get(paste0("screening_design_", intervention)), 
                                       intervention_parameters, cohort_params[n,])  

  # impact
  (detected <- lapply(outputs, function(x) sum(!is.na(x$detection_timing))))
  (symptomdays <- lapply(outputs, function(x) ((time_with_tb(x) %>% filter(scenario == "screening") %>% select(value))/sum(x$TB))[1,]))
  (time <- lapply(outputs, function(x) time_with_tb(x) %>% filter(scenario == "soc") %>% 
    mutate(months = value/30) %>%
    tibble::column_to_rownames('outcome') %>% select(months)))
  (impact <- lapply(outputs, function(x) (time_with_tb(x)  %>% pivot_wider(names_from = "scenario", values_from = "value") %>% 
      mutate(months_averted = (soc - screening)/30)) %>% 
      tibble::column_to_rownames('outcome')  %>% select(months_averted)))
  
#   # costs
#   (cost <- sum(costs(screening_design_aftermath, cohort_screened, covered_screening, 
#                      cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   (cost_earlier <- sum(costs(screening_design_earlier, cohort_screened_earlier, covered_screening_earlier,
#                              cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   (cost_targeted <- sum(costs(screening_design_targeted, cohort_screened_targeted, covered_screening_targeted,
#                               cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   (cost_sputum <- sum(costs(screening_design_sputum, cohort_screened_sputum, covered_screening_sputum,
#                             cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   (cost_moretargeted <- sum(costs(screening_design_moretargeted, cohort_screened_moretargeted, covered_screening_moretargeted,
#                               cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   (cost_alt <- sum(costs(screening_design_aftermath, cohort_screened_alt, covered_screening_alt,
#                          cohort_params$initial_contact_cost_home[n], cohort_params$initial_contact_cost_phone[n], cohort_params$sputum_test_cost[n])))
#   
#   # cost-benefit ratios
#   cost/impact
#   
# 
  for (name in intervention_names)
    results[[name]] [n,] <- c(detected[[name]], symptomdays[[name]], time[[name]]$months[1], time[[name]]$months[2], 
              impact[[name]]$months_averted[1], impact[[name]]$months_averted[2])
              # , 
              # cost,
              # cost/impact$months_averted[1], cost/impact$months_averted[2])
}


# collate results for a table:
# for all dataframe elements of list "results", report the mean and 25th and 75th percentiles of each column

lapply(results, function(x) 
  # mean of each columns of table
  x %>% 
    mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
                                                     inf_reduction = infectious_months_averted/infectious_months_soc*100) %>%
    select(detections, mean_symptom_days, sx_reduction, inf_reduction) %>%
    
    summarise_all(function(y) paste0(round(median(y),0), " (", round(quantile(y, 0.25),0), ", ", round(quantile(y, 0.75),0), ")"))) %>%
  # convert list elements to rows of a single dataframe, with intervention_names as rownames
  bind_rows() %>%
  mutate(intervention = intervention_names) %>%
  # make intervention the first column, and mean_symptom_days the second column
  select(intervention, mean_symptom_days, everything()) %>%
  write_clip() 
  
saveRDS(results, "NTPresults.Rdata")  
readRDS(results, "NTPresults.Rdata")  

quantile_df <- function(x, probs = c(0.25, 0.5, 0.75)) {
  tibble(quantile = probs, value = quantile(x, probs))
}

plotdata <- lapply(results, function(x) 
  # mean of each columns of table
  x %>% 
    mutate(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, 
           inf_reduction = infectious_months_averted/infectious_months_soc*100) %>%
    select(detections, mean_symptom_days, sx_reduction, inf_reduction) %>%
    reframe(across(where(is.double), quantile_df, .unpack = TRUE)) %>%
    # change the first columne name to "quantile" and drop all other columns whose names end with "_quantile"
    mutate(quantile = detections_quantile) %>% 
    select(-ends_with("_quantile")) %>%
        pivot_longer(cols = ends_with("value"), names_to = c(".value","value"), names_pattern = "(.*)_(.*)") %>%
    select(-value)) %>%
    bind_rows() %>%
    mutate(intervention = rep(intervention_names, each=3)) %>% 
  pivot_longer(names_to = "outcome", cols= c(detections, mean_symptom_days, sx_reduction, inf_reduction)) %>% 
  mutate(point = case_when(quantile==0.25 ~ "lci",
                           quantile == 0.5 ~ "median",
                           quantile == 0.75 ~ "uci")) %>% 
  select(-quantile) %>%
  pivot_wider(names_from = point, values_from = value) %>% 
  mutate(intervention = factor(intervention,
                               levels = c("guidelines","counseling", "addearly", "addsputum",
                                          "earlier","limited", "intensified"),
                               labels = c("SOC only", "Add effective counseling", "Add early visit", "Add sputum testing",
                                          "Earlier schedule", "Limit to high risk", "Intensify for high risk")),
         outcome_nicename = case_when(
           outcome == "sx_reduction" ~ "Symptomatic TB time",
           outcome == "inf_reduction" ~ "Infectious TB time",
           TRUE ~ outcome
         ))

library(RColorBrewer)
library(colorspace)
pal <- brewer.pal(name = "YlGnBu", n = 5)[2:5]
fillcols <- c(pal, 
              lighten(pal[1], amount = 0.4, space = "HCL"), 
              darken(pal[1], amount = 0.4, space = "HCL"), darken(pal[4], amount = 0.4, space = "HCL"))

days <- ggplot(plotdata %>% filter(outcome == "mean_symptom_days"), aes(x=intervention, y = median)) + 
  geom_col(aes(fill=intervention)) +
  geom_errorbar(aes(ymin = lci, ymax = uci)) + 
  scale_fill_manual(values = fillcols) +
  ylab("Mean days with symptoms, per recurrence") + 
  theme_cowplot(12) + 
  theme(axis.text.x = element_blank()) +
  labs(fill = "Intervention")

dxs <- ggplot(plotdata %>% filter(outcome == "detections"), aes(x=intervention, y = median)) + 
  geom_col(aes(fill=intervention)) +
  geom_errorbar(aes(ymin = lci, ymax = uci)) + 
  scale_fill_manual(values = fillcols) +
  ylab("Recurrences detected during screening\n(per 100 000 patients") +
  theme_cowplot(12) + 
  theme(axis.text.x = element_blank()) +
  labs(fill = "Intervention")

time <- ggplot(plotdata %>% filter(outcome %in% c("inf_reduction", "sx_reduction")), 
       aes(x=intervention, y = -median/100)) + 
  geom_col(aes(fill=intervention)) +
  scale_fill_manual(values = fillcols) +
  geom_errorbar(aes(ymax = -lci/100, ymin = -uci/100)) + 
  facet_grid(rows = vars(outcome_nicename)) + 
  ylab("Change in time with TB") +
  scale_y_continuous(labels = scales::percent) +
  theme_cowplot(12) + 
  theme(axis.text.x = element_blank()) + #, strip.text = element_text(size = 14)) + 
  labs(fill = "Intervention")

((dxs) | time ) + 
  plot_layout(guides = "collect") +
  theme(legend.position = "left",
        plot.margin = unit(c(1,1,1,1), "cm"))

# 
# 
# 
# 
# # sensitivty analysis
# # for each parameter in cohort_param_ranges, compare results$cost_per_symptomatic_case_averted 
# # between the top 10% and bottom 10% of values for that paramter
# parameters_varied <- names(cohort_param_ranges)[lapply(cohort_param_ranges, length)>1]
# oneways <- array(NA, dim=c(2, length(parameters_varied)))
# dimnames(oneways) = list(c("top", "bottom"), parameters_varied)
# for (p in parameters_varied) {
#   oneways["top", p] <- mean(results_aftermath[cohort_params[p] >= quantile(unlist(cohort_params[p]), 0.8),]$cost_per_symptomatic_month_averted)
#   oneways["bottom", p] <- mean(results_aftermath[cohort_params[p] <= quantile(unlist(cohort_params[p]), 0.2),]$cost_per_symptomatic_month_averted)
# }
# # plot tornado diagram, centered at mean of results$cost_per_symptomatic_month_averted
# oneways %>%
#   as.data.frame() %>%
#   rownames_to_column("top_or_bottom") %>%
#   pivot_longer(cols = -top_or_bottom, names_to = "parameter", values_to = "cost_per_symptomatic_month_averted") %>%
#   ggplot(aes(x = parameter, 
#              col = top_or_bottom)) +
#   geom_errorbar(aes(ymin = mean(results$cost_per_symptomatic_month_averted), ymax=cost_per_symptomatic_month_averted)) +
#   coord_flip() +
#   geom_hline(yintercept = mean(results$cost_per_symptomatic_month_averted), linetype = "dashed") + 
#   # specify text labels in color legend
#   labs(color = "Quantile of parameter values") + 
#   scale_color_manual(labels = c("Lowest 20%", "Highest 20%"), values=c("red","blue"))
#   
#   
# #############
# 
# 
# 
# 
# 
