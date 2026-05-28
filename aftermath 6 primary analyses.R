# Analyze outputs after running create_cohort and run_interventions

#### Filtered and unfiltered results objects ####

if (apply_subclinical_filter) 
{
  results_filtered <- results
  results_unfiltered <- NULL
  } else{
  results_unfiltered <- results
  
  accepted_index <- cohort_features$accepted_subclinical == 1
  
  results_filtered <- lapply(
    results_unfiltered,
    function(x) x[accepted_index, ]
  )
}

results_main <- results_filtered

#### Quick checks - may be redundant ####

mean(cohort_features$accepted_subclinical)

cohort_features %>%
  count(fails_baseline_high, fails_6mo_low, fails_6mo_high) %>%
  mutate(prop = n / sum(n))

cohort_features %>%
  ggplot(aes(
    x = subclinical_baseline_among_micropos,
    y = subclinical_6mo_amongcohort,
    color = factor(accepted_subclinical)
  )) +
  geom_point(alpha = 0.5) +
  geom_vline(
    xintercept = cohort_params$subclinical_baseline_amongTB_max[1],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = cohort_params$subclinical_6m_amongcohort_min[1],
    linetype = "dashed"
  ) +
  geom_hline(
    yintercept = cohort_params$subclinical_6m_amongcohort_max[1],
    linetype = "dashed"
  ) +
  theme_minimal() +
  xlab("Subclinical NAAT+ at treatment completion, among micropositive recurrences") +
  ylab("Subclinical NAAT+ prevalence at 6 months, among full cohort") +
  labs(color = "Accepted")


### Check weibull fits and whether they are driving rejection
ggplot(aes(abs(median_dx_by_540_sim - target_median_dx), accepted_subclinical))

cor(
  cohort_params$objective_value,
  cohort_features$accepted_subclinical
)

# And check distribution of objective values (if tehre's a long tail, should add something to reject that)
hist(cohort_params$objective_value, breaks=50)

#### Evaluate cohort characteristics and natural history results ####

head(cohort_features)
cohort_features[,] %>% summarise_all(median)
cohort_features[,] %>% summarise_all(mean)

sum(cohort_features$accepted_subclinical)
mean(cohort_features$accepted_subclinical)

rbind(
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) median(x, na.rm = TRUE)),
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) quantile(x, 0.025, na.rm = TRUE)),
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) quantile(x, 0.975, na.rm = TRUE))
)

rbind(
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) median(x / 30, na.rm = TRUE)),
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) quantile(x / 30, 0.025, na.rm = TRUE)),
  cohort_features[cohort_features$accepted_subclinical == 1, ] %>%
    summarise_all(function(x) quantile(x / 30, 0.975, na.rm = TRUE))
)

cohort_features %>%
  count(fails_baseline_high, fails_6mo_low, fails_6mo_high) %>%
  mutate(prop = n / sum(n))

#### Figure 3: cascade plot ####

cascade_features[,] %>%
  filter(accepted_subclinical == 1) %>%
  summarise_all(median)

cascade_features$n <- seq_len(nrow(cascade_features))

plotdata <- cascade_features %>%
  filter(accepted_subclinical == 1) %>%
  mutate(
    recur_start = 0,
    recur_end = cumulative_incidence,
    recur_label = "Develops recurrent TB",
    recur_proportion = cumulative_incidence / N_cohort,
    recur_group = "All recurrences",
    
    before_6mo_start = 0,
    before_6mo_end = cumulative_incidence - TB_beyond_6mo,
    before_6mo_label = "Too early (<6mo)",
    before_6mo_proportion = (cumulative_incidence - TB_beyond_6mo) / cumulative_incidence,
    before_6mo_group = "Timing of recurrence",
    
    between618mo_start = cumulative_incidence - TB_beyond_6mo,
    between618mo_end = cumulative_incidence - (TB_beyond_6mo - TB_beyond6_before18),
    between618mo_label = "Between 6-18 months",
    between618mo_proportion = TB_beyond6_before18 / cumulative_incidence,
    between618mo_group = "Timing of recurrence",
    
    after18mo_start = cumulative_incidence - (TB_beyond_6mo - TB_beyond6_before18),
    after18mo_end = cumulative_incidence,
    after18mo_label = "Too late (>18mo)",
    after18mo_proportion = (TB_beyond_6mo - TB_beyond6_before18) / cumulative_incidence,
    after18mo_group = "Timing of recurrence",
    
    asxatvisit_start = between618mo_start,
    asxatvisit_end = between618mo_start + (TB_at_visit - symptomatic_at_visit),
    asxatvisit_label = "Asymptomatic at visit",
    asxatvisit_proportion = (TB_at_visit - symptomatic_at_visit) / TB_beyond6_before18,
    asxatvisit_group = "Status at time of visit",
    
    sxatvisit_start = between618mo_start + (TB_at_visit - symptomatic_at_visit),
    sxatvisit_end = between618mo_start + TB_at_visit,
    sxatvisit_label = "Symptomatic at visit",
    sxatvisit_proportion = symptomatic_at_visit / TB_beyond6_before18,
    sxatvisit_group = "Status at time of visit",
    
    notbatvisit_start = between618mo_start + TB_at_visit,
    notbatvisit_end = between618mo_end,
    notbatvisit_label = "No TB at visit",
    notbatvisit_proportion = (TB_beyond6_before18 - TB_at_visit) / TB_beyond6_before18,
    notbatvisit_group = "Status at time of visit",
    
    undetected_start = sxatvisit_start,
    undetected_end = sxatvisit_start + (symptomatic_at_visit - detected_and_linked),
    undetected_label = "Missed",
    undetected_proportion = (symptomatic_at_visit - detected_and_linked) / symptomatic_at_visit,
    undetected_group = "Outcome of visit",
    
    detectedlinked_start = sxatvisit_start + (symptomatic_at_visit - detected_and_linked),
    detectedlinked_end = sxatvisit_start + symptomatic_at_visit,
    detectedlinked_label = "Detected (early)",
    detectedlinked_proportion = detected_and_linked / symptomatic_at_visit,
    detectedlinked_group = "Outcome of visit",
    
    unaverted_start = detectedlinked_start,
    unaverted_end = detectedlinked_start +
      detected_and_linked * (total_time_of_linked - remaining_time_of_linked) / total_time_of_linked,
    unaverted_label = "Already elapsed",
    unaverted_proportion = (total_time_of_linked - remaining_time_of_linked) / total_time_of_linked,
    unaverted_group = "Time with TB",
    
    averted_start = detectedlinked_start +
      detected_and_linked * (total_time_of_linked - remaining_time_of_linked) / total_time_of_linked,
    averted_end = detectedlinked_start + detected_and_linked,
    averted_label = "Averted",
    averted_proportion = remaining_time_of_linked / total_time_of_linked,
    averted_group = "Time with TB"
  )

plotdata <- plotdata %>%
  pivot_longer(
    cols = c(recur_start:averted_group),
    names_to = c("outcome", ".value"),
    names_pattern = "(.*)_(.*)"
  ) %>%
  select(
    n, outcome, start, end, label, proportion, group,
    cumulative_incidence, TB_beyond6_before18,
    TB_at_visit, symptomatic_at_visit, detected_and_linked
  )

plotdata$outcome <- factor(
  plotdata$outcome,
  levels = c(
    "recur", "before_6mo", "between618mo", "after18mo",
    "asxatvisit", "sxatvisit", "notbatvisit",
    "undetected", "detectedlinked",
    "unaverted", "averted"
  )
)

plotdata$group <- factor(
  plotdata$group,
  levels = c(
    "All recurrences",
    "Timing of recurrence",
    "Status at time of visit",
    "Outcome of visit",
    "Time with TB"
  )
)

plotdata$outcomex <- as.numeric(plotdata$group) + as.numeric(plotdata$outcome)

labelmapping <- plotdata %>%
  count(outcomex, label) %>%
  select(outcomex, label)

pathmapping <- cbind(
  c(
    "recur", "before_6mo", "between618mo", "after18mo",
    "asxatvisit", "sxatvisit", "notbatvisit",
    "undetected", "detectedlinked",
    "unaverted", "averted"
  ),
  c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE)
)

plotdata$path <- as.logical(pathmapping[plotdata$outcome, 2])

axistext <- data.frame(
  text = c("Occurrence", "Timing", "Status at visit", "Detection", "Time with TB"),
  position = c(2, 5, 9, 12.5, 15.5)
)

plotdata <- plotdata %>%
  mutate(
    denominator = case_when(
      group == "All recurrences" ~ N_cohort,
      group == "Timing of recurrence" ~ cumulative_incidence,
      group == "Status at time of visit" ~ TB_beyond6_before18,
      group == "Outcome of visit" ~ symptomatic_at_visit,
      group == "Time with TB" ~ detected_and_linked
    )
  )

cascadefig <- ggplot(
  plotdata %>%
    group_by(outcomex, group) %>%
    summarise(
      start = mean(start),
      end = mean(end),
      label = first(label),
      proportion = median(proportion),
      cumulative_incidence = first(cumulative_incidence),
      TB_beyond6_before18 = mean(TB_beyond6_before18),
      TB_at_visit = mean(TB_at_visit),
      symptomatic_at_visit = mean(symptomatic_at_visit),
      detected_and_linked = mean(detected_and_linked),
      path = ifelse(first(path), "Yes", "No"),
      text = paste0(round(100 * proportion, 0), "%"),
      .groups = "drop"
    ),
  aes(x = outcomex, ymin = start / N_cohort, ymax = end / N_cohort, col = group)
) +
  geom_linerange(size = 10, aes(alpha = factor(path))) +
  geom_label(aes(
    x = outcomex + 0.25,
    y = (start + (end - start) / 2) / N_cohort - 0.003,
    label = text
  )) +
  geom_errorbar(
    data = plotdata %>%
      group_by(outcomex, group) %>%
      mutate(
        endmin = mean(start) / N_cohort +
          mean(denominator) / N_cohort * quantile((end - start) / denominator, 0.025),
        endmax = mean(start) / N_cohort +
          mean(denominator) / N_cohort * quantile((end - start) / denominator, 0.975)
      ),
    aes(ymin = endmin, ymax = endmax),
    col = "black",
    alpha = 0.04,
    width = 0.3
  ) +
  theme_minimal() +
  ylab("Proportion of treatment-completing cohort") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  xlab("Step in recurrent-TB case-finding cascade") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  scale_x_continuous(
    breaks = labelmapping$outcomex,
    labels = labelmapping$label
  ) +
  scale_alpha_discrete(range = c(0.3, 0.9)) +
  guides(
    col = guide_legend(title = "Portion of cascade", override.aes = aes(label = "")),
    alpha = guide_legend(
      title = "On pathway to detection\nthrough screening?",
      reverse = TRUE,
      override.aes = list(colour = "darkgray")
    )
  ) +
  annotate("text", x = axistext$position, y = -0.003, label = axistext$text, fontface = 2)

pdf("Figure 3 - cascade plot.pdf", width = 12, height = 7)
cascadefig
dev.off()

plotdata %>% filter(outcome == "before_6mo") %>% summarise(quantile(proportion, c(0.5, 0.025, 0.975)))
plotdata %>% filter(outcome == "after18mo") %>% summarise(quantile(proportion, c(0.5, 0.025, 0.975)))
plotdata %>% filter(outcome == "sxatvisit") %>% summarise(quantile(proportion, c(0.5, 0.025, 0.975)))
plotdata %>% filter(outcome == "detectedlinked") %>% summarise(quantile(proportion, c(0.5, 0.025, 0.975)))
plotdata %>% filter(outcome == "averted") %>% summarise(quantile(proportion, c(0.5, 0.025, 0.975)))

#### Look at intervention results ####

cbind(
  lapply(results_main, function(x)
    x %>%
      mutate(
        sx_reduction = symptomatic_months_averted / symptomatic_months_soc * 100,
        inf_reduction = infectious_months_averted / infectious_months_soc * 100,
        proportion_detected = detections / recurrences * 100
      ) %>%
      select(proportion_detected, sx_reduction, inf_reduction) %>%
      summarise(across(
        everything(),
        list(percent = function(y) paste0(
          round(median(y, na.rm = TRUE), 0), "% (",
          round(quantile(y, 0.025, na.rm = TRUE), 0), "-",
          round(quantile(y, 0.975, na.rm = TRUE), 0), "%)"
        ))
      ))
  ) %>%
    bind_rows() %>%
    mutate(intervention = intervention_names),
  
  lapply(results_main, function(x)
    x %>%
      mutate(
        proportion_detected_percent = detections / recurrences * 100,
        cost_per_detection = cost / detections,
        cost_per_patient = cost / N_cohort,
        incremental_cost_effectiveness =
          (cost - results_main$guidelines$cost) /
          (detections - results_main$guidelines$detections),
        incremental_cost_effectiveness_vs_3earlier =
          (cost - results_main$earlier_three$cost) /
          (detections - results_main$earlier_three$detections)
      ) %>%
      select(
        recurrences, detections, mean_symptom_days, cost,
        cost_per_patient, cost_per_detection,
        incremental_cost_effectiveness,
        incremental_cost_effectiveness_vs_3earlier
      ) %>%
      summarise(across(
        everything(),
        list(
          int = function(y) paste0(
            "$", round(median(y, na.rm = TRUE), -1), " (",
            round(quantile(y, 0.025, na.rm = TRUE), -1), ", ",
            round(quantile(y, 0.975, na.rm = TRUE), -1), ")"
          ),
          cents = function(y) paste0(
            "$", round(median(y, na.rm = TRUE), 1), " (",
            round(quantile(y, 0.025, na.rm = TRUE), 1), ", ",
            round(quantile(y, 0.975, na.rm = TRUE), 1), ")"
          )
        )
      ))
  ) %>%
    bind_rows() %>%
    mutate(intervention = intervention_names) %>%
    select(-intervention)
) %>%
  select(
    intervention,
    proportion_detected_percent,
    sx_reduction_percent,
    inf_reduction_percent,
    cost_per_patient_cents,
    cost_per_detection_int,
    incremental_cost_effectiveness_int,
    incremental_cost_effectiveness_vs_3earlier_int
  ) %>%
  write_clip()

pairwise_props <-
  lapply(results_main, function(x)
    ((x %>% select(detections, symptomatic_months_averted, infectious_months_averted)) /
       (results_main$guidelines %>% select(detections, symptomatic_months_averted, infectious_months_averted))) %>%
      summarise_all(function(y) c(
        round(median(y, na.rm = TRUE), 2),
        round(quantile(y, 0.025, na.rm = TRUE), 2),
        round(quantile(y, 0.975, na.rm = TRUE), 2)
      ))
  )

pairwise_diffs <-
  lapply(results_main, function(x)
    ((x %>%
        mutate(proportion = detections / recurrences) %>%
        select(proportion, symptomatic_months_averted, infectious_months_averted)) -
       (results_main$guidelines %>%
          mutate(proportion = detections / recurrences) %>%
          select(proportion, symptomatic_months_averted, infectious_months_averted))) %>%
      summarise_all(function(y) c(
        round(median(y, na.rm = TRUE), 2),
        round(quantile(y, 0.025, na.rm = TRUE), 2),
        round(quantile(y, 0.975, na.rm = TRUE), 2)
      ))
  )

########
# Figure 4 
#########


# One more new effort to replace everything below: 
# Rank intervetions in order of incremental cost-effectiveness per case detected, 
# And plot connected cots for incremental cost and increemntal impact. 
# Also show panels for the other outcomes but without the line segments where not sequential. 
# This will require running the following interventions incrementally: 
# one_visit_6 = results$`6m_sx`
# two_visits_36 = results$earlier_two
# three_visits_369 = results$earlier_three
# four_visits_36912
# five_visits_3691215
# six_visits_369121518 = results$frequent
# six_visits_369121518_sputum3 = results$frequent_sputum
# five_visits_3691215_sputum3 # including just in case better than 6 visits



# Initialize list to hold calculations
proportions <- list()

names_to_include <- c("6m_sx", "earlier_two", "earlier_three", "guidelines", "earlier_three_sputum", "four_visits_36912", 
                      "five_visits_3691215", "frequent", "frequent_sputum", "five_visits_sputum", "four_visits_sputum")

# Step 1: Calculate proportions detected for each intervention
for (intervention in names_to_include) {
  df <- results_main[[intervention]]
  # Calculate proportion detected
  df <- df %>%
    mutate(proportion_detected = detections / recurrences,
           symptomatic_averted = symptomatic_months_averted/symptomatic_months_soc,
           infectious_averted = infectious_months_averted/infectious_months_soc,
           sim_number = row_number())
  
  # Store results
  proportions[[intervention]] <- df
}

# Combine all results into one data frame
combined_results <- bind_rows(lapply(proportions, function(x) {
  x %>% select(detections, recurrences, cost, proportion_detected, 
               symptomatic_averted, infectious_averted, sim_number)
}), .id = "intervention")


# Step 2: Plot Cost vs Proportion Detected
# Create a data frame for the plot
plot_data <- combined_results %>%
  select(intervention, cost, proportion_detected, symptomatic_averted, infectious_averted) %>%
  group_by(intervention) %>%
  summarize(
    median_cost_perperson = median(cost/N_cohort, na.rm=T),
    median_proportion_detected = median(proportion_detected, na.rm=T),
    median_symptomatic_averted = median(symptomatic_averted, na.rm=T),
    median_infectious_averted = median(infectious_averted, na.rm=T),
    uci_cost_perperson = quantile(cost/N_cohort, 0.975, na.rm=T),
    lci_cost_perperson = quantile(cost/N_cohort, 0.025, na.rm=T),
    uci_proportion_detected = quantile(proportion_detected, 0.975, na.rm=T),
    lci_proportion_detected = quantile(proportion_detected, 0.025, na.rm=T),
    uci_symptomatic_averted = quantile(symptomatic_averted, 0.975, na.rm=T),
    lci_symptomatic_averted = quantile(symptomatic_averted, 0.025, na.rm=T),
    uci_infectious_averted = quantile(infectious_averted, 0.975, na.rm=T),
    lci_infectious_averted = quantile(infectious_averted, 0.025, na.rm=T),
    .groups = 'drop'
  ) 
plot_data <- plot_data %>%
  pivot_longer(
    cols = c(median_proportion_detected, median_symptomatic_averted, median_infectious_averted, 
             uci_proportion_detected, lci_proportion_detected,
             uci_symptomatic_averted, lci_symptomatic_averted,
             uci_infectious_averted, lci_infectious_averted),
    names_pattern = "(.*)_(.*_.*)",
    names_to = c("stat", "outcome")
  ) %>%
  pivot_wider(
    names_from = stat,
    values_from = value
  ) %>%
  # set outcome as factor with levels in desired order
  mutate(outcome = factor(outcome, levels = c("proportion_detected", "symptomatic_averted", "infectious_averted"))) %>% 
  # rename the interventions
  mutate(
    intervention = case_when(
      intervention == "6m_sx" ~ "One screening (6m)",
      intervention == "earlier_two" ~ "Two screenings (3, 6m)",
      intervention == "earlier_three" ~ "Three screenings (3, 6, 9m)",
      intervention == "guidelines" ~ "Three screenings, guidelines schedule (6, 12, 18m)",
      intervention == "earlier_three_sputum" ~ "Three screenings + 3m micro",
      intervention == "four_visits_36912" ~ "Four screenings (3, 6, 9, 12m)",
      intervention == "four_visits_sputum" ~ "Four screenings + 3m micro",
      intervention == "five_visits_3691215" ~ "Five screenings (3, 6, 9, 12, 15m)",
      intervention == "frequent" ~ "Six screenings (3, 6, 9, 12, 15, 18m)",
      intervention == "five_visits_sputum" ~ "Five screenings + 3m micro",
      intervention == "frequent_sputum" ~ "Six screenings + 3m micro",
      TRUE ~ intervention
    ),
    short_intervention = case_when(
      intervention == "One screening (6m)" ~ "One (6m)",
      intervention == "Two screenings (3, 6m)" ~ "Two (3,6m)",
      intervention == "Three screenings (3, 6, 9m)" ~ "Three (3,6,9m)",
      intervention == "Three screenings, guidelines schedule (6, 12, 18m)" ~ "Guidelines",
      intervention == "Three screenings + 3m micro" ~ "Three + 3m micro",
      intervention == "Four screenings (3, 6, 9, 12m)" ~ "Four (3,6,9,12m)",
      intervention == "Four screenings + 3m micro" ~ "Four + 3m micro",
      intervention == "Five screenings (3, 6, 9, 12, 15m)" ~ "Five (3,6,9,12,15m)",
      intervention == "Five screenings + 3m micro" ~ "Five + 3m micro",
      intervention == "Six screenings (3, 6, 9, 12, 15, 18m)" ~ "Six (3,6,9,12,15,18m)",
      intervention == "Six screenings + 3m micro" ~ "Six + 3m micro",
      TRUE ~ intervention
    )
  )

plot_data <- plot_data %>% 
  mutate(intervention = factor(intervention, levels = plot_data %>% 
                                 arrange(median_cost_perperson) %>% 
                                 pull(intervention) %>% unique()))


# Create the plot
(CEplot <- ggplot(plot_data, aes(y = median_cost_perperson, x = median, label = intervention, col = intervention)) +
    geom_point(size = 4) +  # Points for each intervention
    # geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.2, color = "blue") +
    facet_wrap(.~ outcome, 
               nrow=1,
               labeller = as_labeller(c(
                 proportion_detected = "Recurrences Detected",
                 symptomatic_averted = "Symptomatic Months Averted",
                 infectious_averted = "Infectious Months Averted"
               ))) +
    theme_minimal() +
    xlab("Proportion detected or averted") + 
    ylab("Cost per TB survivor ($US)") +
    # legend on bottom
    theme(legend.position = "bottom")  + 
    # gradient color scale without any very light shades
    scale_color_viridis_d(option = "D") + 
    # add text labels only for last facet, and change hjust to 1 for only the 3rd and 8th label
    geom_text(data = plot_data %>%
                filter(outcome == "proportion_detected",
                       intervention %in%
                         c("Three screenings, guidelines schedule (6, 12, 18m)",
                           "Three screenings + 3m micro")),
              aes(y = median_cost_perperson, x = median, label = short_intervention),
              hjust = 1, vjust = -1, size = 2.5) +
    coord_cartesian(clip = "off") +
    geom_errorbar(aes(xmin = lci, xmax = uci, color = intervention), width = 0.2, alpha = 0.4, lwd=0.3) + 
    geom_errorbar(aes(ymin = lci_cost_perperson, ymax = uci_cost_perperson, color = intervention), width = 0.01, alpha = 0.4, lwd=0.3) +
    geom_text(data = plot_data %>% 
                filter(outcome == "proportion_detected", 
                       !(intervention %in% 
                           c("Three screenings, guidelines schedule (6, 12, 18m)", 
                             "Three screenings + 3m micro"))),
              aes(y = median_cost_perperson, x = median, label = paste0("   ",short_intervention)), 
              hjust = 0,  vjust = 1, size = 2.5) + 
    coord_cartesian(ylim = c(0, 16), xlim=c(0,0.38)))


# Add line segments connecting incrementally cost-effective interventions for each outcome
# To do this, for each outcome, for the ordered list of interventions, estimate incremental cost per incremetnal porportion detected or averted
# Get list of interventions ordered by cost, as plotted along the y axis above 
ordered_interventions <- plot_data %>%
  filter(outcome == "proportion_detected") %>%
  arrange(median_cost_perperson) %>%
  pull(intervention)

# For each of the three outcomes, identify the series of incrementally cost-effective interventions
ce_sequence_list <- list()
for (o in levels(plot_data$outcome))
{
  
  # make a matrix of icers, rows and columns both as ordered_interventions
  icer_matrix <- matrix(NA, nrow = length(ordered_interventions), ncol = length(ordered_interventions))
  for (i in 1:(length(ordered_interventions)-1))
  {
    for (j in (i+1):length(ordered_interventions))
    {
      intervention_i <- plot_data %>% filter(intervention == ordered_interventions[i], outcome == o)
      intervention_j <- plot_data %>% filter(intervention == ordered_interventions[j], outcome == o)
      icer_matrix[i,j] <- (intervention_j$median_cost_perperson - intervention_i$median_cost_perperson) / 
        (intervention_j$median - intervention_i$median)
    }
  }
  # starting with intervention 1, identify the sequence of incrementally cost effective interventions
  icer_matrix[icer_matrix < 0] <- NA
  position <- 1
  ce_sequence <- c(position)
  while (position < length(ordered_interventions)) 
  { position <- which.min(icer_matrix[position, ])
  if (length(position) == 0) break
  ce_sequence <- c(ce_sequence, position)
  }
  ce_sequence_list[[o]] <- ordered_interventions[ce_sequence]
}

# select subset of plot_data for only ce_sequence_list$o for each intervention
segment_data <- data.frame()
for (o in names(ce_sequence_list))
{
  ce_sequence <- ce_sequence_list[[o]]
  segment_data <- rbind(segment_data, plot_data %>% filter(intervention %in% ce_sequence, outcome == o))
}

fig4 <- CEplot + 
  geom_line(data = segment_data, aes(x = median, y = median_cost_perperson), 
            color = "darkgray", size = 0.5) + 
  labs(color = "Intervention") + 
  ggthemes::theme_few() +
  theme(strip.text.x = element_text(size = 12))

pdf("Figure4_CEplot.pdf", width = 12, height = 4.5)
fig4
dev.off()



# calculate icer of 4 screens vs 3 , and of three + sputum vs 3
# but don't use plotdata (which is proportions), instead use raw results
# 4 screens vs 3 screens:
icer_4vs3_detections <- (results_main$four_visits_36912$cost - results_main$earlier_three$cost) / 
  (results_main$four_visits_36912$detections - results_main$earlier_three$detections)
icer_4vs3_symptomatic <- (results_main$four_visits_36912$cost - results_main$earlier_three$cost) / 
  (results_main$four_visits_36912$symptomatic_months_averted - results_main$earlier_three$symptomatic_months_averted)
icer_4vs3_infectious <- pmax(results_main$four_visits_36912$cost - results_main$earlier_three$cost, 0) / 
  pmax(results_main$four_visits_36912$infectious_months_averted - results_main$earlier_three$infectious_months_averted, 0)
icer_sputum_detections <- (results_main$earlier_three_sputum$cost - results_main$earlier_three$cost) / 
  (results_main$earlier_three_sputum$detections - results_main$earlier_three$detections)
icer_sputum_symptomatic <- (results_main$earlier_three_sputum$cost - results_main$earlier_three$cost) / 
  (results_main$earlier_three_sputum$symptomatic_months_averted - results_main$earlier_three$symptomatic_months_averted)
icer_sputum_infectious <- (results_main$earlier_three_sputum$cost - results_main$earlier_three$cost) / 
  (results_main$earlier_three_sputum$infectious_months_averted - results_main$earlier_three$infectious_months_averted)


median(icer_4vs3_detections, na.rm=T)
quantile(icer_4vs3_detections, c(0.025, 0.975), na.rm=T)
median(icer_4vs3_symptomatic, na.rm=T)
quantile(icer_4vs3_symptomatic, c(0.025, 0.975), na.rm=T)
icer
median(icer_4vs3_infectious, na.rm=T)
quantile(icer_4vs3_infectious, c(0.025, 0.975), na.rm=T)
median(icer_sputum_detections, na.rm=T)
quantile(icer_sputum_detections, c(0.025, 0.975), na.rm=T)
median(icer_sputum_symptomatic, na.rm=T)
quantile(icer_sputum_symptomatic, c(0.025, 0.975), na.rm=T)
median(icer_sputum_infectious, na.rm=T)
quantile(icer_sputum_infectious, c(0.025, 0.975), na.rm=T)


# to calculate from results:
#"compared to screening with symptoms at 3, 6, and 9 months after treatment completion, adding a fourth screening visit costs $x more per person, incrementally detects x% of all recurrences, and incrementally averts x% of the time with symptomatic TB and x% of the time with infectious TB that the cohort would have experienced in absence of any screening."
quantile((results_main$four_visits_36912$cost - results_main$earlier_three$cost)/N_cohort, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$four_visits_36912$detections - results_main$earlier_three$detections)/results_main$four_visits_36912$recurrences, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$four_visits_36912$symptomatic_months_averted - results_main$earlier_three$symptomatic_months_averted)/results_main$four_visits_36912$symptomatic_months_soc, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$four_visits_36912$infectious_months_averted - results_main$earlier_three$infectious_months_averted)/results_main$four_visits_36912$infectious_months_soc, c(0.5,0.025,0.975), na.rm=T)

quantile((results_main$earlier_three_sputum$cost - results_main$earlier_three$cost)/N_cohort, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$earlier_three_sputum$detections - results_main$earlier_three$detections)/results_main$earlier_three_sputum$recurrences, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$earlier_three_sputum$symptomatic_months_averted - results_main$earlier_three$symptomatic_months_averted)/results_main$earlier_three_sputum$symptomatic_months_soc, c(0.5,0.025,0.975), na.rm=T)
quantile((results_main$earlier_three_sputum$infectious_months_averted - results_main$earlier_three$infectious_months_averted)/results_main$earlier_three_sputum$infectious_months_soc, c(0.5,0.025,0.975), na.rm=T)


# create a data frame of nice names for all the paramters
param_names <- c(
  incidence_18mo_multiplier = "Cumulative 18-month notification uncertainty multiplier",
  proportion_micro_pos = "Proportion bacteriologically+ at diagnosis",
  # auc = "Predictive accuracy for recurrence (AUC)",
  symptom_duration_meanlog_reported = "Mean log duration of reported symptoms",
  symptom_duration_sdlog_reported = "SD log of duration of reported symptoms",
  reported_fraction_of_true_symptom_duration = "Underestimation of symptom duration",
  programmatic_symptom_duration_factor = "Increase in symptom duration under programmatic conditions",
  proportion_micropos_sputum_first = "Proportion of micropositive recurrences with NAAT+ before symptoms",
  proportion_micropos_subclinical_at_eot = "Proportion of micropositive recurrences NAAT+ at treatment completion",  
  duration_ratio_subclinical_symptomatic= "Relative time asymptomatic vs symptomatic",
  duration_subclinical_cv= "Coefficient of variation in asymptomatic duration",
  subclinical_baseline_amongTB_max = "Maximum proportion with subclinical TB at treatment completion",
  subclinical_6m_amongcohort_min = "Minimum proportion of cohort with subclinical TB at 6 months",
  subclinical_6m_amongcohort_max = "Maximum proportion of cohort with subclinical TB at 6 months",
  coverage_phone =            "Coverage of phone-based screening",
  coverage_home_reduction= "Reduction in coverage for home vs phone screening",
  sensitivity_symptoms_home = "Detection of symptoms if present, home-based screening",
  sensitivity_symptoms_phone_reduction = "Reduction in symptom detection, phone vs home",
  success_sputum_home =           "Sample collection success rate (home)",
  success_sputum_phone_reduction = "Reduction in sample collection success, phone vs home",
  # home_visit_passive_detection_impact = "Impact of home visit on passive detection",
  # intentional_counseling_passive_detection_impact = "Effectiveness of care-seeking-promotion intervention",
  # intentional_counseling_passive_detection_duration = "Duration of effect, care-seeking-promotion intervention",
  # prevention_efficacy= "Effectiveness of prevention intervention",
  initial_contact_cost_home = "Cost of initial home visit",
  initial_contact_cost_home_vs_phone_factor = "Cost of phone contact vs home visit",
  sputum_test_cost = "Cost of NAAT",
  symptom_prevalence_nontb = "Prevalence of non-TB symptoms",
  # prevention_cost = "Cost of prevention intervention",
  # case_fatality = "case fatality",
  initial_contact_cost_phone = "Cost of phone visit",
  coverage_home = "Coverage of home-based screening",
  sensitivity_symptoms_phone = "Detection of symptoms if present, phone-based screening",
  success_sputum_phone= "NAAT completion rate, phone-based screening"
  # incidence_18mo
  # recurrence_shape
  # recurrence_scale
  # recurrence_time_mean
  # recurrence_time_cv
  # probability_dx540_given_recur
  # probability_ever_recur
  # median_dx_by_540_sim
  # p25_dx_by_540_sim
  # p75_dx_by_540_sim
  # prop_dx_le_90_among_dx540_sim
  # prop_dx_le_360_among_dx540_sim
  # target_median_dx
  # target_p25_dx
  # target_p75_dx
  # target_prop_dx_le_90
  # target_prop_dx_le_360
)



# sensitivty analysis
# for each parameter in cohort_param_ranges, compare 
# outcomes of % of symptomatic and infectious TB time averted, comparing
#  (1)  guidelines intervention vs SOC and 
# (2) early, targeted sputum screening vs guidelines intervention. (and perhaps others?)
# And then also include a supplemental figure re: counseling and prevention efficacy/costs:

# guidelines vs soc:
sxarray <- cbind(results$guidelines %>% 
                   reframe(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, ), 
                 as.data.frame(cohort_params))[!is.na(results_main$guidelines$symptomatic_months_averted),]

prcc <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$guidelines$symptomatic_months_averted),]),
  y = (results_main$guidelines %>% 
         reframe(sx_reduction = symptomatic_months_averted/symptomatic_months_soc*100, ))[!is.na(results_main$guidelines$symptomatic_months_averted),],
  rank = TRUE    # Spearman rank correlation (PRCC)
  # nboot = 1000    # optional: bootstrap for CI
)

prcc_inf <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$guidelines$symptomatic_months_averted),]),
  y = (results_main$guidelines %>% 
         reframe(sx_reduction = infectious_months_averted/infectious_months_soc*100, ))[!is.na(results_main$guidelines$symptomatic_months_averted),],
  rank = TRUE    # Spearman rank correlation (PRCC)
  # nboot = 1000    # optional: bootstrap for CI
)



# Extract PRCC values into a dataframe

prcc_results <- data.frame(
  variable = rownames(prcc$PRCC),
  prcc = prcc$PRCC[,1],    # PRCC estimates
  prcc_inf = prcc_inf$PRCC[,1]    # PRCC estimates
  # lower = prcc$PRCC[,2],    # lower CI
  # upper = prcc$PRCC[,3]     # upper CI
) %>%
  
  filter(variable %in% names(param_names)) %>% 
  
  mutate(
    nice_variable = recode(variable, !!!param_names),
    nice_variable = stringr::str_wrap(nice_variable, width = 22),
  ) %>% 
  
  filter(abs(prcc)>0.05 | abs(prcc_inf)>0.05) %>%
  
  pivot_longer(cols = c(prcc, prcc_inf), names_to = "type", values_to = "prcc") %>% 
  mutate(Type = case_when(type=="prcc" ~ "Symptomatic", type=="prcc_inf" ~ "Infectious")) 

### Plot Tornado Diagram ----
prcc_results <- prcc_results %>%
  mutate(nice_variable = reorder(nice_variable, abs(prcc)))

tornado_A <- ggplot(prcc_results, aes(x = nice_variable, y = prcc, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge") +
  # geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  coord_flip() +
  scale_fill_manual(values = c("Symptomatic" = "steelblue4", "Infectious" = "skyblue3")) +
  labs(
    # title = "Tornado Diagram of PRCCs",
    subtitle = "% reduction in time with recurrent TB, vs no screening",
    x = "Model Parameters",
    y = "Partial Rank Correlation Coefficient (PRCC)"
  ) +
  theme_minimal(base_size = 14) +
  # change the legend title
  labs(fill = "Category of recurrent TB time", tag = "A") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text()) +
  # reverse order of items in fill legend
  guides(fill = guide_legend(reverse = TRUE))




# As above, but comparing an alternative intervention to guidelines, in terms of relative change in symptomatic or infectious time averted.

outcomes_sens <- lapply(results_main, function(x) 
  ((x %>% select(detections, symptomatic_months_averted, infectious_months_averted))/ 
     (results_main$guidelines %>% select(detections, symptomatic_months_averted, infectious_months_averted))))

prcc2 <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$guidelines$symptomatic_months_averted),]),
  y = (outcomes_sens$earlier_three[!is.na(results_main$guidelines$symptomatic_months_averted),])$symptomatic_months_averted,
  rank = TRUE
)

prcc2_inf <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$guidelines$symptomatic_months_averted),]),
  y = (outcomes_sens$earlier_three[!is.na(results_main$guidelines$symptomatic_months_averted),])$infectious_months_averted,
  rank = TRUE    # Spearman rank correlation (PRCC)
)

prcc2_sputum <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$earlier_three$symptomatic_months_averted),]),
  y = (outcomes_sens$earlier_three_sputum[!is.na(results_main$earlier_three$symptomatic_months_averted),])$symptomatic_months_averted,
  rank = TRUE
)

prcc2_inf_sputum <- pcc(
  X = as.data.frame((as.matrix(cohort_params))[!is.na(results_main$earlier_three$symptomatic_months_averted),]),
  y = (outcomes_sens$earlier_three_sputum[!is.na(results_main$earlier_three$symptomatic_months_averted),])$infectious_months_averted,
  rank = TRUE    # Spearman rank correlation (PRCC)
  # nboot = 1000    # optional: bootstrap for CI
)

# Extract PRCC values into a dataframe

prcc2_results <- data.frame(
  variable = rownames(prcc2$PRCC),
  prcc = prcc2$PRCC[,1],    # PRCC estimates
  prcc_inf = prcc2_inf$PRCC[,1],
  prcc_sputum = prcc2_sputum$PRCC[,1],
  prcc_inf_sputum = prcc2_inf_sputum$PRCC[,1]
  # lower = prcc$PRCC[,2],    # lower CI
  # upper = prcc$PRCC[,3]     # upper CI
) %>%
  
  filter(variable %in% names(param_names)) %>% 
  
  mutate(
    nice_variable = recode(variable, !!!param_names),
    nice_variable = stringr::str_wrap(nice_variable, width = 22),
  ) %>% 
  
  filter(abs(prcc)>0.1 | abs(prcc_inf)>0.1 | abs(prcc_sputum)>0.1 | abs(prcc_inf_sputum)>0.1) %>%
  pivot_longer(cols = c(prcc, prcc_inf, prcc_sputum, prcc_inf_sputum), 
               names_pattern = "(prcc)(_?inf)?(_?sputum)?", 
               names_to = c("type", "inf", "sputum"),
               values_to = "prcc") %>%
  
  mutate(Type = case_when(inf!="_inf" ~ "Symptomatic", inf=="_inf" ~ "Infectious"),
         Intervention = case_when(sputum!="_sputum" ~ "Earlier screening (vs guidelines)", sputum=="_sputum" ~ "Adding 3m sputum (vs symptoms only)")) 

prcc2_results <- prcc2_results %>%
  mutate(nice_variable = reorder(nice_variable, abs(prcc))) %>%
  # change order of Intervention
  mutate(Intervention = factor(Intervention, levels = c("Earlier screening (vs guidelines)", "Adding 3m sputum (vs symptoms only)")))


tornado_B <- ggplot(prcc2_results, aes(x = nice_variable, y = prcc, fill = Type)) +
  facet_grid(. ~ Intervention) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Symptomatic" = "steelblue4", "Infectious" = "skyblue3")) +
  labs(
    title = "",
    subtitle = "Incremental impact, alternative screening strategies",
    x = "Model Parameters",
    y = "Partial Rank Correlation Coefficient (PRCC)"
  ) +
  theme_minimal(base_size = 14) +
  # change the legend title
  labs(fill = "Category of recurrent TB time", tag = "B") + 
  theme(plot.title = element_text(hjust = 0.5),
        plot.tag = element_text()) +
  # no legend
  theme(legend.position = "none") 



library(gridExtra)
grid.arrange(tornado_A, tornado_B, ncol=1)
pdf("Figure5_tornado_aftermath.pdf", width = 11, height = 10)
grid.arrange(tornado_A, tornado_B, ncol=1)
dev.off()


# # compare deciles of a given parameter. 
# # parameter: 1/symptom_underestimation_factor
# # outcome: symptomatic_months_averted/symptomatic_months_soc
# # scenario: earlier_three_sputum
# # sim numbers considered:
# results_main$earlier_three_sputum %>% 
#   reframe(symptomatic_reduction = symptomatic_months_averted/symptomatic_months_soc*100) %>%
#   bind_cols(cohort_params) %>%
#   group_by(decile = ntile(1/symptom_underestimation_factor, 10)) %>%
#   summarize(median_reduction = median(symptomatic_reduction, na.rm=T),
#             lci_reduction = quantile(symptomatic_reduction, 0.025, na.rm=T),
#             uci_reduction = quantile(symptomatic_reduction, 0.975, na.rm=T))
# 


# Figure 2: Example screening plot for single cohort ----

meanparams <- cohort_params[cohort_features$accepted_subclinical == 1, ] %>%
  summarise(across(where(is.numeric), median, na.rm = TRUE)) %>%
  as.list()
plotcohort <- create_cohort(cohort_params = meanparams)
check_subclinical(plotcohort, cohort_params = meanparams)
intervention_parameters <- list(
  coverage = list("phone" = meanparams$coverage_phone, 
                  "home" = meanparams$coverage_phone * meanparams$coverage_home_reduction),
  sensitivity_symptoms = list("home" = meanparams$sensitivity_symptoms_home, 
                              "phone" = meanparams$sensitivity_symptoms_home * meanparams$sensitivity_symptoms_phone_reduction),
  success_sputum = list("home" = meanparams$success_sputum_home, 
                        "phone" = meanparams$success_sputum_home * meanparams$success_sputum_phone_reduction))


testcohort <- apply_intervention(cohort = plotcohort, 
                                 design = screening_design_guidelines, #screening_design_6m_sx,
                                 intervention_parameters = intervention_parameters, cohort_params = meanparams)
check_subclinical(testcohort, cohort_params = meanparams)
large <- plot_screening(cohort = testcohort,
                        screening_design = screening_design_guidelines, #screening_design_6m_sx, 
                        colorfill = TRUE) + ylim(0,1100)
small <- plot_screening(cohort = testcohort %>% #10% sample of rows
                          sample_n(size = nrow(testcohort)/10),
                        screening_design = screening_design_guidelines, #screening_design_6m_sx, 
                        colorfill = TRUE) + ylim(0,110)
small
large  

library(cowplot)
# arrange as two panels side by side, A large labeled as full cohort and B small labeled as random 20% subset
final_plot <- plot_grid(
  large + theme(legend.position = "none") + labs(tag = "A") + ggtitle("Full Cohort (N=10,000)") +
    theme(plot.tag = element_text(size = 16, face = "bold", hjust = -0.1, vjust = 1.5)),
  small + labs(tag = "B") + ggtitle("Random 10% Subset") +
    theme(plot.tag = element_text(size = 16, face = "bold", hjust = -0.1, vjust = 1.5),
          # no y axis label
          axis.title.y = element_blank()),
  ncol = 2,
  rel_widths = c(1, 1)
)
final_plot
pdf(file = "Figure2_screening_example.pdf", width = 12, height = 8)
final_plot
dev.off()
