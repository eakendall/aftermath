library(tidyverse)

# Simulate a TB cohort and their time to diagnosis

# Data that will eventually come from Aftermath, but for now I'm using a jittered gamma to generate simluated data:
dx_timing_months <- jitter(rgamma(n=50, shape=2, scale = 3), 0.1)
plot(ecdf(dx_timing_months), main="Time to Diagnosis", xlab="Months", ylab="Cumulative Probability")

sx_months_at_dx <- rnbinom(n=50, size=2, mu=12)/4
dx_at_study_visit <- rbinom(n=50, size=1, prob=0.3)
hist(sx_months_at_dx)
# (To test: Is symptom duration shorter for those diagnosed at study visits or through study procedures? # nolint
t.test(sx_months_at_dx[dx_at_study_visit == 1],
       sx_months_at_dx[dx_at_study_visit == 0])

# curve of % of recurrences vs % of total pop you could target.
prediction_curve <- data.frame(
  "target_pop" = c(0.01, 0.05, 0.1, 0.2, 0.4, 1),
  "target_risk_proportion" = c(0.1, 0.25, 0.4, 0.6, 0.8, 1))
ggplot(prediction_curve, aes(x = target_pop, y = target_risk_proportion)) +
  geom_line() +
  xlab("Proportion of population targeted") + 
  ylab("Proportion of recurrences prevented")

# schedule of study visits (probability that they occur on a given day, can use actual intervals)
# will be based on actual timing of completed visits
study_schedule_1 <- rnorm(300, mean = 6 * 30, sd = 14)
study_schedule_2 <- rnorm(300, mean = 12 * 30, sd = 14)
study_schedule_3 <- rnorm(300, mean = 18 * 30, sd = 14)
study_schedule <- c(study_schedule_1, study_schedule_2, study_schedule_3)
hist(study_schedule)

# External data to use in parameter estimation:
proportion_symptomatic <- 0.478 + 0.11
# Will assumes that proportion of time spent asymptomatic is the same in our cohort as among prevalent TB in previously treated people (regardless of time since earlier treatment) #nolint
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0294254 Table 1


# Parameter estimation from above data:
symptoms_mean_months <- mean(sx_months_at_dx[dx_at_study_visit==0])
# symptoms_sd_months <- sd(sx_months_at_dx[dx_at_study_visit==0])
presymptomatic_mean_months <- function(symptoms_mean_months,  prop_sx = proportion_symptomatic) {
  return(symptoms_mean_months * (1 - prop_sx) / prop_sx)
}
presymptomatic_mean_months(symptoms_mean_months)

# Model of processes that generated the dx timing data:
library(distr)

# 1. Onset of TB is a gamma distribution with unknown shape and rate:
tb_onset <- function(shape_onset, scale_onset) {
  Gammad(shape = shape_onset, scale = scale_onset)
}
# 2. After TB onset, onset of symptoms is poisson distributed with mean related to natural mean symptom duration based on prevalence survey symptomatic proportion:  # nolint
symptom_onset_after_tb <- function(presx_months) {
  Pois(lambda = 1 / presx_months)
}

# 3. After symptom onset, routine care seeking and diagnosis also occur at a poisson-distributed rate, with mean 3 mo baesd on external systamatic review data? Or, for now, mean based on study data (will assume but check mean=variance)? # nolint
routine_dx_after_symptoms <- function(sx_months) {
  Pois(lambda = 1 / sx_months)
}
# study visits will be added as competing events that also would lead to dx if they occur first. Maybe with some <1 probability? # nolint

# If we wanted to use MLE, the likelihood of the data given the model would be
# the product of the likelihood of each individual's time to diagnosis given the model. 
# And for each individual, the probability distribution for their timing of diagnosis would be 
# the sum of probabilities for routine and study-visit-based diagnoses.
# Both of these probabilities of diagnosis would depend on having symptoms;
# the latter would also depend on receiving a study visit while symptomatic and before routine diagnosis.
# The density of a sum of distributions is the convolution of the densities of the individual distributions, so the distribution of timing of symptom onset is:
symptom_onset_dist <- function(shape_onset, rate_onset, presx_months) {
  conv <- convpow(tb_onset(shape_onset, rate_onset) +
                    symptom_onset_after_tb(presx_months), 1)
  return(conv)
}

routine_dx_dist <- function(shape_onset, rate_onset, sx_months) {
  conv <- convpow(symptom_onset_dist(shape_onset, rate_onset, presymptomatic_mean_months(sx_months)) +
                  routine_dx_after_symptoms(sx_months), 1)
  return(conv)
}

# study dx occurs if there is a study visit between symptom onset and routine dx



shape_test <- 2; scale_test <- 2; sx_test <- 3
x <- seq(0, 20, 0.1)
# make a dataframe with the densities of all of the above dist functions, at x
dists_test <- data.frame(
  "x" = x,
  "tb_onset" = d(tb_onset(shape_test, scale_test))(x),
  "symptom_onset" = d(symptom_onset_dist(shape_test, scale_test, presymptomatic_mean_months(sx_test)))(x),
  "routine_dx" = d(routine_dx_dist(shape_test, scale_test, sx_test))(x)
  ) %>%
  pivot_longer(cols = c("tb_onset", "symptom_onset", "routine_dx"), names_to = "dist", values_to = "density")
ggplot(dists_test, aes(x = x, y = density, col = dist)) + 
  geom_line() + 
  xlab("Months") + ylab("Density") + 
  geom_histogram(data=data.frame(study_schedule), aes(x=study_schedule/30, y = after_stat(density)), color="black", fill="blue", alpha=0.5)

# Now we can use MLE to estimate the shape and scale of the gamma distribution of TB onset, given the dx_timing_months data.
# (For now, assume mean symptom duration and presymptomatic period are known.)
# The likelihood of the data given the model would be the product of the likelihood of each individual's time to diagnosis given the model.
# And for each individual, the probability distribution for their timing of diagnosis would be the sum of probabilities for routine and study-visit-based diagnoses.
# For a given shape and scale of the gamma distribution of TB onset, we can calculate the likelihood of each individual's time to diagnosis given the model.
likelihood_dx_timing <- function(dx_data = dx_timing_months, shape_onset, rate_onset) {
  n <- length(dx_data)
  onset_timing_months <- rgamma(n=n, shape=shape_onset, rate=rate_onset)
  dx_timing_simulated <- onset_timing_months + 
    rnorm(n, mean=presymptomatic_mean_months(symptoms_mean_months), sd=presymptomatic_sd_months) + 
    rnorm(n, mean=symptoms_mean_months, sd=symptoms_sd_months)
  # compare dx_timing_simulated to dx_timing_months
  return(sum((dx_timing_simulated - dx_timing_months)^2))
}



##### OR ####

# Assume the timing of onset is a gamma distribution with unknown shape and rate, onset_timing_months(shape_onset, rate_onset). 
# And the timing of diagnosis is equal to onset_timing_months + rnorm(1, mean=presymptomatic_mean_months, sd=presymptomatic_sd_months) + rnorm(1, mean=symptoms_mean_months, sd=symptoms_sd_months)
# Fit the resulting simulation to the dx_timing_months data using MLE
dx_timing_sse <- function(dx_data = dx_timing_months, par = list(shape_onset=1, rate_onset=1)) {
    n <- length(dx_data)
    onset_timing_months <- rgamma(n=n, shape=par[[1]], rate=par[[2]])
    dx_timing_simulated <- onset_timing_months + 
        rnorm(n, mean=presymptomatic_mean_months, sd=presymptomatic_sd_months) + 
        rnorm(n, mean=symptoms_mean_months, sd=symptoms_sd_months)
    # compare dx_timing_simulated to dx_timing_months
    return(sum((dx_timing_simulated - dx_timing_months)^2))
}
# find the MLE of shape_onset and rate_onset
fit <- optim(
        c(0.5, 50),
        fn = dx_timing_sse,
        method = "L-BFGS-B",
        lower = c(0, 0),
        upper = c(1000, 1000))

