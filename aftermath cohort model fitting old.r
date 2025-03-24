library(tidyverse)
library(distr)
library(stats4)

## to add: potential for symptom+ sputum- (clinically diagnosable) - can base on % sputum+ at dx, if we know it
#  and different probs of detection at visits by Aftermath arm
# And am I representing cumulative incidecne sensibly?
# And what to do about deaths?

# Simulate a TB cohort and their time to diagnosis

# Data that will eventually come from Aftermath, but for now I'm using a jittered gamma to generate simluated data:
N_participants <- 500
N_events <- 50
dx_timing_days <- jitter(rgamma(n=N_events, shape=2, scale = 50), 0.1)
plot(ecdf(dx_timing_days), main="Time to Diagnosis", xlab="Days", ylab="Cumulative Probability")

sx_days_at_dx <- rnbinom(n=50, size=2, mu=90)
dx_at_study_visit <- rbinom(n=50, size=1, prob=0.3)
hist(sx_days_at_dx)
# (To test: Is symptom duration shorter for those diagnosed at study visits or through study procedures? # nolint
t.test(sx_days_at_dx[dx_at_study_visit == 1],
       sx_days_at_dx[dx_at_study_visit == 0])

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
# Assumes that proportion of time spent asymptomatic is the same in our cohort as among prevalent TB in previously treated people (regardless of time since earlier treatment) #nolint
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0294254 Table 1


# Parameter estimation from above data:
symptoms_mean_days <- mean(sx_days_at_dx[dx_at_study_visit==0])
# symptoms_sd_months <- sd(sx_months_at_dx[dx_at_study_visit==0])
presymptomatic_mean_days <- function(symptoms_mean_days,  prop_sx = proportion_symptomatic) {
  return(symptoms_mean_days * (1 - prop_sx) / prop_sx)
}
presymptomatic_mean_days(symptoms_mean_days)


### Model of processes that generated the dx timing data: ###

# 1. Onset of TB is a gamma distribution with unknown shape and rate:
tb_onset <- function(shape_onset, scale_onset) {
  Gammad(shape = shape_onset, scale = scale_onset)
}
# 2. After TB onset, onset of symptoms is poisson distributed with mean related to natural mean symptom duration based on prevalence survey symptomatic proportion:  # nolint
symptom_onset_after_tb <- function(presx_days) {
  Pois(lambda = presx_days)
}

# 3. After symptom onset, routine care seeking and diagnosis also occur at a poisson-distributed rate, with mean 3 mo baesd on external systamatic review data? Or, for now, mean based on study data (will assume but check mean=variance)? # nolint
routine_dx_after_symptoms <- function(sx_days) {
  Pois(lambda = sx_days)
}
# study visits will be added as competing events that also would lead to dx if they occur first. Maybe with some <1 probability? # nolint

# If we wanted to use MLE, the likelihood of the data given the model would be
# the product of the likelihood of each individual's time to diagnosis given the model. 
# And for each individual, the probability distribution for their timing of diagnosis would be 
# the sum of probabilities for routine and study-visit-based diagnoses.
# Both of these probabilities of diagnosis would depend on having symptoms;
# the latter would also depend on receiving a study visit while symptomatic and before routine diagnosis.
# The density of a sum of distributions is the convolution of the densities of the individual distributions, so the distribution of timing of symptom onset is:
symptom_onset_dist <- function(shape_onset, rate_onset, presx_days) {
  conv <- convpow(tb_onset(shape_onset, rate_onset) +
                    symptom_onset_after_tb(presx_days), 1)
  return(conv)
}

routine_dx_dist <- function(shape_onset, rate_onset, sx_days) {
  conv <- convpow(symptom_onset_dist(shape_onset, rate_onset, presymptomatic_mean_days(sx_days)) +
                  routine_dx_after_symptoms(sx_days), 1)
  return(conv)
}

# study dx occurs if there is a study visit between symptom onset and routine dx


#### plot distributions ####

shape_test <- 3; scale_test <- 40; sx_test <- 70
x <- 0:500
# make a dataframe with the densities of all of the above dist functions, at x
dists_test <- data.frame(
  "x" = x,
  "tb_onset" = d(tb_onset(shape_test, scale_test))(x),
  "symptom_onset" = d(symptom_onset_dist(shape_test, scale_test, presymptomatic_mean_days(sx_test)))(x),
  "routine_dx" = d(routine_dx_dist(shape_test, scale_test, sx_test))(x)
  ) %>%
  pivot_longer(cols = c("tb_onset", "symptom_onset", "routine_dx"), names_to = "dist", values_to = "density")

ggplot(dists_test, aes(x = x, y = density, col = dist)) +
  geom_line() + 
  xlab("Days") + ylab("Density") + 
  geom_histogram(data=data.frame(study_schedule), aes(x=study_schedule, y = after_stat(density)), color="black", fill="blue", alpha=0.5)


#### MLE ####
# Now we can use MLE to estimate the shape and scale of the gamma distribution of TB onset, given the dx_timing_months data.
# (For now, assume mean symptom duration and presymptomatic period are known.)

# The probability of diagnosing an individual on day d is a sum of the probabilities that:
# a) they are routinely diagnosed, or 
# b) they have symptoms, haven't yet been routinely diagnosed, and get a study visit on day d. 
prob_dx_day <- function(day, shape_onset, scale_onset, sx_days,
                          study_schedule_1, study_schedule_2, study_schedule_3) {
  return(
    d(routine_dx_dist(shape_onset, scale_onset, sx_days))(day) + # prob of routine dx on day d 
    d(convpow(symptom_onset_dist(shape_onset, scale_onset, presymptomatic_mean_days(sx_days)) -
            routine_dx_after_symptoms(sx_days), 1))(day/30) * # prob of sx but not yet diagnosed...
      (mean(study_schedule_1 == day) +
       mean(study_schedule_2 == day) +
       mean(study_schedule_3 == day)) # ... and a study visit occurs
  )
}

# #testing
# sapply(150:200, function(x) prob_dx_day(x, shape_test, scale_test, sx_test, 
#                                         study_schedule_1, study_schedule_2, study_schedule_3))

# DEFINE LIKELIHOOD FUNCTION

# The likelihood of the data given the model would be the product of the likelihood of each individual's time to diagnosis given the model.
# This should include those who weren't diagnosed
minus_log_likelihood_dx_timing <- function(shape, scale, # to be estimated
                                     dx_timing = dx_timing_days, # observed data
                                     N_participants =  N_participants,
                                     sx_days = symptoms_mean_days, # known
                                     s1 = study_schedule_1,
                                     s2 = study_schedule_2,
                                     s3 = study_schedule_3) {
  # those who complete follow up without dx (Can modify this to account for LTFU **)
  nodx_ll <- log(1 - sum(sapply(1 : (18 * 30),
                  function(x) prob_dx_day(x, shape, scale, sx_days,
                                          s1, s2, s3)))) *
                              (N_participants - length(dx_timing))
  # those with dx events
  dx_ll <- sum(log(pmax(sapply(dx_timing,
             function(x) prob_dx_day(x, shape, scale, sx_days,
                                     s1, s2, s3)), 1e-100)))
  return(-(nodx_ll + dx_ll))
}

minus_log_likelihood_dx_timing(shape = 2, scale = 36,
                        N_participants =  N_participants,
                        sx_days = symptoms_mean_days,
                        s1 = study_schedule_1,
                        s2 = study_schedule_2,
                        s3 = study_schedule_3)

# FIND BEST FITTIG PARAMETERS
fit <- mle(minuslogl = minus_log_likelihood_dx_timing,
           start = list(shape = 2, scale = 36),
           method = "L-BFGS-B",
           fixed = list("N_participants" = N_participants,
                     "sx_days" = symptoms_mean_days,
                     "s1" = study_schedule_1,
                     "s2" = study_schedule_2,
                     "s3" = study_schedule_3),
           lower = list(shape = 0,
                        scale = 0),
           upper = list(shape = Inf,
                        scale = Inf),
           nobs = N_participants)
