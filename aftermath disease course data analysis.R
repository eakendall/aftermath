library(tidyverse)
library(survival)
library(survminer)
library(Temporal)
library(pROC)
library(conflicted)
library(binom)
conflicts_prefer(dplyr::select, dplyr::filter)

#### Trial Data #####

# Read in aftermath data; will use the dataset Aye/Gayatri prepared for the cohort model fitting

# TBAM_trialdata1512_EK.csv – Data pulled and shared from RedCap by Gayatri.
# variables_trialdata1512.csv – Variables extracted from TBAM_trialdata1512_EK.csv for intervention modeling.
# variables_trialdata1512_codebook.xlsx – A codebook defining the variables in variables_trialdata1512.csv.
# Variables for TBAM intervention modelling.Rmd – The codes used for extracting variables_trialdata1512.csv from TBAM_trialdata1512_EK.csv.

data <- read.csv("../Data June 2025 from Aye/variables_trialdata1806.csv")
dictionary <- readxl::read_excel("../Data Jan 2025 from Aye/variables_trialdata1512_codebook.xlsx")
head(data)
dictionary %>% print(n=100)
# added, not in dictionary:  tb_related_death (based on the verbal autopsy data), death_date
data %>% count(tb_related_death, end_reason)

#### Incident TB ####

# Cumulative TB incidence, overall:
data %>% count(term_reason, end_reason, tb_related_death)

#** Why are there people who completed 18mo but have end_reason of LFU, or term reason of moved out but end_reason of "completion" or in "follow-up"?
#** Assume some portion of deaths with "Probable" or "Unknown" TB relation are recurrences? 
#*There were 11 known TB deaths and 14 definite non-TB deaths.
#*The proportion of the cohort with uknown/probable TB cause of death is 8 /1076 = 0.74% of cohort,
#*so excluding them from recurrences won't change the estimate much, and counting ~4 of them as recurrences would increase it from 
#*sum(data$end_reason == "TB recurrence") / (nrow(data) - 8) = 0.082 to (sum(data$end_reason == "TB recurrence") + 4) / (nrow(data))= 0.086. 
#*But then I'd have to think about timing. 
#*
#*Let's assume that everone who is marked as a TB death is a recurrence, and that if not diagnosed before death they would have been daignosed on their day of death,
#* And let's exclude everyone who died of an unknown cause from the entire analysis. But count the probables as TB recurrences.

data <- data %>% 
  filter(!(end_reason=="death" & tb_related_death == "Unknown")) %>%
  mutate(
    recurrence = case_when(
      end_reason == "TB recurrence" ~ 1,
      end_reason == "death" & tb_related_death == "Yes" ~ 1,
      end_reason == "death" & tb_related_death == "Probable" ~ 1,
      end_reason == "death" & tb_related_death == "No" ~ 0,
      end_reason == "completion" ~ 0,
      end_reason == "LFU" ~ 0,
      TRUE ~ 0),
    end_days = case_when(
      end_reason == "death" & recurrence == 1 ~ interval(txcompl_date, death_date)/days(1),
      TRUE ~ txcompl_endreason_days
    )
  )

data %>% summarise(mean(recurrence)) # 9.0% cumulative incidence (diagnosed during study, not a precisely defined endpoint)
data %>% filter(recurrence==1) %>% select(txcompl_endreason_days) %>% unlist(.) %>% summary(.) # mean ~6, median ~9mo

# Consider a cumulative hazard adjusting for competing events?
# Would assume that death (if not classified as TB recurrence) and LFU outcomes are independent of recurrence risk
# Use survival analysis to estimate cumulative hazard:
survival_dataset <- data %>% select(record_id, term_reason, end_reason, txcompl_endreason_days, recurrence, end_days, ev_sym_durmax) %>% mutate(
  event = case_when(
    recurrence == 1 ~ 1,
    TRUE ~ 0
  ),
  time = case_when(
    end_reason == "TB recurrence" & !is.na(ev_sym_durmax) ~ end_days - 2*ev_sym_durmax,
    recurrence == 1 ~ end_days - 2 * mean(subset(data, recurrence == 1)$ev_sym_durmax, na.rm=T),
    TRUE ~ end_days)
)
# for time of TB recurrence, subtracting off estimated symptom duration (2x reported, or 2x mean reported for deaths with missing symptom durations) to estimate time of symptom onset

# Estimate cumulative hazard (of TB symptom onset):
survival::survfit(Surv(time, event) ~ 1, data = survival_dataset) %>% summary(fun = "cumhaz") # 14% cumulative hazard? but wide uncertainty. 10% at 622 days is probably a better estimate.
# Estimate cumulative hazard to day 540 (i.e. symptom onset by day 540):
survival_dataset_truncated <- survival_dataset; 
survival_dataset_truncated$time[survival_dataset$time > 540] <- 540
survival_dataset_truncated$event[survival_dataset$time > 540] <- 0
survival::survfit(Surv(time, event) ~ 1, data = survival_dataset_truncated %>% mutate()) %>% summary(fun = "cumhaz") # 9.2% (7.4-11.0) cumulative to day 540.

# timing of recurrence, for those diagnosed by 540 days
data %>% filter(recurrence==1, txcompl_endreason_days <= 540) %>% select(txcompl_endreason_days) %>% unlist(.) %>% summary(.) # mean ~6, median ~9mo
223 - 2*17 # mean symptom onset, of those diagnosed by day 540
178 - 2*17 # median symptom onset, of those diagnosed by day 540
data %>% filter(recurrence==1, txcompl_endreason_days <= 540) %>% 
  mutate(sxonset = pmax(txcompl_endreason_days - 2*ev_sym_durmax, 0)) %>% select(sxonset) %>% unlist(.) %>% summary(.) 
# mean 184, median 146 days
# proportion before 90d? 360d?
data %>% filter(recurrence==1, txcompl_endreason_days <= 540) %>% 
  mutate(sxonset = pmax(txcompl_endreason_days - 2*ev_sym_durmax, 0)) %>% summarise(mean(sxonset<=90, na.rm=T)) 
data %>% filter(recurrence==1, txcompl_endreason_days <= 540) %>% 
  mutate(sxonset = pmax(txcompl_endreason_days - 2*ev_sym_durmax, 0)) %>% summarise(mean(sxonset<=360, na.rm=T)) 
# uncertainty in mean 184 + 2*17*c(-1,1)
# from python script, /usr/local/bin/python3 "/Users/emilykendall/Library/CloudStorage/GoogleDrive-emily.kendall@gmail.com/My Drive/Aftermath Golub India R01/Aftermath recurrence intervent
# ions model/aftermath/gamma fitting.py"
# Optimal shape (k): 1.7975235351446937
# Optimal scale (theta): 102.36305626764042
# Fitted truncated mean: 184.0000027704242
# Fitted truncated median: 145.99999588066345


# Cumulative TB incidence, stratified by sputum positivity:
# data %>% filter(end_reason == "TB recurrence") %>% count(ev_micro_reported, ev_micro_test, ev_TBtype)
data <- data %>% mutate(micropos = case_when(ev_micro_test == "Microbiological confirmation" & ev_TBtype == "Pulmonary TB (PTB)" ~ 1, TRUE ~ 0),
                clindx = case_when(ev_micro_test == "Clinical confirmation" | ev_TBtype == "Extra pulmonary TB (EPTB)" ~ 1 , TRUE ~ 0))
data %>% summarise(mean(micropos)) # 4.4% cumulative incidence 
data %>% summarise(mean(clindx)) # 3.8% cumulative incidence (incl 1% micro+ EP TB)
data %>% summarise(mean(clindx | micropos), sum(clindx | micropos))
data %>% filter(recurrence==1) %>% summarise(mean(micropos), n())
# Note: The above classifies all micro+ EPTB with the clinical PTB, as TB that wouldn't be picked up on sputum screening.

# Compare timing of clinical and micro diagnoses:
ggplot(data %>% filter(end_reason == "TB recurrence"), 
       aes(x = txcompl_endreason_days/30, color = as.factor(micropos))) +
  stat_ecdf() +
  xlab("Months to TB recurrence diagnosis") +
  ylab("Density") +
  ggtitle("Timing of TB recurrence diagnosis by diagnostic method") + 
  # set color legend to have title of "TB type" and levels of "Micro+ pulmonary" for 1 and "other" for 0
  scale_color_discrete(name = "TB type", labels = c("Micro+ pulmonary", "Other"))
# Looks pretty similar.

# compare with cumulative incidence curves -- will need to derive each separately and overlay:
survival_dataset <- data %>% 
  select(record_id, term_reason, end_reason, txcompl_endreason_days, micropos, clindx) %>% 
  mutate(
    event = case_when(
      end_reason == "TB recurrence" ~ 1,
      TRUE ~ 0
    ),
    time = txcompl_endreason_days
)
survival_dataset_truncated <- survival_dataset; 
survival_dataset_truncated$time[survival_dataset$time > 540] <- 540
survival_dataset_truncated$event[survival_dataset$time > 540] <- 0

survival_dataset_truncated$eventtype <- case_when(
  survival_dataset_truncated$micropos == 1 & survival_dataset_truncated$event == 1 ~ 1,
  survival_dataset_truncated$clindx == 1 & survival_dataset_truncated$event == 1 ~ 2,
  TRUE ~ 0
)

# # plot only the cumulative hazard:
# micro = survfit(Surv(time, micropos) ~ 1, data = survival_dataset_truncated)
# ggsurvplot(micro)
# clin = survfit(Surv(time, clindx) ~ 1, data = survival_dataset_truncated)
# ggsurvplot(clin)
# # Combine on the same plot
# fit <- list(M = micro, C = clin)
# ggsurvplot_combine(fit, survival_dataset_truncated, xlim = c(0,540), legend.labs = c("Micro+ TB", "clinical"),
#                    ylim = c(0.8, 1))
# # again looks very similar, although really I should censor micros in the clin dataset and vice versa
# 
# 
# # # compare expected diagnosis dates vs aftermath data, derived later
# expected <- as.data.frame(list("days" = rgamma(1e6, shape = 1.78, scale = 103) +
#                                         rlnorm(1e6, meanlog = 2.75, sdlog = 0.75)*2))
# aftermath <- data %>% filter(recurrence==1) %>% select(txcompl_endreason_days)
# 
# 
# ###* Supplemental figure, timing of routine diagnosis *###
# ggplot(data %>% filter(recurrence == 1),
#        aes(x = end_days/30, color = as.factor(micropos))) +
#   stat_ecdf() +
#   xlab("Months to recurrent TB diagnosis") +
#   ylab("Density") +
#   ggtitle("Timing of recurrent TB diagnosis") +
#   # set color legend to have title of "TB type" and levels of "Micro+ pulmonary" for 1 and "other" for 0
#   stat_ecdf(data = expected/30, aes(x = days, color="Simulated")) +
#   scale_color_discrete(name = "TB type",
#                        labels = c("Micro+ pulmonary (Aftermath)", "Other (Aftermath)", "Simulated")) +
#   theme_minimal() +
#   xlim(0,20)
# 
# ggplot(data %>% filter(recurrence == 1),
#        aes(x = end_days/30, color = as.factor(micropos))) +
#   geom_density() +
#   xlab("Months to TB recurrence diagnosis") +
#   ylab("Density") +
#   ggtitle("Timing of TB recurrence diagnosis") +
#   # set color legend to have title of "TB type" and levels of "Micro+ pulmonary" for 1 and "other" for 0
#   scale_color_discrete(name = "TB type",
#                        labels = c("Micro+ pulmonary (Aftermath)", "Other (Aftermath)", "Simulated")) +
#   geom_density(data = expected/30, aes(x = days, color="Simulated")) +
#   theme_minimal()

  


#### Symptoms and Diagnosis Timing ####

# Estimate symptom durations prior to diagnosis:
data %>% count(ev_sympwith_maxdur)
data %>% filter(micropos==1) %>% summarise(mean(ev_sym_durmax), sd(ev_sym_durmax))
data %>% filter(micropos==1) %>% dplyr::select(ev_sym_durmax) %>% summary(.)
# if removing outliers, doesn't change too much:
data %>% filter(micropos==1, ev_sym_durmax <=50) %>% summarise(mean(ev_sym_durmax), sd(ev_sym_durmax))


## different by mode of diagnosis study vs no? Not really, and def not shortened.
data %>% filter(micropos==1) %>% group_by(ev_possible_dxtbam) %>% summarise(mean(ev_sym_durmax))

## different by time since prior treatment?
data %>% filter(end_reason == "TB recurrence") %>% mutate(within60d = txcompl_endreason_days <= 60, 
  within150d = txcompl_endreason_days <= 150) %>%
  group_by(within60d, within150d) %>% summarise(mean(ev_sym_durmax, na.rm=T), sd(ev_sym_durmax, na.rm=T))
data %>% filter(end_reason == "TB recurrence") %>% mutate(quarter = ceiling(txcompl_endreason_days/180)) %>%
  group_by(quarter) %>% summarise(mean(ev_sym_durmax, na.rm=T), sd(ev_sym_durmax, na.rm=T))
# not a big change, except fo rthe very earliest 

## Reasonable to estimate as poisson-distributed? No will use neg binomial, or (changed dec 23) discretized log-normal to avoid zeros
# sxs_days <- rpois(n=1e3, lambda = 17)
# sxs_days_alt <- rnbinom(n=1e3, size = 1.6, mu = 17) # sd =sqrt(mu + mu^2/size) --> when sd = 14 and mu = 17, size = 1.6
# could set size = 3 based on visual comparison of distributions below. corresponds to lower sd of sqrt(17 + 17^2/3) = 10.6 ~ 11.
data %>% filter(end_reason == "TB recurrence", txcompl_endreason_days <= 180) %>% summarise(mean(log(ev_sym_durmax), na.rm=T), sd(log(ev_sym_durmax), na.rm=T))
sxs_days_alt2 <- rlnorm(n=1e3, meanlog = 2.68, sdlog = 0.69)
mean(sxs_days_alt2); sd(sxs_days_alt2)
sd(rlnorm(n=1e3, meanlog = 2.68, sdlog = 0.6))
sd(rlnorm(n=1e3, meanlog = 2.68, sdlog = 0.8))

# from plot below, the log normal seems to fit the data better (and approporately avoids zeros.)
# need to also increase mean and variance with time since diagnosis
# and could vary over time since initial treatmnt completion (t^1/4 for mean and variance, i.e. t^1/2 for sd)
times <- c(30, 120, 300) # illustrative time points. but ctually should define relative to mean time to sx onset (which is ~174d)
meanscale <- (times/174)^0.25
meanscale * 17
sdscale <- (times/174)^sqrt(0.25)
sdscale * 14

#** will need to confirm that simulated data symptoms have the correct mean and SD

ggplot() + #geom_density(aes(x=t), data = data.frame(t = sxs_days)) + 
  # geom_density(aes(x=t), data = data.frame(t = sxs_days_alt), col="blue") + #nbinom
  geom_density(aes(x=t), data = data.frame(t = sxs_days_alt2), col="green") + #lognorm
  geom_density(data = data %>% filter(micropos==1), aes(x = ev_sym_durmax), col="red") + 
  #legend for line color:
  xlab("Symptom duration (days)") +
  ylab("Density") +
  ggtitle("Symptom duration prior to micro+ TB recurrence diagnosis") +
  xlim(0,50) + 
  theme_minimal() + 
  annotate("text", x = 30, y = 0.03, label = "Red: Aftermath data\nBlack: Poisson (mean=17, sd=4.1)\nBlue: Neg binom", hjust=0)
# better fit is wider than the poisson, and a bit narrower than but closer to the nbinom with same sd




# Fit gamma function to timing of micro+ TB recurrence symptom onset:
# assuming symptom onset = diagnosis date - 2 * reported symptom duration (and varying this to get uncertainty range)
data %>% filter(micropos == 1) %>% 
  summarise(mean(txcompl_endreason_days - 2*ev_sym_durmax), sd(txcompl_endreason_days - 2*ev_sym_durmax))
data %>% filter(micropos == 1) %>% 
  summarise(mean(txcompl_endreason_days - 0*ev_sym_durmax), sd(txcompl_endreason_days - 0*ev_sym_durmax))
data %>% filter(micropos == 1) %>% 
  summarise(mean(txcompl_endreason_days - 4*ev_sym_durmax), sd(txcompl_endreason_days - 4*ev_sym_durmax))

# But this SD gives us a really long tail, and a higher median (among those diagnosed by 18 months) than observed in aftermath.
# Aftermath median:
data %>% filter(micropos == 1, txcompl_endreason_days <= 540) %>% 
  summarise(median(txcompl_endreason_days - 2*ev_sym_durmax))

# And SD and mean when limited to those diagnosed by 540 (maybe we should use this 146 for our SD):
data %>% filter(micropos == 1, txcompl_endreason_days <= 540) %>% 
  summarise(mean(txcompl_endreason_days - 2*ev_sym_durmax), sd(txcompl_endreason_days - 2*ev_sym_durmax))
data %>% filter(micropos == 1, txcompl_endreason_days <= 540) %>% 
  summarise(mean(txcompl_endreason_days - 0*ev_sym_durmax), sd(txcompl_endreason_days - 0*ev_sym_durmax))
data %>% filter(micropos == 1, txcompl_endreason_days <= 540) %>% 
  summarise(mean(txcompl_endreason_days - 4*ev_sym_durmax), sd(txcompl_endreason_days - 4*ev_sym_durmax))

# But then I'll need the gamma parameter fitting to look only at those with onset by 540 - 2*17...
# ( I used a python script to estimate gamma parameters here, gamma fitting.py)
# shape 1.6859068628739582
# scale 113.88529634851781
# (Considered updating to include background incidnce rate of 0.4%/y, with updated python script gamma fitting with reinfection incidence.)
# Shape (k): 0.5711
# Scale (theta): 7638.3282
# (but decided the time horizon over which to consider impact complicates things too much - the longer the horizon, the less ACF can reduce recurrence iimpact as a % of total recTB time)

# # * Should I also model a shortening of diagnostic delay for those diagnosed by aftermath (e.g. within 2 weeks after a study visit)? 
# # * No, don't see evidecne that they had shorter symptom durations, or that most reported symptoms that would have led to detection by the study
# # data %>% mutate(studydx = case_when(txcompl_endreason_days - txcompl_m6visit_days %in% 0:14 & !is.na(m6_sym_durmax) ~ 1,
# #                                     txcompl_endreason_days - txcompl_m12visit_days %in% 0:14 & !is.na(m12_sym_durmax) ~ 1,
# #                                     txcompl_endreason_days - txcompl_m18visit_days %in% 0:14 & !is.na(m18_sym_durmax) ~ 1,
# data %>% mutate(studydx = case_when(recurrence == 1 & (txcompl_endreason_days - txcompl_m6visit_days) %in% 0:30 ~ 1,
#                                     recurrence == 1 & (txcompl_endreason_days - txcompl_m12visit_days) %in% 0:30 ~ 1,
#                                     recurrence == 1 & (txcompl_endreason_days - txcompl_m18visit_days) %in% 0:30 ~ 1,                                    
#                                     recurrence == 1 ~ 0,
#                                     TRUE ~ NA)) %>% 
#   group_by(studydx) %>%
#   summarise(median(ev_sym_durmax, na.rm=T), mean(ev_sym_durmax, na.rm=T), sd(ev_sym_durmax, na.rm=T), sum(!is.na(ev_sym_durmax)), n(), mean(!is.na(m6_sym_durmax) | !is.na(m12_sym_durmax) | !is.na(m12_sym_durmax)))


# # old plot, to diagnosis rather than symptom onset time to confirm that gamma looks reasonable:
# 
# scale <- v/m
# shape <- m*m/v
# 
# ggplot() + 
#   geom_density(data = data %>% filter(micropos==1), aes(x = txcompl_endreason_days - 2*ev_sym_durmax)) +
#   geom_density(data = data.frame(t = rgamma(n=1e5, scale = scale, shape = shape)), aes(x=t), col="red")
# # Looks like a reasonable fit. 
# # or we could exclude the diagnoses that were likely made by the study -- but this resulted in even shorter duration. 
# # So we'll stick with the raw estimates abvove.

# And about that possible small uptick in diagnoses at 12 months: there were only 4 (out of 10 diagnoses in that period) that were possibly diagnosed by the study, all in the home arm.
# At six months, 4 out of 6 were in the home arm. At 18 months arm is irrelevant. 
# but most weren't symptomatic at their visit, and reported symptom durations at the even visit were longer than average (not consistent with early detection through aCF)
  

# Numerically solved for the scale and shape parameters of a gamma distribution, such that m and v are the mean and variance when the disribution is truncated at 540.



# Was the rate of diagnosis higher in the visit windows? 
# There are so few Yes's that it's hard to make anything of this, but I don't see much signal.
# data %>% filter(end_reason == "TB recurrence") %>% group_by(floor(txcompl_endreason_days/180)) %>% summarise(mean(ev_possible_dxtbam=="Yes"), n())
data %>% filter(end_reason == "TB recurrence") %>% group_by(floor(txcompl_endreason_days/30)) %>% summarise(mean(ev_possible_dxtbam=="Yes"), n()) %>% print(n=20)
# So, the months with possible study-driven diagnosis are 6-8, 12-14 (and we'll exclude 18+ from all analyses)
data %>% filter(end_reason == "TB recurrence") %>% group_by(floor(txcompl_endreason_days/90)) %>% summarise(sum(end_reason == "TB recurrence"))
# comparing periods 2 and 4 vs 0, 1, 3, and 5, I don't see a difference, apart from maybe very slight inclrease in 4 (12mo fu visit period) 


### Symptom reporting in home vs phone screening
data %>% group_by(arm) %>% 
  pivot_longer(cols = ends_with("sympwith_maxdur"), 
               names_to = "month",
               values_to = "symdur",
               names_pattern = "m(\\d+)_sympwith_maxdur") %>%
  mutate(symptoms = !is.na(symdur)) %>%
  select(symdur, symptoms, arm, month) %>%
  summarise(mean(symptoms))

.0342/.0480

#### Subclinical period ####

# o	Proportion of Symptomatic sputum+ pulmonary TB who become sputum Xpert+ before symptoms
# 	Triangulate from:
  # •	% of prevalent Xpert+ TB that is asymptomatic in cross-sectional surveys 
      # (those with an asymptomatic period are more likely to be prevalent at the time 
      # of a survey than those with immediate symptoms, 
      # but even some who are symptomatic in surveys may have started as asymptomatic; 
      # also symptom-screen-negative doesn’t mean fully asymptomatic) 
      
  # •	1 – (% of symptomatic care-seeking pulmonary TB that is Xpert-negative and culture-psoitive i.e. ~10%): 
      # These definitely have symptoms before becoming xpert+, 
      # so 1 – (this %) sets an upper bound of 90% with a subclinical period. 
      # Could lower that upper bound even more with an estimate of true culture-negative pulmonary TB diagnoses (another 5%  upper bound of 85%?). 
  #  •	Report paper % of diagnoses that were subclinical (specific to recurrent TB) 
      # as a true lower bound on % subclinical at some point. = 28/66 = 42%


  ## So assume that between 42 and 85% have a subclnical period. 
  ## Let's say 3/4 for now. 

# o	Duration of sputum+ before symptom onset: 
  # 	Define mean relative to estimated duration of symptoms (above): 
  # •	Proportion of sputum+ recurrent TB time spent with symptoms: 
      # 0.588 (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0294254 Table 1) --
        # but these 381 recurrent TB include 190 already on treatment.
            # Of those, there were 99 with “symptoms, past TB, on ATT, Xray abnormal” 
            # and 43 with “Past TB, on ATT, xray abnormal [i.e. no symptoms]” (Fig 1) 
            # leaving 48 with symptom status unspecified. 
        # Of the *untreated* recurrent TB, somewhere between (157 -43 ) /191 = 60%  and (157 -43-48 ) /191 =  35% was symptom-screen negative. 
        # That's consistent with the ~50% asymptomatic desribed among overall preavlence. 
        
        ## So let's assume mean subclinical duration of 1x the symptom duration. 
        # That's short! (17d)
        # But we could assume symptoms (as in positive symptom screen) had lasted longer than people reported...


# 	Model distribution of durations as truncated (at prior tx completion) poisson distribution, 
# with variance chosen to require little truncation and explored in sensitivity analysis. 


# How many should be cx+ subclinical from the end of their initial treatment? 
# I looked at the CTriumph/TBDM data for patients who were culture positive at the end of treatment, and I found 4 out of ~1300 with cultures near that time point. Presuming that those 4 would have become recurrences later if they hadn’t been cultured.
# But actually I did this wrong because the dataset shared with me is only those followed for recurrence??
reportdata <- readxl::read_xlsx("../../TBDM CTRIUMPH recurrent TB data/TBDM_CT_Dataset_9Jun.xlsx", sheet = 1, col_names = T, trim_ws = T, 
                          col_types = c("text", "text", "date", "text","text",
                                        "numeric","text","numeric", "text", "numeric",
                                        "text", "numeric", "text", "text","text",
                                        "text","text", "numeric", "text", "numeric",
                                        "numeric","numeric", "numeric","numeric", "text",
                                        "numeric","text","text","date","text",
                                        "text","text","numeric","text","numeric",
                                        "numeric","numeric","numeric", "text","text",
                                        "date","date", "date","date", "date"),
                          na = c("NA","", "Missing"))
reportdata %>% ungroup() %>% count(`External ID`)
reportdata %>% filter(Event %in% c("EOT", "Month 6", "Week 24", "Month 5", "Week 20")) %>% group_by(`External ID`) %>%
  count(cxpos = `LJ Result` == "Positive" | `MGIT Result` == "Positive", Event) %>% filter(cxpos)
reportdata %>% filter(Event %in% c("EOT", "Month 6", "Week 24")) %>% group_by(`External ID`) %>%
  count(cxpos = `LJ Result` == "Positive" | `MGIT Result` == "Positive") %>% filter(cxpos) %>% summarise(n())
reportdata %>% group_by(`External ID`) %>% count() %>% nrow()
reportdata %>%  filter(Foll_Days > 0)  %>% count(`External ID`) %>% nrow()
reportdata %>% filter(Foll_Days > 0)  %>% filter(`Clinical/subclinical` %in% c("Clinical","Subclinical")) %>% count(`External ID`) %>% nrow()

# So, 4/861 is 0.46% of the cohort, and 4/(68+4) = 5.5%  of recurrences (not same people becaeu they were detcted and retreated).
# and 9/865 is 1% of the cohort and 9/(68+9) = 12% of all recurrences were already subclinical sputum+ at the end of treatment.
# These are positive cultures, so the proportion detectable by Xpert screening at the end of treatment would be fewer,
 # say 3-7 of 865 or 0.3-0.8% of the cohort, and 4-9% of recurrences.

# And how many should be subclinical in a slice at 6m post-treatment?
reportdata %>% count(Event) %>% print(n=100)
reportdata %>% filter(Foll_Days > 0) %>% group_by(`External ID`) %>% count(n())
reportdata  %>% filter(Event %in% c("Week 48", "Month 12"), `Visit Date`== `Outcome Date`) %>% group_by(`External ID`) %>% 
  filter(`Clinical/subclinical` %in% c("Clinical","Subclinical")) %>% ungroup() %>% count(`Clinical/subclinical`, Event)
reportdata  %>% filter(Event %in% c("Week 72", "Month 18"), `Visit Date`== `Outcome Date`) %>% group_by(`External ID`) %>% 
  filter(`Clinical/subclinical` %in% c("Clinical","Subclinical")) %>% ungroup() %>% count(`Clinical/subclinical`, Event)
# So in a cohort of 861 followed for recurrence, there were 11 subclinical at the 12mo visit where eveyone got sputum. (Some not Xpert+ though.)
binom.agresti.coull(11, 861)[c("lower", "upper")] * 0.7 # 0.7 is to estimate those that would be Xpert+


#### Symptoms underreporting ####

# ~3x difference bewteen durations reported in clinical cohorts (https://bmcpublichealth.biomedcentral.com/articles/10.1186/s12889-019-7026-4)
(81 + 29.5 + 7.9 )/30 ; (70 + 26 + 7)/30; (92+33+9)/30
# vs in preavlence surveys (https://bmcmedicine-biomedcentral-com.proxy1.library.jhu.edu/articles/10.1186/s12916-021-02128-9)
total_duration <- c(36, 22, 22, 18, 19, 21, 22, 16, 13, 13, 9)
asx_duration <- c(14, 14, 13, 9, 9, 6, 6, 5, 5, 5, 4)
sx_duration <- total_duration - asx_duration
summary(sx_duration)

# Could also consier (as more indirect support for at least 2x underreporting) the ~2x difference in peavlence between any cough vs 2 weeks cough in Stuck et al 2024, 
# and the ~2x difference in duration between symptoms+ and care-seeking in Ku et al. 


####### Aftermath Costs #########

# Top-down estimates (from Shreeya's union slides):
# Home visit arm, $ 6322 for 539 HHs ($11.73 per HH for 3 visits --> $3.91 per home visit) 
# Telephonic arm, $ 3132 for 537 HHs, but that includes the same 18m home visit as home arm.
# Difference = ($(537 * 11.73) - 3132)/2/537 = $2.9 as incremental cost per home vs phone visit --> 
# cost $1.01 per phone vist.
# (Note, ignoring costs of contact screening as negligible and/or somethign that would occur in our modeled intervention too)

# Bottom-up: The cost per follow-up visit at 6 and 12 months in the 
# home-based arm was $5.49 and $4.86, respectively, 
#  compared to $0.42 and $0.40 in the phone arm. 


