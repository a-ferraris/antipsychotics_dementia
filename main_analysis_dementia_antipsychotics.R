###############################################################################
# Project Title: Antipsychotic Use and Mortality in Patients with Dementia
# Script Author: AF
# Contact: aferra@uw.edu
# Institutions: Hospital Italiano de Buenos Aires (HIBA); University of Washington (UW)
# Version: February 2025
#
# Objective: Following PLOS Mental Health policies, we published a version of the code used 
# that  is intended to show the key features of the analysis implemented in the manuscript. 
# This code should be run usingthe "long_mock" and "baseline_mock" files.
# The results will approximate but  not exactly reproduce the findings of the manuscript 
#
# Study Summary:
# This script implements the full analytic pipeline for estimating the 
# time-varying association between antipsychotic initiation and all-cause 
# mortality among patients with dementia. The workflow includes:
#   - Preparation of baseline and longitudinal datasets
#   - Multiple imputation of missing baseline covariates using MICE
#   - Construction of discrete 3-month time intervals for time-varying covariates
#   - Identification and splitting of intervals at the exact date of 
#     antipsychotic initiation
#   - Creation of counting-process style survival data (tstart, tstop)
#   - Fitting of time-varying Cox proportional hazards models
#   - Pooling of model estimates across imputed datasets
#
# Time Scale:
#   - Analyses use days since cohort entry as the underlying time scale.
#
# Exposure Definition:
#   - Antipsychotic initiation is defined by the first prescription date.
#   - Exposure is modeled as a time-varying "ever exposed" indicator.
#
# Outcome Definition:
#   - All-cause mortality during follow-up.
#   - Death is coded at the final interval for each individual.
#
# Follow-up:
#   - Default follow-up window is 5 years (365 * 5 days), with sensitivity 
#     analyses possible by modifying the 'follow_up' parameter.
#
# Software and Package Versions:
# Analyses were conducted in R using the following packages and versions:
#   - tidyverse (2.0.0)
#   - reshape2 (1.4.4)
#   - survival (3.7-0)
#   - rms (6.8-2)
#   - mice (3.19.0)
#   - survminer (0.4.9)
#   - purrr (1.0.2)
#   - splines (4.4.1)
#   - lubridate (1.9.3)
#
# Reproducibility:
#   - Users should update file paths to match their local environment.
#   - For reproducible environments, consider using renv or sessionInfo().
###############################################################################

# code starts. 

if (!require("tidyverse")) install.packages("tidyverse"); library(tidyverse)
if (!require("reshape2")) install.packages("reshape2"); library(reshape2)
if (!require("survival")) install.packages("survival"); library(survival)
if (!require("rms")) install.packages("rms"); library(rms)
if (!require("mice")) install.packages("mice"); library(mice)
if (!require("survminer")) install.packages("survminer"); library(survminer)
if (!require("purrr")) install.packages("purrr"); library(purrr)
if (!require("splines")) install.packages("splines"); library(splines)
if (!require("lubridate")) install.packages("lubridate"); library(lubridate)

# mice parameters
n_imputations <- 10
n_iterations <- 50

follow_up <- 365*5 # change this for sensitivity analyses 

# load data
data_dir <- "path/to/your/data"   # <-- users edit this

long <- read.csv(file.path(data_dir, "long_mock.csv"))
baseline_data <- read.csv(file.path(data_dir, "baseline_mock.csv"))

####################################
# important definitions for MICE 
analysis_vars <-  c("id", "time",
                    "bdi_score", "mmse_score", "education", "age", "female", 
                           "parkinson", "dementia_lewy", "tobacco", "epilepsy", "asthma", 
                           "hypertension", "ckd", "diabetes", "insomnia", "depression", 
                           "anxiety", "peripheric_vascular", "cardiac_failure", "coronary", 
                           "opioids", "antidepressants", "benzodiazepines", "z_drugs", 
                           "antiepileptics")

baseline_vars <-  c("id", "time",
                    "cancer", "diabetes","copd", "insomnia",  
                    "opioids", "antidepressants", "benzodiazepines", "z_drugs", 
                    "antiepileptics")

# Identify variables to exclude from imputation
exclude_vars <- c("id", "time","death_date", "date",  
                  "disaffiliation_date", "drug", "milligrams",  "first_prescription_date", 
                  "time_to_prescription", "death", "time_to_event", "date_end_fup", 
                  "time_baseline", "time_end_fup", "time_start_ap")

# Impute in wide format
baseline_data <- merge(baseline_data, 
                       long[, baseline_vars], # adding variables only to MICE 
                       by = c("id", "time"), 
                       all.x = TRUE
                       )

# Step 1: Create the predictor matrix
pred_matrix <- make.predictorMatrix(baseline_data)

# Step 2: Set columns to 0 for variables ot exclude from mice process

# Prevent these variables from being imputed
pred_matrix[exclude_vars, ] <- 0

# Prevent these variables from being used to impute others
pred_matrix[, exclude_vars] <- 0

method_vec <- rep("rf", ncol(baseline_data))
names(method_vec) <- names(baseline_data)
method_vec[exclude_vars] <- ""  # do not impute these

# Convert appropriate variables to factors for mice
factor_vars <- setdiff(names(baseline_data), c("id", "bdi_score", "mmse_score", "age", 
                                               exclude_vars))
baseline_data[factor_vars] <- lapply(baseline_data[factor_vars], as.factor)

# transforming to date format ----
for(date in c("death_date", "date", "disaffiliation_date",
              "first_prescription_date", "date_end_fup")){
  baseline_data[[date]] <- as.Date(baseline_data[[date]], 
                                   format = "%Y-%m-%d")
}

# Run MICE
imputed_data <- mice(baseline_data, 
                     seed = 87593921, 
                     m = n_imputations, 
                     maxit = n_iterations, 
                     predictorMatrix = pred_matrix,
                     method = method_vec, 
                     print = TRUE) 

#### merge with long ----
# variables to keep before left join: those that are time-fixed
vars_keep <- c(
  "id", "time", "death_date",
  "date", "bdi_score", "education",
  "age", "female", "mmse_score",
  "disaffiliation_date", "drug", "milligrams",
  "first_prescription_date", "time_to_prescription", "death",
  "time_to_event", "date_end_fup", "time_baseline",
  "time_end_fup", "time_start_ap"
)

# Create a list of imputed datasets
imputed_list <- mice::complete(imputed_data, action = "all")

# Function to reshape and process each dataset
process_imputed <- function(df) {
  df <- df[, vars_keep]
  df <- df %>% # time, time_baseline, time_end_fup are in months since cohort entry
    mutate(month_seq = map2(time_baseline, time_end_fup, ~ seq(.x, .y, by = 3))) %>%
    
    unnest(month_seq) %>%
    
    mutate(time = month_seq) %>%
    
    left_join(long[!names(long) %in% c( "death_date", "age", 
                                        "female", "disaffiliation_date", "death" )], 
              by = c("id", "time")) %>% 
    # create start and end of discrete time-unit periods dates.
    # In this case, we selected 3 months as the time interval to update variable's values. 
    mutate(
      date_start = as.Date(date %m+% months(month_seq - time_baseline), 
                           format = "%Y-%m-%d"), 
      date_end = as.Date(date %m+% months(month_seq - time_baseline + 3), 
                         format = "%Y-%m-%d")
    )   %>%

    group_by(id) %>%
    mutate(
      date_start = ifelse(test = min(time) == time, 
                          yes = as.Date(date, format = "%Y-%m-%d"), 
                          no = as.Date(date_start, format = "%Y-%m-%d")
                          ),
      
      date_end = ifelse(test = max(time) == time,  
                        yes = as.Date(date_end_fup, format = "%Y-%m-%d"), 
                        no = as.Date(date_end, format = "%Y-%m-%d")
                        ),
      
      # flaggin start of AP
      ap_starts = ifelse(test = date_start <= first_prescription_date & 
                           date_end >= first_prescription_date &
                           !is.na(first_prescription_date), 
                         yes = 1, 
                         no = 0),
      # start of antipsychotic exposure interval - ever exposed
      antipsychotic_exposure = cumsum(ap_starts), # once exposed, ever exposed
      
      # reshaping tstart/tstop to allow for delayed entry methods in Cox's model
      tstart = as.numeric(difftime(as.Date(date_start, 
                                           format = "%Y-%m-%d"), 
                                   as.Date(date, 
                                           format = "%Y-%m-%d"), 
                                   units = "days")),
      
      tstop = as.numeric(difftime(as.Date(date_end, format = "%Y-%m-%d"), 
                                  as.Date(date, format = "%Y-%m-%d"), 
                                  units = "days"))
    )
  
  # Duplicate rows where AP starts, to create one exposed, one unexposed
  rows_ap_starts <- df %>% filter(ap_starts == 1) %>%
    mutate(tstart = as.numeric(difftime(
                    as.Date(first_prescription_date, 
                            format = "%Y-%m-%d"), 
                    as.Date(date, 
                            format = "%Y-%m-%d"), 
                    units = "days")))
  # prepare to merge with main dataset
  df <- df %>%
    mutate(
      antipsychotic_exposure = ifelse(test = ap_starts == 1, # if the row is the one that starts
                                      yes = 0,  # send the initial one to unexposed
                                      no = antipsychotic_exposure),
      
      tstop = ifelse(test = ap_starts == 1, 
                     yes = as.numeric(difftime(
                            as.Date(first_prescription_date, format = "%Y-%m-%d"), 
                            as.Date(date, format = "%Y-%m-%d"), 
                            units = "days")), 
                     tstop)
    )
  
  df <- bind_rows(df, rows_ap_starts)

  # Final flags
  df <- df %>%
    group_by(id) %>%
    # death_fup = 1 only in the final interval if death_date is present (death during follow-up)
    mutate(death_fup = ifelse(test = max(time) == time & !is.na(death_date), 
                              yes = 1, 
                              no = 0)) %>%
    ungroup()
  
  df <- df %>% mutate(time_interval = tstop - tstart)
  
  if (any(df$time_interval < 0)) warning("Negative intervals detected before filtering.")
  df <- df %>% filter(time_interval >= 0)
  
  df <- df[order(df$id, df$time), ]
  return(df)
}

long_imputed_list <- lapply(imputed_list, process_imputed)

# create function for Cox PH model
fit_main_model <- function(data) {
  surv_tvc_ap <- Surv(time =  data$tstart, 
                      time2 = data$tstop, 
                      event  = data$death_fup)
  
  model <- coxph(surv_tvc_ap ~ 
                   antipsychotic_exposure +
                   bdi_score + mmse_score + education + age + 
                   female + parkinson + dementia_lewy + tobacco +  
                   epilepsy + asthma +
                   hypertension + ckd + diabetes + insomnia +
                   anxiety + depression + peripheric_vascular + 
                   cardiac_failure +
                   coronary + copd + cancer +
                   opioids + benzodiazepines + z_drugs + 
                   antidepressants + antiepileptics 
                   , 
                data = data)
  return(model)
}

models <- lapply(1:n_imputations, 
                 function(i) fit_main_model(long_imputed_list[[i]]))

# Pooling Cox models across imputed datasets using Rubin's rules
pooled_model <- pool(models) # some warnings may appear referring to dates having inconsistent format 
                              # that is expected due to variation introduced in the creation of the dataset. 

# Print and save as cvs the results
main_model <- summary(pooled_model)
main_model

# end of R script

