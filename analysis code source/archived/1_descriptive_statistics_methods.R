
# This script include three methods to get descriptive statistics

# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Method 1. summarytools::view(summarytools::dfSummary())
# Method 2. table1::table1()
# Method 3. group_by() --> summarise()

# Example-----
## Step1. load data----
data <- read_rds(fs::path(dir_data, "example_teenlabs_data.rds"))

## Step2. Define variables----
# covariates
covars <- c("race_binary", 
            "smoke_0", 
            "age_0", 
            "sex", 
            "parents_income_0")

# categorical outcomes(two categories or more)
cat_outcomes <- c("nafld_nash_mul_0", "nafld_di_0","steatgrd_mul_0", "bhepa_0",
                   "dyslipid_di_0")

# continuous outcomes
cont_outcomes <- c("alt_0", "ast_0")

# plasma pfas
plasma_pfas <- colnames(data)[grepl("imputed",colnames(data))]
## or
# plasma_pfas <- c("pfda_untargeted_plasma_0_imputed",
#                  "pf_hx_s_untargeted_plasma_0_imputed",
#                  "pf_hp_s_untargeted_plasma_0_imputed",
#                  "pfna_untargeted_plasma_0_imputed",
#                  "pfoa_untargeted_plasma_0_imputed",
#                  "pfos_untargeted_plasma_0_imputed")

# Step 3. Descriptive statistics-----
## Method 1. summarytools::view(summarytools::dfSummary())
summarytools::view(summarytools::dfSummary(data))


## Method 2. table1::table1()
stats <- table1(~., data = data %>% dplyr::select(-key))

stats

## Method 3. pivot_longer --> group_by() --> summarise()
## Take pfas as example

# long format data for plasma pfas
data_l <- data %>% 
  pivot_longer(
  names_to = "variables",
  values_to = "value",
  cols = c(all_of(plasma_pfas)))

# descriptive statistics
stats <- data_l %>%
  drop_na(value) %>%
  group_by(variables) %>%
  dplyr::summarise(
    geometric_mean = fungm(value),
    percentile_50 = qntle_fxn(value, .50),
    percentile_75 = qntle_fxn(value, .75),
    percentile_90 = qntle_fxn(value, .9)) %>%
  ungroup()

# Rename columns
stats_summary <- stats %>%
  dplyr::rename(
    "PFAS" = variables,
    "Geometric Mean [95% CI]" = geometric_mean,
    "50th percentile" = percentile_50,
    "75th percentile" = percentile_75,
    "90th percentile" = percentile_90)

# save result
# write_csv(stats_summary, fs::path(dir_report,
#                                   "1_0_descriptive_stats.csv"))
