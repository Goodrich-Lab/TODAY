# This script is to get categorical exposures.

# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Example-----
## Step1. load data----
data <- read_rds(fs::path(dir_data, "example_teenlabs_data.rds"))

## Step2. Define variables----

# plasma pfas
plasma_pfas <- colnames(data)[grepl("imputed",colnames(data))]
## or
# plasma_pfas <- c("pfda_untargeted_plasma_0_imputed",
#                  "pf_hx_s_untargeted_plasma_0_imputed",
#                  "pf_hp_s_untargeted_plasma_0_imputed",
#                  "pfna_untargeted_plasma_0_imputed",
#                  "pfoa_untargeted_plasma_0_imputed",
#                  "pfos_untargeted_plasma_0_imputed")


#Step3. Adding the categorical exposures into the data set----
data <- data %>%
  mutate_at(.vars = vars(all_of(plasma_pfas)), 
            .funs = list(tertile = ~ntile(.,  3), 
                         median = ~if_else(. > median(., na.rm = TRUE), 
                                           2,
                                           1),
                         quantile = ~as.integer(cut(., quantile(., c(0, 1/4,1/2,3/4,1), 
                                                                na.rm = TRUE))),
                         quintile = ~as.integer(cut(., quantile(., c(0, 1/5,2/5,3/5,4/5,1), 
                                                                na.rm = TRUE))),
                         sextile = ~as.integer(cut(., quantile(., c(0, 1/6,2/6,3/6,4/6,5/6,1), 
                                                               na.rm = TRUE))),
                         septile = ~as.integer(cut(., quantile(., c(0, 1/7,2/7,3/7,4/7,5/7,6/7,1), 
                                                               na.rm = TRUE))),
                         octile = ~as.integer(cut(., quantile(., c(0, 1/8,2/8,3/8,4/8,5/8,6/8,7/8,1), 
                                                              na.rm = TRUE)))
            ))

## Step4. Adding categorical outcome to the data-----
# creating new dichotomous ALT variable
data <- data %>% tidylog::drop_na(alt_0) %>%
  mutate(alt_di_0 = case_when(
    alt_0 > 25.8 & sex =="male" ~ 'high',
    alt_0 > 22.1 & sex =="female" ~ 'high',
    TRUE ~ 'low'
  ), 
  alt_di_0 = factor(alt_di_0, levels = c("low", "high"), 
                    ordered = TRUE)
  )

# extracted the related column names and stored into a variable-----
pfas_names_categorical <- colnames(data)[
  grepl("median|tertile|quantile|quintile|sextile|septile|octile",
        colnames(data))]

#Step5. Getting the intervals for each tiled data-----
data_interval <- data %>%
  mutate_at(.vars = vars(all_of(plasma_pfas)), 
            .funs = list(tertile = ~cut(., quantile(., c(0, 1/3, 2/3, 1), 
                                                    na.rm = TRUE)), 
                         median = ~cut(., quantile(., c(0, 1/2, 1), 
                                                   na.rm = TRUE)),
                         quantile = ~cut(., quantile(., c(0, 1/4,1/2,3/4,1), 
                                                     na.rm = TRUE)),
                         quintile = ~cut(., quantile(., c(0, 1/5,2/5,3/5,4/5,1), 
                                                     na.rm = TRUE)),
                         sextile = ~ cut(., quantile(., c(0, 1/6,2/6,3/6,4/6,5/6,1), 
                                                     na.rm = TRUE)),
                         septile = ~ cut(., quantile(., c(0, 1/7,2/7,3/7,4/7,5/7,6/7,1), 
                                                     na.rm = TRUE)),
                         octile  = ~ cut(., quantile(., c(0, 1/8,2/8,3/8,4/8,5/8,6/8,7/8,1), 
                                                     na.rm = TRUE))
            ))

#Step6. saving the intervals information to a file-----
interval <- data_interval %>% 
  dplyr::select(all_of(pfas_names_categorical)) %>%
  map(.,table) %>% 
  map(.,as.data.frame) %>% 
  bind_rows(.id = "pfas") %>%
  dplyr::rename(interval = Var1,
                count = Freq)

n <- length(plasma_pfas)

num <- c(rep(c(1,2,3),n),rep(c(1,2),n),
         rep(c(1,2,3,4),n), rep(c(1,2,3,4,5),n),
         rep(c(1,2,3,4,5,6),n), rep(c(1,2,3,4,5,6,7),n),
         rep(c(1,2,3,4,5,6,7,8),n))

interval_df<- interval %>% mutate(level = num,
                                  pfas = str_c(pfas, level, ""))

# write_csv(interval_df,fs::path(dir_report,
#                                "Categorical_targeted_plasma_pfas_intervals.csv"))

