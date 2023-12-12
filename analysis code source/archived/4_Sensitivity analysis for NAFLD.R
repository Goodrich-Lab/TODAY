# This script is to do sensitivity analysis for NAFLD
# control for dichotomous site
# control for alb

# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Load data
data <- read_rds(fs::path(dir_project_data, "tl_analysis_ready_data.rds"))

data <- data %>% mutate(site_bi = ifelse(site == "CIN",
                                         "CIN", "NON CIN"))
# NAFLD related outcomes----
cat_outcomes <- c("nafld_nash_mul_0", "nafld_di_0",
                  "steato_mul_0", "steatgrd_mul_0", "bhepa_0",
                  "fibrostg_di_0","lob_di_0","nash_mul_0")

cat_outcomes_di <- c("nafld_di_0","fibrostg_di_0",
                     "lob_di_0")

cat_outcomes_mul <- c("nafld_nash_mul_0", "steato_mul_0", 
                      "steatgrd_mul_0", "bhepa_0", "nash_mul_0")

cont_outcomes <- c("alt_0", "ast_0", "ggt_0")

# Control for site -----
covars <- c("bmi_0", 
            "race_binary", 
            "site_bi", 
            "age_0", 
            "sex", 
            "parents_income_0")


# exposures: plasma pfas----
plasma_pfas <- colnames(data)[grep("0_imputed|_targeted_plasma", colnames(data))]

# Model1: Logistic regression: cat_outcomes_di ~ plasma pfas + covars----
## Create regression model formula----
eo_combinations1 <- list(
  x = plasma_pfas,
  y = cat_outcomes_di) %>%
  cross_df()

models1 <- eo_combinations1 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

## Run logistic regression----
models1$output <- map(models1$formula, 
                      ~glm(as.formula(.),
                           data = data,
                           na.action = na.exclude,
                           family = "binomial") %>%
                        tidy(., conf.int = TRUE)
)

## Create final dataframe of all results----
mod_output_df1 <- models1 %>%
  unnest(output) %>%
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 

# Model2: Ordinal Logistic regression: cat_outcomes_mul ~ plasma pfas + covars----
## Create regression model formula----
eo_combinations2 <- list(
  x = plasma_pfas,
  y = cat_outcomes_mul) %>%
  cross_df()

models2 <- eo_combinations2 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

## Run ordinal logistic regression----
models2$output <- map(models2$formula, 
                      ~polr(as.formula(.),
                            data = data,
                            na.action = na.exclude,
                            Hess = TRUE) %>%
                        tidy(., conf.int = TRUE))

## Create final data frame of all results----
mod_output_df2 <- models2 %>%
  unnest(output) %>%   
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 

# Model3: Linear regression: cont_outcomes ~ plasma pfas + covars----
## Create group mean centered values for the continuous outcomes

## Create regression model formula----
eo_combinations3 <- list(
  x = plasma_pfas,
  y = cont_outcomes) %>%
  cross_df()

models3 <- eo_combinations3 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c("scale(", y,")","~","log2(", x,")", "+", covars))

## Run linear regression----
models3$output <- map(models3$formula, 
                      ~lm(as.formula(.),
                          data = data,
                          na.action = na.exclude) %>%
                        tidy(., conf.int = TRUE)
)

## Create final data frame of all results----
mod_output_df3 <- models3 %>%
  unnest(output) %>% 
  mutate(matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 
## Merge all the results----
mod_output_df <- mod_output_df1 %>% 
  mutate(method = "logistic regression") %>%
  bind_rows(mod_output_df2 %>% 
              mutate(method = "ordinal logistic regression")%>% 
              dplyr::select(-coef_type)) %>%
  bind_rows(mod_output_df3 %>% 
              mutate(method = "linear regression")) %>%
  dplyr::select(method, everything())

## Saving results----
write_csv(mod_output_df,fs::path(dir_report,
                                 "4_adjusted_univariate_associations_NAFLD_with_site.csv"))

# Coefficient plot------
## untargeted pfas and categorical outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "untargeted" & method != "linear regression") %>%
  ggplot(aes(x = pfas,y = odds_ratio)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = exp_conf_low,
                    ymax = exp_conf_high),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab("Odds Ratio (95% CI)") +
  xlab("Untargeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_untargeted_pfas_cat_outcomes_with_site.png"),
       width = 7, height = 5.5, dpi = 300)

## targeted pfas and categorical outcomes---- 
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "targeted" & method != "linear regression") %>%
  ggplot(aes(x = pfas,y = odds_ratio)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = exp_conf_low,
                    ymax = exp_conf_high),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab("Odds Ratio (95% CI)") +
  xlab("targeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_targeted_pfas_cat_outcomes_with_site.png"),
       width = 7, height = 5.5, dpi = 300)

## untargeted pfas and continuous outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "untargeted" & method == "linear regression") %>%
  ggplot(aes(x = pfas,y = estimate)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Beta (95% CI)") +
  xlab("Untargeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_untargeted_pfas_cont_outcomes_with_site.png"),
       width = 7, height = 3, dpi = 300)

## targeted pfas and continuous outcomes---- 
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "targeted" & method == "linear regression") %>%
  ggplot(aes(x = pfas,y = estimate)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Beta (95% CI)") +
  xlab("targeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity", "4_adjusted_NAFLD_targeted_pfas_cont_outcomes_with_site.png"),
       width = 7, height = 5, dpi = 300)


# Control for site -----
covars <- c("bmi_0", 
            "race_binary", 
            "alb_0", 
            "age_0", 
            "sex", 
            "parents_income_0")

data$alb_0 <- scale(data$alb_0)

# exposures: plasma pfas----
plasma_pfas <- colnames(data)[grep("0_imputed|_targeted_plasma", colnames(data))]

# Model1: Logistic regression: cat_outcomes_di ~ plasma pfas + covars----
## Create regression model formula----
eo_combinations1 <- list(
  x = plasma_pfas,
  y = cat_outcomes_di) %>%
  cross_df()

models1 <- eo_combinations1 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

## Run logistic regression----
models1$output <- map(models1$formula, 
                      ~glm(as.formula(.),
                           data = data,
                           na.action = na.exclude,
                           family = "binomial") %>%
                        tidy(., conf.int = TRUE)
)

## Create final dataframe of all results----
mod_output_df1 <- models1 %>%
  unnest(output) %>%
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 

# Model2: Ordinal Logistic regression: cat_outcomes_mul ~ plasma pfas + covars----
## Create regression model formula----
eo_combinations2 <- list(
  x = plasma_pfas,
  y = cat_outcomes_mul) %>%
  cross_df()

models2 <- eo_combinations2 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

## Run ordinal logistic regression----
models2$output <- map(models2$formula, 
                      ~polr(as.formula(.),
                            data = data,
                            na.action = na.exclude,
                            Hess = TRUE) %>%
                        tidy(., conf.int = TRUE))

## Create final data frame of all results----
mod_output_df2 <- models2 %>%
  unnest(output) %>%   
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 

# Model3: Linear regression: cont_outcomes ~ plasma pfas + covars----
## Create group mean centered values for the continuous outcomes

## Create regression model formula----
eo_combinations3 <- list(
  x = plasma_pfas,
  y = cont_outcomes) %>%
  cross_df()

models3 <- eo_combinations3 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c("scale(", y,")","~","log2(", x,")", "+", covars))

## Run linear regression----
models3$output <- map(models3$formula, 
                      ~lm(as.formula(.),
                          data = data,
                          na.action = na.exclude) %>%
                        tidy(., conf.int = TRUE)
)

## Create final data frame of all results----
mod_output_df3 <- models3 %>%
  unnest(output) %>% 
  mutate(matrices = ifelse(grepl("_untargeted_", x), 
                           "untargeted", 
                           "targeted")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  clean_names() 
## Merge all the results----
mod_output_df <- mod_output_df1 %>% 
  mutate(method = "logistic regression") %>%
  bind_rows(mod_output_df2 %>% 
              mutate(method = "ordinal logistic regression")%>% 
              dplyr::select(-coef_type)) %>%
  bind_rows(mod_output_df3 %>% 
              mutate(method = "linear regression")) %>%
  dplyr::select(method, everything())

## Saving results----
write_csv(mod_output_df,fs::path(dir_report,
                                 "Sensitivity",
                                 "4_adjusted_univariate_associations_NAFLD_with_alb.csv"))

# Coefficient plot------
## untargeted pfas and categorical outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "untargeted" & method != "linear regression") %>%
  ggplot(aes(x = pfas,y = odds_ratio)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = exp_conf_low,
                    ymax = exp_conf_high),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab("Odds Ratio (95% CI)") +
  xlab("Untargeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_untargeted_pfas_cat_outcomes_with_alb.png"),
       width = 7, height = 5.5, dpi = 300)

## targeted pfas and categorical outcomes---- 
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "targeted" & method != "linear regression") %>%
  ggplot(aes(x = pfas,y = odds_ratio)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = exp_conf_low,
                    ymax = exp_conf_high),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  ylab("Odds Ratio (95% CI)") +
  xlab("targeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_targeted_pfas_cat_outcomes_with_alb.png"),
       width = 7, height = 5.5, dpi = 300)

## untargeted pfas and continuous outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "untargeted" & method == "linear regression") %>%
  ggplot(aes(x = pfas,y = estimate)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Beta (95% CI)") +
  xlab("Untargeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_untargeted_pfas_cont_outcomes_with_alb.png"),
       width = 7, height = 3, dpi = 300)

## targeted pfas and continuous outcomes---- 
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(matrices == "targeted" & method == "linear regression") %>%
  ggplot(aes(x = pfas,y = estimate)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Beta (95% CI)") +
  xlab("targeted Plasma PFAS") + 
  theme_bw() +
  theme(text = element_text(size = 10))+
  theme(legend.title = element_blank()) +
  coord_flip() +
  facet_wrap(.~ y,scales = "free_x")

ggsave(fs::path(dir_figure, "Sensitivity","4_adjusted_NAFLD_targeted_pfas_cont_outcomes_with_alb.png"),
       width = 7, height = 5, dpi = 300)




rm(cat_outcomes, cat_outcomes_di, 
   cat_outcomes_mul, cont_outcomes,
   models1, models2, models3, 
   eo_combinations1, eo_combinations2, eo_combinations3,
   mod_output_df1, mod_output_df2, mod_output_df3,
   mod_output_df, covars, plasma_pfas)
