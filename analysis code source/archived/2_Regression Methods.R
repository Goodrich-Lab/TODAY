# This script include regression models for different type of outcomes

# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Overview-----
# Model 1. Binary Logistic Regression--- Binary outcome
# Model 2. Ordinal Logistic Regression --- ordinal outcome
# Model 3. Linear regression--- continuous outcome
# Model 4. Multinomial Logistic Regression --- multi categorical outcome
# Creating Coefficient plots

# Example-----
## Step1. load data----
data <- read_rds(fs::path(dir_data, "example_teenlabs_data.rds"))

## Step 3. Model2: Ordinal Logistic regression: cat_outcomes_mul ~ plasma pfas + covars----
# log2 transfored for the exposures/plasma pfas

### Create regression model formula----
eo_combinations2 <- list(
  x = plasma_pfas,
  y = cat_outcomes_mul) %>%
  cross_df()

models2 <- eo_combinations2 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

### Run ordinal logistic regression----
models2$output <- map(models2$formula, 
                     ~polr(as.formula(.),
                          data = data,
                          na.action = na.exclude,
                          Hess = TRUE) %>%
                       tidy(., conf.int = TRUE))

### Create final data frame of all results----
mod_output_df2 <- models2 %>%
  unnest(output) %>%   
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low)) %>%
  clean_names() 

## Step3. Model3: Linear regression: cont_outcomes ~ plasma pfas + covars----
# group mean centered values for the continuous outcomes
# log2 transfored for the exposures/plasma pfas

### Create regression model formula----
eo_combinations3 <- list(
  x = plasma_pfas,
  y = cont_outcomes) %>%
  cross_df()

models3 <- eo_combinations3 %>%
  mutate(covars = str_c(covars, collapse = "+"),
    formula = str_c("scale(", y,")","~","log2(", x,")", "+", covars))

### Run linear regression----
models3$output <- map(models3$formula, 
                      ~lm(as.formula(.),
                            data = data,
                            na.action = na.exclude) %>%
                        tidy(., conf.int = TRUE)
)

### Create final data frame of all results----
mod_output_df3 <- models3 %>%
  unnest(output) %>% 
  clean_names() 

## Step3. Model4: Multinomial Logistic regression: cat_outcomes_mul ~ plasma pfas + covars----
### Create regression model formula----
eo_combinations4 <- list(
  x = plasma_pfas,
  y = cat_outcomes_mul) %>%
  cross_df()

models4 <- eo_combinations4 %>%
  mutate(covars = str_c(covars, collapse = "+"),
         formula = str_c(y,"~","log2(", x,")", "+", covars))

### Run multinomial logistic regression----
models4$output <- map(models4$formula, 
                      ~nnet::multinom(as.formula(.),
                                      data = data,
                                      na.action = na.exclude) %>%
                        tidy(., conf.int = TRUE))

### Create final data frame of all results----
mod_output_df4 <- models4 %>%
  unnest(output) %>%   
  mutate(odds_ratio = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low)) %>%
  clean_names() 

## Step4. Merge all the results from four models----
mod_output_df <- mod_output_df1 %>% 
  mutate(method = "logistic regression") %>%
  bind_rows(mod_output_df2 %>% 
              mutate(method = "ordinal logistic regression")%>% 
              dplyr::select(-coef_type)) %>%
  bind_rows(mod_output_df3 %>% 
              mutate(method = "linear regression")) %>%
  bind_rows(mod_output_df4 %>% 
              mutate(method = "multinomial logsitic regression")) %>%
  separate(x, 
           into = c("pfas", NA), 
           sep = "_targeted|_untargeted",
           remove = FALSE) %>%
  dplyr::select(method, everything())
  

## Step5. Saving results----
# write_csv(mod_output_df,fs::path(dir_report,
#                                  "1_Regression Result.csv"))
  
## Step6. Coefficient plots------
### Example1.Binary categorical outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
    filter(method == "logistic regression") %>%
    ggplot(aes(x = pfas,y = odds_ratio)) +
    geom_point(size = 1) +
    coord_flip() +
    geom_errorbar(aes(ymin = exp_conf_low,
                      ymax = exp_conf_high),
                  width = 0) +
    geom_hline(yintercept = 1, linetype = 2) +
    ylab("Odds Ratio (95% CI)") +
    xlab("Plasma PFAS") + 
    theme_bw() +
    theme(text = element_text(size = 10))+
    theme(legend.title = element_blank()) +
    coord_flip() +
    facet_wrap(.~ y,scales = "free_x")
### Example2.Continuous categorical outcomes----
mod_output_df[grep("log2", mod_output_df$term),] %>%
  filter(method == "linear regression") %>%
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

##Step7. Saving coefficient plots-----
# ggsave(fs::path(dir_figure,"1_coefficient_plot.png"),
#        width = 7, height = 3, dpi = 300)


