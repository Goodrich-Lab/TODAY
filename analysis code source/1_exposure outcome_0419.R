# Exposure outcome associations
library(ggbeeswarm)
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))

# All potential outcomes: 
# MIC (dichot): Development to mic
# DAYSTOMIC (continuous): Days to microalbumineria
# MAC (dichot): macroalbuminera
# DAYSTOMAC (continuous): Days to MAC
# RAPID (dichot): 
# DAYSTORAPID (continuous):  
# HYP0 (dichot):
# DAYSTOHYP (continuous):
outcome <- c("time_to_glyc_scld", 
             "codi", 
             "si_1_ins0",
             "hb_a1c", 
             "mic", 
             "daystomic", 
             # "mac",
             "daystomac", 
             "rapid",
             "daystorapid", 
             "hyp0", 
             "daystohyp")
## MAC only has 0 in that, so I exclude it
# Descriptive plots ----
# Days to glycemic failure, per group

## Summary stats of reuslt
data = full_data %>% 
  dplyr::select(all_of(outcome)) %>%
  mutate(mic = as.factor(mic),
         rapid = as.factor(rapid),
         hyp0 = as.factor(hyp0))

table1::table1(~., data = data)

full_data <- full_data %>% 
  mutate_at(.vars = vars(outcome[c(1:4, 6,7,9, 11)]),
            .funs = ~scale(.))

# Other outcomes_HW -----
# this is a function to run owas
owas_mul <- function(outcome, family){
  # filter the full data to only include cases
  exwas_lm_case <- full_data |>
    tidylog::filter(case1_control0 == 1) |>
    # log2 transform the pfas variables
    mutate(across(all_of(pfas_names), log2)) |>
    # run the owas function from the epiomics package
    epiomics::owas(var = outcome,
                   var_exposure_or_outcome = "outcome",
                   covars = c("sex", "agebase", "tx"),
                   omics = pfas_names, 
                   family = family,
                   conf_int = TRUE)
  # if the family is gaussian, then calculate the confidence intervals for the beta estimates
  if(family == "gaussian"){
    exwas_lm_case <- exwas_lm_case %>% 
      mutate(beta_ci = jag2::effest_ci(estimate, conf_low, conf_high)) %>%
      mutate(case_type = "case") %>%
      select(case_type, everything())
    # if the family is binomial, then calculate the confidence intervals for the odds ratios
  }else{
    exwas_lm_case <- exwas_lm_case %>%
      mutate(odds_ratio_ci = jag2::effest_ci(exp(estimate), exp(conf_low), exp(conf_high)))%>%
      mutate(case_type = "case") %>%
      select(case_type, everything())
  }
  # filter out the pfas variables that are not included in the analysis
  exwas_lm_case <- exwas_lm_case %>% 
    filter(!feature_name %in% 
             c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
    # create a new variable that is the name of the pfas variable
    mutate(pfas_name = case_when(grepl("pfos", feature_name) ~ "PFOS",
                                 grepl("pfoa", feature_name) ~ "PFOA",
                                 grepl("pfda", feature_name) ~ "PFDA",
                                 grepl("pfhpa", feature_name) ~ "PFHpA",
                                 grepl("pfhps", feature_name) ~ "PFHpS",
                                 grepl("pfhxs", feature_name) ~ "PFHxS",
                                 grepl("pfna", feature_name) ~ "PFNA",
                                 grepl("pfuna", feature_name) ~ "PFUnA",
                                 grepl("nmefosaa", feature_name) ~ "NMeFOSAA"))
  
  # filter the full data to only include controls
  exwas_lm_control <- full_data |>
    tidylog::filter(case1_control0 == 0) |>
    # log2 transform the pfas variables
    mutate(across(all_of(pfas_names), log2)) |>
    # run the owas function from the epiomics package
    epiomics::owas(var = outcome,
                   var_exposure_or_outcome = "outcome",
                   covars = c("sex", "agebase", "tx"),
                   omics = pfas_names, 
                   family = family,
                   conf_int = TRUE) 
  # if the family is gaussian, then calculate the confidence intervals for the beta estimates
  if(family == "gaussian"){
    exwas_lm_control <- exwas_lm_control %>% 
      mutate(beta_ci = jag2::effest_ci(estimate, conf_low, conf_high)) %>%
      mutate(case_type = "control") %>%
      select(case_type, everything())
    # if the family is binomial, then calculate the confidence intervals for the odds ratios
  }else{
    exwas_lm_control <- exwas_lm_control %>%
      mutate(odds_ratio_ci = jag2::effest_ci(exp(estimate), exp(conf_low), exp(conf_high)))%>%
      mutate(case_type = "control") %>%
      select(case_type, everything())
  }
  
  # filter out the pfas variables that are not included in the analysis
  exwas_lm_control <- exwas_lm_control %>% 
    filter(!feature_name %in% 
             c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
    # create a new variable that is the name of the pfas variable
    mutate(pfas_name = case_when(grepl("pfos", feature_name) ~ "PFOS",
                                 grepl("pfoa", feature_name) ~ "PFOA",
                                 grepl("pfda", feature_name) ~ "PFDA",
                                 grepl("pfhpa", feature_name) ~ "PFHpA",
                                 grepl("pfhps", feature_name) ~ "PFHpS",
                                 grepl("pfhxs", feature_name) ~ "PFHxS",
                                 grepl("pfna", feature_name) ~ "PFNA",
                                 grepl("pfuna", feature_name) ~ "PFUnA",
                                 grepl("nmefosaa", feature_name) ~ "NMeFOSAA"))
  
  # bind the case and control data together
  result <- exwas_lm_case %>% bind_rows(exwas_lm_control)
  return(result)
}



# Categorical outcomes
mic_result    <- owas_mul("mic",   "binomial")
rapid_result  <- owas_mul("rapid", "binomial")
hyp0_result   <- owas_mul("hyp0",  "binomial")
result_cat <- mic_result  %>% bind_rows(rapid_result) %>%
  bind_rows(hyp0_result) 


# Continuous outcomes
time_to_glyc_scld_result    <- owas_mul("time_to_glyc_scld", "gaussian")
codi_result                 <- owas_mul("codi",              "gaussian")
si_1_ins0_result            <- owas_mul("si_1_ins0",         "gaussian")
hb_a1c_result               <- owas_mul("hb_a1c",            "gaussian")
daystomic_result   <- owas_mul("daystomic",   "gaussian")
daystomac_result   <- owas_mul("daystomac",   "gaussian")
daystorapid_result <- owas_mul("daystorapid", "gaussian")
daystohyp_result   <- owas_mul("daystohyp",   "gaussian")

result_cont<- daystomic_result %>%
  bind_rows(time_to_glyc_scld_result)  |>
  bind_rows(codi_result)  |>
  bind_rows(si_1_ins0_result)  |>
  bind_rows(hb_a1c_result)  |>
  bind_rows(daystomac_result) %>% 
  bind_rows(daystorapid_result) %>% 
  bind_rows(daystohyp_result)

# Combine continuous and categorical
result <- result_cat %>% bind_rows(result_cont)

write_csv(result, 
          fs::path(dir_results,
                   "Stratified case control analysis multiple outcomes.csv"))

levels <- c("PFUnA", "PFHpS","PFHpA","NMeFOSAA","PFDA",
            "PFNA","PFHxS", "PFOA", "PFOS")
## coefficient plot for rapid progression
(p1 <- result_cont %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid( ~ var_name + case_type, scales = "free") +
  ylab("Coefficient (95% CI)") +
  # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip())
table(result_cont$p_value < 0.1)


(p1 <- result_cont %>%
    filter(p_value < 0.1) |>
    ggplot(aes(x = factor(pfas_name, levels = levels), 
               y = estimate,
               color = threshold)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = conf_low,
                      ymax = conf_high),
                  width = 0, position = position_dodge()) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid( case_type + var_name ~ ., scales = "free") +
    ylab("Coefficient (95% CI)") +
    # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
    theme(text = element_text(size = 10),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill="white"), 
          strip.background = element_rect(fill = "white"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black")))# +
    coord_flip())


ggsave(fs::path(dir_figure, "coef_plot_cont_outcome.png"), bg = "white",
       width = 12, height = 3, dpi = 300)

p2 <- result_cat %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = exp(estimate))) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = exp(conf_low),
                    ymax = exp(conf_high)),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid( ~ var_name + case_type, scales = "free") +
  ylab("Odds Ratio (95% CI)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 

mic_result <- owas_mul("mic", "binomial")
daystomic_result <- owas_mul("daystomic", "gaussian")
daystomac_result <- owas_mul("daystomac", "gaussian")
rapid_result <- owas_mul("rapid", "binomial")
daystorapid_result <- owas_mul("daystorapid", "gaussian")
hyp0_result <- owas_mul("hyp0", "binomial")
daystohyp_result <- owas_mul("daystohyp", "gaussian")

result_cat <- mic_result  %>% bind_rows(rapid_result) %>%
  bind_rows(hyp0_result) 

result_cont<- daystomic_result %>%
  bind_rows(daystomac_result) %>% 
  bind_rows(daystorapid_result) %>% 
  bind_rows(daystohyp_result)

result <- result_cat %>% bind_rows(result_cont)

write_csv(result, 
          fs::path(dir_results,
                   "Stratified case control analysis multiple outcomes.csv"))

levels <- c("PFUnA", "PFHpS","PFHpA","NMeFOSAA","PFDA","PFNA","PFHxS", "PFOA", "PFOS" 
)
## coefficient plot for rapid progression
p1 <- result_cont %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid( ~ var_name + case_type, scales = "free") +
  ylab("Coefficient (95% CI)") +
  # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip()

ggsave(fs::path(dir_figure, "coef_plot_cont_outcome.png"), bg = "white",
       width = 12, height = 3, dpi = 300)

p2 <- result_cat %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = exp(estimate))) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = exp(conf_low),
                    ymax = exp(conf_high)),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid( ~ var_name + case_type, scales = "free") +
  ylab("Odds Ratio (95% CI)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 


# cowplot::plot_grid(p1, p2, labels = c("A. Rapid Progression Y-T2D","B.Slow Progression Y-T2D"),
#                    label_size = 12 , vjust = 0.5, hjust = 0) + 
#   theme(plot.margin = unit(c(0.7, 0.5, 0.5, 0.5), "cm"))
# 
# ggsave(fs::path(dir_figure, "coef_plot.png"), bg = "white",
#        width = 6, height = 3, dpi = 300)
# Volcano plot
epiomics::volcano_owas(result_cont %>% filter(case_type == "case"))
epiomics::volcano_owas(result_cont %>% filter(case_type == "control"))

epiomics::volcano_owas(result_cat %>% filter(case_type == "case"))
epiomics::volcano_owas(result_cat %>% filter(case_type == "control"))

# # Plot PFDA
# full_data |>
#   # tidylog::filter(case1_control0 == 1) |>
#   ggplot(aes(x = pfas_pfda, 
#              y = monthstoglyc)) +
#   geom_point() + 
#   facet_wrap(~case1_control0, scales =  "free_y")
# 


# Analysis with whole group adjusted for "case1_control0", "sex", "agebase", "tx"--------
data <- full_data %>% 
  mutate_at(.vars = vars(outcome[c(2,3,5,7)], pfas_names),
            .funs = ~scale(.))

cont_outcome_result <- model_output(pfas_names,
                                    outcomes = outcome[c(2,3,5,7)],
                                    covars = c("case1_control0", "sex", "agebase", "tx") ,
                                    outcome_family = "gaussian",
                                    data = data)
cont_result <- cont_outcome_result %>% 
  filter(term == "exposure_concentration") %>%
  select(-term) %>%
  filter(!exposures %in% 
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", exposures) ~ "PFOS",
                               grepl("pfoa", exposures) ~ "PFOA",
                               grepl("pfda", exposures) ~ "PFDA",
                               grepl("pfhpa", exposures) ~ "PFHpA",
                               grepl("pfhps", exposures) ~ "PFHpS",
                               grepl("pfhxs", exposures) ~ "PFHxS",
                               grepl("pfna", exposures) ~ "PFNA",
                               grepl("pfuna", exposures) ~ "PFUnA",
                               grepl("nmefosaa", exposures) ~ "NMeFOSAA"))

cat_outcome_result <- model_output(pfas_names,
                                   outcomes = outcome[c(1, 4, 6)],
                                   covars = c("case1_control0", "sex", "agebase", "tx"),
                                   outcome_family = "binomial",
                                   data = data)
cat_result <- cat_outcome_result %>% 
  filter(term == "exposure_concentration") %>%
  select(-term) %>%
  filter(!exposures %in% 
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", exposures) ~ "PFOS",
                               grepl("pfoa", exposures) ~ "PFOA",
                               grepl("pfda", exposures) ~ "PFDA",
                               grepl("pfhpa", exposures) ~ "PFHpA",
                               grepl("pfhps", exposures) ~ "PFHpS",
                               grepl("pfhxs", exposures) ~ "PFHxS",
                               grepl("pfna", exposures) ~ "PFNA",
                               grepl("pfuna", exposures) ~ "PFUnA",
                               grepl("nmefosaa", exposures) ~ "NMeFOSAA"))


cont_outcome_result1 <- model_output(pfas_names,
                                     outcomes = outcome[c(2,3,5,7)],
                                     covars = c("case1_control0", "sex", "agebase", "tx", "est_creat_clear") ,
                                     outcome_family = "gaussian",
                                     data = data) 
cont_result1 <- cont_outcome_result1 %>% 
  filter(term == "exposure_concentration") %>%
  select(-term)%>%
  filter(!exposures %in% 
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", exposures) ~ "PFOS",
                               grepl("pfoa", exposures) ~ "PFOA",
                               grepl("pfda", exposures) ~ "PFDA",
                               grepl("pfhpa", exposures) ~ "PFHpA",
                               grepl("pfhps", exposures) ~ "PFHpS",
                               grepl("pfhxs", exposures) ~ "PFHxS",
                               grepl("pfna", exposures) ~ "PFNA",
                               grepl("pfuna", exposures) ~ "PFUnA",
                               grepl("nmefosaa", exposures) ~ "NMeFOSAA"))

cat_outcome_result1 <- model_output(pfas_names,
                                    outcomes = outcome[c(1, 4, 6)],
                                    covars = c("case1_control0", "sex", 
                                               "agebase", "tx",
                                               "est_creat_clear"),
                                    outcome_family = "binomial",
                                    data = data)

cat_result1 <- cat_outcome_result1 %>% 
  filter(term == "exposure_concentration") %>%
  select(-term)%>%
  filter(!exposures %in% 
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", exposures) ~ "PFOS",
                               grepl("pfoa", exposures) ~ "PFOA",
                               grepl("pfda", exposures) ~ "PFDA",
                               grepl("pfhpa", exposures) ~ "PFHpA",
                               grepl("pfhps", exposures) ~ "PFHpS",
                               grepl("pfhxs", exposures) ~ "PFHxS",
                               grepl("pfna", exposures) ~ "PFNA",
                               grepl("pfuna", exposures) ~ "PFUnA",
                               grepl("nmefosaa", exposures) ~ "NMeFOSAA"))

levels <- c("PFUnA", "PFHpS","PFHpA","NMeFOSAA","PFDA","PFNA","PFHxS", "PFOA", "PFOS" 
)
## coefficient plot
p1 <- cont_result %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid( ~ outcome, scales = "free") +
  ylab("Coefficient (95% CI)") +
  # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip()

ggsave(fs::path(dir_figure, "coef_plot_cont_outcome_full_group.png"), bg = "white",
       width = 8, height = 3, dpi = 300)

p2 <- cat_result %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = exp(estimate))) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = exp(conf.low),
                    ymax = exp(conf.high)),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid( ~ outcome, scales = "free") +
  ylab("Odds Ratio (95% CI)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 

ggsave(fs::path(dir_figure, "coef_plot_cat_outcome_full_group.png"), bg = "white",
       width =  8, height = 3, dpi = 300)

p3 <- cont_result1 %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid( ~ outcome, scales = "free") +
  ylab("Coefficient (95% CI)") +
  # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip()

ggsave(fs::path(dir_figure, "coef_plot_cont_outcome_full_group_with_extra_covar.png"), bg = "white",
       width =  8, height = 3, dpi = 300)

p4 <- cat_result1 %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = exp(estimate))) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = exp(conf.low),
                    ymax = exp(conf.high)),
                width = 0) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_grid( ~ outcome, scales = "free") +
  ylab("Odds Ratio (95% CI)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 

ggsave(fs::path(dir_figure, "coef_plot_cat_outcome_full_group_with_extra_covar.png"), bg = "white",
       width = 8, height = 3, dpi = 300)



# QG-comp, stratified by cases vs. controls -----
pfas_qgcomp <- pfas_names[!(pfas_names %in% c("pfas_br_pfos",
                                              "pfas_l_pfos",
                                              "pfas_l_pfoa"))]


library(qgcomp)

qgcomp_fun_cont <- function(outcome){
  (qgcomp_case <- qgcomp(as.formula(paste0(outcome, "~",
                                           paste0(paste(pfas_qgcomp, collapse = "+"),
                                                  "+sex+agebase+tx"))), 
                         expnms = pfas_qgcomp,
                         family = gaussian(),
                         q=5,
                         data = full_data |>
                           tidylog::filter(case1_control0 == 1)))
  
  
  (qgcomp_cont <- qgcomp(as.formula(paste0(outcome, "~",
                                           paste0(paste(pfas_qgcomp, collapse = "+"),
                                                  "+sex+agebase+tx"))),
                         expnms = pfas_qgcomp,
                         family = gaussian(),
                         q=2,
                         data = full_data |>
                           tidylog::filter(case1_control0 == 0)))
  
  return(list(case = qgcomp_case, control = qgcomp_cont))
}


qgcomp_fun_cat <- function(outcome){
  (qgcomp_case <- qgcomp(as.formula(paste0(outcome, "~",
                                           paste0(paste(pfas_qgcomp, collapse = "+"),
                                                  "+sex+agebase+tx"))), 
                         expnms = pfas_qgcomp,
                         family = binomial(),
                         q=NULL,
                         data = full_data |>
                           tidylog::filter(case1_control0 == 1)))
  
  
  (qgcomp_cont <- qgcomp(as.formula(paste0(outcome, "~",
                                           paste0(paste(pfas_qgcomp, collapse = "+"),
                                                  "+sex+agebase+tx"))),
                         expnms = pfas_qgcomp,
                         family = binomial(),
                         q=NULL,
                         data = full_data |>
                           tidylog::filter(case1_control0 == 0)))
  
  return(list(case = qgcomp_case, control = qgcomp_cont))
}

daystomic <- qgcomp_fun_cont(outcome[2])
daystomac <- qgcomp_fun_cont(outcome[3])
daystorapid <- qgcomp_fun_cont(outcome[5])
daystohyp <- qgcomp_fun_cont(outcome[7])
rapid <- qgcomp_fun_cat(outcome[4]) # doesn't work
mic <- qgcomp_fun_cat(outcome[1]) #
hyp0 <- qgcomp_fun_cat(outcome[6])

# "sample_id"      "case1_control0" "matched_cc"     "glyc"           "daystoglyc"    # "agebase"        "sex"            "dxtime"         "tx"             "codi"          
# "si_1_ins0"      "bmi"            "hb_a1c"         "log_trig"       "hdl"           
# "monthstoglyc"      "pfas_br_pfos"   "pfas_l_pfos"    "pfas_l_pfoa"    "pfas_nmefosaa" 































# PFDA and outcomes ----
# library(lmerTest)
# 
# library(CMAverse)
# 
# cmdag(outcome = "monthstoglyc", exposure = "PFDA", mediator = c("M1", "M2"), 
#       basec = c("C1", "C2"), postc = NULL, node = FALSE, text_col = "black")
# 
# 
# 
# cmtest_results <- cmest(
#   data = full_data |> tidylog::filter(case1_control0 == 1),
#   model = "rb", 
#   outcome = "monthstoglyc", 
#   exposure = "pfas_pfda",
#   mediator = c("bmi"), 
#   basec = c("sex", "agebase"), 
#   EMint = TRUE,
#   mreg = list("linear"), 
#   yreg = "linear",
#   astar = 0.1, a = 0.5, 
#   # mval = list(median(full_data$bmi)),
#   estimation = "imputation", 
#   # multimp = TRUE,
#   inference = "bootstrap", 
#   nboot = 20)
# 
# 
# cmtest_results
