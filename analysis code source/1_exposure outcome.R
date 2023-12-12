# Exposure outcome associations
library(ggbeeswarm)
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))

pfas_names_reduced <- pfas_names[-c(1:3)]


# All potential outcomes: 
# MIC (dichot): Development to mic
# DAYSTOMIC (continuous): Days to microalbumineria
# MAC (dichot): macroalbuminera
# DAYSTOMAC (continuous): Days to MAC
# RAPID (dichot): 
# DAYSTORAPID (continuous):  
# HYP0 (dichot):
# DAYSTOHYP (continuous):

ggplot(full_data, aes(x = as.factor(case1_control0), 
                      y = time_to_glyc_scld)) +
  geom_beeswarm()
  

# PFAS boys vs. girls
ggplot(full_data, aes(x = as.factor(case1_control0), 
                      y = pfas_pfhps)) +
  geom_beeswarm() +
  facet_wrap(~sex)

# # Conditional logistic regression -----
# exwas <- full_data |>
#   mutate(across(all_of(pfas_names), log2)) |>
#   # tidylog::filter(sample_id != "3533017006-6107123") |>
#   epiomics::owas_clogit(cc_status = "case1_control0",
#                         cc_set = "matched_cc",
#                         covars = c("agebase", "tx"),
#                         omics = pfas_names,
#                         conf_int = TRUE) |>
#   mutate(
#     or = exp(estimate),
#     lcl_or = exp(conf_low),
#     ucl_or = exp(conf_high),
#     or_ci = jag2::effest_ci(or, lcl_or, ucl_or))
# 
# # Plot
# epiomics::volcano_owas(exwas)
# 
# 
# # Plot PFOS in cases vs controls
# ggplot(full_data, 
#        aes(x = as.factor(case1_control0), 
#            y = pfas_pfhps)) +
#   geom_beeswarm()
# 
# 
# # Table:
# pfas_gm <- full_data |>
#   select(case1_control0, all_of(pfas_names)) |>
#   pivot_longer(cols = all_of(pfas_names), names_to = "feature_name") |>
#   group_by(case1_control0, feature_name) |> 
#   summarise(geomean_sd = jag2::fungsd(value, n.digits = 2, na.rm = TRUE)) |>
#   pivot_wider(values_from = c(geomean_sd), 
#               names_from = case1_control0, names_prefix = "gm_sd_")
# 
# 
# # Full results 
# full_results <- left_join(pfas_gm, 
#                           exwas) |>
#   arrange(desc(gm_sd_0)) |>
#   mutate(feature_name = str_remove(feature_name, "pfas_")) |>
#   select(feature_name, gm_sd_1,gm_sd_0,  or_ci, p_value)
# 
# write_csv(full_results, fs::path(dir_results, "PFAS outcome associations.csv"))
# 
# 
# ## Sensitivity: Girls only ----
# exwas_girls <- full_data |>
#   mutate(across(all_of(pfas_names), log2)) |> 
#   tidylog::filter(sex == "female") |>
#   epiomics::owas_clogit(cc_status = "case1_control0",
#                         cc_set = "matched_cc", 
#                         covars = c("agebase", "tx"),
#                         omics = pfas_names, 
#                         conf_int = TRUE)
# 
# # Volcano plot
# epiomics::volcano_owas(exwas_girls)


# Linear regression, stratified by cases vs. controls -----
# (exwas_lm_case <- full_data |>
#    tidylog::filter(case1_control0 == 1) |>
#    mutate(across(all_of(pfas_names), log2)) |>
#    epiomics::owas(var = "monthstoglyc",
#                   var_exposure_or_outcome = "outcome",
#                   covars = c("sex", "agebase", "tx"),
#                   omics = pfas_names, 
#                   conf_int = TRUE) |>
#    mutate(beta_ci_case = jag2::effest_ci(estimate, conf_low, conf_high),
#           p_value_case = p_value))
# 
# exwas_lm_control <- full_data |>
#   tidylog::filter(case1_control0 == 0) |>
#   mutate(across(all_of(pfas_names), log2)) |>
#   epiomics::owas(var = "monthstoglyc",
#                  var_exposure_or_outcome = "outcome",
#                  covars = c("sex", "agebase", "tx"),
#                  omics = pfas_names, 
#                  conf_int = TRUE) |>
#   mutate(beta_ci_cont = jag2::effest_ci(estimate, conf_low, conf_high), 
#          p_value_cont = p_value) ; filter(exwas_lm_control, 
#                                           estimate == min(estimate)) 


# Coefficient plot_HW
## restrict to 9 pfas
exwas_lm_case_plot <- exwas_lm_case %>%
  filter(!feature_name %in%
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", feature_name) ~ "PFOS",
                               grepl("pfoa", feature_name) ~ "PFOA",
                               grepl("pfda", feature_name) ~ "PFDA",
                               grepl("pfhpa", feature_name) ~ "PFHpA",
                               grepl("pfhps", feature_name) ~ "PFHpS",
                               grepl("pfhxs", feature_name) ~ "PFHxS",
                               grepl("pfna", feature_name) ~ "PFNA",
                               grepl("pfuna", feature_name) ~ "PFUnA",
                               grepl("nmefosaa", feature_name) ~ "NMeFOSAA"))


exwas_lm_control_plot <- exwas_lm_control %>%
  filter(!feature_name %in%
           c("pfas_br_pfos", "pfas_l_pfoa", "pfas_l_pfos")) %>%
  mutate(pfas_name = case_when(grepl("pfos", feature_name) ~ "PFOS",
                               grepl("pfoa", feature_name) ~ "PFOA",
                               grepl("pfda", feature_name) ~ "PFDA",
                               grepl("pfhpa", feature_name) ~ "PFHpA",
                               grepl("pfhps", feature_name) ~ "PFHpS",
                               grepl("pfhxs", feature_name) ~ "PFHxS",
                               grepl("pfna", feature_name) ~ "PFNA",
                               grepl("pfuna", feature_name) ~ "PFUnA",
                               grepl("nmefosaa", feature_name) ~ "NMeFOSAA"))

levels <- c("PFUnA", "PFHpS","PFHpA","NMeFOSAA","PFDA","PFNA","PFHxS", "PFOA", "PFOS" 
               )
## coefficient plot for rapid progression
p1 <- exwas_lm_case_plot %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Coefficient (95% CI)") +
  # ggtitle("Coefficient Plot (Binary Logistic Regression)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 


p2 <- exwas_lm_control_plot %>%
  ggplot(aes(x = factor(pfas_name, levels = levels),y = estimate)) +
  geom_point(size = 1) +
  coord_flip() +
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Coefficient (95% CI)") +
  theme(text = element_text(size = 10),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        axis.line.x = element_line(color = "black"),
        axis.line.y = element_line(color = "black")) +
  coord_flip() 

cowplot::plot_grid(p1, p2, labels = c("A. Rapid Progression Y-T2D","B.Slow Progression Y-T2D"),
                   label_size = 12 , vjust = 0.5, hjust = 0) + 
  theme(plot.margin = unit(c(0.7, 0.5, 0.5, 0.5), "cm"))

ggsave(fs::path(dir_figure, "coef_plot.png"), bg = "white",
       width = 6, height = 3, dpi = 300)
# Volcano plot
epiomics::volcano_owas(exwas_lm_case)
epiomics::volcano_owas(exwas_lm_control)

# Plot PFDA
full_data |>
  # tidylog::filter(case1_control0 == 1) |>
  ggplot(aes(x = pfas_pfda, 
             y = monthstoglyc)) +
  geom_point() + 
  facet_wrap(~case1_control0, scales =  "free_y")


## Table of results ----
# Table:
full_results_stratified <- exwas_lm_case |> 
  select(feature_name, beta_ci_case, p_value_case) |>
  full_join(exwas_lm_control |>
              select(feature_name, beta_ci_cont, p_value_cont) ) |>
  mutate(feature_name = str_remove(feature_name, "pfas_"), 
         across(contains("p_value"), ~round(., 2))) |>
  filter(str_detect(feature_name, "l_", negate = TRUE), 
         str_detect(feature_name, "br_", negate = TRUE))

write_csv(full_results_stratified, 
          fs::path(dir_results,
                   "Stratified case control analysis.csv"))


# QG-comp, stratified by cases vs. controls -----
pfas_qgcomp <- pfas_names[!(pfas_names %in% c("pfas_br_pfos",
                                              "pfas_l_pfos",
                                              "pfas_l_pfoa"))]


library(qgcomp)
(qgcomp_case <- qgcomp(as.formula(paste0("monthstoglyc~",
                                     paste(pfas_qgcomp, collapse = "+"),
                                     "+sex+agebase+tx")), 
                   expnms = pfas_qgcomp,
                   q=5,
                   data = full_data |>
                     tidylog::filter(case1_control0 == 1)))

plot(qgcomp_case)

(qgcomp_cont <- qgcomp(as.formula(paste0("monthstoglyc~",
                                     paste(pfas_qgcomp, collapse = "+"),
                                     "+sex+agebase+tx")),
                   expnms = pfas_qgcomp,
                   q=2,
                   data = full_data |>
                     tidylog::filter(case1_control0 == 0)))

# "sample_id"      "case1_control0" "matched_cc"     "glyc"           "daystoglyc"    # "agebase"        "sex"            "dxtime"         "tx"             "codi"          
# "si_1_ins0"      "bmi"            "hb_a1c"         "log_trig"       "hdl"           
# "monthstoglyc"      "pfas_br_pfos"   "pfas_l_pfos"    "pfas_l_pfoa"    "pfas_nmefosaa" 

# Volcano plot
epiomics::volcano_owas(exwas_lm_case)
epiomics::volcano_owas(exwas_lm_control)

# Plot PFDA
full_data |>
  # tidylog::filter(case1_control0 == 1) |>
  ggplot(aes(x = pfas_pfda, 
             y = monthstoglyc)) +
  geom_point() + 
  facet_wrap(~case1_control0, scales =  "free_y")


## Table of results ----
# Table:
full_results_stratified <- exwas_lm_case |> 
  select(feature_name, beta_ci_case, p_value_case) |>
  full_join(exwas_lm_control |>
              select(feature_name, beta_ci_cont, p_value_cont) ) |>
  mutate(feature_name = str_remove(feature_name, "pfas_"), 
         across(contains("p_value"), ~round(., 2))) |>
  filter(str_detect(feature_name, "l_", negate = TRUE), 
         str_detect(feature_name, "br_", negate = TRUE))

write_csv(full_results_stratified, 
          fs::path(dir_results,
                   "Stratified case control analysis.csv"))


# Cross Sectional ----
pfas_names
outcome_names_xs <- c("u_alb_creat", "uacid", "est_creat_clear")

exwas <- full_data |>
  mutate(across(all_of(pfas_names), \(x) scale(x) |> as.numeric()),
         across(all_of(outcome_names_xs), \(x) scale(x) |> as.numeric())) |>
  # tidylog::filter(sample_id != "3533017006-6107123") |>
  epiomics::owas(var = outcome_names_xs, 
                 omics = pfas_names, 
                 covars = c("agebase", "sex"),
                 var_exposure_or_outcome = "exposure", 
                 conf_int = TRUE) |>
  mutate(
    or = exp(estimate),
    lcl_or = exp(conf_low),
    ucl_or = exp(conf_high),
    or_ci = jag2::effest_ci(or, lcl_or, ucl_or))

# Plot
epiomics::volcano_owas(exwas, annotate_ftrs = FALSE)


## Cross sectional, non-linear ----
library(mgcv)
nrow <- 3
ncol <- 3
par(mfrow = c(nrow, ncol), mar = c(4, 4, 2, 1))
full_data$est_creat_clear
full_data$uacid


# Loop over the dependent variables in pafs_names
for (paf in pfas_names_reduced) {
  # Construct the formula
  formula_string <- paste("log2(", paf, ")", "~ s(log2(uacid), bs = \"tp\") + agebase + sex")
  model_formula <- as.formula(formula_string)
  
  # Fit the GAM model
  gam_model <- gam(model_formula, data = full_data)
  
  # Plot the association
  plot(gam_model, scheme = 1, main = paf)
}



#
#
 




















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
