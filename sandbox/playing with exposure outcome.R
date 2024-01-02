


# 1) PFAS Burden Score --------
# PFAS names
c("pfas_pfna",
  "pfas_pfda",
  "pfas_pfhpa",
  "pfas_pfoa",
  "pfas_pfos",
  "pfas_pfhps",
  "pfas_pfhxs",
  "pfas_nmefosaa",
  "pfas_pfuna")   

## A. Create categorical PFAS which are needed for the score ------
# This has all been moved to "!load clean data" 
# data_analysis <- full_data %>%
#   mutate_at(.vars = vars(all_of(pfas_names_all[c(1,4,5,7,8)])), 
#             .funs = list(
#               quartile = ~as.integer(cut(., quantile(., probs = seq(0, 1, 0.25)), 
#                                        include.lowest = TRUE
#               )))) %>%
#   mutate(pfas_pfda_detected = ifelse(pfas_pfda<0.05|pfas_pfda==0.2/sqrt(2), 1, 2)) %>%
#   mutate_at(.vars = vars(c("pfas_pfhpa","pfas_pfuna")),
#             .funs = list(detected = ~ifelse(.<0.05|.==0.1/sqrt(2), 1, 2))) %>%
#   tidylog::mutate_at(.vars = vars(c("pfas_pfhps", "pfas_pfhpa")), #, 
#             .funs = list(dichotomous = ~cut(., 
#                                         quantile(., probs = seq(0, 1, .5)), 
#                                         include.lowest = TRUE) %>% 
#                            as.integer()))
# 
# # Select just the categorical PFAS
# burden_score_pfas <- data_analysis %>% 
#   dplyr::select(contains("detected"), 
#                 contains("quartile"), 
#                 contains("dichotomous"), 
#                 -pfas_pfhpa_dichotomous) %>% 
#   as.data.frame()
# # Select PFCAs
# burden_score_pfcas <- burden_score_pfas |> 
#   dplyr::select(pfas_pfhpa_detected, pfas_pfoa_quartile, pfas_pfna_quartile, pfas_pfda_detected, pfas_pfuna_detected)
# 
# # Select PFSAs
# burden_score_pfsas <- burden_score_pfas |> 
#   dplyr::select(pfas_pfos_quartile, pfas_pfhxs_quartile, pfas_pfhps_dichotomous)
# 
# ## B. Calculate Scores ------
# eap.pfas <- ltm::factor.scores(grm(burden_score_pfas), 
#                                method="EAP", 
#                                resp.patterns = burden_score_pfas)$score.dat$z1
# 
# eap.pfcas <- ltm::factor.scores(grm(burden_score_pfcas), 
#                                 method="EAP", 
#                                 resp.patterns = burden_score_pfcas)$score.dat$z1
# eap.pfsas <- ltm::factor.scores(grm(burden_score_pfsas), 
#                                 method="EAP", 
#                                 resp.patterns = burden_score_pfsas)$score.dat$z1
# 
# # Plot scores
# plot(grm(burden_score_pfcas), type = "IIC",legend = TRUE,
#      xlab = "PFAS Burden", main = "", cx = "topright", cex = 0.6)
# 
# # Add burden scores back into data
# data_analysis1 <- data_analysis %>% 
#   bind_cols(score = eap.pfas) %>%  
#   bind_cols(score_pfcas = eap.pfcas) %>% 
#   bind_cols(score_pfsas = eap.pfsas) 
# 
# ## C. Create Categorical Scores  ------
# data <- data_analysis1 %>% 
#   mutate_at(.vars = vars(c(all_of(pfas_names_all), score)),
#             .funs = list(median = ~ifelse(.<median(.),"0", "1")
#             )) %>%
#   mutate(score_tertile  = cut(score, quantile(score, probs = seq(0, 1, 1/3)), include.lowest = TRUE) %>% as.integer(), #%>% as.character(),
#          score_quartile = cut(score, quantile(score, probs = seq(0, 1, 1/4)), include.lowest = TRUE) %>% as.integer(), #%>% as.character(),
#          score_quintile = cut(score, quantile(score, probs = seq(0, 1, 1/5)), include.lowest = TRUE) %>% as.integer()  #%>% as.character()
#   ) %>%
#   mutate_at(.vars = vars(all_of(pfas_names_all)), .funs = ~  scale(.) %>% as.numeric(.)) # #############################################################
# 
# ## D. Create dummy vars, scale outcomes -----
# data <- data |>
#   fastDummies::dummy_cols(select_columns = c('sex'),
#                           remove_selected_columns = FALSE, 
#                           remove_most_frequent_dummy = TRUE) |>
#   mutate_at(.vars = vars(c(outcome_glu)), 
#             .funs = ~scale(.) %>% as.numeric(.))


# 2. Run Models ------------  
## A. Set up for analysis ------
# Get the name of all PFAS exposure variables that were just created 
ind_vars <- data %>% 
  dplyr::select(#contains("median"),
    #contains("tertile"),
    #contains("quartile"),
    contains("score"), 
    contains("_sum"),
    all_of(pfas_names_all)) |> 
  colnames() %>%  
  tibble(exposure = .)

# Get the name of all outcome variables created 
dep_vars <- tibble(time = c(#"yrstohyp2", "yrstomic2",   #"yrstoglyc2",
  "daystohyp", "daystomic"),   #"daystoglyc",
  event = c(#"hyp", "mic",    #"glyc", 
    "hyp", "mic"))   #"glyc", 

# Get dataframe of all exposure outcome combinations
eo_comb <- list(pfas = ind_vars$exposure, event = dep_vars$event) %>%
  cross_df() %>% 
  left_join(dep_vars, by = "event", relationship = "many-to-many")

# B. Set covars ---------
covars2 <- c("sex_male", "agebase", "serum_creat", "dxtime")
# original covars:  "case1_control0", "sex_male", "agebase", "serum_creat"   
# colnames(data[,1:30])


# run_models <- function(covars2){
# Get the formula for all models
models <- eo_comb %>%
  mutate(covar = paste0("sex_male+",str_c(covars2, collapse = "+")),
         formula = str_c("Surv(", time , ",", event, ")", "~", pfas, "+", covar))

## c. Run the models -----
models$output <- map(models$formula,
                     ~coxph(as.formula(.), data = data) %>%
                       tidy(., conf.int = TRUE))

# Clean up results
pfas_survival_models <- models %>%
  unnest(output) %>%
  tidylog::filter(grepl("score", term)| grepl("pfas_", term)| grepl("_sum", term)&!grepl("quintile",term)) %>%
  mutate(HR = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         sig = ifelse(p.value < 0.05, "Sig.", "Not Sig.")) %>%
  rename_pfas() %>%
  mutate(pfas = ifelse(pfas %in% levels, pfas, term))
# Note: 1 unit increase in pfas resulting in xx-fold increase in hazard.


# Plot results
(p <- pfas_survival_models %>%
    ggplot(aes(x = term,y = estimate, color = sig)) +
    geom_point(size = 1) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width = 0) +
    geom_hline(yintercept = 0, linetype = 2) +
    facet_grid( ~ time, scales = "free") +
    ylab("Log HR (95% CI)") +
    theme(text = element_text(size = 10),
          axis.title.y = element_blank(),
          panel.background = element_rect(fill="white"),
          strip.background = element_rect(fill = "white"),
          axis.line.x = element_line(color = "black"),
          axis.line.y = element_line(color = "black"),
          legend.position = "none") +
    coord_flip() +
    scale_color_manual(values = c("grey", "red")))

# }
pfas_survival_models |> filter(str_detect(term, "pfcas")) |> dplyr::select(pfas, event, p.value)

