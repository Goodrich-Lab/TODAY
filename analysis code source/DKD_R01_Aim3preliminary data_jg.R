# sandbox
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))
library(gplots)
library(ggrepel)
library(CMAverse)


# Names of Key Variables ----
pfas_min <-  c("pfas_pfna", "pfas_pfda", "pfas_pfhpa", 'pfas_pfoa', 
               "pfas_pfos", "pfas_pfhps", "pfas_pfhxs")

pfas_ca <- c("pfas_pfna", "pfas_pfda", "pfas_pfhpa", 'pfas_pfoa')
pfas_sa <- c("pfas_pfos", "pfas_pfhps", "pfas_pfhxs")

outcome_names <- as_tibble(matrix(c("seq.2836.68" ,  "NGAL",  "biomarker",
                                    "seq.9021.1"  ,  "KIM-1", "biomarker",
                                    "seq.5661.15" ,  "IL18",  "biomarker",
                                    "seq.11516.7" ,  "FABPL", "biomarker",
                                    "seq.17138.8" ,  "a-GST", "biomarker",
                                    "seq.12446.49",  "a-GST", "biomarker", 
                                    "seq.15509.2" ,  "NAG",   "biomarker",
                                    "daystomic"   ,  "daystomic"   , "clinical",
                                    "daystomac"   ,  "daystomac"   , "clinical",
                                    "daystorapid" ,  "daystorapid" , "clinical",
                                    "daystohyp"   ,  "daystohyp"   , "clinical"), 
                                  ncol=3, byrow = TRUE)) |>
  janitor::clean_names() |>
  rename(colname =  v1, 
         gene_outcome_name = v2, 
         class = v3)

# Log transform PFAS
data_pfas_lg2_tmp <- full_data |>
  mutate(across(all_of(pfas_names), log), 
         across(all_of(c(pfas_names, met_names, prot_names, outcome_names$colname[1:7])), 
                \(x) scale(x) %>% as.numeric)) 

# Remove Empty metabolite names
data_pfas_lg2 <- data_pfas_lg2_tmp |> 
  janitor::remove_constant()

met_names_removed <- colnames(data_pfas_lg2_tmp)[!colnames(data_pfas_lg2_tmp) %in% colnames(data_pfas_lg2)]

met_names_red <- setdiff(met_names, met_names_removed)


# 1: exposure --> PFAS associations ----
e_m_res <- epiomics::owas(df = data_pfas_lg2, 
                          omics = c("uacid", met_names_red), 
                          var = pfas_min,  
                          covars = c("sex", "est_creat_clear"), 
                          var_exposure_or_outcome = "outcome") |>
  rename(pfas = var_name, 
         exposure = feature_name, 
         em_est = estimate, 
         em_p = p_value) |> 
  select(pfas, exposure, em_est, em_p)

# 2: PFAS --> outcome associations ----
m_o_res <- epiomics::owas(df = data_pfas_lg2, 
                          omics = outcome_names$colname, 
                          var = pfas_min,  
                          covars = c("sex", "est_creat_clear", "tx"), 
                          var_exposure_or_outcome = "exposure") |>
  rename(pfas = var_name, 
         outcome = feature_name, 
         mo_est = estimate, 
         mo_p = p_value) |> 
  select(pfas, outcome, mo_est, mo_p) #|>
# filter(mo_p < 0.1); table(m_o_res$pfas)


# 3: exposure --> outcome associations ----
e_o_res <- epiomics::owas(df = data_pfas_lg2, 
                          omics = outcome_names$colname, 
                          var = c("uacid", met_names_red),  
                          covars = c("sex", "est_creat_clear", "tx"), 
                          var_exposure_or_outcome = "exposure") |>
  rename(exposure = var_name, 
         outcome = feature_name, 
         eo_est = estimate, 
         eo_p = p_value) |> 
  select(exposure, outcome, eo_est, eo_p)

# combine datasets 
full_res <- left_join(e_m_res, m_o_res, 
                      by= "pfas", 
                      relationship = "many-to-many") |>
  tidylog::left_join(e_o_res) |> 
  tidylog::left_join(outcome_names, by = c("outcome" = "colname"))

# multiply alpha times beta, calculate significance
# significance threshold: 
th = 0.2
full_res <- full_res |> 
  mutate(mediated_effect = em_est * mo_est,
         increased_risk_mo = if_else(class == "clinical", 
                                     mo_est < 0, 
                                     mo_est > 0), 
         increased_risk_eo = if_else(class == "clinical", 
                                     eo_est < 0, 
                                     eo_est > 0), 
         sig = case_when(em_p < th & mo_p < th ~ "Both", 
                         em_p < th | mo_p < th ~ "One", 
                         TRUE ~ "N.S.")) 

increased_risk <- full_res |>
  tidylog::filter(increased_risk_mo == TRUE,
                  increased_risk_eo == TRUE,
                  mediated_effect > 0 & class == "biomarker" |
                    mediated_effect < 0 & class == "clinical",
                  sig == "Both") |>
  droplevels()

table(increased_risk$pfas, increased_risk$gene_outcome_name)

# Uric acid
uacid_only <- full_res |> 
  filter(exposure == "uacid", 
         outcome == "seq.15509.2")

# Causal Mediation Analysis ----
### NAG -----
exposure_name = "uacid"
mediator_name = c("pfas_pfna", "pfas_pfhxs")

cmest_out_m.2 <- cmest(data = data_pfas_lg2,
                       model = "gformula",
                       outcome = "seq.15509.2",
                       exposure = exposure_name,
                       mediator = mediator_name, # note: can be multiple mediators
                       basec = c("sex", "est_creat_clear"), #"tx"
                       yreg  = "linear",
                       mreg  = list("linear",  "linear"),
                       mval  = list(quantile(data_pfas_lg2[[mediator_name[1]]], .2),
                                    quantile(data_pfas_lg2[[mediator_name[2]]], .2)
                       ),
                       a     = quantile(data_pfas_lg2[[exposure_name]], .25),
                       astar = quantile(data_pfas_lg2[[exposure_name]], .75),
                       EMint = TRUE, full = TRUE)
# Results
summary(cmest_out_m.2)
# summary(cmest_out_m.8)

ggcmest(cmest_out_m.2) +
  theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = .5)) +
  coord_flip()
# 
# 
### days to mic ----- 
cmest_out_surv <- cmest(data = data_pfas_lg2,
                        model = "gformula",
                        outcome = "daystohyp",
                        exposure = exposure_name,
                        mediator = mediator_name, # note: can be multiple mediators
                        basec = c("sex", "est_creat_clear", "seq.15509.2"), #
                        yreg  = "coxph",
                        mreg  = list("linear",  "linear"),
                        mval  = list(quantile(data_pfas_lg2[[mediator_name[1]]], .2),
                                     quantile(data_pfas_lg2[[mediator_name[2]]], .2)
                        ),
                        a     = quantile(data_pfas_lg2[[exposure_name]], .1),
                        astar = quantile(data_pfas_lg2[[exposure_name]], .9),
                        EMint = TRUE, full = TRUE)

summary(cmest_out_surv)

# Process results -----
# Results, NAG
cmest_effect_nag <- cbind(cmest_out_m.2$effect.pe, 
                          cmest_out_m.2$effect.ci.low,
                          cmest_out_m.2$effect.ci.high,
                          cmest_out_m.2$effect.pval) |>
  as_tibble(rownames = "effect") |> 
  mutate(outcome = "NAG")

# Results
cmest_effect_surv <- cbind(log(cmest_out_surv$effect.pe), 
                           log(cmest_out_surv$effect.ci.low),
                           log(cmest_out_surv$effect.ci.high),
                           cmest_out_surv$effect.pval) |>
  as_tibble(rownames = "effect") |> 
  mutate(outcome = "daystomic")

cmest_effect <- bind_rows(cmest_effect_nag, cmest_effect_surv)

# Rename columns
colnames(cmest_effect) <- c("effect", "pe", "ci.low", "ci.high", "pval", 'outcome')

# subset only key results
cm_df <- cmest_effect |>
  filter(effect %in% c("te", "tnde", "pnie", 
                       "Rte", "Rtnde", "Rpnie")) %>%
  mutate(effect = case_when(str_detect(effect, "te")  ~ "Total Effect", 
                            str_detect(effect, "nde") ~ "Direct Effect", 
                            str_detect(effect, "nie") ~ "Indirect Effect") |> 
           fct_relevel("Total Effect", 
                       "Indirect Effect",
                       "Direct Effect"))

# Plot Results -----
# NAG
(plot_nag <- ggplot(cm_df |> filter(outcome == "NAG"),
                    aes(x = pe, xmin = ci.low,
                        xmax = ci.high, y= effect)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50")  +
    geom_pointrange() +
    xlab("Change in NAG (SD)") +
    theme(axis.line.x = element_line(color = "black", size = 0.8),
          axis.title.y=element_blank(),
          # axis.title.x= element_blank(),
          axis.line.y = element_line(color = "black", size = 0.8),
          axis.text = element_text(size =10),
          plot.margin = margin(10, 10, 10, 10)) + 
  geom_text(x = .25,y= 1.65, label = "*", size = 11))

# microbl
(plot_mic <- ggplot(cm_df |> 
                      filter(outcome == "daystomic"),
                    aes(x = pe, xmin = ci.low, 
                        xmax = ci.high, y= effect)) + 
    geom_vline(xintercept = 0, 
               linetype = 2, color = "grey50")  +
    geom_pointrange() + 
    xlab("log(HR) DKD") + 
    theme(axis.line.x = element_line(color = "black", size = 0.8),
          axis.title.y=element_blank(),
          axis.title.x= element_text(hjust = 1),
          axis.line.y = element_line(color = "black", size = 0.8),
          axis.text = element_text(size =10),
          plot.margin = margin(10, 10, 10, 10)) )


# theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
# ylim(c(0, 10)) +
# coord_flip()

(plotout <- cowplot::plot_grid(plot_nag, #NULL,
                               plot_mic, #NULL,
                               ncol = 1, align = "v"))

ggsave(plotout, filename = fs::path(dir_figure, "Figure 7 Mediation plot.jpg"),
       width = 3, height = 3, dpi=600, bg = "white")
