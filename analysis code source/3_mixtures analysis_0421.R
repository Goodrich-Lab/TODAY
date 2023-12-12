# Exposure outcome associations
library(ggbeeswarm)
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))

library(epiomics)

pfas_exposure <- pfas_names[str_detect(pfas_names, "_br", negate = TRUE)]
pfas_exposure <- pfas_exposure[str_detect(pfas_exposure, "_l", negate = TRUE)]


pfas_names_all <- c("pfas_pfos", 
                    "pfas_pfoa", 
                    "pfas_pfhxs", 
                    "pfas_pfna", 
                    "pfas_pfda", 
                    "pfas_nmefosaa", 
                    "pfas_pfhps", 
                    "pfas_pfhpa", 
                    "pfas_pfuna")



pfas_names_reduced <- c("pfas_pfos", 
                        "pfas_pfoa", 
                        "pfas_pfhxs", 
                        "pfas_pfna", 
                        "pfas_pfda", 
                        "pfas_nmefosaa")

# outcome_names <- c("time_to_glyc_scld", 
#                    "codi", 
#                    "si_1_ins0",
#                    "hb_a1c") 
outcome_name<- c("mic", "daystomic", 
             # "mac",
             "daystomac", 
             "rapid", "daystorapid", "hyp0", "daystohyp")

# outcome_names_full <- c("Time to glycemic failure", 
#                         "β-cell function (C-Peptide oDI)", 
#                         "Ins. Sensitivity (1/fast. ins)",
#                         "Hyperglycemia (HbA1c)") 

scld_dat <- full_data |>
  mutate(across(all_of(outcome_names[c(2,3,5,7)]), 
                ~scale(.) |> as.numeric()))

# qgcomp --------
qgcomp_res_cont <- owas_qgcomp(df = scld_dat, 
                          expnms = pfas_names_all, 
                          omics = outcome_names[c(2,3,5,7)], 
                          covars = c("case1_control0",
                                     "tx",
                                     "agebase",
                                     "sex"),
                          confidence_level = 0.90,
                          q = 3)

qgcomp_res_cat <- owas_qgcomp(df = scld_dat, 
                               expnms = pfas_names_all, 
                               omics = outcome_names[c(1,4,6)], 
                               covars = c("case1_control0",
                                          "tx",
                                          "agebase",
                                          "sex"),
                               confidence_level = 0.90,
                               q = NULL)

qgcomp_res <- qgcomp_res_cont %>% bind_rows(qgcomp_res_cat)
# qgcomp_res$feature_names  <- outcome_names_full

# qgcomp_res <- qgcomp::qgcomp(time_to_glyc_scld ~ pfas_pfos+ 
#                                pfas_pfoa+ pfas_pfhxs+ 
#                                pfas_pfna+ pfas_pfda+ pfas_nmefosaa+
#                                case1_control0 + 
#                                dxtime+
#                                agebase,
#                              expnms = pfas_names_reduced, 
#                              q = 3,
#                              data = full_data)
# plot(qgcomp_res)
reorder_outcome_name <- outcome_names
# Plot all outcomes 
(plotout <- qgcomp_res |>
  ggplot(aes(y = psi, 
             x = fct_relevel(feature, 
                             reorder_outcome_name) |> fct_rev())) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lcl_psi,
                    ymax = ucl_psi),
                width = 0, size = 1) +
  geom_hline(yintercept = 0, linetype = 2) +
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_text(size = 14, vjust = .5)) +
  ylab("PFAS Mixture") + 
  coord_flip())

ggsave(plotout, 
        filename = fs::path(dir_figure, "QGComp pfas mixture TODAY_0421.jpeg"),
       height = 2, width = 4.5, bg = "white")

# plot HbA1c ----------------------
# (ttgf_a <- qgcomp_res |> 
#    filter(feature == "hb_a1c") |>
#    ggplot(aes(y = psi, x = feature)) + 
#    geom_point(size = 3) + 
#    geom_errorbar(aes(ymin = lcl_psi, 
#                      ymax = ucl_psi), 
#                  width = 0, size = 1) +
#    geom_hline(yintercept = 0, linetype = 2) +
#    theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), 
#          axis.title = element_blank()) +
#    ylim(c(-.15, 1.3)))
# 
# (ttgf_b <- qgcomp_res |> 
#     filter(feature == "hb_a1c") |>
#     pivot_longer(cols = c(coef_pfas_pfos:coef_pfas_pfuna), 
#                  names_to = "pfas") |>
#     mutate(pfas = str_remove(pfas, "coef_pfas_") |> toupper(), 
#            pfas = fct_reorder(as.factor(pfas), -1*value)) |>
#     ggplot(aes(x = pfas, y = value)) + 
#     geom_bar(stat = "identity") + 
#     geom_hline(yintercept = 0, linetype = 2) + 
#     theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1), 
#           axis.title = element_blank()) +
#     ylim(c(-.15, 1.3))) 
# 
# 
# cowplot::plot_grid(NULL, NULL,  
#                    ttgf_a, ttgf_b, 
#                    align = "h",
#                    nrow = 2,
#                    rel_widths = c(.3, .75), 
#                    rel_heights = c(.05, 1),
#                    labels = c("A. Mixture", 
#                               "B. PFAS Contributions"))

# WQS -------------------------
library(gWQS)
wqs_result <- gWQS::gwqs(hb_a1c ~ wqs + 
                           case1_control0 + 
                           dxtime + 
                           agebase,
                         mix_name = pfas_names_all, 
                         data = full_data %>%
                           select(all_of(c("hb_a1c", 
                                           "case1_control0", 
                                           "dxtime",
                                           "agebase",
                                           pfas_names_all))),
                         # family=gaussian,
                         family=gaussian,
                         seed = 1001)

gwqs_barplot(wqs_result)
summary(wqs_result)

full_data$si_1_ins0





# BKMR --------------------------------------------------------------------
library(bkmr)


# Fit model
fitkm <- kmbayes(y = y, 
                 Z = Z,
                 X = X,  
                 iter = 10000, 
                 verbose = FALSE, varsel = TRUE)

fitkm

# Get predicted response univar
pred.resp.univar <- PredictorResponseUnivar(fit = fitkm)
ggplot(pred.resp.univar, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se)) + 
  geom_smooth(stat = identity) + 
  facet_wrap(~ variable) +
  ylab(h(z))


# Get Risks overall 
risks.overall <- OverallRiskSummaries(fit = fitkm, 
                                      y = y, Z = Z, X = X, 
                                      qs = seq(0.15, 0.85, by = 0.05), 
                                      q.fixed = 0.15, method = exact)



ggplot(risks.overall, 
       aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd)) + 
  geom_pointrange() + 
  geom_hline(yintercept =  0, color = red)


