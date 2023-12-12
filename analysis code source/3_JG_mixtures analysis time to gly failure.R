# Exposure outcome associations
library(ggbeeswarm)
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))

library(epiomics)



# full_data |>
#   group_by(tx) |> 
#   summarise(mean(daystomac), sd(daystomac))
table(full_data$tx)
full_data$pfas_pfos

full_data |>
  filter(tx != 2) %>%
  glm(hyp ~ tx*pfas_pfna, 
      data = ., family = "binomial") |>
  tidy() |> mutate(or = exp(estimate))


# Set exposures and outcomes ----
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

## Outcomes ----
outcome_names <- c("time_to_glyc_scld", 
                   "codi", 
                   "si_1_ins0",
                   "hb_a1c", 
                   # "daystomic", 
                   "est_creat_clear",
                   "daystomac")
# "uacid", 
# "u_alb_creat", 
# "serum_creat"

## Outcome names ----

rename_outcomes <- function(x){
  nms <- case_when(
    x=="time_to_glyc_scld" ~ "Time to glycemic failure", 
    x=="codi" ~ "β-cell function (C-Peptide oDI)",
    x=="si_1_ins0" ~  "Ins. Sensitivity (1/fast. ins)",
    x=="hb_a1c" ~ "Hyperglycemia (HbA1c)", 
    x=="daystomic" ~ "Time to microalbuminuria", 
    x=="est_creat_clear" ~ "Estimated Creatinine Clearance",
    x=="daystomac" ~ "Time to macroalbuminuria",
    x=="uacid" ~ "Uric Acid",
    x=="u_alb_creat" ~ "u_alb_creat",
    x=="serum_creat" ~ "Serum Creatanine")
}


## Scale data ----
scld_dat <- full_data |>
  mutate(across(all_of(outcome_names), 
                ~scale(.) |> as.numeric()))


# All single exposure outcome associations ------
(single_eo_res <- owas(df = scld_dat, 
                       var = outcome_names, 
                       var_exposure_or_outcome = "outcome",
                       omics = pfas_names_all, 
                       covars = c("case1_control0",
                                  "dxtime",
                                  "agebase", 
                                  "tx", 
                                  "sex"
                                  # "est_creat_clear"
                       ),
                       confidence_level = 0.90) |>
   mutate(sig = if_else(p_value < 0.1, "*", ""))) |>
  filter(p_value < 0.05) 

# qgcomp --------
## Glucose outcomes ----
(qgcomp_glu <- owas_qgcomp(df = scld_dat, 
                           expnms = pfas_names_all, 
                           omics = outcome_names[1:4], 
                           covars = c("case1_control0",
                                      "dxtime",
                                      "agebase", 
                                      "sex"
                                      # "tx"
                           ),
                           confidence_level = 0.90,
                           q = 3) |>
   mutate(outcome_type = "pancreas")) |>
  select(feature, psi, p_value)

## Kidney outcomes ----
(qgcomp_kid1 <- owas_qgcomp(df = scld_dat, 
                            expnms = pfas_names_all, omics = outcome_names[6], 
                            covars = c("dxtime","agebase","tx","est_creat_clear"),
                            confidence_level = 0.90,q = 3) |>
   mutate(outcome_type = "kidney")) |> 
  select(feature, psi, p_value) 


## Kidney outcomes ----
(qgcomp_kid2 <- owas_qgcomp(df = scld_dat, 
                            expnms = pfas_names_all, omics = outcome_names[5], 
                            covars = c("dxtime","agebase","tx"),
                            confidence_level = 0.90,q = 3) |>
   mutate(outcome_type = "kidney")) |> 
  select(feature, psi, p_value)

# Combine results
qgcomp_res <- bind_rows(qgcomp_glu, qgcomp_kid1, qgcomp_kid2) |>
  mutate(feature_names=rename_outcomes(feature))



# Check if creatanine clearance is driven by high levels
# qgcomp_res <- full_data |>
#   mutate(
#     high_creat_clear = case_when(
#       sex == "male" & est_creat_clear > 150 ~ 1, 
#       sex == "female" & est_creat_clear > 130 ~ 1, 
#       TRUE ~ 0)) %>%
#   qgcomp::qgcomp(high_creat_clear ~ pfas_pfos+
#                    pfas_pfoa+ pfas_pfhxs+
#                    pfas_pfna+ pfas_pfda+ pfas_nmefosaa+
#                    # case1_control0 +
#                    dxtime+
#                    tx+
#                    agebase,
#                  expnms = pfas_names_reduced,
#                  q = 3,
#                  family = "binomial",
#                  data = .)


reorder_outcome_name <- c("Hyperglycemia (HbA1c)", 
                          "β-cell function (C-Peptide oDI)",  
                          "Ins. Sensitivity (1/fast. ins)",
                          "Time to glycemic failure", 
                          # "Time to microalbuminuria", 
                          "Estimated Creatinine Clearance",
                          "Time to macroalbuminuria") 
# Get color scale
colors <- RColorBrewer::brewer.pal(n=3, name = "Dark2")[1:2]

# Reorder outcome name and set color scale
qgcomp_res_fin <- qgcomp_res |>
  mutate(feature_names = fct_relevel(feature_names, 
                                     reorder_outcome_name) |> 
           fct_rev(), 
         color = if_else(outcome_type == "kidney", 
                         colors[1], colors[2])) |>
  arrange(feature_names)

# Plot all outcomes 
(plotout <- qgcomp_res_fin |>
    ggplot(aes(y = psi, 
               x = feature_names, 
               color = color)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lcl_psi,
                      ymax = ucl_psi),
                  width = 0, size = .75) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = .5) +
    theme(axis.title.y = element_blank(), 
          axis.text.y = element_text(size = 14, vjust = .5)) +
    scale_color_manual(values = colors) + 
    ylab("PFAS Mixture") + 
    coord_flip() + 
    theme(legend.position = "None", 
          axis.text.y = element_text(color = qgcomp_res_fin$color, 
                                     face = "bold")))


ggsave(plotout,
       filename = fs::path(dir_figure, 
                           "QGComp pfas mixture TODAY_0425.jpeg"),
       height = 2.25, width = 5, bg = "white")

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






# Survival analysis- DKD --------
# Notes- the only PFAS consistently associated with macroalbuminuria is pfas_pfna 
library(survival)
library(ggfortify)
full_data$m
data <- full_data |> 
  tidylog::select(-all_of(prot_names), -contains("met_")) |>
  fastDummies::dummy_cols(select_columns = c('tx', 'sex'),
                          remove_selected_columns = TRUE, 
                          remove_most_frequent_dummy = TRUE) |>
  mutate(yrstohyp = (dxtime*12 + daystohyp)/365, 
         yrstohyp2 = yrstohyp + dxtime/12, 
         yrstomic = (daystomic)/365, 
         yrstomic2 = (yrstomic+dxtime/12),
         # if_else(yrstomac > 8.5, 8.5, yrstomac), 
         mic2 = if_else(yrstomic > 8.5, 0, 1), 
         hyp2 = if_else(yrstohyp > 8.5, 0, 1), 
         across(all_of(pfas_names), scale))

data$hyp2

hist(data$daystohyp)

## no boot ---- 
pfas_all <- pfas_names[-c(1:3)]
pfas_min <-  c("pfas_pfna", "pfas_pfda", "pfas_pfhpa", 'pfas_pfoa', 
               "pfas_pfos", "pfas_pfhps", "pfas_pfhxs")

pfas_ca <- c("pfas_pfna", "pfas_pfda", "pfas_pfhpa", 'pfas_pfoa')
pfas_sa <- c("pfas_pfos", "pfas_pfhps", "pfas_pfhxs")

### hyperperfusion ----
qc.survfit <- qgcomp::qgcomp.cox.noboot(
  survival::Surv(yrstohyp2, hyp) ~ ., 
  expnms = pfas_min,
  data = data |> 
    select(all_of(pfas_min), 
           u_alb_creat,
           case1_control0,
           # dxtime, 
           yrstohyp2, hyp),
  q = 2)
qc.survfit


### microalbuminuria -----
qc.survfit <- qgcomp::qgcomp.cox.noboot(
  survival::Surv(yrstomic2, mic) ~ ., 
  expnms = pfas_ca,
  data = data |> 
    select(all_of(pfas_ca), 
           # tx_1, tx_2,
           sex_male,
           # uacid,
           # est_creat_clear,
           u_alb_creat,
           dxtime,
           # codi,
           case1_control0,
           yrstomic2, 
           mic),
  q = 2); qc.survfit
plot(qc.survfit)
str(qc.survfit)
weights <- c(qc.survfit$pos.weights*qc.survfit$pos.size, 
             -1*qc.survfit$neg.weights*qc.survfit$neg.size)
data.frame(coef = weights) |>
  rownames_to_column("pfas")

#
table(data$hyp0, data$hyp)
#

## Bootstrap ----
qc.survfit2_nomale <- qgcomp::qgcomp.cox.boot(
  survival::Surv(yrstohyp2, hyp) ~ ., 
  expnms = pfas_min,
  data = data |>
    select(all_of(pfas_min), 
           u_alb_creat,
           # sex_male,
           case1_control0,
           yrstohyp2, hyp),
  # daystomac, mac),
  q = 2, B=5*32, MCsize=20000, parallel=TRUE, parplan=TRUE)

qc.survfit2_nomale
plot(qc.survfit2_nomale)
# qc.survfit2 # with male as covar
# p <- plot(qc.survfit2, suppressprint = TRUE) 
p <- plot(qc.survfit2_nomale, suppressprint = TRUE, linewidth = 10) 
(out <- p + scale_color_manual(values = c("00FFFF", "black")) + 
    scale_linetype_manual(values = c(0,1,1)) + 
    geom_hline(yintercept = 0.5) + 
    labs(y = "DKD Free (%)") +
    theme(legend.position = "none"))

# Save plot
ggsave(out, 
       filename = fs::path(dir_results, "Jesse_K01_KM_curve_PFAS_dkd.jpeg"), 
       width = 1.5, 
       height = 1.5, dpi = 600)


View(data |> select(daystofailure, hyp0, hyp))



## Macroalbuminuria ----
### PFNA ----
# qc.survfit2 <- qgcomp::qgcomp.cox.boot(
#   survival::Surv(daystomac, mac) ~  pfas_pfna + #pfas_pfna + pfas_pfos +
#     tx_1 + tx_2 + sex_male + case1_control0 + dxtime,
#   expnms = c( 'pfas_pfna'),
#   data = data,# |> filter(case1_control0 == 1),
#   q = 4, B=32, MCsize=1000, parallel=TRUE, parplan=TRUE)
# plot(qc.survfit2)
# survival::cox.zph(qc.survfit2$fit)
# qc.survfit2

View(data |> select(daystomac))
### PFAS MIN ----
qc.survfit_mac <- qgcomp::qgcomp.cox.boot(
  # survival::Surv(daystomic, mic) ~ ., 
  survival::Surv(yrstomac2, mac2) ~  .,
  expnms = pfas_min,
  data = data |> 
    select(all_of(pfas_min), 
           sex_male, 
           # case1_control0, 
           yrstomac2, mac2),
           # daystomic, mic),
  q = 2, B=320, MCsize=20000, parallel=TRUE, parplan=TRUE)
qc.survfit_mac
plot(qc.survfit_mac)
survival::cox.zph(qc.survfit_mac$fit)
qc.survfit_mac

data$y
# p <- plot(qc.survfit2, suppressprint = TRUE) 
p <- plot(qc.survfit_mac, suppressprint = TRUE, linewidth = 10) 
(mac_plot <- p + scale_color_manual(values = c("black", "00FFFFFF")) +  
    # scale_linewidth_manual(values = 5) +
    # geom_hline(yintercept = 0.5, linetype = 2, color = "grey50") + 
    scale_linetype_manual(values = c(0,1,1)) + 
    coord_cartesian(xlim = c(0, 8)) +
    labs(y = "DKD Free (%)") +
    theme(legend.position = "none"))

# Save plot
ggsave(mac_plot, 
       filename = fs::path(dir_results, "Jesse_K01_KM_curve_PFAS_MAC.jpeg"), 
       width = 1.5, 
       height = 1.5, dpi = 600)


# Other------
full_data[,1:20]
sfit <- survfit(Surv(daystofailure, hyp)~sex_male, data=data)

autoplot(qc.survfit2_nomale$msmfit$y )

hist(data$daystofailure)
View(data |> select(daystofailure, hyp0, hyp))
