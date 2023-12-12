
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))
library(gplots)
library(ggrepel)
library(qgcomp)

pfcas <- pfas_names[c(4,5,7,9,11)] 
pfsas <- pfas_names[c(6,8,12)] 

# Outcome Names -----
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
         across(all_of(c(pfas_names, met_names, prot_names, outcome_names$colname)), 
                \(x) scale(x) %>% as.numeric)) 


# Individual exposure outcome associations -------
tubular_biomarkers <- epiomics::owas(df = data_pfas_lg2_tmp, 
               var = outcome_names$colname, 
               omics = pfas_names, 
               covars = c("sex", "tx", "uacid"), 
               var_exposure_or_outcome = "outcome", 
               conf_int = TRUE) |> 
  left_join(outcome_names, by = c("var_name" = "colname"))


# Individual outcome analyses ----
qgcomp_mod <- function(prot_name, pfas_names_in_mix, q=3){
  qgcomp(as.formula(paste0(prot_name, "~", 
                           paste0(pfas_names_in_mix, collapse = "+"), #pfas_names[-c(1:3)]
                           "+ sex + tx + uacid")),
         expnms = pfas_names_in_mix,
         data = data_pfas_lg2_tmp, q = q)
}

# NGAL ----
# Carboxcilic acids
qgcomp_mod("seq.2836.68", pfcas)
# Sulfonic acids
qgcomp_mod("seq.2836.68", pfsas)
# all pfas
qgcomp_mod("seq.2836.68", pfas_names_in_mix = pfas_names[-c(1:3)], q = 4)

# NAG seq.15509.2 -----------------------------------------------------------
# Carboxcilic acids
qgcomp_mod("seq.15509.2", pfcas)
# Sulfonic acids
qgcomp_mod("seq.15509.2", pfsas)
# all pfas
qgcomp_mod("seq.15509.2", pfas_names_in_mix = pfas_names[-c(1:4, 6)], q = 4)

# Add mixture effect to data
(fit <- qgcomp_mod("seq.15509.2", pfas_names_in_mix = pfas_names[-c(1:4, 6)], q = 4))

tubular_biomarkers1 <- tubular_biomarkers %>% 
  add_row(var_name = "seq.15509.2", 
          feature_name = "Mixture",
          estimate = fit$psi,
          conf_low = fit$ci.coef[2],
          conf_high = fit$ci.coef[4]
          )

df <- tubular_biomarkers %>% filter (var_name == "seq.15509.2") %>%
  filter(!feature_name %in% c("pfas_nmefosaa","pfas_br_pfos","pfas_l_pfoa","pfas_l_pfos")) 

# Plot
tubular_biomarkers1 %>% 
  filter (var_name == "seq.15509.2") %>%
  filter(!feature_name %in% c("pfas_nmefosaa","pfas_br_pfos","pfas_l_pfoa","pfas_l_pfos", "pfas_pfhps", "pfas_pfhpa")) %>%
  mutate(feature_name = str_replace_all(feature_name, "pfas_", "") %>% toupper()) %>%
  mutate(feature_name = case_when(feature_name == "PFUNA"~"PFUnA",
                                  feature_name == "PFHXS"~"PFHxS",
                                  feature_name == "MIXTURE"~"Mixture",
                                  TRUE ~ feature_name)) %>%
  ggplot(aes(x = factor(feature_name, 
                        levels = c("PFDA", "PFOS", "PFOA",
                                       "PFNA", "PFUnA", "PFHxS","Mixture")),
                        y = estimate)) +
  geom_point(size = 1) +
  coord_flip()+
  geom_errorbar(aes(ymin = conf_low,
                    ymax = conf_high),
                width = 0) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey") +
  # ylab("Odds Ratio (95% CI)") +
  
  theme(panel.background = element_rect(fill="white"), 
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        axis.title.x=element_blank(),
        axis.line.x = element_line(color = "black", size = 0.8),
        axis.title.y=element_blank(),
        axis.line.y = element_line(color = "black", size = 0.8),
        axis.text = element_text(size =15),
        legend.title = element_blank(),
        plot.margin = margin(10, 10, 10, 10))+
  coord_flip()

ggsave(filename = fs::path(dir_figure, "Figure 3. panel A coef plot.jpg"),
       width = 2.5, height = 2.5, dpi=600, bg = "white")

# KIM-1 seq.9021.1 --------------------------------------------------------
# Carboxcilic acids
qgcomp_mod("seq.9021.1", pfcas)
# Sulfonic acids
qgcomp_mod("seq.9021.1", pfsas)
# all pfas
qgcomp_mod("seq.9021.1", 
           pfas_names_in_mix = pfas_names[-c(1:3)])


# IL-18 (seq.5661.15) ----
# Carboxcilic acids
qgcomp_mod("seq.5661.15", pfcas)
# Sulfonic acids
qgcomp_mod("seq.5661.15", pfsas)
# all pfas
qgcomp_mod("seq.5661.15", pfas_names_in_mix = pfas_names[-c(1:3)])


#  seq.11516.7 ----
# Carboxcilic acids
qgcomp_mod("seq.11516.7", pfcas)
# Sulfonic acids
qgcomp_mod("seq.11516.7", pfsas)
# all pfas
qgcomp_mod("seq.11516.7", pfas_names_in_mix = pfas_names[-c(1:3)])


# α-GST seq.17138.8 -----------------------------------------------------------
# Carboxcilic acids
qgcomp_mod("seq.17138.8", pfcas)
# Sulfonic acids
qgcomp_mod("seq.17138.8", pfsas)
# all pfas
qgcomp_mod("seq.17138.8", pfas_names_in_mix = pfas_names[-c(1:3)])


# α-GST seq.12446.49 -----------------------------------------------------------
# Carboxcilic acids
qgcomp_mod("seq.12446.49", pfcas)
# Sulfonic acids
qgcomp_mod("seq.12446.49", pfsas, q = 3)
# all pfas
qgcomp_mod("seq.12446.49", pfas_names_in_mix = pfas_names[-c(1:3)])



# R2 for POWER calculations -------------

mod1 <- lm(scale(seq.15509.2) ~ scale(pfas_pfhxs) + sex + tx + uacid, 
   data = data_pfas_lg2_tmp) |> glance()

mod2 <- lm(scale(seq.15509.2) ~ sex + tx + uacid, 
                   data = data_pfas_lg2_tmp) |> glance()

mod1$r.squared-mod2$r.squared


