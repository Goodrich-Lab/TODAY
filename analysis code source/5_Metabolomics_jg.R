# sandbox
source(here::here("!libraries.R"))
source(here::here("!directories.R"))
source(here::here("!load_clean_data.R"))
library(gplots)
library(ggrepel)

pfas_names_reduced <- pfas_names[-c(1:3)]

hist(full_data$daystomic)
# Scale proteins
full_data <- full_data |> 
  mutate(across(all_of(c(prot_names, met_names)), scale))

View(data_frame(c(met_names, prot_names)))

# Single Exposures, pfas as outcome ----
# owas
resout <- epiomics::owas(df = full_data, 
                         var = pfas_names, 
                         omics = c(met_names, prot_names), #,  
                         covars = c("agebase", "sex", 
                                    "est_creat_clear",
                                    "uacid"), 
                         var_exposure_or_outcome = "outcome")

# Mean Estimates
mean_estimates <- resout |> 
  group_by(feature_name) |>
  summarise(mean_est = mean(estimate), 
            min = min(estimate), 
            max = max(estimate))

max_per_compound <- resout |> 
  group_by(var_name) |>
  tidylog::filter(abs(estimate) == max(abs(estimate)))
  
# See which PFAS are associated with daystohyp
exp_out_res <- full_data |> 
  mutate(across(all_of(pfas_names), scale)) |>
  # tidylog::filter(sample_id != "3533017006-6107123") |>
  epiomics::owas(var = "daystohyp", 
                 covars = c("agebase", "tx"),
                 omics = pfas_names,
                 var_exposure_or_outcome = "outcome",
                 conf_int = TRUE)






# "Exposures":
# met_homovanillic_acid
# met_histidine
# met_hydroxyphenyllactic_acid

# Omics layers:
# L1: Carboxylic acids
pfas_pfda 
pfas_pfhpa
pfas_pfna
pfas_pfuna
pfas_pfoa

# L2: Sulfonic Acids
pfas_pfhps
pfas_pfhxs
pfas_pfos
  
# Outcomes:
# mic   
# daystomic      
# daystomac                               
# rapid          
# daystorapid  


View(data.frame(colnames(full_data)))

epiomics::volcano_owas(resout, annotate_ftrs = FALSE)

plot(resout$p_value, resout$test_statistic)
resoutfin <- tidylog::left_join(resout, prot_metadata,
                                by = c("feature_name" = "AptName")) |>
  clean_names() |>
  select(target_full_name, target, everything()) |>
  tidylog::filter(#p_value < 0.05,
                  # str_detect(var_name, "_l_", negate = TRUE), 
                  # str_detect(var_name, "_br_", negate = TRUE)
                  ) |>
  as_tibble()

# Pivot wider
resout_w <- resoutfin |> 
  mutate(var_name = str_remove(var_name, "pfas_") |> 
           toupper()) |>
  select(feature_name, target_full_name, 
         target:entrez_gene_symbol, var_name, test_statistic) |>
  pivot_wider(id_cols = feature_name:entrez_gene_symbol,
              names_from = var_name,
              values_from = test_statistic)
 # write_csv(resout_w, fs::path(dir_results, "proteins_to_pfas_test_stats.csv"))

 ## Heamap of results ----
# Filter dataset to only solute carrier proteins:
slc_prot <- resout_w |>
  tidylog::filter(str_detect(entrez_gene_symbol, "SLC") |
                  str_detect(tolower(target_full_name), "transport")) |>
  mutate(any_pfas_sig = feature_name %in% 
           resoutfin$feature_name[resoutfin$p_value<0.05])

# write_csv(slc_prot, fs::path(dir_results, "slc_proteins_teststat.csv"))

slc_mat <- slc_prot |>
  column_to_rownames("target_full_name") |>
  select(BR_PFOS:PFUNA) |>
  mutate(across(everything(), \(x) replace_na(x, 0))) |>
  as.matrix() |> t()

heatmap.2(slc_mat,
          col = greenred(100),
          trace = "none", )
#
RColorBrewer::display.brewer.all()




# mixtures, qgcomp ----
resout <- epiomics::owas_qgcomp(df = full_data, 
                                expnms = pfas_names[-c(1:4)],
                                omics = prot_names, 
                                covars = c("agebase", "sex", 
                                           "est_creat_clear",
                                           "uacid"), 
                                q = 3)

table(abs(resout$psi>1))
# Res
df <- tidylog::left_join(resout, prot_metadata,
                         by = c("feature" = "AptName")) |>
  clean_names() |>
  select(target_full_name, target, everything()) 
# Sig
marginally_sig <- df |> tidylog::filter(p_value < 0.001)

# Visualize results -------
df$ftr_names_label <- ifelse(df$p_value < 0.001, 
                             as.character(df$target), "")

# Make plot
(main_plot <- ggplot(df, 
                     aes(x = psi,
                         y = -log10(p_value), 
                         alpha =  p_value<0.001)) + 
    geom_point(aes(color = p_value<0.001)) + 
    geom_hline(yintercept = -log10(0.001), 
               linetype = 2, color = "grey70") + 
    scale_alpha_manual(values = c(.25, 1)) + 
    scale_color_manual(values = c("grey50", "darkred")) + 
    ylim(c(0,4)) + 
    ylab("-log10 p-value") + 
    theme(legend.position = "none")) + 
    geom_text_repel(aes(label = ftr_names_label),
                    min.segment.length = 0,  force = 5)


# Feature Annotation
ggsave(main_plot, 
       filename = fs::path(dir_results, "Jesse_K01_PFAS Proteomics.jpeg"), 
       width = 2.5, 
       height = 2.5, dpi = 600)
