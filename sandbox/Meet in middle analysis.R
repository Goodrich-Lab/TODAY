## Meet in the middle analysis
(analysis_pfas_names <- c(pfas_names_all, "score_pfsas", "score_pfcas")) 
# JG removed overall PFAS burden score, becuase the heatmap looks better without it


# 1. Exposure-Mediator regressions ----------------
result_em <- epiomics::owas(df = data_scaled, 
                            var = "pfas_pfna", #analysis_pfas_names,
                            omics = omic_names, #prot_names,
                            covars = covars,
                            var_exposure_or_outcome = "exposure",
                            family = "gaussian",
                            confidence_level = 0.95, 
                            conf_int = TRUE)

# Calculate significance based on effective number of tests (42)
result_em <- result_em |> 
  mutate(omic_layer = if_else(feature_name %in% prot_names, "Proteomics", "Metabolomics")) |> 
  group_by(var_name, omic_layer) |>
  mutate(sig_efnum = if_else(p_value < (0.05/42), "Sig", "Not sig"), 
         adjusted_pval = p.adjust(p_value, method = "fdr"),
         sig_fdr = if_else(adjusted_pval < 0.2, "Sig", "Not sig")) |> 
  ungroup()

# Check seq.8032.23 (MCCD1)

## filter significant proteins only -----
# get name of all sig proteins
omic_names_sig <- result_em |> 
  dplyr::filter(p_value < 0.05)
omic_names_sig <- unique(omic_names_sig$feature_name)

# Select key columns, filter to significant omics only
result_em_sig <- result_em %>% 
  tidylog::filter(feature_name %in% omic_names_sig) |> 
  dplyr::select(omic_layer, var_name, feature_name, test_statistic,
                estimate, conf_low, conf_high, p_value, sig_efnum) |> 
  rename(estimate_em = estimate, 
         conf_low_em =  conf_low, 
         conf_high_em = conf_high,
         p_value_em = p_value, 
         test_statistic_em = test_statistic) 

length(unique(result_em_sig$feature_name))


# 2. Mediator-outcome regressions ----------------
ind_vars <- tibble(exposure = omic_names)
dep_vars <- tibble(time = c("daystomic"), event = c("mic"))

# Get exposure and outcome combinations
mo_comb <- list(omic = ind_vars$exposure, event = dep_vars$event) %>%
  cross_df() %>% 
  left_join(dep_vars, by = "event")

# Create combination of exposures and outcomes
mo_res_models <- mo_comb %>%
  mutate(covar = str_c(covars, collapse = "+"),
         formula = str_c("Surv(", time, ",", event, ")", "~", omic, "+", covar))

# Run all models
mo_res_models$output <- map(mo_res_models$formula,
                            ~coxph(as.formula(.),
                                   data = data_scaled) %>%
                              tidy(., conf.int = TRUE))

# Modify output
mo_res_models <- mo_res_models %>%
  unnest(output) %>%
  filter(grepl("seq", term)|grepl("met",term))%>%
  mutate(HR = exp(estimate),
         exp_conf_high = exp(conf.high),
         exp_conf_low = exp(conf.low),
         sig = ifelse(p.value < 0.05, "Sig.", "Not Sig.")) 

## Filter sig results only ----
# Select key columns, filter significant only
result_mo_sig <- mo_res_models %>% 
  tidylog::filter(p.value < 0.05) %>%
  dplyr::select(omic, estimate, p.value) |> 
  rename(feature_name = omic, 
         estimate_mo = estimate, 
         p.value_mo = p.value)


# 3. Combine and analyze meet in middle --------------------
# Combine em and mo
mim_res <- tidylog::inner_join(result_em_sig, 
                               result_mo_sig, 
                               by = "feature_name")

# examine meet in middle proteins
mim_res <- mim_res |> 
  mutate(emxmo = estimate_em*estimate_mo, 
         max_p = if_else(p_value_em < p.value_mo, p.value_mo, p_value_em))

# Create cleaned PFAS name variable
mim_res <- mim_res |> 
  rename(pfas = var_name) |>
  rename_pfas() |>
  rename(pfas_name = pfas) 


## Combine with proteomics metadata ------
mim_res2 <- mim_res %>% 
  tidylog::left_join(prot_metadata, by = c("feature_name" = "AptName")) |> 
  tidylog::filter(Organism == "Human") 

## Filter to features only overlapping in same direction (didn't work)
# same_dir_results <- mim_res2 |> 
#   filter(emxmo > 0, p_value_em < 0.05)
# 
# length(unique(mim_res2$feature_name))
# length(unique(same_dir_results$feature_name))
# 
# mim_res2 <- mim_res2 %>% 
#   tidylog::filter(feature_name %in% same_dir_results$feature_name)

# Examine results
table(mim_res2$highest_expression, mim_res2$second_highest_expression)
table(mim_res2$category)
table(mim_res2$omic_layer)
length(unique(mim_res2$feature_name))


# 4. Heatmap ---------------
library(ComplexHeatmap)
library(circlize)

## reformat data for heatmap ----
# Select correct variables, filter out the pfas burden score
ind_pfas_prot <-  mim_res2 |> 
  tidylog::select(EntrezGeneSymbol, test_statistic_em,
                  pfas_name, estimate_em, p_value_em)


length(unique(mim_res2$EntrezGeneSymbol))

table(is.na(mim_res2$EntrezGeneID))
length(unique(mim_res2$feature_name))

### Pivot effect estimates wider to create data matrix -----
# effect estimates
pfas_prot_ee_w <- ind_pfas_prot |> 
  dplyr::select(-p_value_em, -estimate_em) |>
  pivot_wider(names_from = pfas_name,  
              values_from = test_statistic_em, 
              id_cols = EntrezGeneSymbol) |>
  column_to_rownames("EntrezGeneSymbol")

# p-values
pfas_prot_p_w <- ind_pfas_prot |> 
  dplyr::select(-estimate_em) |>
  pivot_wider(names_from = pfas_name,  
              values_from = p_value_em, 
              id_cols = EntrezGeneSymbol) |>
  column_to_rownames("EntrezGeneSymbol")


# split individual PFAS and the PFAS burden score
ind_pfas_ee_prot_w   <- pfas_prot_ee_w |> dplyr::select(-contains("Score")) |> as.matrix()
pfas_score_ee_prot_w <- pfas_prot_ee_w |> dplyr::select( contains("Score")) |> as.matrix()
ind_pfas_p_prot_w    <- pfas_prot_p_w  |> dplyr::select(-contains("Score")) |> as.matrix()
pfas_score_p_prot_w  <- pfas_prot_p_w  |> dplyr::select( contains("Score")) |> as.matrix()

### Create protein clusters ----
row_dend = hclust(dist(ind_pfas_ee_prot_w), method="complete") # row clustering

# Determine optimal number of clusters 
scale <- scale(ind_pfas_ee_prot_w)
factoextra::fviz_nbclust(scale, kmeans, method = "gap_stat")
# optimal number of clusters: k = 4

# Set color scale
max_val <- max(abs(min(ind_pfas_ee_prot_w)), max(ind_pfas_ee_prot_w))
max_val_round <- ceiling(max_val * 4) / 4
col_range <- c(-1*max_val_round, 0, max_val_round)
col_fun = circlize::colorRamp2(col_range, c("blue", "white", "red"))
lgd = Legend(col_fun = col_fun, title = "Effect estimate")

# Heatmap
set.seed(121)
# Individual PFAS and proteins 
pfas_prot_hm = Heatmap(ind_pfas_ee_prot_w, 
                        name = "Effect estimate",
                        column_title_gp = gpar(fontsize = 26),
                        column_title = "Single PFAS analysis",
                        col = col_fun,
                        row_km = 3,
                        show_heatmap_legend = FALSE,
                        show_column_names = TRUE, 
                        show_row_names = TRUE,
                        column_names_gp = gpar(fontsize = 8),
                        row_names_gp = gpar(fontsize = 8),
                        column_names_rot = 45,
                        cell_fun = function(j, i, x, y, w, h, fill){
                          if(ind_pfas_p_prot_w[i, j] < 0.05/42) {
                            grid.text("**", x, y, vjust = .75)
                          } else if(ind_pfas_p_prot_w[i, j] < 0.05) {
                            grid.text("*", x, y, vjust = .75)
                          }
                        })

# Add PFAS burden score column
(pfas_prot_hm_fin <- pfas_prot_hm +
    Heatmap(pfas_score_ee_prot_w, 
            name = "PFAS Burden",
            show_heatmap_legend = FALSE ,
            col = colorRamp2(c(-4, 0, 4), c( "blue","white", "red")),
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 8),
            column_names_gp = gpar(fontsize = 8), 
            column_names_rot = 45,
            heatmap_legend_param = list(title = "Test statistic", 
                                        title_gp = gpar(fontsize = 20),
                                        legend_direction = "horizontal",
                                        legend_width = unit(5, "cm"), 
                                        title_position = "lefttop"),
            cell_fun = function(j, i, x, y, w, h, fill){
              if(pfas_score_p_prot_w[i, j] < 0.05/42) {
                grid.text("**", x, y, vjust = .75)
              } else if(pfas_score_p_prot_w[i, j] < 0.05) {
                grid.text("*", x, y, vjust = .75)
              }
            }))


# 5. Heatmap of Protein Tissue Expression ----
# Select correct variables, remove duplicates with group_by and slicehead
mim_res_for_outcome <-  mim_res2 |> 
  group_by(EntrezGeneSymbol) |> 
  tidylog::slice_head() |>
  ungroup() |> 
  arrange(match(EntrezGeneSymbol, rownames(pfas_score_ee_prot_w))) |> 
  column_to_rownames("EntrezGeneSymbol") 


# Convert to matrix for heatmap 
prot_tissue_expression <- mim_res_for_outcome |>
  tidylog::select(contains("ratio"),
                  -high_second_high_ratio)  |> 
  as.matrix()

# Reorder the dataframe to match the other heat maps
prot_tissue_expression_mat <- prot_tissue_expression[rownames(pfas_score_ee_prot_w),]

# Descriptives on tissue expression
high_expressed_tissues <- prot_tissue_expression |> 
  as_tibble() |> 
  summarise(across(where(is.numeric), ~max(., na.rm = TRUE))) |>
  pivot_longer(cols = everything()) |>
  filter(value > 0.4, !(name %in% c("ratio_placenta", "ratio_testis")))

# Filter to only columns which have a value > 0.5
prot_tissue_fin <- prot_tissue_expression_mat[,high_expressed_tissues$name] 


# Heatmap with tissue protein expression  -------
(pfas_tissue_expression_prot_hm <- pfas_prot_hm_fin +
    Heatmap(prot_tissue_fin, 
            name = "Tissue Expression",
            show_heatmap_legend = FALSE,
            col = circlize::colorRamp2(c(0,.5,1), c("white", "darkgreen", "black")),
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 8),
            # Add percent expression in the kidney: 
            right_annotation = rowAnnotation(kidney = anno_barplot(prot_tissue_expression_mat[,"ratio_kidney"])) ,
            column_names_gp = gpar(fontsize = 8),  
            column_names_rot = 45))

# Heatmap with protein to outcome association ----------- 
(pfas_tissue_expression_prot_hm <- pfas_prot_hm_fin +
    Heatmap(as.matrix(mim_res_for_outcome %>% dplyr::select(estimate_mo)), 
            name = "Tissue Expression",
            show_heatmap_legend = FALSE,
            col = circlize::colorRamp2(c(-1.9,0,1.9), c("blue", "white", "red")),
            show_row_names = TRUE,
            row_names_gp = gpar(fontsize = 8),
            # Add percent expression in the kidney: 
            right_annotation = rowAnnotation(kidney = anno_barplot(prot_tissue_expression_mat[,"ratio_kidney"])) ,
            column_names_gp = gpar(fontsize = 8),  
            column_names_rot = 45))
