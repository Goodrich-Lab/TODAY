# Playing around with qgcomp

qgcomp <- owas_qgcomp(df = data, 
                      expnms = pfas_names_all, 
                      omics = c(outcome_glu, outcome_biomaker), 
                      covars = c("case1_control0", "sex", "agebase", "serum_creat"),
                      confidence_level = 0.95,
                      q = 3) |>
  dplyr::select(feature, psi, lcl_psi, ucl_psi, p_value) %>%
  mutate(type = ifelse(feature %in% outcome_glu, "clinical", "biomarker")) %>%
  mutate(feature = rename_outcomes(feature))


reorder_outcome_name <- c("Hyperglycemia (HbA1c)", 
                          "β-cell function (C-Peptide oDI)",  
                          "Ins. Sensitivity (1/fast. ins)",
                          "NGAL","KIM-1", "IL18",
                          "FABPL","a-GST", "a-GST1","NAG") 
# Get color scale
colors <- RColorBrewer::brewer.pal(n=3, name = "Dark2")[1:2]

# Reorder outcome name and set color scale
qgcomp_res <- qgcomp %>%
  mutate(feature = factor(feature,
                          levels = reorder_outcome_name) %>%
           fct_rev(),
         color = if_else(type == "biomarker", 
                         colors[1], colors[2])) |>
  arrange(feature)

# Plot all outcomes 
(plotout <- qgcomp_res |>
    ggplot(aes(y = psi, 
               x = feature,
               color = color)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lcl_psi,
                      ymax = ucl_psi),
                  width = 0) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = .5) +
    theme() +
    scale_color_manual(values = colors) +
    ylab("PFAS Mixture (95% CI)") + 
    coord_flip() + 
    theme(legend.position = "None",
          axis.text.y = element_text(color = qgcomp_res$color, 
                                     face = "bold"), 
          axis.title.y = element_blank()))
  