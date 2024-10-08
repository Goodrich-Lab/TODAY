# Sandbox analysis
source(here::here("!libraries.r"))
source(here::here("!directories.r"))
source(here::here("!load_clean_data.R"))
library(igraph)

# 1. ComptoxAI graph generated in python ----
# 35 protein graph
g <- read_graph(fs::path(dir_results,
                         "ComptoxAI",
                         "PFAS_prot_dkd_expanded_020924.graphml"),
                format = "graphml")

# g <- read_graph(fs::path(dir_results,
#                          "ComptoxAI",
#                          "PFAS_prot_dkd_all_prot_032724.graphml"),
#                 format = "graphml")


# Calculate graph data 
node_metadata <- tibble(
  geneSymbol = V(g)$id, 
  eigen_centrality = igraph::eigen_centrality(g)$vector,
  authority = igraph::authority_score(g)$vector,
  harmonic_centrality = harmonic_centrality(g, mode = "out"),
  hub_score = hub_score(g)$vector, 
  page_rank = page_rank(g)$vector,
  strength = strength(g),
  # power_centrality = power_centrality(g),
  # subgraph_centrality = subgraph_centrality(g),
  # alpha_centrality = alpha_centrality(g),
  betweenness = igraph::betweenness(g),
  degree = igraph::degree(g),
  eccentricity = eccentricity(g),
  coreness = coreness(g),
  closeness = closeness(g)
)
# plot(node_metadata$harmonic_centrality_out, node_metadata$harmonic_centrality_in)


# 2. Proteomics meet in middle results ----
# Meet in middle results:
meet_in_middle_effects <- read.csv(fs::path(dir_results, "med_res_df_020924.csv"))
# Results of all mediation analyses
pfas_proteomics_res <- read.csv(fs::path(dir_results,
                                         "pairwise_mediation_result.csv"))

# Combine protein metadata and proteomics results
# pfas_proteomics_all <- meta_pro |> 
#   # tidylog::select(-c(category:ratio_salivary_gland)) |>
#   tidylog::select(-c(SeqId:EntrezGeneID),
#                   -c(Organism:ratio_salivary_gland)) |>
#   tidylog::left_join(pfas_proteomics_res, by = c("EntrezGeneSymbol","AptName")) 

# For now, since we didnt include an interaction, we can exclude the .T and 
# .C vars from the mediation results
pfas_proteomics_all <- pfas_proteomics_res |>
  tidylog::select(-contains(".T"), 
                  -contains(".C"), 
                  # -contains(".pval"),
                  -contains(".hi"), 
                  -contains(".lo")) |>
  rename_all(~str_remove(., ".avg"))

## a. Filter to select only 35 sig from meet in middle -------
pfas_proteomics_35 <- pfas_proteomics_all |>
  tidylog::filter(AptName %in% meet_in_middle_effects$feature_name)

hist(pfas_proteomics_35$PMed)

# Clean proteomics gene names to merge better
pfas_proteomics_35 <- pfas_proteomics_35 %>%
  # pfas_proteomics_all <- pfas_proteomics_all %>%
  dplyr::mutate(
    EntrezGeneSymbol = case_when(EntrezGeneSymbol == "FGF8_3" ~ "FGF8", 
                                 EntrezGeneSymbol == "HSPA1A_2" ~ "HSPA1A",
                                 EntrezGeneSymbol == "CLEC4G_2" ~ "CLEC4G",
                                 EntrezGeneSymbol == "HEY1_2" ~ "HEY1",
                                 EntrezGeneSymbol == "GPNMB_4" ~ "GPNMB",
                                 EntrezGeneSymbol == "RSPO2_1" ~ "RSPO2",
                                 EntrezGeneSymbol == "IL27|EBI3" ~ "IL27",
                                 EntrezGeneSymbol == "GDF11|MSTN" ~ "GDF11", 
                                 TRUE ~ EntrezGeneSymbol)) |>
  # Remove the "_#" pattern from EntrezGeneSymbol
  mutate(EntrezGeneSymbol = gsub("_\\d+", "", EntrezGeneSymbol)) %>%
  group_by(EntrezGeneSymbol) |> 
  tidylog::filter(abs(PMed) == max(abs(PMed))) |> 
  slice_head() |> 
  ungroup() #|>
# # Handle the "|" in EntrezGeneSymbol
# mutate(geneSymbol = strsplit(geneSymbol, "\\|")) %>%
# unnest_longer(geneSymbol)


## b. Combine data from proteomics and from comptoxai ----
ppw <- tidylog::inner_join(node_metadata, 
                           pfas_proteomics_35, #pfas_proteomics_all
                           by = c("geneSymbol" = "EntrezGeneSymbol"))

# Determine top 30% of effects:
ppw2 <- ppw |> 
  mutate(
         # top_30_hub = ifelse(hub_score > quantile(hub_score, 0.66), 1, 0), 
         top_30_hc = ifelse(harmonic_centrality > quantile(harmonic_centrality, 0.66), 1, 0),
         top_30_prop_med = ifelse(PMed > quantile(PMed, 0.66), 1, 0))

ppw_top_30_pct <- ppw2 |> 
  tidylog::filter(top_30_hc == 1, 
                  top_30_prop_med == 1)
table(ppw2$top_30_hc, ppw2$top_30_prop_med)


## c. Filter to only positive percent mediated ---------------------------------
ppw_pos <- ppw2 |> 
  tidylog::filter(PMed > 0) |>
  mutate(ACME_abs = abs(ACME))

# Rank the data values
ppw_rank <- ppw_pos |> 
  tidylog::select(-geneSymbol, -top_30_hc, -top_30_prop_med) |> #, -c(AptName:ensembl_gene_id)
  tidylog::mutate_all(~rank(.)) %>%
  bind_cols(ppw_pos[,c(1,23,24)], .) 

# Add color groups
ppw_rank1 <- ppw_rank %>% 
  mutate(Group = case_when(top_30_hc == 1&top_30_prop_med == 1 ~ "#f54278",
                                 top_30_hc == 1 ~ "#fc8c03",
                                 top_30_prop_med == 1 ~ "#4287f5",
                                 top_30_hc != 1&top_30_prop_med != 1 ~ "white"),
         )


# Get name of all mediation variables of interest
mediation_res_colnames <- ppw_rank |> 
  dplyr::select(-nobs, -geneSymbol, -AptName, 
                -all_of(colnames(node_metadata[,-1])), 
                -contains(".pval")) |>
  colnames()

# Plot relationships between graph characteristics and mediation results
# (cor_plot <- ggplot(ppw_rank, aes(x = hub_score, y = Rpnie, label = geneSymbol)) +
(cor_plot <- ppw_rank1 %>% ggplot(aes(x = harmonic_centrality, 
                                     y = PMed, 
                                     label = geneSymbol,
                                     color = Group)) +
    ggrepel::geom_label_repel(force = 1, box.padding = .5, min.segment.length = 0) +
    # scale_fill_manual(values = c("#4287f5","#f54278","#fc8c03", "white"))+
    scale_color_manual(values = c("#4287f5","#f54278","#fc8c03", "grey10"))+
    # annotate("rect", xmin = -Inf, xmax = 15, ymin = -Inf, ymax = 15, fill = "lightgreen", alpha = 0.3) +
    # annotate("rect", xmin = 15, xmax = Inf, ymin = -Inf, ymax = 15, fill = "lightblue", alpha = 0.3) +
    # annotate("rect", xmin = -Inf, xmax = 15, ymin = 15, ymax = Inf, fill = "lightyellow", alpha = 0.3) +
    # annotate("rect", xmin = 15, xmax = Inf, ymin = 15, ymax = Inf, fill = "lightpink", alpha = 0.3) +
    geom_point(size = 2, shape = 21, color = "black", 
               fill = "grey50") + 
    # geom_point(aes(size = pm), shape = 21, color = "black", fill = "grey50") +
    # geom_point(aes(size = pm_quartile, fill = pm_quartile), shape = 21, color = "black") +
    scale_size(name = "Percent\nMediated\n(Rank)") +
    xlab("Harmonic Centrality from ComptoxAI (Rank)") +
    ylab("Proportion Mediated\nfrom TODAY Study (Rank)") + 
  xlim(c(2, 32))  +
    theme(legend.position = "none",
          legend.title = element_blank()
          ) 
    # +
    # guides(fill=guide_legend(ncol=2))
) 

quantile(ppw_rank$harmonic_centrality, 2/3)


ggsave(cor_plot, 
       filename = fs::path(dir_figure, "Fig3b_Centrality vs PM V2.jpg"), 
       height = 6, width = 6)

