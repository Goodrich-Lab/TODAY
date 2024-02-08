
source(here::here("!libraries.r"))
source(here::here("!directories.r"))
library(igraph)

# 1. ComptoxAI graph generated in python ----
g <- read.graph("~/Documents/comptoxai/comptox_ai/PFAS_prot_dkd_expanded_020524.graphml",
                format = "graphml")


# Calculate graph data 
node_metadata <- data.frame(
  geneSymbol = V(g)$id, 
  eigen_centrality = igraph::eigen_centrality(g)$vector,
  authority = igraph::authority_score(g)$vector, 
  harmonic_centrality = harmonic_centrality(g),
  hub_score = hub_score(g)$vector,
  page_rank = page_rank(g)$vector,
  power_centrality = power_centrality(g),
  strength = strength(g),
  subgraph_centrality = subgraph_centrality(g),
  betweenness = igraph::betweenness(g),
  degree = igraph::degree(g), 
  alpha_centrality = alpha_centrality(g), 
  eccentricity = eccentricity(g),
  coreness = coreness(g),
  closeness = closeness(g))
plot(node_metadata$authority)

# 2. Proteomics meet in middle results ----
pfas_proteomics <- read.csv(fs::path(dir_results, "med_res_df.csv"))

# Clean proteomics gene names to merge better
pfas_proteomics <- pfas_proteomics |> 
  dplyr::mutate(
    geneSymbol = case_when(EntrezGeneSymbol == "FGF8_3" ~ "FGF8", 
                           EntrezGeneSymbol == "HSPA1A_2" ~ "HSPA1A",
                           EntrezGeneSymbol == "CLEC4G_2" ~ "CLEC4G",
                           EntrezGeneSymbol == "HEY1_2" ~ "HEY1",
                           EntrezGeneSymbol == "GPNMB_4" ~ "GPNMB",
                           EntrezGeneSymbol == "RSPO2_1" ~ "RSPO2",
                           EntrezGeneSymbol == "IL27|EBI3" ~ "IL27",
                           EntrezGeneSymbol == "GDF11|MSTN" ~ "GDF11", 
                           TRUE ~ EntrezGeneSymbol)) 


## a. Combine data from proteomics and from comptoxai ----
pfas_prot_network <- tidylog::inner_join(node_metadata, 
                                         pfas_proteomics,
                                         by = "geneSymbol")

# Remove RSPO2, CHRNA5, FSCN1 because effect estimates are not possible
pfas_prot_network <- pfas_prot_network |>
  tidylog::filter(!(geneSymbol %in% c ("RSPO2", "CHRNA5", "FSCN1")))

# pivot data wider on mediation effect estimates
ppw <- pfas_prot_network |> 
  pivot_wider(names_from = Effect, 
              values_from = pe, 
              id_cols = c(geneSymbol:closeness)) 

# Rank the data values
ppw_rank <- ppw |> 
  tidylog::select(-geneSymbol) |> 
  tidylog::mutate_all(rank)

## b. correlation plots ----
(cor_matrix <- cor(ppw_rank, method = "spearman"))
corrplot::corrplot(cor_matrix, method = 'ellipse', type = 'upper')
# PerformanceAnalytics::chart.Correlation(ppw_rank, 
#                                         histogram=TRUE, pch=19,
#                                         method = "pearson")

# Test all associations 
resout <- epiomics::owas(
  df = ppw_rank, 
  var = colnames(node_metadata[,-1]), 
  omics = setdiff(colnames(ppw_rank),
                  colnames(node_metadata[,-1])), 
  var_exposure_or_outcome = "exposure")

cor(ppw_rank$hub_score, ppw_rank$Rpnie)
cor.test(ppw_rank$hub_score, ppw_rank$Rpnie)
ppw_rank$geneSymbol <- ppw$geneSymbol

# Plot relationships between graph characteristics and mediation results
ggplot(ppw_rank, aes(x = hub_score, y = Rpnie, label = geneSymbol)) +
  geom_point(aes(size = pm), shape = 21, color = "purple", fill = "white") + 
  ggrepel::geom_text_repel() +
  xlim(c(2, 32))
# 
# ggplot(pfas_prot_network, aes(x = authority, y = pe)) + 
#   geom_point() + 
#   facet_wrap(~Effect, scales = "free")

jag2::reformat_names("Nikos Stratakis,1 Augusto Anguita-Ruiz,1 Lorenzo Fabbri,1  Léa Maitre,1 Juan R. González,1  Sandra Andrusaityte,2 Xavier Basagaña,1 Eva Borràs,3 Hector C. Keun,4 Lida Chatzi,5 David Conti,5 Jesse Goodrich,5 Regina Grazuleviciene,2 Line Småstuen Haug,6 Barbara Heude,7 Rosemary McEachan,8 Mark Nieuwenhuijsen,1 Theano Roumeliotaki,9 Eduard Sabidó,3 Rémy Slama,10 Cathrine Thomsen,6 Jose Urquiza,1 Marina Vafeiadi,9 John Wright,8 Mariona Bustamante,1 Martine Vrijheid")
