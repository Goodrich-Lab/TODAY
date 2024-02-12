
source(here::here("!libraries.r"))
source(here::here("!directories.r"))
library(igraph)

# 1. ComptoxAI graph generated in python ----
g <- read.graph(fs::path(dir_results, 
                         "ComptoxAI",
                         "PFAS_prot_dkd_expanded_020924.graphml"),
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
pfas_proteomics <- read.csv(fs::path(dir_results, "med_res_df_020924.csv"))

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
length(unique(pfas_prot_network$geneSymbol))
# Remove RSPO2, CHRNA5, FSCN1 because effect estimates are not possible
# pfas_prot_network <- pfas_prot_network |>
#   tidylog::filter(!(geneSymbol %in% c ("CHRNA5", "FSCN1")))

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
(cor_plot <- ggplot(ppw_rank, aes(x = hub_score, y = Rpnie, label = geneSymbol)) +
    geom_hline(yintercept = 15, linetype = 3) + 
    geom_vline(xintercept = 15, linetype = 3) + 
    ggrepel::geom_label_repel(force = 1, box.padding = .5, min.segment.length = 0) +
  geom_point(aes(size = pm), shape = 21, color = "black", fill = "grey50") + 
  scale_size(name = "Percent\nMediated\n(Ranking)") + 
    xlab("Hub Score (Ranking)") + 
    ylab("Pure Natural Indirect Effect (Rank)") +
  xlim(c(2, 32)))
 
# ggplot(pfas_prot_network, aes(x = authority, y = pe)) + 
#   geom_point() + 
#   facet_wrap(~Effect, scales = "free")



# 3. Plot Graph ------


V(g)$node_color <- case_when(
  V(g)$type == "Disease"                ~ "purple",  
  V(g)$type == "mediating-protein"      ~ "Light Blue",
  V(g)$type == "non-identified protein" ~ "grey50", 
  V(g)$type == "PFAS"                   ~ "red")
  

str(g)


plot.igraph(g, layout = layout_with_fr(g), 
            vertex.color=V(g)$node_color, 
            vertex.frame.color="#555555", 
            vertex.label=V(g)$id, vertex.label.color="black")
            edge.
            vertex.label.cex=.7) 
                    vertex.size = 6, vertex.label = NA, edge.curved=.1)


# With NetworkD3


g_d3 <- igraph_to_networkD3(g, group = members)

# Create force directed network plot
forceNetwork(Links = karate_d3$links, Nodes = karate_d3$nodes, 
             Source = 'source', Target = 'target', 
             NodeID = 'name', Group = 'group')


# With GGraph
library(ggraph)
