# lucidus with multi-omics In Parallel Sankey Diagram ----


# ----------- reorder LUCID model by specifying reference cluster ------------
# note: only works for K = 2 in each omic layer
# reference = c(1,1,1,2)
# lucidus_fit <- fit
reorder_lucidM <- function(lucidus_fit,
                           reference = NULL) {
  if(is.null(reference)) {
    warning("no reference specified, return the original model")
    return(lucidus_fit)
  }
  
  n_omic <- length(reference)
  
  # reorder beta
  GtoX <- lucidus_fit$res_Beta$Beta
  lucidus_fit$res_Beta$Beta <- lapply(1:n_omic, function(i) {
    (-1)^(reference[i] - 1) * GtoX[[i]] # if reference = 1, no changes; 
    # if reference = 2, flip the reference and negate the estimates
  })
  # reorder mu
  XtoZ <- lucidus_fit$res_Mu_Sigma$Mu
  lucidus_fit$res_Mu_Sigma$Mu <- lapply(1:n_omic, function(i) {
    x <- c(1, 2) # order of clusters
    if(reference[i] == 2) {
      x <- c(2, 1)
      -1*XtoZ[[i]][, x]
    } else{
      XtoZ[[i]][, x]
    }
  }) 
  # reorder gamma
  XtoY <- lucidus_fit$res_Delta$Delta$mu
  XtoY[1] <- XtoY[1] + sum(XtoY[-1] * (reference - 1)) # reference level using the new reference
  XtoY[-1] <- (-1)^(reference - 1) * XtoY[-1] # if reference = 2, flip the estimates
  lucidus_fit$res_Delta$Delta$mu <- XtoY
  
  # return the object using the new reference
  return(lucidus_fit)
}

# Example: 
# fit_new <- reorder_lucidM(lucidus_fit = fit,
#                           reference = c(1, 1))


# Colors ----
col_pal <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")

# Color pallet for sankey
sankey_colors <- matrix(c("exposure", col_pal[6],
                          "lc1",      col_pal[1],
                          "lc2",      col_pal[2],
                          "lc3",      col_pal[3],
                          "lc4",      col_pal[4],
                          "layer1",   col_pal[1],
                          "layer2",   col_pal[2],
                          "layer3",   col_pal[3],
                          "layer4",   col_pal[4],
                          "outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red", 
                          "neg_clus_to_out", "#e4e5f2"), 
                        byrow = TRUE, nrow = 14) |> 
  as_tibble() |>
  rename("domain" = V1, 
         "range" = V2)


# sankey diagram function ----
# only works for K= 2
# lucidus_fit <- fit_new
plot_lucidM <- function(lucidus_fit, 
                        sankey_colors, 
                        text_size = 15){
  
  K          <- lucidus_fit$K
  dimG       <- lucidus_fit$res_Beta$Beta[[1]] |> ncol()-1
  n_layers   <- length(lucidus_fit$res_Beta$Beta)
  # Get list of mu 
  mu_lst <- purrr::map(lucidus_fit$res_Mu_Sigma$Mu, 
                       ~as_tibble(.x, rownames = "name")) 
  names(mu_lst) <- paste0("layer", c(1:n_layers))
  dimZ       <- purrr::map(mu_lst, ncol) |> as.numeric()-1
  n_features <- purrr::map(mu_lst, nrow) |> as.numeric()
  names(n_features) <- paste0("layer", c(1:n_layers))
  names_features <- bind_rows(mu_lst, .id = "color_group") |> 
    rowwise() |>
    mutate(sum = sum(abs(V1)+abs(V2))) |>
    group_by(color_group) |> arrange(-sum) |>ungroup() |> 
    mutate(rnum = row_number()) |>
    group_by(name) |> slice_head() |> ungroup() |>
    arrange(color_group, rnum) |>
    select(name, color_group)
  
  # Values for g --> x association
  valueGtoX <- c(lapply(lucidus_fit$res_Beta$Beta, 
                        function(x)(x[-1])) |>
                   unlist(), 
                 rep(0, dimG*n_layers))
  
  # For Cluster 2 (which needs effect estimates): 
  valueGtoX_c1 <- do.call(rbind, lucidus_fit$res_Beta$Beta)[,-1] |>
    as_tibble() |>
    dplyr::mutate(layer = str_c("(Layer ", row_number(), ")"),
                  cluster = "Cluster 2") 
  
  # For cluster 1 (ref. cluster, effect est = 0):
  valueGtoX_c2 <- valueGtoX_c1 |>
    mutate(across(where(is.numeric), ~0), 
           cluster = "Cluster 1")
  
  # combine, pivot longer, and create source and target columns
  GtoX <- bind_rows(valueGtoX_c1, valueGtoX_c2) |>
    mutate(target = str_c(cluster, layer, sep = " ")) |>
    pivot_longer(cols = setdiff(colnames(valueGtoX_c1), 
                                c("layer", "cluster")), 
                 names_to = "source", values_to = "value") |>
    mutate(color_group = as.factor(value < 0), 
           value = abs(value)) |>
    select(source, target, value, color_group) |>
    as.data.frame()
  
  
  valueXtoZ <- c(lapply(lucidus_fit$res_Mu_Sigma$Mu, 
                        function(x)x[, 1]) |> 
                   unlist(), 
                 lapply(lucidus_fit$res_Mu_Sigma$Mu, 
                        function(x)x[, 2]) |> 
                   unlist())
  
  valueXtoY <- c(rep(0, n_layers), 
                 # rep(lucidus_fit$res_Delta$Delta$mu[1] / n_layers, n_layers),
                 lucidus_fit$res_Delta$Delta$mu[-1])
  

  # X to Z data (only works for 2 clusters)
  XtoZ <- data.frame(source =  c(rep("Cluster 1 (Layer 1)", n_features[1]),
                                 rep("Cluster 1 (Layer 2)", n_features[2]),
                                 rep("Cluster 2 (Layer 1)", n_features[1]),
                                 rep("Cluster 2 (Layer 2)", n_features[2])), 
                     target = rep(c(lapply(lucidus_fit$res_Mu_Sigma$Mu,
                                           rownames) |> unlist()),
                                  K[1]), 
                     value = abs(valueXtoZ), 
                     color_group = as.factor(valueXtoZ > 0))
  
  
  XtoY <- data.frame(source = rep("Outcome", 2*n_layers), # originally target
                     target = c("Cluster 1 (Layer 1)", # Originally source
                                "Cluster 1 (Layer 2)",
                                "Cluster 2 (Layer 1)",
                                "Cluster 2 (Layer 2)"), 
                     value = abs(valueXtoY), 
                     color_group = as.factor(valueXtoY > 0))
  
  # Create Sankey diagram
  # Links ----
  links <- rbind(GtoX, XtoZ, XtoY) |>
    mutate(
      # Group: one of exposure, clusters, or outcomes 
      # (doesn't include Z.order by desired order)
      source_group = case_when(
        str_detect(source, "Cluster") ~ "2_Cluster", 
        source == "Outcome" ~ "3_outcome", 
        TRUE ~ "1_exposure"), 
      # Source Omics Layer: lc1-lc4 (for omics layers), outcome, or other 
      source_layer = case_when(
        str_detect(source, "Layer 1") ~ "lc1", 
        str_detect(source, "Layer 2") ~ "lc2", 
        source == "Outcome" ~ str_sub(target, start = -3, end = -2), 
        TRUE ~ "exposure"), 
      # Source group_ for color (one of: exposure, : lc1-lc4 (for omics layers), outcome, or other 
      color_group_node = if_else(source == "Outcome", 
                                 "outcome", 
                                 source_layer)) |>
    group_by(source_group) |>
    arrange(source_layer, .by_group = TRUE) |>
    ungroup() |>
    select(source, target, value, color_group, color_group_node)
  
  
  # Nodes
  nodes <- links |>
    select(source, color_group_node) |>
    mutate(rownum = row_number()) |>
    rename(name = source, 
           color_group = color_group_node) |>
    group_by(name) |>
    slice_head() |>
    ungroup() |>
    arrange(rownum) |>
    select(-rownum) |>
    bind_rows(names_features)
  
  
  links$IDsource <- match(links$source, nodes$name) - 1
  links$IDtarget <- match(links$target, nodes$name) - 1
  # Change group to allow changing color of links to outcome
  links <- links |>
    mutate(color_group = case_when(
      source != "Outcome" ~ as.character(color_group), 
      color_group == "TRUE"  ~ "pos_clus_to_out", 
      color_group == "FALSE" ~ "neg_clus_to_out")) 
  
  
  color_scale <- sankey_colors #data.frame(domain = c("exposure", "outcome",
  #                                      "lc1", "lc2", "lc3", "lc4",
  #                                      "layer1", "layer2", "layer3", "layer4",
  #                                      "TRUE", "FALSE", 
  #                                      "pos_clus_to_out", "neg_clus_to_out"), 
  #                           range = sankey_colors)
  
  p <- sankeyNetwork(Links = as.data.frame(links), 
                     Nodes = as.data.frame(nodes),
                     Source = "IDsource", 
                     Target = "IDtarget", 
                     Value = "value", 
                     NodeID = "name", 
                     colourScale = JS(
                       sprintf("d3.scaleOrdinal()\n.domain(%s)\n.range(%s)\n", 
                               jsonlite::toJSON(color_scale$domain), 
                               jsonlite::toJSON(color_scale$range))), 
                     LinkGroup = "color_group", 
                     NodeGroup = "color_group", 
                     sinksRight = FALSE, 
                     fontSize = text_size, 
                     iterations = 0)
  p
}





