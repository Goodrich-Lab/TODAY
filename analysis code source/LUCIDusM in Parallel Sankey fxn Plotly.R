# lucidus_fit <- fit
# n_z_ftrs_to_plot <- c(7, 7, 7)

plot_lucid_in_parallel_plotly<- function(lucidus_fit,
                                         sankey_colors,
                                         text_size = 10, 
                                         n_z_ftrs_to_plot = NULL){
  # Get number of clusters, layers, etc.
  K <- lucidus_fit$K
  dimG <- lucidus_fit$res_Beta$Beta[[1]] |> ncol()-1
  n_layers   <- length(lucidus_fit$res_Beta$Beta)
  
  # Get top omics features based on effect size
  if(!is.null(n_z_ftrs_to_plot)){
    top_ftrs <- vector("list", n_layers)
    for(i in seq_along(top_ftrs)){
      top_ftrs[[i]] <- lucidus_fit$res_Mu_Sigma$Mu[[i]] %>%
        as_tibble(rownames = "name") %>%
        # rownames_to_column("name") %>%
        mutate(effect_size = abs(V1) + abs(V2)) %>%
        arrange(desc(effect_size))
      top_ftr_nms <- top_ftrs[[i]]$name[1:n_z_ftrs_to_plot[i]]
      lucidus_fit$res_Mu_Sigma$Mu[[i]] <- 
        lucidus_fit$res_Mu_Sigma$Mu[[i]][
          rownames(lucidus_fit$res_Mu_Sigma$Mu[[i]]) %in% top_ftr_nms, ]
    }
  }
  
  mu_lst <- purrr::map(lucidus_fit$res_Mu_Sigma$Mu, 
                       ~as_tibble(.x, rownames = "name"))
  names(mu_lst) <- paste0("layer", c(1:n_layers))
  dimZ <- purrr::map(mu_lst, ncol) |> as.numeric()-1
  n_features <- purrr::map(mu_lst, nrow) |> as.numeric()
  names(n_features) <- paste0("layer", c(1:n_layers))
  # Names of features and set order of omics features
  names_features <- bind_rows(mu_lst, .id = "color_group") |> 
    rowwise() |>
    mutate(sum = sum(abs(V1)+abs(V2)), 
           pos_c2 = if_else(V2>0, "pos", "neg")) |>
    group_by(color_group, pos_c2) |> arrange(-sum, .by_group = TRUE) |> ungroup() |> 
    mutate(rnum = row_number()) |>
    group_by(name) |> slice_head() |> ungroup() |>
    arrange(color_group, rnum) |>
    dplyr::select(name, color_group)
  
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
    mutate(color_group = as.factor(value > 0), 
           value = abs(value)) |>
    dplyr::select(source, target, value, color_group) |>
    as.data.frame()
  
  valueXtoZ <- c(lapply(lucidus_fit$res_Mu_Sigma$Mu, 
                        function(x)x[, 1]) |> 
                   unlist(), 
                 lapply(lucidus_fit$res_Mu_Sigma$Mu, 
                        function(x)x[, 2]) |> 
                   unlist())
  
  valueXtoY <- c(rep(0, n_layers), 
                 # rep(lucidus_fit$res_Delta$Delta$mu[1] / n_layers, n_layers),
                 lucidus_fit$res_Delta$fit$coefficients[-1])
  
  # n features in each layer
  XtoZ <- data.frame(source = c(rep("Cluster 1 (Layer 1)", n_features[1]),
                                rep("Cluster 1 (Layer 2)", n_features[2]),
                                # rep("Cluster 1 (Layer 3)", n_features[3]),
                                # rep("Cluster 1 (Layer 4)", n_features[4]),
                                rep("Cluster 2 (Layer 1)", n_features[1]),
                                rep("Cluster 2 (Layer 2)", n_features[2])
                                # rep("Cluster 2 (Layer 3)", n_features[3]), 
                                # rep("Cluster 2 (Layer 4)", n_features[4])
  ), 
  target = rep(c(lapply(lucidus_fit$res_Mu_Sigma$Mu,
                        rownames) |> unlist()),
               K[1]), 
  value = abs(valueXtoZ), 
  color_group = as.factor(valueXtoZ > 0))
  
  # To change the outcome from left to right hand side, flip source and target
  XtoY <- data.frame(target = rep("Outcome", 2*n_layers), 
                     source = c("Cluster 1 (Layer 1)", 
                                "Cluster 1 (Layer 2)",
                                # "Cluster 1 (Layer 3)",
                                # "Cluster 1 (Layer 4)",
                                "Cluster 2 (Layer 1)",
                                "Cluster 2 (Layer 2)"
                                # "Cluster 2 (Layer 3)" ,
                                # "Cluster 2 (Layer 4)"
                     ), 
                     value = abs(valueXtoY), 
                     color_group = as.factor(valueXtoY > 0))
  
  # create Sankey diagram
  # Create Links ----
  links <- rbind(GtoX, XtoZ, XtoY) |>
    mutate(
      # Group: one of exposure, clusters, or outcomes 
      # (doesn't include Z.order by desired order)
      source_group = case_when(
        str_detect(source, "Cluster") ~ "2_Cluster", 
        # source == "Outcome" ~ "3_outcome", # removed when moving outcome to right
        TRUE ~ "1_exposure"), 
      # Source Omics Layer: lc1-lc4 (for omics layers), outcome, or other 
      source_layer = case_when(
        str_detect(source, "Layer 1") ~ "lc1", 
        str_detect(source, "Layer 2") ~ "lc2", 
        str_detect(source, "Layer 3") ~ "lc3", 
        str_detect(source, "Layer 4") ~ "lc4",
        # source == "Outcome" ~ str_sub(target, start = -3, end = -2),  # removed when moving outcome to right
        TRUE ~ "exposure"), 
      # Source group_ for color (one of: exposure, : lc1-lc4 (for omics layers), outcome, or other 
      color_group_node = if_else(source == "Outcome", 
                                 "Outcome", 
                                 source_layer)) |>
    group_by(source_group) |>
    arrange(source_layer, .by_group = TRUE) |>
    ungroup() |>
    dplyr::select(source, target, value, color_group, color_group_node)
  
  
  # Create Nodes ----
  nodes <- links |>
    dplyr::select(source, color_group_node) |>
    mutate(rownum = row_number()) |>
    rename(name = source, 
           color_group = color_group_node) |>
    # Add outcome (only if outcome is on right side)
    bind_rows(data.frame(name = "Outcome", color_group = "Outcome")) |>
    group_by(name) |>
    slice_head() |>
    ungroup() |>
    arrange(rownum) |>
    dplyr::select(-rownum) |>
    # Add feature names
    bind_rows(names_features) |> 
    mutate(id = row_number()-1) |>
    tidylog::left_join(sankey_colors, by = c( "color_group"= "domain"))
  
  # Join links and nodes for color names -----
  links <- links %>%
    tidylog::left_join(nodes %>% 
                         dplyr::select(id, name), 
                       by = c("source" = "name")) %>%
    rename(source_id = id) %>% 
    dplyr::select(source_id, everything()) %>%
    tidylog::left_join(nodes %>% 
                         dplyr::select(id, name), 
                       by = c("target" = "name")) %>%
    rename(target_id = id) %>% 
    dplyr::select(source_id,source, target_id,everything()) 
  
    
  # Manually change colors ----
  links <- links  |>
    mutate(
      link_color = case_when(
        # Ref link color
        value == 0 ~     "#f3f6f4", 
        # Outcome
        str_detect(target, "Outcome") &  color_group == TRUE  ~  "red",
        str_detect(target, "Outcome") &  color_group == FALSE  ~  "#d3d3e3",
        # Methylation 
        str_detect(source, "Layer 2") &  color_group == TRUE  ~  "#bf9000",
        str_detect(source, "Layer 2") &  color_group == FALSE ~  "#ffd966",
        # Transcriptome
        str_detect(source, "Layer 1") &  color_group == TRUE  ~  "#38761d",
        str_detect(source, "Layer 1") &  color_group == FALSE ~  "#b6d7a8",
        # # mirna
        # str_detect(source, "Layer 3") &  color_group == TRUE  ~  "#a64d79",
        # str_detect(source, "Layer 3") &  color_group == FALSE ~  "#ead1dc",
        # # mirna
        # str_detect(source, "Layer 4") &  color_group == TRUE  ~  "#fee2ef",
        # str_detect(source, "Layer 4") &  color_group == FALSE ~  "#d90166",
        
        links$color_group == FALSE ~ "#d3d3e3", # Negative association
        links$color_group == TRUE ~  "red"))
  
  ## change node names
  nodes <- nodes %>%
    mutate(name = case_when(name == "value" ~ "<b>Met</b>",
                            name == "Cluster 1 (Layer 1)" ~ "<b>Carboxylic acids\nProfile 0</b>",
                            name == "Cluster 2 (Layer 1)" ~ "<b>Carboxylic acids\nProfile 1</b>",
                            name == "Cluster 1 (Layer 2)" ~ "<b>Sulfonic Acids\nProfile 0</b>",
                            name == "Cluster 2 (Layer 2)" ~ "<b>Sulfonic Acids\nProfile 1</b>",
                            # name == "Cluster 1 (Layer 3)" ~ "<b>Lipoxygenase\nProfile 0</b>",
                            # name == "Cluster 2 (Layer 3)" ~ "<b>Lipoxygenase\nProfile 1</b>",
                            # name == "Cluster 1 (Layer 4)" ~ "<b>Parent Compound\nProfile 0</b>",
                            # name == "Cluster 2 (Layer 4)" ~ "<b>Parent Compound\nProfile 1</b>",
                            TRUE ~ name),
           x = case_when(
             name == "Met" ~ 0,
             str_detect(name, "acids") ~ 1/3, 
             str_detect(color_group, "layer")~  2/3, 
             str_detect(name, "Outcome") ~ 2.5/3))
  
  (fig <- plot_ly(
    type = "sankey",
    orientation = "h",
    arrangement = "snap",
    node = list(
      label = nodes$name,
      color = nodes$range,
      pad = 15,
      thickness = 20,
      line = list(
        color = "black",
        width = 0.5
      ),
      x = nodes$x
    ),
    
    link = list(
      source = links$source_id,
      target = links$target_id,
      value =  links$value+.00000000000000000000001,
      label = links$source,
      color = links$link_color
    )
  )
  )
  
  fig <- fig %>% layout(
    font = list(
      size = text_size
    ),
    xaxis = list(showgrid = F, zeroline = F),
    yaxis = list(showgrid = F, zeroline = F)
  )
  
  fig
}

