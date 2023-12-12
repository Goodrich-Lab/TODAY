# ! Load clean data 

full_data <- read_rds(fs::path(dir_data, "clean_data", "TODAY_complete_data_0417.RDS"))

prot_metadata <- read_csv(fs::path(dir_data, "raw_data", "protein_analytes.csv"))

# Set var names
pfas_names <- full_data |> tidylog::select(contains("pfas_")) |> colnames()
met_names <- full_data |> tidylog::select(contains("met_")) |> colnames()
prot_names <- full_data |> tidylog::select(contains("seq")) |> colnames()




# sankey colors
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
                          "Outcome",  col_pal[8],
                          "TRUE",     "#6372e0", # Blue
                          "FALSE",    "#d1d4ff", # Light grey
                          "pos_clus_to_out", "red", 
                          "neg_clus_to_out", "#e4e5f2"), 
                        byrow = TRUE, nrow = 14) |> 
  as_tibble(.name_repair = "universal") |>
  janitor::clean_names() |>
  rename("domain" = x1, 
         "range" = x2)
