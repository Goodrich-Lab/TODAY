# ! Load clean data 
full_data <- read_rds(fs::path(dir_data, "clean_data", "TODAY_complete_data_0417.RDS"))

# Set var names
pfas_names <- full_data |> tidylog::select(contains("pfas_")) |> colnames()
met_names <- full_data |> tidylog::select(contains("met_")) |> colnames()
prot_names <- full_data |> tidylog::select(contains("seq")) |> colnames()


# Get informative metabolites/protein names
omic_names <- c(prot_names, met_names)
var_uni_value <- full_data %>% 
  pivot_longer(cols = all_of(omic_names)) %>%
  group_by(name) %>%
  summarise(n = length(unique(value))) %>%
  filter(n == 1)

omic_names <- full_data %>% 
  dplyr::select(all_of(omic_names), 
                -all_of(var_uni_value$name)) |> 
  colnames()


# Log transform and scale Metabolites and Proteins
full_data <- full_data |>
  mutate_at(.vars = vars(contains("seq.")),
            .funs = ~log(.) %>% scale(.) %>% as.numeric(.)) %>% # jg added log transform
  mutate_at(.vars = vars(contains("met_")),
            .funs = ~log(.) %>% scale(.) %>% as.numeric(.)) # jg added log transform


# Proteomics metadata -------------------------------
# proteomics metadata specific to this project
prot_metadata1 <- read_csv(fs::path(dir_data, "raw_data", "protein_analytes.csv")) |>
  dplyr::select(-c(Dilution:Dilution2)) 

prot_names_reduced <- prot_metadata1 |> 
  tidylog::filter(Organism == "Human")
prot_names_human <- prot_names_reduced$AptName

# read in gene expression by tissue information from Fagerberg et al
gene_exp_by_tissue <- read_csv(
  fs::path(dir_data,
           "protein_expression_by_organ",
           "Protein Atlas Fagerberg et al",
           "Fagerberg et al MCP 2014_cleaned.csv"))


# # Determine which tissue has highest expression of the genes --------
# Define a function to extract the highest and second highest values and their column names
extract_values <- function(row) {
  sorted_indices <- order(row, decreasing = TRUE)
  highest_value <- row[sorted_indices[1]]
  second_highest_value <- row[sorted_indices[2]]
  highest_column <- names(gene_exp_by_tissue)[c(4:30)][sorted_indices[1]]
  second_highest_column <- names(gene_exp_by_tissue)[c(4:30)][sorted_indices[2]]
  return(c(highest_column, highest_value, second_highest_column, second_highest_value))
}

# Apply the function to each row
extracted_values <- t(apply(gene_exp_by_tissue[,c(4:30)], 1, extract_values))

# Add the results to the dataframe
gene_exp_by_tissue$highest_expression <- extracted_values[, 1]
gene_exp_by_tissue$highest_expression_value <- as.numeric(extracted_values[, 2])
gene_exp_by_tissue$second_highest_expression <- extracted_values[, 3]
gene_exp_by_tissue$second_highest_expression_value <- as.numeric(extracted_values[, 4])
gene_exp_by_tissue <- gene_exp_by_tissue |> 
  mutate(high_second_high_ratio = highest_expression_value/second_highest_expression_value)



# Merge metadata and gene expression by tissue information
prot_metadata2 <- prot_metadata1 |> 
  filter(!is.na(EntrezGeneID)) |>
  tidylog::left_join(gene_exp_by_tissue, by = c("EntrezGeneID", "EntrezGeneSymbol")) #|>
  # dplyr::select(-c(colon:salivary_gland))

# Fix issues with missing Entrez gene name
prot_metadata2 <- prot_metadata2 |> 
  mutate(EntrezGeneSymbol = case_when(
    EntrezGeneID == "648791" ~ "PPP1R3G", # Missing Entrez gene name
    AptName == "seq.21813.171" ~ "MGAT5_1", # Same gene, different Apt Name 
    AptName == "seq.21768.9"   ~ "MGAT5_2", # Same gene, different Apt Name 
    AptName == "seq.16768.3" ~ "JUP_1", # Same gene, different Apt Name 
    AptName == "seq.23007.8" ~ "JUP_2", # Same gene, different Apt Name 
    EntrezGeneID == "647087"  ~ "STMP1",
    TRUE ~ EntrezGeneSymbol))

# Calculate ratio of expression in kidney versus other tissues. 
# Select key variables
prot_tissue_expression <- prot_metadata2 %>%
  dplyr::select(AptName, colon:salivary_gland) %>% 
  group_by(AptName) %>% 
  slice_head() %>%
  ungroup()

# Calculate sum of expression across tissues
prot_tissue_expression$sum_expression <-  prot_tissue_expression |> 
  dplyr::select(colon:salivary_gland) %>% 
  rowSums(na.rm = TRUE)

# Calculate the ratios
prot_tissue_expression_ratio <- prot_tissue_expression %>%
  mutate(across(where(is.numeric), ~(.)/sum_expression, 
                .names = "ratio_{.col}")) |> 
  dplyr::select(-c(colon:salivary_gland), -contains("sum_expression"))

# Add the ratios to the metadata
prot_metadata <- prot_metadata2 |> 
  tidylog::full_join(prot_tissue_expression_ratio) |>
  tidylog::filter(Organism == "Human") |> 
  dplyr::select(-c(colon:salivary_gland))

