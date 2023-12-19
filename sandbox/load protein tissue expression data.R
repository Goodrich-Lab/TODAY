# Load 
source(here::here("!directories.R"))

# data from https://doi.org/10.1074/mcp.M113.035600 
prot_exp_by_tissue <- readxl::read_excel(
  fs::path(dir_data, 
           "protein_expression_by_organ", 
           "Protein Atlas Fagerberg et al",
           "Fagerberg et al MCP 2014.xlsx")) |> 
  janitor::clean_names()
colnames(prot_exp_by_tissue)



library("EnsDb.Hsapiens.v86")
edb <- EnsDb.Hsapiens.v86


gene_names_symbol <- mapIds(edb, 
       keys = prot_exp_by_tissue$ensembl_gene_id,
       column = "SYMBOL", 
       keytype = "GENEID", multiVals = "list")

gene_id_symbol <- mapIds(edb, 
                            keys = prot_exp_by_tissue$ensembl_gene_id,
                            column = "ENTREZID", 
                            keytype = "GENEID", multiVals = "list")




out <- lapply(gene_id_symbol, FUN = function(x){length(x)}) |> unlist()
table(out)

outdf <- data.frame(out) |> rownames_to_column("egid")



# Single Cell Data from the Human Protein Atlas 

dat <- read_tsv(fs::path(dir_data, 
                  "protein_expression_by_organ", 
                  "rna_single_cell_type_tissue.tsv")) 


dat <- read_rds(fs::path(dir_data, 
                   "protein_expression_by_organ", 
                   "rna_single_cell_type_tissue.rds")) |>
  janitor::clean_names()


head(dat)
          
mccd1 <- dat |> 
  filter(gene_name == "MCCD1")

