## g prior analysis

source(here::here("4_1_g_prior_sel_logistic.R"))

data_imputed_recode <- readRDS(fs::path(dir_merged_data, "tl_imputed_recode_without_log.RDS"))

test_plasma <- data_imputed_recode %>% drop_na(nafld_1_di,
                              all_of(pfas_names_imputed_1[grepl("plasma",pfas_names_imputed_1)])) %>%
  mutate(nafld = ifelse(nafld_1_di == "NAFLD", 1, 0))
test_liver <- data_imputed_recode %>% drop_na(nafld_1_di,
                              all_of(pfas_names_imputed_1[grepl("liver",pfas_names_imputed_1)])) %>%
  mutate(nafld = ifelse(nafld_1_di == "NAFLD", 1, 0))

X_plasma <- test_plasma %>% dplyr::select(all_of(pfas_names_imputed_1[grepl("plasma",pfas_names_imputed_1)])) %>%
  scale()
X_liver <- test_liver %>% dplyr::select(all_of(pfas_names_imputed_1[grepl("liver",pfas_names_imputed_1)])) %>%
  scale()

Y_plasma <- test_plasma$nafld

Y_liver <- test_liver$nafld

U_plasma <- test_plasma %>% dplyr::select(all_of(covars_names_1)) %>%
  sapply(function(x) ifelse(is.na(x),NaN,x)) %>%
  fastDummies::dummy_cols(remove_selected_columns = TRUE,
                          remove_first_dummy = TRUE)

U_liver <- test_plasma %>% dplyr::select(all_of(covars_names_1)) %>%
  sapply(function(x) ifelse(is.na(x),NaN,x)) %>%
  fastDummies::dummy_cols(remove_selected_columns = TRUE,
                          remove_first_dummy = TRUE)

result_plasma <- g.prior.sel.logistic(X_plasma,Y_plasma,U_plasma)
result_liver <- g.prior.sel.logistic(X_liver,Y_liver,U_liver)

write_xlsx(result_plasma, fs::path(dir_report, "mixture_model_plasma_pfas.xlsx"))
write_xlsx(result_liver, fs::path(dir_report, "mixture_model_liver_pfas.xlsx"))


