## this script is to recode exposures to dichotomous variables like low(0) and high(1) 
## and then redo the g prior analysis

source(here::here("4_1_g_prior_sel_logistic.R"))

data1 <- data %>% 
  dplyr::mutate(
    row = row_number()
  )
data1_l <- data1 %>% 
  labelled::remove_labels() %>%
  pivot_longer(
  names_to = "pfas",
  values_to = "value",
  cols = all_of(all_pfas_names_1)
)

data1_w <- data1_l %>% 
  group_by(pfas) %>%
  dplyr::mutate(
    m = median(value, na.rm = TRUE),
    new_value = ifelse(value <= m, "low","high")
  ) %>% 
  dplyr::select(-m, -value) %>% 
  pivot_wider(
    names_from = pfas,
    values_from = new_value
  )


test_plasma <- data1_w %>% drop_na(nafld_type_1,
                                               all_of(all_pfas_names_1[grepl("plasma",all_pfas_names_1)])) %>%
  mutate(nafld = ifelse(nafld_type_1 == "No NAFLD", 0, 1))
test_liver <- data1_w %>% drop_na(nafld_type_1,
                                              all_of(all_pfas_names_1[grepl("liver",all_pfas_names_1)])) %>%
  mutate(nafld = ifelse(nafld_type_1 == "No NAFLD", 0, 1))

X_plasma <- test_plasma %>% dplyr::select(all_of(all_pfas_names_1[grepl("plasma",all_pfas_names_1)])) %>%
  sapply(function(x) ifelse(is.na(x),NaN,x)) %>%
  fastDummies::dummy_cols(remove_selected_columns = TRUE,
                          remove_first_dummy = TRUE)

X_liver <- test_liver %>% dplyr::select(all_of(all_pfas_names_1[grepl("liver",all_pfas_names_1)])) %>%
  sapply(function(x) ifelse(is.na(x),NaN,x)) %>%
  fastDummies::dummy_cols(remove_selected_columns = TRUE,
                          remove_first_dummy = TRUE)


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

write_xlsx(result_plasma, fs::path(dir_report, "mixture_model_plasma_pfas_dichotomous.xlsx"))
write_xlsx(result_liver, fs::path(dir_report, "mixture_model_liver_pfas_dichotomous.xlsx"))
