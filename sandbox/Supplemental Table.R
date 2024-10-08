est <- read_csv(fs::path(dir_results, 
                          "meet_in_middle_res_sig_020924.csv")) %>%
  mutate(`PFNA β[95%CI]` = paste(round(estimate_em, 2), "[", 
                               round(conf_low_em, 2), ",",
                               round(conf_high_em, 2),"]"),
         `PFNA P-Value` = round(p_value_em, 2),
         `Albuminuria β[95%CI]` = paste(round(estimate_mo, 2), "[", 
                                          round(conf_low_mo, 2), ",",
                                          round(conf_high_mo, 2),"]"),
         `Albuminuria P-Value` = round(p.value_mo, 2))
  

med <- read_csv(fs::path(dir_results, "med_res_df_020924.csv"))

med_w <- med %>% filter(Effect %in% c("Rtnde", "Rpnie", "Rte")) %>%
  dplyr::select(feature_name, EntrezGeneSymbol,Effect, pe, ci_low, ci_high) %>%
  pivot_wider(
              names_from = "Effect",
              values_from = c(pe, ci_low, ci_high)) %>%
  mutate(`Rtnde[95%CI]` = paste(round(pe_Rtnde, 2), "[", 
                                round(ci_low_Rtnde, 2), ",",
                                round(ci_high_Rtnde, 2),"]"),
         `Rpnie[95%CI]` = paste(round(pe_Rpnie, 2), "[", 
                                round(ci_low_Rpnie, 2), ",",
                                round(ci_high_Rpnie, 2),"]"),
         `Total Effect[95%CI]` = paste(round(pe_Rte, 2), "[", 
                                round(ci_low_Rte, 2), ",",
                                round(ci_high_Rte, 2),"]"))

df <- est %>% 
  dplyr::select(EntrezGeneSymbol, contains("PFNA"), contains("Albuminuria")) %>%
  tidylog::left_join(med_w %>% 
                       dplyr::select(EntrezGeneSymbol, contains("[95%CI]"))) %>%
  rename(`Protein name` = EntrezGeneSymbol)


writexl::write_xlsx(df, fs::path(dir_results, "Supplemental Table 2.xlsx"))                   
