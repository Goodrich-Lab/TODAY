# ###Loading Variables
# 
# # outcome_names: all the interested longitudinal outcomes
# # liver_enzymes_outcome_names: all longitudinal liver enzymes outcomes
# # cat_outcome_names: all categorical outcomes
# # outcome_names_1: baseline outcome
# # liver_enzymes_outcome_names_1: baseline liver enzymes outcome
# # cat_outcome_names: baseline categorical outcomes
# # plasma_pfas_names: all longitudinal plasma pfas name
# # liver_pfas_names: all longitudinal liver pfas name
# # all_pfas_names: all pfas name
# 
# outcome_names<- data %>% 
#   dplyr::select(matches("ckd|nafld_type|nash|fibrostg|fibrospe|steatgrd|lob|alt|ast|ggt|alp")) %>% 
#   colnames()
# 
# liver_enzymes_outcome_names <- outcome_names[grepl("alt|ast|ggt|alp",outcome_names)]
# 
# #the _0$ indicates the end ... $ means the end
# cat_outcome_names<- outcome_names[!outcome_names %in% liver_enzymes_outcome_names]
# ## notes: ckd, lob, nash and steatosis only have baseline outcomes.
# outcome_names_1 <- outcome_names[grepl("_0$",outcome_names)]
# 
# liver_enzymes_outcome_names_1 <- liver_enzymes_outcome_names[grepl("_0$",liver_enzymes_outcome_names)]
# 
# cat_outcome_names_1_raw <- cat_outcome_names[grepl("_0$",cat_outcome_names)]
# 
# di_outcomes_1<- c("ckd_1_di", "nafld_1_di", "fibrostg_1_di","lob_1_di")
# 
# di_outcomes_nash_nafld_1 <- c("ckd_1_di", "nafld_nash", "fibrostg_1_di","lob_1_di")
# 
# mul_outcomes_1<- c("nafld_nash_1_mul","nash_1_mul", "steatgrd_1_mul")
# 
# cat_outcome_names_1 <- c(di_outcomes_1, mul_outcomes_1)
# 
# #this is fine as is --> this will select all plasma and liver 
# plasma_pfas_names<-
#   data %>% dplyr::select(contains("_plasma"),
#                          -contains("run_order"), -contains("batch_number"), -contains("niddk")) %>% 
#   colnames()
# liver_pfas_names<- 
#   data %>% dplyr::select(matches("_liver$"),
#                          -contains("identifier_targeted_liver"), -contains("weightofbiopsy_targeted_liver")) %>% 
#   colnames()
# 
# #vector of all pfas names (column names of the pfas)
# all_pfas_names <- c(plasma_pfas_names,liver_pfas_names)
# 
# all_pfas_names_1 <- all_pfas_names[grepl("_0$|_targeted_plasma|liver",all_pfas_names)]
# 
# #age is created in 0_1_restrict to variables in analysis.R (it is year not month)
# covars_names_1 <- c("age", "sex", "race_binary", "parents_income_0", "site_bi")
