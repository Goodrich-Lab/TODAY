# 
# 
# ## clean pfas name
# clean_names<- function(df){
#   df <- df %>% 
#     mutate(pfas_names = 
#              case_when(
#                pfas_names=="8_2FTS"   ~ "8:2FTS", 
#                pfas_names=="6_2FTS"   ~ "6:2FTS",
#                pfas_names=="NMeFOSAA" ~ "N-MeFOSAA",
#                TRUE ~ pfas_names))
#   return(df)
# }
# ## rename varibales in the plot.
# rename_variables_plot <- function(df){
#   df <- df %>% mutate(
#     pfas_name = case_when(
#       grepl("8_2FTS", pfas_name, ignore.case = TRUE) ~ "8:2FTS",
#       grepl("pfbs", pfas_name, ignore.case = TRUE) ~ "PFBS",
#       grepl("pfds", pfas_name, ignore.case = TRUE) ~ "PFDS",
#       grepl("pfhpa", pfas_name, ignore.case = TRUE) ~ "PFHpA",
#       grepl("pfpes", pfas_name, ignore.case = TRUE) ~ "PFPeS",
#       grepl("pfuda", pfas_name, ignore.case = TRUE) ~ "PFUdA",
#       grepl("pfda", pfas_name, ignore.case = TRUE) ~ "PFDA",
#       grepl("pfhps", pfas_name, ignore.case = TRUE) ~ "PFHpS",
#       grepl("pfhxs", pfas_name, ignore.case = TRUE) ~ "PFHxS",
#       grepl("pfna", pfas_name, ignore.case = TRUE) ~ "PFNA",
#       grepl("pfoa", pfas_name, ignore.case = TRUE) ~ "PFOA",
#       grepl("^pfos", pfas_name, ignore.case = TRUE) ~ "PFOS",
#       grepl("6_2fts", pfas_name, ignore.case = TRUE) ~ "6:2FTS",
#       grepl("brpfos", pfas_name, ignore.case = TRUE) ~ "brPFOS",
#       grepl("nmefosaa", pfas_name, ignore.case = TRUE) ~ "NMeFOSAA",
#       grepl("npfos", pfas_name, ignore.case = TRUE) ~ "nPFOS"
#     ),
#     output_name = case_when(
#       grepl("nafld_1_di", output_name, ignore.case = TRUE) ~ "NAFLD(NAFLD or No NAFLD)",
#       grepl("fibrostg", output_name, ignore.case = TRUE) ~ "Fibrosis(Yes or No)",
#       grepl("lob", output_name, ignore.case = TRUE) ~ "Lobullar Inflammation(Yes or No)",
#       grepl("nafld_nash", output_name, ignore.case = TRUE) ~ "NAFLD(No NAFLD, NAFLD not NASH, NASH)",
#       grepl("steatgrd", output_name, ignore.case = TRUE) ~ "Steatosis(<5%, 5-33%, >33%)",
#       grepl("nash", output_name, ignore.case = TRUE) ~ "NASH(0, 1, 2, 3+)"
#     )
#   ) 
#   return(df)
# }
# 
# 
# 
# ## Recode outcome -- this creates the categorical variables for many of the outcomes
# 
# recode_outcome <- function(df){
#   # new columns added, ckd_0_dichotomous, nash_0_recode,
#   df <- df %>% mutate(
#     # recode ckd to two categories
#     ckd_0_di = 
#       ifelse(ckd_0!="Normal", "Not Normal", "Normal"),
#     # recode nafld_0 to three categories
#     nafld_nash_0_mul =
#       case_when(
#         nafld_type_0 == "No NAFLD" ~ "No NAFLD",
#         nafld_type_0 == "NAFLD not NASH" ~ "NAFLD not NASH",
#         nafld_type_0 == "Borderline NASH" ~ "NASH",
#         nafld_type_0 == "Definite NASH" ~ "NASH"
#       ),
#     # recode nash_0 to four categories
#     nash_0_mul = 
#       case_when(
#         nash_0 == 0 ~ "0",
#         nash_0 == 1 ~ "1",
#         nash_0 == 2 ~ "2",
#         nash_0 >= 3 ~ "3+",
#       ),
#     # recode fibrostg_0 to dichonomous
#     fibrostg_0_di =
#       ifelse(fibrostg_0=="None", "No", "Yes"),
#     steatgrd_0_mul = 
#       case_when(
#         steatgrd_0 =="0%"|steatgrd_0 =="<5%" ~ "<5%",
#         steatgrd_0 =="5-33%" ~ "5-33%",
#         steatgrd_0 =="33-67%"|steatgrd_0==">67%" ~ ">33%"
#       ),
#     lob_0_di =
#       ifelse(lob_0=="None", "No", "Yes"),
#   )
#   ## reorder the factor outcome
#   df$nafld_nash_0_mul <- df$nafld_nash_0_mul %>%
#     factor(levels = c("No NAFLD", "NAFLD not NASH", "NASH"),ordered = TRUE) 
#   # NASH
#   df$nash_0_mul <- df$nash_0_mul %>%
#     factor(levels = c("0", "1", "2", "3+"),ordered = TRUE) 
#   # CKD
#   df$ckd_0_di <- df$ckd_0_di %>%
#     factor(levels = c("Normal", "Not Normal"),ordered = FALSE)
#   
#   #Fibrostg
#   df$fibrostg_0_di <- df$fibrostg_0_di %>%
#     factor(levels = c("No", "Yes"), ordered = FALSE)
#   #steatgrd_0
#   df$steatgrd_0_mul <- df$steatgrd_0_mul%>%
#     factor(levels = c("<5%", "5-33%", ">33%" ), ordered = TRUE)
#   #lob
#   df$lob_0_di <- df$lob_0_di%>%
#     factor(levels = c("No", "Yes"), ordered = FALSE)
#   #adding two category nafld
#   df <- df %>% 
#     mutate(nafld_0_di = ifelse(nafld_nash_0_mul == "No NAFLD", "No NAFLD", "NAFLD"))
#   df$nafld_0_di <- df$nafld_0_di %>%
#     factor(levels = c("No NAFLD", "NAFLD"), ordered = FALSE)
#   df <- df %>% 
#     dplyr::select(ckd_0,ckd_0_di,
#                   nafld_type_0, nafld_nash_0_mul, nafld_0_di,
#                   nash_0, nash_0_mul,
#                   fibrostg_0, fibrostg_0_di,
#                   steatgrd_0, steatgrd_0_mul,
#                   lob_0, lob_0_di,
#                   everything())
#   
#   return (df)
# }
# 
# data_working<- recode_outcome(data_working)
# 
# 
# 
# # analysis function 
# 
# #model_ouput = crude model no covariates includes LR and OLR 
# 
# model_output <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       # covars = str_c(cov_names,collapse = "+"),
#       formula = str_c(y,"~","`", x,"`"))
#   # Run ordinal logistic regression models
#   if(grepl("nafld_nash|nash|steatgrd",models$y[1])){
#     models$output <- map(models$formula,
#                          ~polr(as.formula(.),
#                                data =
#                                  df%>%
#                                  droplevels(),
#                                Hess=TRUE,
#                                na.action = na.exclude
#                                # ,
#                                # method = "loglog"
#                          ) %>%
#                            tidy(., conf.int = TRUE)%>% 
#                            na.omit())
#     
#   } else {
#     models$output <- map(models$formula,
#                          ~glm(as.formula(.),
#                               data =
#                                 df%>%
#                                 droplevels(),
#                               na.action = na.exclude,
#                               family = binomial
#                          ) %>%
#                            tidy(., conf.int = TRUE)%>% 
#                            na.omit())
#   }
#   return(models)
# }
# 
# #adjusted for covariates LR and OLR model 
# model_output_adjusted <- function(eo_combinations, df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","`", x,"`", "+", covars))
#   # Run ordinal logistic regression models
#   if(grepl("nafld_nash|nash|steatgrd",models$y[1])){
#     models$output <- map(models$formula,
#                          ~polr(as.formula(.),
#                                data =
#                                  df%>%
#                                  droplevels(),
#                                Hess=TRUE,
#                                na.action = na.exclude
#                                # ,
#                                # method = "loglog"
#                          ) %>%
#                            tidy(., conf.int = TRUE, na.rm = TRUE)%>% 
#                            na.omit())
#     
#   } else {
#     models$output <- map(models$formula,
#                          ~glm(as.formula(.),
#                               data =
#                                 df%>%
#                                 droplevels(),
#                               na.action = na.exclude,
#                               family = binomial
#                          ) %>%
#                            tidy(., conf.int = TRUE, na.rm = TRUE)%>% 
#                            na.omit())
#   }
#   return(models)
# }
# 
# #adjusted LR and OLR without site model
# model_output_adjusted_without_site <- function(eo_combinations, df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1[-5],collapse = "+"),
#       formula = str_c(y,"~","`", x,"`", "+", covars))
#   # Run ordinal logistic regression models
#   if(grepl("nafld_nash|nash|steatgrd",models$y[1])){
#     models$output <- map(models$formula,
#                          ~polr(as.formula(.),
#                                data =
#                                  df%>%
#                                  droplevels(),
#                                Hess=TRUE,
#                                na.action = na.exclude
#                          ) %>%
#                            tidy(., conf.int = TRUE))
#     
#   } else {
#     models$output <- map(models$formula,
#                          ~glm(as.formula(.),
#                               data =
#                                 df%>%
#                                 droplevels(),
#                               na.action = na.exclude,
#                               family = binomial
#                          ) %>%
#                            tidy(., conf.int = TRUE))
#   }
#   return(models)
# }
# 
# #verify what this does and if we need it .... 
# #Hongxu please verify 
# olr_model_mine <- function(one_pfas_name, one_outcome_name,df){
#   eo_combinations <- list(x = one_pfas_name,
#                           y = one_outcome_name) %>%
#     cross_df()
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","`", x,"`", "+", covars))
#   if(grepl("nafld_nash|nash|steatgrd",models$y[1])){
#     models$output <- polr(as.formula(models$formula),
#                           data =
#                             df%>%
#                             droplevels(),
#                           Hess=TRUE,
#                           na.action = na.exclude
#     ) %>%
#       tidy_mine(.,one_pfas_name, conf.int = TRUE)%>% 
#       na.omit()
#   } else {
#     models$output <- glm(as.formula(models$formula),
#                          data =
#                            df%>%
#                            droplevels(),
#                          na.action = na.exclude,
#                          family = binomial
#     ) %>%
#       tidy_mine(.,one_pfas_name,conf.int = TRUE)%>% 
#       na.omit()
#   }
#   return(models)
# }
# 
# 
# #used by current script "p-trend" -- adjusted logistic regression model 
# model_output_tertile <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","`", x,"`", "+", covars))
#   # Run ordinal logistic regression models
#   models$output <- map(models$formula,
#                        ~glm(as.formula(.),
#                             data = df,
#                             na.action = na.exclude,
#                             family = binomial
#                        ) %>%
#                          tidy(., conf.int = TRUE)%>% 
#                          na.omit())
#   return(models)
# }
# 
# #adjusted gam model that is log2 transformed exposure 
# model_output_gam_log <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","s(","log2(","`", x,"`" , ")",")", "+", covars))
#   # Run GAM models
#   models$output <- map(models$formula,
#                        ~mgcv::gam(as.formula(.),
#                                   data = df,
#                                   na.action = na.exclude,
#                                   family = binomial
#                        )%>%
#                          tidy(., conf.int = FALSE))
#   
#   models$mod <- map(models$formula,
#                     ~mgcv::gam(as.formula(.),
#                                data = df,
#                                na.action = na.exclude,
#                                family = binomial
#                     ))
#   return(models)
# }
# 
# #GAM model with covariates no log transformation includes spline (logistic)
# model_output_gam <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","s(",  "`", x,"`",")", "+", covars))
#   # Run GAM models
#   models$output <- map(models$formula,
#                        ~mgcv::gam(as.formula(.),
#                                   data = df,
#                                   na.action = na.exclude,
#                                   family = binomial
#                        )%>%
#                          tidy(., conf.int = FALSE))
#   
#   models$mod <- map(models$formula,
#                     ~mgcv::gam(as.formula(.),
#                                data = df,
#                                na.action = na.exclude,
#                                family = binomial
#                     ))
#   return(models)
# }
# 
# #GAM model with covariates no log -- no spline (logistic)
# model_output_gam_linear <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~", "`", x,"`", "+", covars))
#   # Run GAM models
#   models$output <- map(models$formula,
#                        ~glm(as.formula(.),
#                             data = df,
#                             na.action = na.exclude,
#                             family = binomial
#                        )%>%
#                          tidy(., conf.int = TRUE))
#   
#   models$mod <- map(models$formula,
#                     ~mgcv::gam(as.formula(.),
#                                data = df,
#                                na.action = na.exclude,
#                                family = binomial
#                     ))
#   return(models)
# }
# 
# #log transformed GAM model (adjusted) no spline
# model_output_gam_linear_log <- function(eo_combinations,df){
#   # Create regression model formula
#   models <- eo_combinations %>%
#     mutate(
#       covars = str_c(covars_names_1,collapse = "+"),
#       formula = str_c(y,"~","log2(", "`", x,"`", ")", "+", covars))
#   # Run GAM models
#   models$output <- map(models$formula,
#                        ~glm(as.formula(.),
#                             data = df,
#                             na.action = na.exclude,
#                             family = binomial
#                        )%>%
#                          tidy(., conf.int = TRUE))
#   
#   models$mod <- map(models$formula,
#                     ~mgcv::gam(as.formula(.),
#                                data = df,
#                                na.action = na.exclude,
#                                family = binomial
#                     ))
#   return(models)
# }
# 
