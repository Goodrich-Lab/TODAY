rename_pfas <- function(df){
  df <- df %>%
    mutate(pfas = case_when(grepl("pfos", pfas) ~ "PFOS",
                            grepl("pfoa", pfas) ~ "PFOA",
                            grepl("pfda", pfas) ~ "PFDA",
                            grepl("pfhpa", pfas) ~ "PFHpA",
                            grepl("pfhps", pfas) ~ "PFHpS",
                            grepl("pfhxs", pfas) ~ "PFHxS",
                            grepl("pfna", pfas) ~ "PFNA",
                            grepl("pfuna", pfas) ~ "PFUnA",
                            grepl("nmefosaa", pfas) ~ "NMeFOSAA",
                            grepl("score_pfsas", pfas) ~ "PFSA Burden Score",
                            grepl("score_pfcas", pfas) ~ "PFCA Burden Score",
                            grepl("score", pfas) ~ "PFAS Burden Score",
                            TRUE ~ pfas))
  return(df)
}

rename_outcomes <- function(x){
  nms <- case_when(
    x=="time_to_glyc_scld" ~ "Time to glycemic failure", 
    x=="codi" ~ "β-cell function (C-Peptide oDI)",
    x=="si_1_ins0" ~  "Ins. Sensitivity (1/fast. ins)",
    x=="hb_a1c" ~ "Hyperglycemia (HbA1c)", 
    x=="daystomic" ~ "Time to microalbuminuria", 
    x== "daystorapid" ~ "Time to rapid",
    x== "daystohyp" ~ "Time to hyperfiltration",
    x=="est_creat_clear" ~ "Estimated Creatinine Clearance",
    x=="daystomac" ~ "Time to macroalbuminuria",
    x=="uacid" ~ "Uric Acid",
    x=="u_alb_creat" ~ "u_alb_creat",
    x=="serum_creat" ~ "Serum Creatanine",
    x== "seq.2836.68" ~"NGAL", 
    x=="seq.9021.1" ~ "KIM-1", 
    x == "seq.5661.15" ~ "IL18",
    x == "seq.11516.7" ~ "FABPL",
    x == "seq.17138.8" ~ "a-GST", 
    x == "seq.12446.49" ~ "a-GST1",
    x == "seq.15509.2" ~ "NAG")
}
