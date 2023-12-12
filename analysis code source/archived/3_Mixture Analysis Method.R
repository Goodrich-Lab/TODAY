# This script include mixture analysis using qgcomp()

# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Example-----
## Step1. load data----
data <- read_rds(fs::path(dir_data, "example_teenlabs_data.rds"))

## Step2. Define variables----
# covariates
covars <- c("race_binary", 
            "smoke_0", 
            "age_0", 
            "sex", 
            "parents_income_0")

# plasma pfas
plasma_pfas <- colnames(data)[grepl("imputed",colnames(data))]
## or
# plasma_pfas <- c("pfda_untargeted_plasma_0_imputed",
#                  "pf_hx_s_untargeted_plasma_0_imputed",
#                  "pf_hp_s_untargeted_plasma_0_imputed",
#                  "pfna_untargeted_plasma_0_imputed",
#                  "pfoa_untargeted_plasma_0_imputed",
#                  "pfos_untargeted_plasma_0_imputed")

## Step3. Mixture analysis using qgcomp()

# qgcomp-----
qgcomp_output <- qgcomp.noboot(as.formula(paste("nafld_di_0~",
                                                str_c(plasma_pfas, collapse = "+") ,
                                                "+",
                                                str_c(covars, collapse = "+"))), # Covariates - same as MLR model
                               expnms=plasma_pfas, # Name of Mixtures
                               data=data,
                               family=binomial(), # family=binomial() for logistic
                               q = 3)

## Step4. View Models --------------

## weighted plot
png(fs::path(dir_figure,
             "3_weighted_plot_gcomp_nafld.png"),
    width = 4000, height = 2000, res = 350)
plot(qgcomp_output)
dev.off()


## coefficient plot
coef <- data.frame(qgcomp_output$fit$coefficients) %>% 
  rownames_to_column("Predictor") %>%
  filter(grepl("targeted_plasma", Predictor)) %>%
  rename('Log Odds Ratio' = qgcomp_output.fit.coefficients) %>%
  mutate(direction = ifelse(`Log Odds Ratio` > 0, "pos", "neg"),
         direction = factor(direction, levels = c("pos", "neg")))

# change the name of the predictor
coef1 <- coef %>% mutate(Predictor = case_when(
  grepl("pfda", Predictor) ~ "PFDA",
  grepl("pfos", Predictor) ~ "PFOS",
  grepl("pfna", Predictor) ~ "PFNA",
  grepl("pf_hx_s", Predictor) ~"PFHxS",
  grepl("pfoa", Predictor) ~ "PFOA",
  grepl("pf_hp_s", Predictor) ~ "PFHpS"
))

p <- ggplot(coef1,aes(Predictor, `Log Odds Ratio`,
                      fill = direction)) + 
  geom_bar(stat = "identity", show.legend = FALSE) + 
  coord_flip() +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  geom_text(aes(label = round(`Log Odds Ratio`,2)), hjust= -0.05) +
  theme(axis.text = element_text(size = 15),
        text = element_text(size = 15),
        axis.title.y = element_blank())


ggsave(p, 
       filename = fs::path(dir_figure,
                           "3_coef_plot_qgcomp.jpeg"),
       width = 7,
       height = 4)

