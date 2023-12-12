# This script is the comparison of targeted and untargeted plasma PFAS

#	Part A: Spearman correlation plot
#	Part B: Bland altman plots 
#	Part C: scatter plots 


# Load libraries and directories.
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Load data-----
data <- read_rds(fs::path(dir_project_data, "tl_analysis_ready_data.rds"))

# exposures: plasma pfas----
plasma_pfas <- colnames(data)[grep("0_imputed|_targeted_plasma", colnames(data))]

# subset plasma pfas data-----
data_plasma_pfas <- data %>% dplyr::select(plasma_pfas) %>% 
  rename_at(vars(everything()), 
            ~ sub("_plasma$|_plasma_0_imputed","",.x))

# Part A: Spearman correlation plot------

cor_matrix <- cor(data_plasma_pfas, use = "complete.obs")

# cor_matrix1 <- cor_matrix[str_detect(colnames(cor_matrix), "untargeted"),
#                           str_detect(rownames(cor_matrix), "^pfos_targeted|pfda_targeted|pf_hx_s_targeted|pf_hp_s_targeted|pfna_targeted|pfoa_targeted")]

png(fs::path(dir_figure, "1_2_correlation_plot.png"),
    width = 2000, height = 2000, res = 350)
corrplot(cor_matrix, 
         type="lower", 
         order="original", 
         method = "circle", 
         col=brewer.pal(n=8, name="Blues"), 
         addCoef.col = "black",
         tl.cex=0.5, number.cex=0.3)
dev.off()


#	Part B: Bland altman plots -----

data_plasma_pfas_w <- data %>% dplyr::select(key, plasma_pfas) %>% 
  rename_at(vars(everything()), 
            ~ sub("_plasma$|_plasma_0_imputed","",.x))

data_plasma_pfas_l_targeted <- data_plasma_pfas_w %>% 
  dplyr::select(key, contains("_targeted")) %>% 
  dplyr::select(key, matches("pfda_targeted|pf_hx_s_targeted|pf_hp_s_targeted|pfna_targeted|pfoa_targeted"),
                total_pfos_targeted)%>%
  pivot_longer(cols = contains("_targeted"),
               names_to = "pfas_name_targeted",
               values_to = "targeted")%>%
  separate(col = "pfas_name_targeted",
           into = c("pfas_name",NA),
           sep = "_targeted",
           remove = FALSE
  ) %>%
  mutate(pfas_name = ifelse(pfas_name == "total_pfos", "pfos", pfas_name))

data_plasma_pfas_l_untargeted <- 
  data_plasma_pfas_w %>% 
  dplyr::select(key, contains("_untargeted")) %>% 
  pivot_longer(cols = contains("_untargeted"),
               names_to = "pfas_name_untargeted",
               values_to = "untargeted") %>%
  separate(col = "pfas_name_untargeted",
           into = c("pfas_name",NA),
           sep = "_untargeted",
           remove = FALSE
  )

data_plot <- data_plasma_pfas_l_targeted %>% 
  tidylog::left_join(data_plasma_pfas_l_untargeted, 
                     by = c("key","pfas_name")) %>%
  dplyr::select(key, "pfas_name", "targeted", "untargeted") %>%
  rowwise() %>%
  mutate(avg = mean(c(targeted,untargeted)),
         diff = targeted - untargeted) %>% 
  drop_na() %>% group_by(pfas_name) %>%
  mutate(mean_diff = mean(diff),
         lower = mean_diff - 1.96*sd(diff),
         upper = mean_diff + 1.96*sd(diff)) %>%
  ungroup()

#create Bland-Altman plot
plots <- data_plot %>% group_by(pfas_name) %>%
  do(
      plots = ggplot(.,aes(x = avg, y = diff)) +
      geom_point(size=2) +
      geom_hline(yintercept = .$mean_diff[1]) +
      geom_hline(yintercept = .$lower[1],
                 color = "red", linetype="dashed") +
      geom_hline(yintercept = .$upper[1] ,
                 color = "red", linetype="dashed") +
      ggtitle(str_c("Bland-Altman Plot of ", .$pfas_name)) +
      ylab("Difference") + 
      xlab("Average") +
      theme(text = element_text(size = 6)))

png(fs::path(dir_figure, "1_2_Bland_Altman plot"),
    width = 2000, height = 1500, res = 350)
merged_plot <- grid.arrange(plots$plots[[1]],
                            plots$plots[[2]],
                            plots$plots[[3]],
                            plots$plots[[4]],
                            plots$plots[[5]],
                            plots$plots[[6]],
                            ncol=3,
                            top = textGrob("Bland Altman plots",
                                           gp=gpar(fontsize=15,
                                                   font=3,
                                                   cex=1,
                                                   col="darksalmon"),
                                           vjust=0.4))
dev.off()

#	Part C: scatter plots -----

scatter_plot <- data_plot %>% ggplot(aes(x=untargeted, y=targeted))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  labs(x="Untargeted Plasma", y="Targeted Plasma")+
  geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free") +
  theme_cowplot()

png(fs::path(dir_figure, "1_2_scatter_plot_plasma"),
    width = 3000, height = 2000, res = 350)
scatter_plot
dev.off()


