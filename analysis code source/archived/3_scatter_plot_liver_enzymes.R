
# Load libraries and directories----
source(here::here("!libraries.R"))
source(here::here("!directories.R"))

# Load data-----
data <- read_rds(fs::path(dir_project_data, "tl_analysis_ready_data.rds"))

data_no_nafld_nash <- data %>% 
  filter(nafld_di_0 == "No NAFLD") 

data_nafld_nash <- data %>% 
  filter(nafld_di_0 == "NAFLD")

## exposures: plasma pfas----
plasma_pfas <- colnames(data)[grep("0_imputed|_targeted_plasma", colnames(data))]

##outcomes: liver enzymes
outcomes <- c("alt_0", "ast_0", "ggt_0")

# Full data---
# subset plasma pfas data
data_plasma_pfas <- data %>% dplyr::select(key,plasma_pfas,outcomes) %>% 
  rename_at(vars(everything()),
            ~ sub("_plasma$|_plasma_0_imputed","",.x))
## Long format----
data_plasma_pfas_l <-  data_plasma_pfas %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_scaled <- data_plasma_pfas_w %>%
  mutate_at(.vars = vars(contains("targeted")),
            .funs = scale) %>% 
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_log <- data_plasma_pfas_w %>% 
  mutate_at(.vars = vars(contains("targeted")),
            .funs = log2) %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

## Scatter plot-----
### alt, raw pfas-----
scatter_plot_raw_alt <- data_plasma_pfas_l %>% 
  group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_raw.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_alt
dev.off()

### ast, raw pfas----
scatter_plot_raw_ast <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =   4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_raw.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ast
dev.off()

### ggt, raw pfas-----
scatter_plot_raw_ggt <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =   4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_raw.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ggt
dev.off()


### alt, group mean centered pfas-----
scatter_plot_scaled_alt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_scaled.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_alt

dev.off()


### ast, group mean centered pfas----
scatter_plot_scaled_ast <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_scaled.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_ast

dev.off()

### ggt, group mean centered pfas----
scatter_plot_scaled_ggt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_scaled.png"),
    width = 3000, height = 2000, res = 350)

scatter_plot_scaled_ggt

dev.off()

### alt, log2 pfas-----
scatter_plot_log_alt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_lg2.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_alt
dev.off()

### ast, log2 pfas----
scatter_plot_log_ast <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_lg2.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ast

dev.off()

### ggt, log2 pfas----
scatter_plot_log_ggt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_lg2.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ggt

dev.off()

# NAFLD/NASH data-----
## Long format data-----

## nafld or nash data---- 

data_plasma_pfas_w <- data_nafld_nash %>% dplyr::select(key, plasma_pfas, cont_outcomes) %>% 
  rename_at(vars(everything()),
            ~ sub("_plasma$|_plasma_0_imputed","",.x))

data_plasma_pfas_l <- data_plasma_pfas_w %>% 
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_scaled <- data_plasma_pfas_w %>%
  mutate_at(.vars = vars(contains("targeted")),
            .funs = scale) %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_log <- data_plasma_pfas_w %>% 
  mutate_at(.vars = vars(contains("targeted")),
            .funs = log2) %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

##	scatter plots -----
### alt, raw pfas-----
scatter_plot_raw_alt <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_raw_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_alt

dev.off()

###ast, raw pfas----
scatter_plot_raw_ast <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_raw_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ast

dev.off()

### ggt, raw pfas----
scatter_plot_raw_ggt <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_raw_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ggt

dev.off()

### alt, group mean centered pfas-----
scatter_plot_scaled_alt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_scaled_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_alt

dev.off()

### ast, group mean centered pfas----
scatter_plot_scaled_ast <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_scaled_NAFLD.png"),
    width = 3000, height = 2000, res = 350)

scatter_plot_scaled_ast

dev.off()


### ggt, group mean centered pfas----
scatter_plot_scaled_ggt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_scaled_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_ggt

dev.off()

### alt, log2 pfas-----
scatter_plot_log_alt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_lg2_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_alt

dev.off()

### ast, log2 pfas----
scatter_plot_log_ast <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_lg2_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ast

dev.off()

### ggt, log2 pfas----
scatter_plot_log_ggt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_lg2_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ggt

dev.off()

# Non NAFLD/NASH data-----
## Long format data-----
data_plasma_pfas_w <- data_no_nafld_nash %>% dplyr::select(key, plasma_pfas, cont_outcomes) %>% 
  rename_at(vars(everything()),
            ~ sub("_plasma$|_plasma_0_imputed","",.x))

data_plasma_pfas_l <- data_plasma_pfas_w %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_scaled <- data_plasma_pfas_w %>%
  mutate_at(.vars = vars(contains("targeted")),
            .funs = scale) %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")

data_plasma_pfas_l_log <- data_plasma_pfas_w %>%
  mutate_at(.vars = vars(contains("targeted")),
            .funs = log2) %>%
  pivot_longer(cols = contains("targeted"),
               names_to = "pfas_name",
               values_to = "pfas_value")


##	scatter plots -----

### alt, raw pfas-----
scatter_plot_raw_alt <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_raw_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_alt

dev.off()


###ast, raw pfas----
scatter_plot_raw_ast <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_raw_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ast

dev.off()


### ggt, raw pfas----
scatter_plot_raw_ggt <- data_plasma_pfas_l %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_raw_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_raw_ggt
dev.off()


### alt, group mean centered pfas-----
scatter_plot_scaled_alt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm, 
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_scaled_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_alt

dev.off()

### ast, group mean centered pfas----
scatter_plot_scaled_ast <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_scaled_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_scaled_ast

dev.off()

### ggt, group mean centered pfas----
scatter_plot_scaled_ggt <- data_plasma_pfas_l_scaled %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Group mean centered PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_scaled_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)

scatter_plot_scaled_ggt

dev.off()

### alt, log2 pfas-----
scatter_plot_log_alt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(alt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline ALT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_alt_lg2_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_alt
dev.off()


### ast, log2 pfas----
scatter_plot_log_ast <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ast_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred", 
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline AST")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ast_lg2_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ast

dev.off()

### ggt, log2 pfas----
scatter_plot_log_ggt <- data_plasma_pfas_l_log %>% group_by(pfas_name) %>%
  ggplot(aes(x= pfas_value, y=scale(ggt_0)))+
  geom_point(shape=18, color="blue") +
  geom_smooth(method = lm,
              linetype="dashed",
              color="darkred",
              fill="blue") +
  stat_regline_equation(label.y =  4, aes(label= ..rr.label..))+
  labs(x="Log2 PFAS value", y="Baseline GGT")+
  # geom_abline(slope = 1, intercept = 0, color = "grey")+
  facet_wrap(.~pfas_name, scales = "free_x") +
  theme_cowplot()

png(fs::path(dir_figure, "3_scatter_plot_pfas_ggt_lg2_NO_NAFLD.png"),
    width = 3000, height = 2000, res = 350)
scatter_plot_log_ggt
dev.off()

rm(scatter_plot_log_ggt, scatter_plot_log_ast, scatter_plot_log_alt,
   scatter_plot_scaled_ggt, scatter_plot_scaled_ast, scatter_plot_scaled_alt,
   scatter_plot_raw_ggt, scatter_plot_raw_ast, scatter_plot_raw_alt,
   data_plasma_pfas_l_log, data_plasma_pfas_l_scaled, data_plasma_pfas_l,
   data_plasma_pfas, data, data_no_nafld_nash, data_nafld_nash, plasma_pfas,
   outcomes)
