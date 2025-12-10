# ---- Load dependencies and script arguments ----
library(rstatix)
library(pROC)
library(yaml)
library(openxlsx)
library(ggplot2)
library(here)
source(here("src", "R", "pipeline_utils.R"))

model_colours <- c("SG" = "black",
                   "ICS" ="#79706E",
                   "proteomics" = "#4E79A7",
                   "ISO" = "#59A14F",
                   "freq" ="#B07AA1", 
                   "PIIL246" ="#F28E2B",
                   "US" ="#E15759", 
                   "LPS" ="#9D7660"
)
currentDate <- Sys.Date() 

args <- commandArgs(trailingOnly = TRUE)


cfg  <- yaml::yaml.load_file(args[1])
cfg$data$save_dir          <- here(cfg$data$save_dir)
cfg$data$result_dir        <- here(cfg$data$result_dir)
cfg$data$processed_data    <- here(cfg$data$processed_data)
cfg$data$sample_metadata   <- here(cfg$data$sample_metadata)


dir.create(cfg$data$result_dir,recursive = TRUE)




# ---- Load data ----

files <- list.files(cfg$data$save_dir, recursive = F)
files <- files[str_detect(files, ".rds")]
predictions <- list()
coefs <- list()
for(i in 1:length(files)) {
  
  i <- here(cfg$data$save_dir, files[i])
  dat <- readRDS(i)
  predictions <- append(predictions, list(dat$pred_result))
  coefs <-  append(coefs, list(dat$coef_result))
  
}


predictions <- predictions %>%
  bind_rows() %>% 
  mutate(pair = paste(pmin(inner_id, outer_id), 
                      pmax(inner_id, outer_id),
                      sep = "_")) %>%
  mutate(model = case_when(model %in% c("glmnet_classification", "glmnet_regression") ~ "Lasso",
                           model %in% c("ranger_classification", "ranger_regression") ~ "RF"))


length(unique(predictions$pair))

if (cfg$data$data_config == "Patient") {
  predictions <- predictions %>% 
    distinct(pair, patientID, model, stimulation, .keep_all = TRUE) %>%
    mutate(obs = patientID)
    
} else if (cfg$data$data_config == "Sample") {
  predictions <- predictions %>% 
    distinct(pair, patientID, timepoint, model, stimulation, .keep_all = TRUE) %>%
    mutate(obs = paste0(patientID, "_", timepoint))

}

coefs <- coefs %>% 
  bind_rows() %>%
  mutate(pair = paste(pmin(inner_id, outer_id, na.rm = TRUE), 
                      pmax(inner_id, outer_id, na.rm = TRUE),
                      sep = "_")) %>%
  mutate(model = case_when(model %in% c("glmnet_classification", "glmnet_regression") ~ "Lasso",
                           model %in% c("ranger_classification", "ranger_regression") ~ "RF")) %>%
  distinct(pair, term, model, stimulation, .keep_all = TRUE)











# ---- Model predictive performance ----

if (cfg$training$mode == "classification") {
  
  
  
  calculateRoc <- predictions %>%
    mutate(U = paste0(stimulation, "_", model)) %>%
    droplevels() %>%
    
    
    split(.$U) %>%
    purrr::map(., function(x) {
      stimulation_i <- unique(x$stimulation)
      model_i <- unique(x$model)
      U_i <- unique(x$U)
      
      r <- roc(x$response, x$.pred_1, ci = TRUE, direction="<")
      r.stats <- data.frame(r$ci) %>%
        t() %>%
        data.frame() %>%
        dplyr::rename(ci.lower = X1,
                      auc = X2,
                      ci.upper = X3) %>%
        mutate(label = paste0("AUC = ", round(auc, 2), " (", round(ci.lower, 2), " - ", round(ci.upper, 2), ")"))
      r.stats <- r.stats %>%
        mutate(label = paste0("AUC=", round(r.stats$auc, 2), " (", round(r.stats$ci.lower, 2), "-", round(r.stats$ci.upper, 2), ")")) %>%
        mutate(stimulation = stimulation_i,
               model = model_i,
               U = U_i)
      
      r <- smooth(r, method = "density", bw = 0.001)
      r.data <- data.frame(x = 1- r$specificities,
                           y = r$sensitivities) %>%
        mutate(stimulation = stimulation_i,
               model = model_i,
               U = U_i)
      
      result <- list(r.stats = r.stats,
                     r.data = r.data)
      return(result)
      
    }) 
  
  r.data <- list()
  r.stats <- list()
  i <- 1
  for(i in 1:length(calculateRoc)) {
    r.data[[i]] <- calculateRoc[[i]]$r.data
    r.stats[[i]] <-  calculateRoc[[i]]$r.stats
  }
  r.stats <- r.stats %>%
    bind_rows() %>%
    arrange(auc) %>%
    mutate(U = factor(U, levels = unique(.$U))) %>%
    mutate(x = 0.3,
           y = seq(0, 0.56, by = 0.04)) %>%
    mutate(label = paste0(U, " | ", label))
  
  SG_ROC <- r.data %>%
    bind_rows() %>%
    mutate(SG_flag = if_else(str_detect(stimulation, "^SG"), "SG", "alpha")) %>%
    ggplot(aes(x = x, y = y, colour = stimulation)) +
    geom_abline(intercept = 0, slope = 1, linetype = "solid", size = 0.5, alpha = 1) +
    geom_line(aes(group = U, linetype = model, alpha = SG_flag), size = 0.5) +
    theme_bw() +
    scale_colour_manual(values = model_colours) +
    scale_linetype_manual(values = c("Lasso" = "solid", "RF" = "dashed")) +
    scale_alpha_manual(values = c("alpha" = 0.5, "SG" = 1)) +
    geom_text(data = r.stats, aes(label = label), size = 5 / (14/5), hjust = 0) +
    scale_x_continuous("False Positive Rate (1-Specificity)") +
    scale_y_continuous("True Positive Rate (Sensitivity)") +
    coord_equal() +
    theme(axis.title.x = element_text(size = 6.5),
          axis.title.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6.5),
          axis.text.y = element_text(size = 6.5),
          plot.title = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 6.5),
          strip.text.y = element_text(size = 6.5, angle = 0),
          strip.background = element_blank(),
          legend.key.width = unit(1.5, "cm"),
          legend.key.height = unit(0.5, "cm"),
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.2),
          panel.border = element_rect(fill=NA, colour = "black", size=0.2),
          plot.subtitle = element_text(size = 6.5))
  SG_ROC
  
  ggsave(plot = SG_ROC,
         file = paste("SG_ROC", ".pdf", sep=""),
         path = cfg$data$result_dir,
         width = 10,
         height = 10,
         units = "cm",
         device = "pdf",
         useDingbats = FALSE)
  
  predictions <- predictions %>%
    mutate(response = factor(response, levels = c(0, 1))) %>%
    group_by(obs, response, model, stimulation) %>%
    dplyr::summarise(pred = mean(.pred_1)) %>%
    ungroup
  
  p.value <- predictions %>%
    filter(stimulation == "SG") %>%
    wilcox_test(pred ~ response)
  p.value <- p.value %>%
    mutate(p = paste0("p.value = ", p))
  
  
  SG_preds <- predictions %>%
    filter(stimulation == "SG") %>%
    ggplot(aes(x = response, y = pred)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.15, height = 0, shape = 21, size = 2, colour = "darkgrey", aes(fill = response)) +
    theme_bw() +
    scale_fill_manual(values = c("1" = "#F8931D", "0" = "#034E7B")) +
    theme(axis.title.x = element_text(size = 6.5),
          axis.title.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6.5),
          axis.text.y = element_text(size = 6.5),
          plot.title = element_text(size = 6.5),
          plot.subtitle = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 6.5),
          strip.text.y = element_text(size = 6.5, angle = 0),
          strip.background = element_blank(),
          legend.key.width = unit(0.2, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.2),
          panel.border = element_rect(fill=NA, colour = "black", size=0.2),
          plot.caption = element_text(size = 4.5)) +
    scale_x_discrete(labels = c("Term", "PreTerm")) +
    geom_text(data = p.value, aes(x = 1.5, y = 1.1, label = p), size = 6.5 / (14/5)) +
    scale_y_continuous(breaks = seq(0,1,by= 0.2))+
    labs(y = "prediction", x = "")
  SG_preds
  ggsave(plot = SG_preds,
         file = paste("SG_preds", ".pdf", sep=""),
         path = cfg$data$result_dir,
         width = 4,
         height = 6,
         units = "cm",
         device = "pdf",
         useDingbats = FALSE)
  
  r.stats <- r.stats %>%
    dplyr::select(-c(x, y, U))
  
  
  

  
  
} else if (cfg$training$mode == "regression") {
  
  predictions <- predictions %>%
    mutate(U = paste0(stimulation, "_", model)) %>%
    droplevels() %>%
    group_by(obs, response, model, stimulation, U) %>%
    dplyr::summarise(pred = median(.pred)) %>%
    ungroup 
 
  
  r.stats <- predictions %>%
    split(.$U) %>%
    purrr::map(., function(x) {
      stimulation_i <- unique(x$stimulation)
      model_i <- unique(x$model)
      U_i <- unique(x$U)
      
      r.stats <- x %>%
        dplyr::summarise(
          n        = n(),
          spearman = cor(response, pred, method = "spearman", use = "pairwise.complete.obs"),
          rmse     = sqrt(mean((pred - response)^2, na.rm = TRUE)),
          p_value  = suppressWarnings(cor.test(response, pred, method = "spearman", exact = FALSE)$p.value)
        ) %>%
        mutate(stimulation = stimulation_i,
               model = model_i,
               U = U_i)
      
      
      return(r.stats)
      
    }) %>%
    bind_rows() %>%
    arrange(desc(spearman))
  
  p.value <- r.stats %>%
    filter(stimulation == "SG") %>%
    pivot_longer(cols = c(n, spearman, rmse, p_value), names_to = "metric") %>%
    mutate(metric = factor(metric, levels = c("n", "spearman", "p_value", "rmse"))) %>%
    arrange(desc(metric)) %>%
    mutate(label = paste0(metric, " = ", round(value, 3))) %>%
    mutate(x = 38,
           y = seq(0, 0.03, by = 0.01) * 40 + 32)
    
    
  
  
  SG_preds <- predictions %>%
    filter(stimulation == "SG") %>%
    ggplot(aes(x = response, y = pred)) +
    geom_point(shape = 21, size = 2, colour = "black", fill = "darkgrey") +
    theme_bw() +
    theme(axis.title.x = element_text(size = 6.5),
          axis.title.y = element_text(size = 6.5),
          axis.text.x = element_text(size = 6.5),
          axis.text.y = element_text(size = 6.5),
          plot.title = element_text(size = 6.5),
          plot.subtitle = element_text(size = 6.5),
          legend.text = element_text(size = 6.5),
          legend.title = element_text(size = 6.5),
          strip.text.y = element_text(size = 6.5, angle = 0),
          strip.background = element_blank(),
          legend.key.width = unit(0.2, "cm"),
          legend.key.height = unit(0.2, "cm"),
          legend.position = "none",
          panel.grid = element_blank(),
          axis.ticks = element_line(colour = "black", size = 0.2),
          panel.border = element_rect(fill=NA, colour = "black", size=0.2),
          plot.caption = element_text(size = 4.5)) +
    coord_cartesian(ylim = c(32, 40)) +
    coord_equal() +
    scale_x_continuous(breaks = c(32, 34, 36, 38, 40)) +
    scale_y_continuous(breaks = c(32, 34, 36, 38, 40)) +
    geom_smooth(method = "lm") +
    geom_text(data = p.value, aes(x = x, y = y, label = label), size = 6.5 / (14/5)) +
    labs(y = "prediction", x = "")
  SG_preds
  
  ggsave(plot = SG_preds,
         file = paste("SG_preds", ".pdf", sep=""),
         path = cfg$data$result_dir,
         width = 6,
         height = 6,
         units = "cm",
         device = "pdf",
         useDingbats = FALSE)
    
    
  
}



# ---- SG Model feature attribution ---- 

SG_coefs <- coefs %>%
  filter(stimulation == "SG") %>%
 
  mutate(selected = case_when(estimate != 0 ~ 1,
                              TRUE ~ 0)) %>%
  group_by(term, model, stimulation) %>%
  dplyr::summarise(estimate = median(estimate),
                   n_selected = sum(selected),
                   n_iterations = n()) %>%
  mutate(freq_selected = round(n_selected/n_iterations, 3)) %>%
  arrange(desc(estimate)) %>%
  mutate(stimulation = str_split(term, pattern = "_", simplify = T)[,1],
         model = str_split(term, pattern = "_", simplify = T)[,2]) %>%
  mutate(model = case_when(model %in% c("glmnet") ~ "Lasso",
                           model %in% c("ranger") ~ "RF")) %>%
  mutate(U = paste0(stimulation, "_", model))

SG_coefs_plot <- SG_coefs %>%
  mutate(U = factor(U, levels = unique(.$U))) %>%
  ggplot(aes(x = U, y = estimate, fill = stimulation)) +
  geom_col() +
  scale_fill_manual(values = model_colours) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw() +
  labs(y = "median coefficient", x = "") +
  theme(axis.text.x = element_text(size = 6.5, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_blank(),
        strip.background = element_rect(colour="white", fill="white"), 
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 6.5), 
        axis.title.y = element_text(size = 6.5),
        strip.text.x = element_text(size = 6.5),
        legend.text = element_text(size = 6.5),
        legend.title = element_text(size = 6.5),
        legend.position = "none",
        plot.caption = element_text(size = 6.5),
        plot.title = element_text(size = 5),
        axis.ticks = element_line(colour = "black", size = 0.2),
        panel.border = element_rect(fill=NA, colour = "black", size=0.2),
        legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.2, "cm"))

ggsave(plot = SG_coefs_plot,
       file = paste("SG_coefs", ".pdf", sep=""),
       path = cfg$data$result_dir,
       width = 5,
       height = 7.25,
       units = "cm",
       device = "pdf",
       useDingbats = FALSE)



# ---- Write model metrics ---- 

fname <- "model_metrics.xlsx"
r.stats %>%
  write.xlsx(., file = here(cfg$data$result_dir, fname),  rowNames = FALSE)
  
fname <- "model_predictions.xlsx"
predictions %>%
  write.xlsx(., file = here(cfg$data$result_dir, fname), rowNames = FALSE)

fname <- "sg_coefs.xlsx"
SG_coefs %>%
  write.xlsx(., file = here(cfg$data$result_dir, fname), rowNames = FALSE)




# ---- Submodel feature attribution ---- 



selected_models <- SG_coefs %>%
  filter(estimate > 0) %>%
  pull(U)
total_hits <- coefs %>%
  mutate(U = paste0(stimulation, "_", model)) %>%
  
  
  filter(U %in% selected_models) %>%
  dplyr::select(-c(inner_id, outer_id, penalty)) %>%
  mutate(term = str_remove(term, pattern = "raw\\+\\+")) %>%
  split(.$U) %>%
  purrr::map(., function(x) {
    
    x <- x %>%
      pivot_wider(names_from = pair,
                  values_from = estimate) %>%
      pivot_longer(cols = -c(term, model, stimulation, U),
                   names_to = "pair", 
                   values_to = "estimate"
      ) %>%
      # ---- Calculate median coefficient / importance scores and the frequency of selection ----
      mutate(estimate = case_when(is.na(estimate) == T ~ 0,
                                  TRUE ~ estimate)) %>%
      mutate(selected = case_when(estimate != 0 ~ 1,
                                  TRUE ~ 0)) %>%
      group_by(term, model, stimulation, U) %>%
      dplyr::summarise(estimate = median(estimate),
                       n_selected = sum(selected),
                       n_iterations = n()) %>%
      ungroup %>%
      mutate(freq_selected = n_selected/n_iterations) %>%
      mutate(freq_selected = round(freq_selected, 5)) %>%
      arrange(desc(n_selected)) %>%
      filter(abs(estimate) > 0) %>%
      data.frame() 
    return(x)
    
  }) %>%
  bind_rows() %>%
  mutate(importance = abs(estimate))
maxes <- total_hits %>%
  dplyr::group_by(U) %>%
  dplyr::summarise(max_importance = max(importance))



metadata <- read.csv(file = cfg$data$sample_metadata) %>%
  dplyr::select(patientID, GA_sampling, GA_delivery) %>%
  unique()
totalData <- readRDS(file = cfg$data$processed_data) %>%
  left_join(., metadata, by = c("patientID", "GA_sampling")) %>%
  mutate(variable = str_split(variable, pattern = "__", simplify = T)[,1]) %>%
  dplyr::select(patientID, group, GA_delivery, timepoint, GA_sampling, stimulation, variable, value) 

totalData <- totalData %>%
  filter(timepoint %in% cfg$data$timepoints) 


if (cfg$data$data_config == "Patient") {
  totalData <- totalData %>%
    mutate(variable = paste0(variable, "__", stimulation, "_", timepoint)) %>%
    dplyr::select(-timepoint)
} else if (cfg$data$data_config == "Sample") {
  totalData <- totalData %>%
    mutate(variable = paste0(variable, "__", stimulation))
}

if (cfg$training$mode == "classification") {
  univariable_results <- totalData %>%
    filter(variable %in% total_hits$term) %>%
    mutate(response = factor(group, levels = c(0, 1))) %>%
    split(.$variable) %>%
    purrr::map(., function(x) {
      stats <- x %>%
        wilcox_test(value ~ response) %>%
        mutate(term = unique(x$variable))
      median_1 <- x %>% filter(response == "1") %>% pull(value) %>% median(.)
      median_0 <- x %>% filter(response == "0") %>% pull(value) %>% median(.)
      sign_diff <- sign(median_1 - median_0)
      effectsize <- x %>%
        wilcox_effsize(value ~ response, ci = FALSE)
      stats <- stats %>%
        mutate(r = effectsize$effsize * sign_diff)
      return(stats)
    }) %>%
    bind_rows() %>%
    arrange(p) %>%
    dplyr::select(-c(.y., statistic))
  
  
  feature_importance <- total_hits %>%
    left_join(., maxes, by = "U") %>%
    mutate(relative_importance = importance / max_importance) %>%
    dplyr::select(-max_importance) %>%
    left_join(., univariable_results, by = "term") %>%
    mutate(index = relative_importance * -log10(p)) %>%
    arrange(desc(index)) %>%
    rownames_to_column("rank") %>%
    select(rank, term, model, stimulation, estimate, relative_importance, n_selected, n_iterations, freq_selected, group1, group2, n1, n2, p, r, index)

  
  } else if (cfg$training$mode == "regression") {
    
    
  univariable_results <- totalData %>%
    filter(variable %in% total_hits$term) %>%
    dplyr::rename(response = GA_delivery) %>%
    split(.$variable) %>%
    purrr::map(., function(x) {
      stats <- x %>%
        dplyr::summarise(
          n        = n(),
          r = cor(response, value, method = "spearman", use = "complete.obs"),
          p = suppressWarnings(cor.test(response, value, method = "spearman", exact = FALSE)$p.value)
        ) %>%
        mutate(term = unique(x$variable))
      return(stats)
    }) %>%
    bind_rows() %>%
    arrange(p) 
    
  feature_importance <- total_hits %>%
    left_join(., maxes, by = "U") %>%
    mutate(relative_importance = importance / max_importance) %>%
    dplyr::select(-max_importance) %>%
    left_join(., univariable_results, by = "term") %>%
    mutate(index = relative_importance * -log10(p)) %>%
    arrange(desc(index)) %>%
    rownames_to_column("rank") %>%
    select(rank, term, model, stimulation, estimate, relative_importance, n_selected, n_iterations, freq_selected, n, p, r, index)

}





fname <- "feature_importance.xlsx"
feature_importance %>%
  write.xlsx(., file = here(cfg$data$result_dir, fname), rowNames = FALSE)


# ---- Confounder analyses ---- 

confounders <- read.csv(here(cfg$data$sample_metadata)) %>%
  select(patientID, Age, Parity, Gravidity,  BMI, Infant_Sex, History_priorPTB) %>%
  unique() %>%
  mutate(patientID = as.character(patientID))
selected_models <- c(selected_models, "SG_Lasso")
model_preds <- predictions %>%
  mutate(U = paste0(stimulation, "_", model)) %>%
  filter(U %in% selected_models) %>%
  mutate(patientID = str_split(obs, pattern = "_", simplify = T)[,1]) %>%
  group_by(obs, patientID, response, model, stimulation, U) %>%
  dplyr::summarise(pred = median(pred)) %>%
  ungroup() %>%
  select( obs,patientID,response, U,  pred) %>%
  pivot_wider(names_from = "U", values_from = "pred")
confounder.table <- model_preds %>%
  full_join(., confounders, by = "patientID") %>%
  pivot_longer(cols = selected_models, names_to = "model", values_to = "pred") %>%
  select(patientID, response, model, pred, Age, BMI, Gravidity, Parity, History_priorPTB, Infant_Sex ) %>%
  split(.$model) %>%
  purrr::map(., function(x) {
    model.i <- unique(x$model)
    
    x <- x %>%
      na.omit()
    if(cfg$training$mode == "regression") {
      fit.tmp <- lm(response ~ pred + Age + BMI + Gravidity + Parity + History_priorPTB + Infant_Sex, data = x)
    } else if (cfg$training$mode == "classification") {
      fit.tmp <- glm(response ~ pred + Age + BMI + Gravidity + Parity + History_priorPTB + Infant_Sex, data = x, family = binomial(link = "logit"))
    }
    summ.tmp <- summary(fit.tmp)
    result <- summ.tmp$coefficients %>%
      data.frame() %>%
      mutate(model = model.i) %>%
      rownames_to_column("param") %>%
      filter(param != "(Intercept)")
    return(result)
    
    
  }) %>%
  bind_rows() 

fname <- "confounder_table.xlsx"
confounder.table %>%
  write.xlsx(., file = here(cfg$data$result_dir, fname), rowNames = FALSE)




# ---- Write a sentinel file for snakemake to track progression ----
done_file <- here(cfg$data$save_dir, "MODEL_ANALYSIS_DONE.txt")
writeLines(sprintf("Finished at %s", Sys.time()), done_file)




