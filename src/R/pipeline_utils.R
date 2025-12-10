library(tidymodels) 
library(stringr)
library(healthyR.ai)  
library(vip)
library(here)


`%out%` <- function(a,b) ! a %in% b 

LOO_splitter <- function(id, data) {
  
  data_train <- data %>%
    filter(patientID %out% id)
  
  data_test <- data %>%
    filter(patientID %in% id)
  
  ret <- list(data_train = data_train,
              data_test = data_test)
  
  return(ret)
}






fitModel <- function(trainingdata, 
                     weighted = TRUE, 
                     metadata_cols = c("patientID", "timepoint"),
                     repeats = 5, 
                     folds = 20,
                     seed = 123,
                     mode = c("classification", "regression"),
                     model_type = c("rf", "lasso", "sg")) {
  
  set.seed(seed)

  
  if (mode == "classification") {
    
    missing_impute <- function(rec) step_impute_median(rec, all_predictors())
    step_winsor <- function(rec) step_hai_winsorized_truncate(rec, all_predictors(), fraction = 0.05)
    step_rm_ <- function(rec) step_rm(rec, starts_with("raw++"))
    zv_filter <- function(rec) step_zv(rec, all_predictors(), group = "response") 
    step_center_ <- function(rec) step_center(rec, all_predictors())
    step_scale_ <- function(rec) step_scale(rec, all_predictors())
    

    
    if (model_type == "rf") {
      
      model_spec <- 
        rand_forest(trees = 1000) %>%
        set_engine("ranger", importance = "impurity", seed = seed, num.threads = 1) %>%
        set_mode(mode)

      
    } else if (model_type == "lasso") {
      model_spec <- 
        logistic_reg(mode = mode, penalty = tune(), mixture = 1) %>%
        set_engine("glmnet") 


    } else if (model_type == "sg") {
      model_spec <-
        logistic_reg(mode = mode, penalty = tune(), mixture = 1) %>%
        set_engine("glmnet", lower.limits = 0)
      missing_impute <- function(rec) step_impute_mean(rec, all_predictors())
      step_winsor <- function(rec) rec
      step_rm_ <- function(rec) rec
      zv_filter <- function(rec) rec
      step_center_ <- function(rec) rec
      step_scale_ <- function(rec) rec
      
      
      
      
    }
    
    
    train.y <-  trainingdata$response
    fraction_0 <- rep(1 - sum(train.y == 0) / length(train.y), sum(train.y == 0))
    fraction_1 <- rep(1 - sum(train.y == 1) / length(train.y), sum(train.y == 1))
    weights <- numeric(length(train.y))
    if (weighted == TRUE) {
      weights[train.y == 0] <- fraction_0
      weights[train.y == 1] <- fraction_1
    } else {
      weights <- rep(1, length(train.y))
    }  
    metric <- metric_set(roc_auc)
    metric_ <- "roc_auc"
    
    
    
  } else if (mode == "regression") {
    
    missing_impute <- function(rec) step_impute_median(rec, all_predictors())
    step_winsor <- function(rec) step_hai_winsorized_truncate(rec, all_predictors(), fraction = 0.05)
    step_rm_ <- function(rec) step_rm(rec, starts_with("raw++"))
    zv_filter <- function(rec) step_zv(rec, all_predictors())
    step_center_ <- function(rec) step_center(rec, all_predictors())
    step_scale_ <- function(rec) step_scale(rec, all_predictors())
    
    
    if (model_type == "rf") {
      
      model_spec <- 
        rand_forest(trees = 1000) %>%
        set_engine("ranger", importance = "impurity", seed = seed, num.threads = 1) %>%
        set_mode(mode)
      

    } else if (model_type == "lasso") {
      
      model_spec <-       
        linear_reg(mode = mode, penalty = tune(), mixture = 1) %>%
        set_engine("glmnet") 
      
      
      
    } else if (model_type == "sg") {
      model_spec <-
        linear_reg(mode = mode, penalty = tune(), mixture = 1) %>%
        set_engine("glmnet", lower.limits = 0)
      missing_impute <- function(rec) step_impute_mean(rec, all_predictors())
      step_winsor <- function(rec) rec
      step_rm_ <- function(rec) rec
      zv_filter <- function(rec) rec
      step_center_ <- function(rec) rec
      step_scale_ <- function(rec) rec
      
    }
    
    
    train.y <-  trainingdata$response
    weights <- rep(1, length(train.y))
    metric <- metric_set(rmse)
    metric_ <- "rmse"
    
    
  }

  trainingdata <-  trainingdata %>%
    mutate(
      case_wts = weights,
      case_wts = importance_weights(case_wts)
    )
  rec <- recipe(response ~ ., data =  trainingdata)
  metadata_cols <- metadata_cols[metadata_cols %in% colnames(trainingdata)]
  for(i in metadata_cols) {
    rec <- rec %>%
      update_role(i, new_role = i)
  }
  rec <- rec %>%
    missing_impute %>%
    step_winsor %>%
    step_rm_ %>%
    zv_filter %>%
    step_center_ %>%
    step_scale_
  
  
 
  
  
  wfset <-
    workflow_set(
      preproc = list(preprocess = rec),
      models = list(model_spec = model_spec),
      case_weights = case_wts 
    )
  
  
  grid_folds <- group_vfold_cv(
    trainingdata,
    group   = patientID,
    v       = folds,
    repeats = repeats,
    strata  = response
  )
  
  grid_ctrl <- control_grid( 
    verbose = FALSE,
    save_pred = TRUE,
    save_workflow = FALSE,
    parallel_over = "resamples",
    allow_par = FALSE
  )
  
  
  tune_results <-
    wfset %>%
    workflow_map(
      "tune_grid",
      seed = seed,
      resamples = grid_folds,
      grid = 100,
      control = grid_ctrl,
      verbose = FALSE,
      metrics = metric
    )
  
  tuned <- tune_results %>%
    extract_workflow_set_result("preprocess_model_spec") %>% 
    select_best(metric = metric_)
  model <- tune_results %>%
    extract_workflow("preprocess_model_spec") %>%
    finalize_workflow(tuned) %>%
    fit(data = trainingdata)
  
  return(model)
  
}


predictModel <- function(testingdata, model, metadata_cols = c("patientID", "timepoint"), top25 = FALSE) {
  
  metadata_idx <- which(colnames(testingdata) %in% c("response", metadata_cols))
  
  model_spec <- workflows::extract_spec_parsnip(model)
  model_settings = paste0(model_spec$engine, "_", model_spec$mode)

  if (model_spec$mode == "classification") {
    pred <- cbind(testingdata[,metadata_idx], predict(model, testingdata, type = "prob")[,2]) %>%
      mutate(model = model_settings) 
    
  } else if (model_spec$mode == "regression") {
    pred <- cbind(testingdata[,metadata_idx], predict(model, testingdata)) %>%
      mutate(model = model_settings) 
  }
  
  
  if (model_spec$engine == "glmnet") {
    feature <- model %>%
      pull_workflow_fit() %>%
      tidy() %>%
      mutate(model = model_settings,
             term = str_remove(term, "winsorized_truncate_"))
  
    } else if (model_spec$engine == "ranger") {
    feature <- model  %>%
      pull_workflow_fit() %>%
      vip(num_features = model$fit$fit$fit$num.independent.variables)
    
    feature <- feature$data %>%
      mutate(Variable = str_remove(Variable, "winsorized_truncate_")) %>%
      mutate(model = model_settings,
             penalty = NA) %>%
      dplyr::rename(estimate = Importance,
                    term = Variable)
    
    if(top25 ==T) {
      feature <- feature %>%
        arrange(desc(estimate)) %>%
        dplyr::slice(1:25)
      } 
    
    
  }

  
  result <- list(prediction = pred,
                 feature = feature)
  return(result)
}


process_i_submodel <- function(iteration, totalData, save_dir = "cfg$data$save_dir") {
  
  stimulation_name <- str_split(iteration, pattern = "_", simplify = T)[,1]
  outer_id <- str_split(iteration, pattern = "_", simplify = T)[,2]
  file = paste0(str_remove_all(currentDate, "-"),"_", cfg$data$data_config, "_", paste0(cfg$data$timepoints, collapse = ""), "_", stimulation_name, "_", outer_id, ".rds")
  file = here(save_dir, file)
  
  
  outer_split <- totalData %>%
    filter(stimulation == stimulation_name) %>%
    LOO_splitter(outer_id, .) %>%
    purrr::map(., function(x) {
      x <- x %>%
        dplyr::select(-stimulation) %>%
        pivot_wider(names_from = variable)
      return(x)
    })

  
  inner_ids <- unique(outer_split$data_train$patientID)
  
  pred_result <- list()
  coef_result <- list()
  for(inner_id in inner_ids){
    print(inner_id)
    inner_split <- LOO_splitter(inner_id, outer_split$data_train)
    
    data_test <- bind_rows(inner_split$data_test, outer_split$data_test)
    
    
    LassoTrained <- fitModel(trainingdata = inner_split$data_train,
                             mode = cfg$training$mode,
                             model_type = "lasso")
    LassoResult <- predictModel(testingdata = data_test, 
                                model = LassoTrained)
    
    pred_ <- LassoResult$prediction %>%
      mutate(inner_id = inner_id,
             outer_id = outer_id,
             stimulation = stimulation_name)
    coef_ <- LassoResult$feature %>%
      filter(term != "(Intercept)") %>%
      mutate(inner_id = inner_id,
             outer_id = outer_id,
             stimulation = stimulation_name)
    
    pred_result <- append(pred_result, list(pred_))
    coef_result <-  append(coef_result, list(coef_))
    
    
    RFTrained <- fitModel(inner_split$data_train,
                          mode = cfg$training$mode,
                          model_type = "rf")
    RFresult <- predictModel(data_test, RFTrained, top25 = TRUE)
    data_train_filt <- inner_split$data_train %>%
      dplyr::select(any_of(c("patientID", "response", "timepoint")), RFresult$feature$term)
    data_test_filt <- data_test %>%
      select(any_of(c("patientID", "response", "timepoint")), RFresult$feature$term)
    
    RFTrained <- fitModel(data_train_filt,
                          mode = cfg$training$mode,
                          model_type = "rf")
    RFresult <- predictModel(data_test_filt, RFTrained, top25 = FALSE)
    
    
    pred_ <- RFresult$prediction %>%
      mutate(inner_id = inner_id,
             outer_id = outer_id,
             stimulation = stimulation_name)
    coef_ <- RFresult$feature %>%
      mutate(inner_id = inner_id,
             outer_id = outer_id,
             stimulation = stimulation_name)
    
    pred_result <- append(pred_result, list(pred_))
    coef_result <-  append(coef_result, list(coef_))
    
    
    
  } 
  pred_result <- pred_result %>%
    bind_rows()
  coef_result <- coef_result %>%
    bind_rows()
  
  ret <- list(pred_result = pred_result,
              coef_result = coef_result)
  
  saveRDS(ret, file = file)
  return(ret)
  
  
}



process_i_sg <- function(iteration, predictions, save_dir = "cfg$data$save_dir") {
  
  stimulation_name <- "SG"
  file = paste0(str_remove_all(currentDate, "-"),"_", cfg$data$data_config, "_", paste0(cfg$data$timepoints, collapse = ""), "_", stimulation_name, "_", iteration, ".rds")
  file = here(save_dir, file)
  
  predictions_fold <- predictions %>%
    filter(outer_id == iteration)
  
  outer_split <- predictions_fold %>%
    LOO_splitter(iteration, .) 

  
  LassoTrained <- fitModel(trainingdata = outer_split$data_train,
                           metadata_cols = c("patientID", "timepoint", "inner_id", "outer_id"),
                           mode = cfg$training$mode,
                           model_type = "sg")
  
  LassoResult <- predictModel(testingdata = outer_split$data_test, 
                              metadata_cols = c("patientID", "timepoint", "inner_id", "outer_id"),
                              model = LassoTrained)
  
  pred_result <- LassoResult$prediction %>%
    mutate(stimulation = "SG")
  coef_result <- LassoResult$feature %>%
    filter(term != "(Intercept)") %>%
    mutate(inner_id = NA,
           outer_id = iteration,
           stimulation = "SG")
  
  ret <- list(pred_result = pred_result,
              coef_result = coef_result)
  
  saveRDS(ret, file = file)
  return(ret)
}

