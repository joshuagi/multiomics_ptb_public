# ---- Load dependencies and script arguments ----
library(yaml)
library(parallel)
library(doParallel)
library(foreach)
library(here)
source(here("src", "R", "pipeline_utils.R"))

currentDate <- Sys.Date() 

args <- commandArgs(trailingOnly = TRUE)


cfg  <- yaml::yaml.load_file(args[1])
cfg$data$save_dir          <- here(cfg$data$save_dir)
cfg$data$result_dir        <- here(cfg$data$result_dir)
cfg$data$processed_data    <- here(cfg$data$processed_data)
cfg$data$sample_metadata   <- here(cfg$data$sample_metadata)





# ---- Load data ----
print('loading data')
files <- list.files(cfg$data$save_dir, recursive = F)
files <- files[str_detect(files, ".rds")]
i<-1
predictions <- list()
for(i in 1:length(files)) {
  
  i <- here(cfg$data$save_dir, files[i])
  dat <- readRDS(i)$pred_result
  predictions <- append(predictions, list(dat))
  
}
predictions <- predictions %>%
  bind_rows() %>%
  filter(stimulation != "SG")


if (cfg$training$mode == "classification") {
  predictions  <- predictions  %>%
    mutate(response = factor(response, levels = c(0,1)))
} else if (cfg$training$mode == "regression") {
  predictions <- predictions  %>%
    dplyr::rename(.pred_1 = .pred)

  
}

if (cfg$data$data_config == "Patient") {
  predictions <- predictions %>%
    mutate(model = paste0(stimulation, "_", model)) %>%
    pivot_wider(
      id_cols   = c(patientID, response, inner_id, outer_id),
      names_from  = model,
      values_from = .pred_1
    )
} else if (cfg$data$data_config == "Sample") {
  predictions <- predictions %>%
    mutate(model = paste0(stimulation, "_", model)) %>%
    pivot_wider(
      id_cols   = c(patientID, response, timepoint, inner_id, outer_id),
      names_from  = model,
      values_from = .pred_1
    )
  
}



# ---- Initialise iterations  ----


all_iterations <-  predictions %>% pull(outer_id) %>% unique

files <- list.files(cfg$data$save_dir, pattern = "_SG_[0-9]+\\.rds$", full.names = FALSE) %>% str_remove(., pattern = ".rds")
m <- str_split(files, pattern = "_", simplify = TRUE)
done_iterations <- m[,5]
iterations <- setdiff(all_iterations, done_iterations)

message(sprintf(
  "%d remaining iterations",
  length(iterations)
))

if (length(iterations) == 0L) {
  message("Nothing to do, all iterations completed.")
  done_file <- here(cfg$data$save_dir, "SG_DONE.txt")
  writeLines(sprintf("Finished at %s", Sys.time()), done_file)
  quit(save = "no")
}









# ---- Train model  ----

print("registering parallel backend")

threads_requested   <- as.integer(args[2])
max_cores <- max(1L, detectCores() - 1L)
ncores    <- min(threads_requested, max_cores)

if (ncores < threads_requested) {
  message(sprintf("Requested %d cores, capping at %d (available-1)", threads_requested, ncores))
} else {
  message(sprintf("Using %d cores", ncores))
}


try(stopCluster(cl), silent = TRUE); closeAllConnections()

cl <- makeCluster(ncores, type = "FORK")
registerDoParallel(cl)

pkgs <- c("dplyr","stringr","readr","tibble",
          "tidymodels","glmnet","ranger","vip","healthyR.ai")

print('running analyses')

pipeline_result <- foreach(
  i = iterations,
  .packages = pkgs,
  .multicombine = TRUE,
  .inorder = FALSE,
  .verbose = TRUE
) %dopar% {
  result <- process_i_sg(i, predictions, save_dir = cfg$data$save_dir)
  return(result)
}
stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()
closeAllConnections()


# ---- Write a sentinel file for snakemake to track progression ----
done_file <- here(cfg$data$save_dir, "SG_DONE.txt")
writeLines(sprintf("Finished at %s", Sys.time()), done_file)
