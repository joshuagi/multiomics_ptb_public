# ---- Load dependencies and script arguments ----
library(yaml)
library(parallel)
library(doParallel)
library(foreach)
library(here)
Sys.setenv(NO_COLOR = "1", CLICOLOR_FORCE = "0", R_CLI_NUM_COLORS = "0")
options(crayon.enabled = FALSE, cli.num_colors = 0)
source(here("src", "R", "pipeline_utils.R"))

currentDate <- Sys.Date() 

args <- commandArgs(trailingOnly = TRUE)

cfg  <- yaml::yaml.load_file(args[1])
cfg$data$save_dir          <- here(cfg$data$save_dir)
cfg$data$result_dir        <- here(cfg$data$result_dir)
cfg$data$processed_data    <- here(cfg$data$processed_data)
cfg$data$sample_metadata   <- here(cfg$data$sample_metadata)


dir.create(cfg$data$save_dir,recursive = TRUE)




# ---- Load data ----

print('loading data')
metadata <- read.csv(file = cfg$data$sample_metadata) %>%
  dplyr::select(patientID, GA_sampling, GA_delivery) %>%
  unique()
totalData <- readRDS(file = cfg$data$processed_data) %>%
  left_join(., metadata, by = c("patientID", "GA_sampling")) %>%
  mutate(variable = str_split(variable, pattern = "__", simplify = T)[,1]) %>%
  dplyr::select(patientID, group, GA_delivery, timepoint, GA_sampling, stimulation, variable, value) %>%
  mutate(variable = paste0("raw++", variable))


totalData <- totalData %>%
  filter(timepoint %in% cfg$data$timepoints) %>%
  select(-GA_sampling)

if (cfg$training$mode == "classification") {
  totalData <- totalData %>%
    select(-GA_delivery) %>%
    mutate(group = factor(group, levels = c(0,1))) %>%
    dplyr::rename(response = group)
} else if (cfg$training$mode == "regression") {
  totalData <- totalData %>%
    select(-group) %>%
    dplyr::rename(response = GA_delivery)
}



if (cfg$data$data_config == "Patient") {
  totalData <- totalData %>%
    mutate(variable = paste0(variable, "__", stimulation, "_", timepoint)) %>%
    dplyr::select(-timepoint)
  
  pairs <- totalData %>% distinct(stimulation, variable)
  patients <- totalData %>% distinct(patientID, response)
  totalData <- patients %>%
    crossing(pairs) %>%
    left_join(.,
              totalData %>%
                dplyr::select(-response),
              by = c("stimulation", "patientID", "variable")) %>%
    dplyr::select(patientID, response, stimulation, any_of("timepoint"), variable, value)

  } else if (cfg$data$data_config == "Sample") {
    
    totalData <- totalData %>%
      mutate(variable = paste0(variable, "__", stimulation))

    pairs <- totalData %>% distinct(stimulation, variable) 
    patients <- totalData %>% distinct(patientID, response, timepoint)
    totalData <- patients %>%
      crossing(pairs) %>%                               
      left_join(.,
              totalData %>%
                dplyr::select(-response),
              by = c("stimulation", "patientID", "variable", "timepoint")) %>%
      dplyr::select(patientID, response, stimulation, any_of("timepoint"), variable, value)
  
  
}










# ---- Initialise iterations  ----

all_iterations <- totalData %>% mutate(iteration = paste0(stimulation, "_", patientID)) %>% pull(iteration) %>% unique


files <- list.files(cfg$data$save_dir, pattern = "\\.rds$", full.names = FALSE) %>% str_remove(., pattern = ".rds")
m <- str_split(files, pattern = "_", simplify = TRUE)
done_iterations <- paste0(m[,4], "_", m[,5])
iterations <- setdiff(all_iterations, done_iterations)

message(sprintf(
  "%d remaining iterations",
  length(iterations)
))

if (length(iterations) == 0L) {
  message("Nothing to do, all iterations completed.")
  done_file <- here(cfg$data$save_dir, "SUBMODELS_DONE.txt")
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

Sys.setenv(OMP_NUM_THREADS=1, OPENBLAS_NUM_THREADS=1, MKL_NUM_THREADS=1,
           BLIS_NUM_THREADS=1, VECLIB_MAXIMUM_THREADS=1)

cl <- makeCluster(ncores, type = "FORK")
registerDoParallel(cl)

lockBinding("totalData", .GlobalEnv)

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
    result <- process_i_submodel(i, totalData, save_dir = cfg$data$save_dir)
    return(result)
    }
stopCluster(cl)
stopImplicitCluster()
registerDoSEQ()
closeAllConnections()


# ---- Write a sentinel file for snakemake to track progression ----
done_file <- here(cfg$data$save_dir, "SUBMODELS_DONE.txt")
writeLines(sprintf("Finished at %s", Sys.time()), done_file)
