### Parse simulation output and get summary RDS files.

rm(list = ls(all.names = TRUE))

library(parallel)

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

cutoffs <- c("nu0.50", "nu0.65", "nu0.80", "nu0.95")



for (cutoff in cutoffs) {

  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]

  POMS_func.based_output <- readRDS(file = paste("POMS_out/ref.based_POMS_sim_func.based_1000reps_cutoff_", cutoff, ".rds", sep = ""))

  regress_func.based_output <- readRDS(file = paste("regress_out/ref.based_regress_sim_func.based_1000reps_cutoff_", cutoff, ".rds", sep = ""))

  wilcoxon.musicc_func.based_output <- readRDS(file = paste("wilcoxon.musicc_out/ref.based_wilcoxon.musicc_sim_func.based_500reps_cutoff_", cutoff, ".rds", sep = ""))

  func.based_summary <- simulation_summaries(POMS_sims = POMS_func.based_output[1:500],
                                             regress_sims = regress_func.based_output[1:500],
                                             alt_tool_sims = wilcoxon.musicc_func.based_output,
                                             focal_func_present = TRUE,
                                             func_table = ref_func_cutoff,
                                             sig_cutoffs = 0.05,
                                             num_cores = 20)
  
  saveRDS(object = func.based_summary,
          file = paste("simulation_summaries/func.based_cutoff_", cutoff, ".rds", sep = ""))

}

