### Parse simulation output and get summary RDS files.

rm(list = ls(all.names = TRUE))

library(parallel)

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

cutoffs <- c("nu0.50", "nu0.65", "nu0.80", "nu0.95")



for (cutoff in cutoffs) {

  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]
  
  func_rand_output <- readRDS(file = paste("POMS_out/ref.based_POMS_sim_rand_func_1000reps_cutoff_", cutoff, ".rds", sep = ""))

  func_rand_wilcoxon.musicc <- readRDS(file = paste("wilcoxon.musicc_out/ref.based_wilcoxon.musicc_sim_rand_func_500reps_cutoff_", cutoff, ".rds", sep = ""))
  
  func_rand_summary_wilcoxon.musicc <- simulation_summaries(POMS_sims = func_rand_output[1:500],
                                                            alt_tool_sims = func_rand_wilcoxon.musicc,
                                                            focal_func_present = TRUE,
                                                            func_table = ref_func_cutoff,
                                                            num_cores = 20)
  
  saveRDS(object = func_rand_summary_wilcoxon.musicc,
          file = paste("simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc_cutoff_", cutoff, ".rds", sep = ""))

}



for (cutoff in cutoffs) {
  
  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]
  
  taxa_rand_output <- readRDS(file = paste("POMS_out/ref.based_POMS_sim_rand_taxa_1000reps_cutoff_", cutoff, ".rds", sep = ""))
  
  taxa_rand_wilcoxon.musicc <- readRDS(file = paste("wilcoxon.musicc_out/ref.based_wilcoxon.musicc_sim_rand_taxa_500reps_cutoff_", cutoff, ".rds", sep = ""))
  
  taxa_rand_summary_wilcoxon.musicc <- simulation_summaries(POMS_sims = taxa_rand_output[1:500],
                                                            alt_tool_sims = taxa_rand_wilcoxon.musicc,
                                                            focal_func_present = TRUE,
                                                            func_table = ref_func_cutoff,
                                                            num_cores = 20)
  
  saveRDS(object = taxa_rand_summary_wilcoxon.musicc,
          file = paste("simulation_summaries/taxa_rand_summary_POMS_wilcoxon.musicc_cutoff_", cutoff, ".rds", sep = ""))
  
}


