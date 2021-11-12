### Parse simulation output and get summary RDS files.

rm(list = ls(all.names = TRUE))

library(parallel)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

func_subset <- readRDS("MAG.based_prepped_func.rds")

unperturbed_POMS_output <- readRDS(file = "MAG.based_POMS_sim_unperturbed_1000reps.rds")
unperturbed_alt.tools_output <- readRDS(file = "MAG.based_wilcoxon.musicc_wilcoxon.relab_limma.voom_unperturbed_1000reps.rds")

for (rep_i in 1:1000) {
  rep_rds_in <- readRDS(paste("MAG.based_unperturbed_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  unperturbed_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  unperturbed_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}


taxa_rand_output <- readRDS(file = "MAG.based_POMS_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.musicc <- readRDS(file = "MAG.based_wilcoxon.musicc_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.relab_limma_voom <- readRDS(file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_taxa_1000reps.rds")

taxa_rand_alt.tools_output <- list()

for (rep_i in 1:1000) {
  
  taxa_rand_alt.tools_output[[rep_i]] <- list()
  
  taxa_rand_alt.tools_output[[rep_i]]$func <- taxa_rand_wilcoxon.musicc[[rep_i]]$func
  
  taxa_rand_alt.tools_output[[rep_i]]$wilcoxon.musicc <- taxa_rand_wilcoxon.musicc[[rep_i]]$wilcoxon.musicc
  taxa_rand_alt.tools_output[[rep_i]]$wilcoxon.relab <- taxa_rand_wilcoxon.relab_limma_voom[[rep_i]]$wilcoxon.relab
  taxa_rand_alt.tools_output[[rep_i]]$limma.voom <- taxa_rand_wilcoxon.relab_limma_voom[[rep_i]]$limma.voom

  
  rep_rds_in <- readRDS(paste("MAG.based_taxa_sim_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  taxa_rand_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  taxa_rand_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}



func_rand_output <- readRDS(file = "MAG.based_POMS_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.musicc <- readRDS(file = "MAG.based_wilcoxon.musicc_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.relab_limma_voom <- readRDS(file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_func_1000reps.rds")

func_rand_alt.tools_output <- list()

for (rep_i in 1:1000) {
  
  func_rand_alt.tools_output[[rep_i]] <- list()
  
  func_rand_alt.tools_output[[rep_i]]$func <- func_rand_wilcoxon.musicc[[rep_i]]$func
  
  func_rand_alt.tools_output[[rep_i]]$wilcoxon.musicc <- func_rand_wilcoxon.musicc[[rep_i]]$wilcoxon.musicc
  func_rand_alt.tools_output[[rep_i]]$wilcoxon.relab <- func_rand_wilcoxon.relab_limma_voom[[rep_i]]$wilcoxon.relab
  func_rand_alt.tools_output[[rep_i]]$limma.voom <- func_rand_wilcoxon.relab_limma_voom[[rep_i]]$limma.voom
  
  
  rep_rds_in <- readRDS(paste("MAG.based_func_sim_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  func_rand_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  func_rand_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}



unperturbed_summary_vs_alt.tools <- simulation_summaries(POMS_sims = unperturbed_POMS_output,
                                                         alt_tool_sims = unperturbed_alt.tools_output,
                                                                        focal_func_present = FALSE,
                                                                        func_table = func_subset,
                                                                        num_cores = 20)


taxa_rand_summary_vs_alt.tools <- simulation_summaries(POMS_sims = taxa_rand_output,
                                                       alt_tool_sims = taxa_rand_alt.tools_output,
                                                       focal_func_present = TRUE,
                                                       func_table = func_subset,
                                                       num_cores = 20)



func_rand_summary_vs_alt.tools <- simulation_summaries(POMS_sims = func_rand_output,
                                                       alt_tool_sims = func_rand_alt.tools_output,
                                                       focal_func_present = TRUE,
                                                       func_table = func_subset,
                                                       num_cores = 20)


saveRDS(object = unperturbed_summary_vs_alt.tools,
        file = "simulation_summaries/unperturbed_summary.rds")

saveRDS(object = taxa_rand_summary_vs_alt.tools,
        file = "simulation_summaries/taxa_rand_summary.rds")

saveRDS(object = func_rand_summary_vs_alt.tools,
        file = "simulation_summaries/func_rand_summary.rds")

