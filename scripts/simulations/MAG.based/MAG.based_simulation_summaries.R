### Parse simulation output and get summary RDS files.

rm(list = ls(all.names = TRUE))

library(parallel)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

func_subset <- readRDS("MAG.based_prepped_func.rds")

unperturbed_POMS_output <- readRDS(file = "POMS_out/MAG.based_POMS_sim_unperturbed_1000reps.rds")
unperturbed_regress_output <- readRDS(file = "regression_out/MAG.based_regression_sim_unperturbed_1000reps.rds")
unperturbed_alt.tools_output <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.musicc_wilcoxon.relab_limma.voom_unperturbed_1000reps.rds")

for (rep_i in 1:1000) {
  rep_rds_in <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_unperturbed_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  unperturbed_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  unperturbed_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}



taxa_rand_POMS_output <- readRDS(file = "POMS_out/MAG.based_POMS_sim_rand_taxa_1000reps.rds")
taxa_rand_regress_output <- readRDS(file = "regression_out/MAG.based_regression_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.musicc <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.musicc_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.relab_limma_voom <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.relab_limma.voom_sim_rand_taxa_1000reps.rds")

taxa_rand_alt.tools_output <- list()

for (rep_i in 1:1000) {
  
  taxa_rand_regress_output[[rep_i]]$func <- taxa_rand_POMS_output[[rep_i]]$func
  
  taxa_rand_alt.tools_output[[rep_i]] <- list()
  
  taxa_rand_alt.tools_output[[rep_i]]$func <- taxa_rand_wilcoxon.musicc[[rep_i]]$func
  
  taxa_rand_alt.tools_output[[rep_i]]$wilcoxon.musicc <- taxa_rand_wilcoxon.musicc[[rep_i]]$wilcoxon.musicc
  taxa_rand_alt.tools_output[[rep_i]]$wilcoxon.relab <- taxa_rand_wilcoxon.relab_limma_voom[[rep_i]]$wilcoxon.relab
  taxa_rand_alt.tools_output[[rep_i]]$limma.voom <- taxa_rand_wilcoxon.relab_limma_voom[[rep_i]]$limma.voom

  
  rep_rds_in <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_taxa_sim_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  taxa_rand_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  taxa_rand_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}



func_rand_POMS_output <- readRDS(file = "POMS_out/MAG.based_POMS_sim_rand_func_1000reps.rds")
func_rand_regress_output <- readRDS(file = "regression_out/MAG.based_regression_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.musicc <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.musicc_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.relab_limma_voom <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.relab_limma.voom_sim_rand_func_1000reps.rds")

func_rand_alt.tools_output <- list()

for (rep_i in 1:1000) {
  
  func_rand_regress_output[[rep_i]]$func <- func_rand_POMS_output[[rep_i]]$func
  
  func_rand_alt.tools_output[[rep_i]] <- list()
  
  func_rand_alt.tools_output[[rep_i]]$func <- func_rand_wilcoxon.musicc[[rep_i]]$func
  
  func_rand_alt.tools_output[[rep_i]]$wilcoxon.musicc <- func_rand_wilcoxon.musicc[[rep_i]]$wilcoxon.musicc
  func_rand_alt.tools_output[[rep_i]]$wilcoxon.relab <- func_rand_wilcoxon.relab_limma_voom[[rep_i]]$wilcoxon.relab
  func_rand_alt.tools_output[[rep_i]]$limma.voom <- func_rand_wilcoxon.relab_limma_voom[[rep_i]]$limma.voom
  
  
  rep_rds_in <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_func_sim_aldex2_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  func_rand_alt.tools_output[[rep_i]]$aldex2 <- rep_rds_in$aldex2
  func_rand_alt.tools_output[[rep_i]]$deseq2 <- rep_rds_in$deseq2 
}


clade.based_rand_POMS_output <- readRDS(file = "POMS_out/MAG.based_POMS_sim_rand_clade.based_693reps.rds")
clade.based_rand_regress_output <- readRDS(file = "regression_out/MAG.based_regression_sim_rand_clade.based_693reps.rds")
clade.based_rand_wilcoxon.relab_wilcoxon.musicc_limma_voom <- readRDS(file = "DA_tool_out/MAG.based_wilcoxon.relab_wilcoxon.musicc_limma.voom_sim_rand_clade.based_693reps.rds")

clade.based_rand_alt.tools_output <- list()

for (rep_i in 1:693) {
  
  clade.based_rand_alt.tools_output[[rep_i]] <- list()
  
  clade.based_rand_alt.tools_output[[rep_i]]$wilcoxon.musicc <- clade.based_rand_wilcoxon.relab_wilcoxon.musicc_limma_voom[[rep_i]]$wilcoxon.musicc
  clade.based_rand_alt.tools_output[[rep_i]]$wilcoxon.relab <- clade.based_rand_wilcoxon.relab_wilcoxon.musicc_limma_voom[[rep_i]]$wilcoxon.relab
  clade.based_rand_alt.tools_output[[rep_i]]$limma.voom <- clade.based_rand_wilcoxon.relab_wilcoxon.musicc_limma_voom[[rep_i]]$limma.voom
  
  
  # aldex2_rep_rds_in <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_clade.based_sim_aldex2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  # deseq2_rep_rds_in <- readRDS(paste("DA_tool_out/aldex2_deseq2/MAG.based_clade.based_sim_deseq2_rds/rep", as.character(rep_i), ".rds", sep = ""))
  # clade.based_rand_alt.tools_output[[rep_i]]$aldex2 <- aldex2_rep_rds_in$aldex2
  # clade.based_rand_alt.tools_output[[rep_i]]$deseq2 <- deseq2_rep_rds_in$deseq2 
}


unperturbed_summary <- simulation_summaries(POMS_sims = unperturbed_POMS_output,
                                            regress_sims = unperturbed_regress_output,
                                            alt_tool_sims = unperturbed_alt.tools_output,
                                            focal_func_present = FALSE,
                                            func_table = func_subset,
                                            num_cores = 20,
                                            sig_cutoffs = c(0.05))


taxa.based_summary <- simulation_summaries(POMS_sims = taxa_rand_POMS_output,
                                           regress_sims = taxa_rand_regress_output,
                                          alt_tool_sims = taxa_rand_alt.tools_output,
                                          focal_func_present = TRUE,
                                          func_table = func_subset,
                                          num_cores = 20,
                                          sig_cutoffs = c(0.05))



func.based_summary <- simulation_summaries(POMS_sims = func_rand_POMS_output,
                                           regress_sims = func_rand_regress_output,
                                           alt_tool_sims = func_rand_alt.tools_output,
                                           focal_func_present = TRUE,
                                           func_table = func_subset,
                                           num_cores = 20,
                                           sig_cutoffs = c(0.05))

clade.based_summary <- simulation_summaries(POMS_sims = clade.based_rand_POMS_output,
                                           regress_sims = clade.based_rand_regress_output,
                                           alt_tool_sims = clade.based_rand_alt.tools_output,
                                           focal_func_present = FALSE,
                                           func_table = func_subset,
                                           num_cores = 20,
                                           sig_cutoffs = c(0.05))

saveRDS(object = unperturbed_summary,
        file = "simulation_summaries/unperturbed_summary.rds")

saveRDS(object = taxa.based_summary,
        file = "simulation_summaries/taxa.based_summary.rds")

saveRDS(object = func.based_summary,
        file = "simulation_summaries/func.based_summary.rds")

saveRDS(object = clade.based_summary,
        file = "simulation_summaries/clade.based_summary.rds")

