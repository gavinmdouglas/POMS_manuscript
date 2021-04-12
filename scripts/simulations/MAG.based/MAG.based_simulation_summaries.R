### Parse simulation output and get summary RDS files.

rm(list = ls(all.names = TRUE))

library(parallel)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")
source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

func_subset <- readRDS("MAG.based_prepped_func.rds")

unperturbed_rand_output <- readRDS(file = "MAG.based_POMS_sim_unperturbed_1000reps.rds")
unperturbed_rand_wilcoxon <- readRDS(file = "MAG.based_wilcoxon.musicc_wilcoxon.relab_sim_unperturbed_100reps.rds")

taxa_rand_output <- readRDS(file = "MAG.based_POMS_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.musicc <- readRDS(file = "MAG.based_wilcoxon.musicc_sim_rand_taxa_1000reps.rds")
taxa_rand_wilcoxon.relab_limma_voom <- readRDS(file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_taxa_1000reps.rds")
taxa_rand_aldex2_deseq2 <- readRDS(file = "MAG.based_aldex2_deseq2_sim_rand_taxa_10reps.rds")

func_rand_output <- readRDS(file = "MAG.based_POMS_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.musicc <- readRDS(file = "MAG.based_wilcoxon.musicc_sim_rand_func_1000reps.rds")
func_rand_wilcoxon.relab_limma_voom <- readRDS(file = "MAG.based_wilcoxon.relab_limma.voom_sim_rand_func_1000reps.rds")
func_rand_aldex2_deseq2 <- readRDS(file = "MAG.based_aldex2_deseq2_sim_rand_func_10reps.rds")

unperturbed_summary <- simulation_summaries(POMS_sims = unperturbed_rand_output[1:100],
                                            alt_tool_sims = unperturbed_rand_wilcoxon,
                                            focal_func_present = FALSE,
                                            func_table = func_subset,
                                            num_cores = 20)


taxa_rand_summary_wilcoxon.musicc <- simulation_summaries(POMS_sims = taxa_rand_output,
                                                          alt_tool_sims = taxa_rand_wilcoxon.musicc,
                                                          focal_func_present = TRUE,
                                                          func_table = func_subset,
                                                          num_cores = 20)

taxa_rand_summary_wilcoxon.relab_limma_voom <- simulation_summaries(POMS_sims = taxa_rand_output,
                                                                    alt_tool_sims = taxa_rand_wilcoxon.relab_limma_voom,
                                                                    focal_func_present = TRUE,
                                                                    func_table = func_subset,
                                                                    num_cores = 20)

taxa_rand_summary_aldex2_deseq2 <- simulation_summaries(POMS_sims = taxa_rand_output[1:10],
                                                        alt_tool_sims = taxa_rand_aldex2_deseq2,
                                                        focal_func_present = TRUE,
                                                        func_table = func_subset,
                                                        num_cores = 20)


func_rand_summary_wilcoxon.musicc <- simulation_summaries(POMS_sims = func_rand_output,
                                                          alt_tool_sims = func_rand_wilcoxon.musicc,
                                                          focal_func_present = TRUE,
                                                          func_table = func_subset,
                                                          num_cores = 20)

func_rand_summary_wilcoxon.relab_limma_voom <- simulation_summaries(POMS_sims = func_rand_output,
                                                                    alt_tool_sims = func_rand_wilcoxon.relab_limma_voom,
                                                                    focal_func_present = TRUE,
                                                                    func_table = func_subset,
                                                                    num_cores = 20)

func_rand_summary_aldex2_deseq2 <- simulation_summaries(POMS_sims = func_rand_output[1:10],
                                                        alt_tool_sims = func_rand_aldex2_deseq2,
                                                        focal_func_present = TRUE,
                                                        func_table = func_subset,
                                                        num_cores = 20)

saveRDS(object = unperturbed_summary, file = "simulation_summaries/unperturbed_summary.rds")

saveRDS(object = taxa_rand_summary_wilcoxon.musicc, file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.musicc.rds")
saveRDS(object = taxa_rand_summary_wilcoxon.relab_limma_voom, file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.relab_limma_voom.rds")
saveRDS(object = taxa_rand_summary_aldex2_deseq2, file = "simulation_summaries/taxa_rand_summary_POMS_aldex2_deseq2.rds")

saveRDS(object = func_rand_summary_wilcoxon.musicc, file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")
saveRDS(object = func_rand_summary_wilcoxon.relab_limma_voom, file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.relab_limma_voom.rds")
saveRDS(object = func_rand_summary_aldex2_deseq2, file = "simulation_summaries/func_rand_summary_POMS_aldex2_deseq2.rds")
