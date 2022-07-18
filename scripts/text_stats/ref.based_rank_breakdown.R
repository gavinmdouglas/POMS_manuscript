rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)

setwd("/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

nu0.50_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.50.rds")
nu0.65_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.65.rds")
nu0.80_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.80.rds")
nu0.95_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.95.rds")

median(nu0.50_sim_func_summary$POMS_rank_0.05, na.rm = TRUE)
median(nu0.50_sim_func_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)

median(nu0.65_sim_func_summary$POMS_rank_0.05, na.rm = TRUE)
median(nu0.65_sim_func_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)

median(nu0.80_sim_func_summary$POMS_rank_0.05, na.rm = TRUE)
median(nu0.80_sim_func_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)

median(nu0.95_sim_func_summary$POMS_rank_0.05, na.rm = TRUE)
median(nu0.95_sim_func_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)


cor.test(nu0.95_sim_func_summary$wilcoxon.musicc_rank_0.05, nu0.95_sim_func_summary$num_focal_pos_mags)
