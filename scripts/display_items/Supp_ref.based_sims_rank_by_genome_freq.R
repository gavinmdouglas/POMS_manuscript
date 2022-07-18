rm(list = ls(all.names = TRUE))

library(ggplot2)
library(reshape2)

setwd("/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

nu0.50_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.50.rds")
nu0.65_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.65.rds")
nu0.80_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.80.rds")
nu0.95_sim_func_summary <- readRDS(file = "simulation_summaries/func.based_cutoff_nu0.95.rds")

nu0.50_sim_func_summary_ranks <- nu0.50_sim_func_summary[, c("num_focal_pos_mags", 
                                                             "POMS_rank_0.05",
                                                             "regress_specificity_rank_0.05",
                                                             "regress_sig_taxa_rank_0.05",
                                                             "wilcoxon.musicc_rank_0.05")]
nu0.50_sim_func_summary_ranks$cutoff <- "nu=0.50"
nu0.50_sim_func_summary_ranks_melt <- melt(nu0.50_sim_func_summary_ranks, id.vars = c("cutoff", "num_focal_pos_mags"))


nu0.65_sim_func_summary_ranks <- nu0.65_sim_func_summary[, c("num_focal_pos_mags", 
                                                             "POMS_rank_0.05",
                                                             "regress_specificity_rank_0.05",
                                                             "regress_sig_taxa_rank_0.05",
                                                             "wilcoxon.musicc_rank_0.05")]
nu0.65_sim_func_summary_ranks$cutoff <- "nu=0.65"
nu0.65_sim_func_summary_ranks_melt <- melt(nu0.65_sim_func_summary_ranks, id.vars = c("cutoff", "num_focal_pos_mags"))


nu0.80_sim_func_summary_ranks <- nu0.80_sim_func_summary[, c("num_focal_pos_mags", 
                                                             "POMS_rank_0.05",
                                                             "regress_specificity_rank_0.05",
                                                             "regress_sig_taxa_rank_0.05",
                                                             "wilcoxon.musicc_rank_0.05")]
nu0.80_sim_func_summary_ranks$cutoff <- "nu=0.80"
nu0.80_sim_func_summary_ranks_melt <- melt(nu0.80_sim_func_summary_ranks, id.vars = c("cutoff", "num_focal_pos_mags"))


nu0.95_sim_func_summary_ranks <- nu0.95_sim_func_summary[, c("num_focal_pos_mags", 
                                                             "POMS_rank_0.05",
                                                             "regress_specificity_rank_0.05",
                                                             "regress_sig_taxa_rank_0.05",
                                                             "wilcoxon.musicc_rank_0.05")]
nu0.95_sim_func_summary_ranks$cutoff <- "nu=0.95"
nu0.95_sim_func_summary_ranks_melt <- melt(nu0.95_sim_func_summary_ranks, id.vars = c("cutoff", "num_focal_pos_mags"))


sim_func_ranks <- do.call(rbind, list(nu0.50_sim_func_summary_ranks_melt, nu0.65_sim_func_summary_ranks_melt,
                                      nu0.80_sim_func_summary_ranks_melt, nu0.95_sim_func_summary_ranks_melt))

sim_func_ranks$variable <- as.character(sim_func_ranks$variable)

sim_func_ranks[which(sim_func_ranks$variable == "POMS_rank_0.05"), "variable"] <- "POMS"
sim_func_ranks[which(sim_func_ranks$variable == "regress_specificity_rank_0.05"), "variable"] <- "Phylo. regress.\n(specificity)"
sim_func_ranks[which(sim_func_ranks$variable == "regress_sig_taxa_rank_0.05"), "variable"] <- "Phylo. regress.\n(sig. taxa)"
sim_func_ranks[which(sim_func_ranks$variable == "wilcoxon.musicc_rank_0.05"), "variable"] <- "Wilcoxon test\n(corrected)"

sim_func_ranks$variable <- factor(sim_func_ranks$variable, 
                          levels = c("POMS",
                                     "Phylo. regress.\n(specificity)",
                                     "Phylo. regress.\n(sig. taxa)",
                                     "Wilcoxon test\n(corrected)"))

ref.based_ranks_and_genome_nums <- ggplot(sim_func_ranks, aes(y = num_focal_pos_mags, x = value, colour = variable)) +
                                                          geom_point(size = 2) +
                                                          theme_bw() +
                                                          xlab("") +
                                                          scale_colour_manual(name = "Method", values = c("#1f78b4", "blue", "blue4", "mediumorchid")) +
                                                          ylab("No. genomes encoding focal gene") +
                                                          xlab("Rank of focal gene in output") +
                                                          facet_wrap(cutoff ~ variable) +
                                                          xlim(0, 2600) + ylim(0, 500)

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_ref.based_ranks_and_genome_nums.pdf",
       plot = ref.based_ranks_and_genome_nums,
       device = "pdf",
       width = 12,
       height = 8,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_ref.based_ranks_and_genome_nums.png",
       plot = ref.based_ranks_and_genome_nums,
       device = "png",
       width = 12,
       height = 8,
       dpi = 300)

