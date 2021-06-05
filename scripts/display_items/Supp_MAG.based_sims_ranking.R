rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")

MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG <- ggplot(func_rand_summary_POMS_wilcoxon.musicc,
                                                         aes(y = num_focal_pos_mags,
                                                             x = wilcoxon.musicc_rank_0.05)) +
                                                         geom_point(colour = "springgreen4", size = 2) +
                                                         theme_bw() +
                                                         xlab("Focal gene ranking") +
                                                         ylab("No. MAGs encoding focal gene") +
                                                         theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                                         ggtitle("Wilcoxon test") +
                                                         ylim(0, 1600) +
                                                         xlim(0, 4000)

MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG_marginal <- ggMarginal(MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG, type = "histogram", size = 10)

MAG_sim_func.based_ranking_POMS_vs_num.MAG <- ggplot(func_rand_summary_POMS_wilcoxon.musicc,
                                                         aes(y = num_focal_pos_mags,
                                                             x = POMS_rank_0.05)) +
                                                         geom_point(colour = "steelblue4", size = 2) +
                                                         theme_bw() +
                                                         xlab("Focal gene ranking") +
                                                         ylab("No. MAGs encoding focal gene") +
                                                         theme(plot.title = element_text(hjust = 0.5, vjust = -1)) +
                                                         ggtitle("POMS") +
                                                         ylim(0, 1600) +
                                                         xlim(0, 4000)

MAG_sim_func.based_ranking_POMS_vs_num.MAG_marginal <- ggMarginal(MAG_sim_func.based_ranking_POMS_vs_num.MAG, type = "histogram", size = 10)


MAG.based_POMS_vs_wilcoxon.musicc_ranking_marginal_combined <- plot_grid(MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG_marginal,
                                                                         MAG_sim_func.based_ranking_POMS_vs_num.MAG_marginal,
                                                                          nrow = 1, ncol = 2, labels = c('a', 'b'))

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking.pdf",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_ranking_marginal_combined,
       device = "pdf",
       width = 8,
       height = 4,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking.png",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_ranking_marginal_combined,
       device = "png",
       width = 8,
       height = 4,
       dpi = 300)


# Summary statistics
median(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
mean(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
sd(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rel_rank_0.05, na.rm = TRUE)

median(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05, na.rm = TRUE)
mean(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05, na.rm = TRUE)
sd(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05, na.rm = TRUE)

# Percent in top 10
length(which(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05 < 10)) / length(which(!is.na(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05))) * 100
length(which(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05 < 10)) / length(which(!is.na(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05))) * 100


# No. NAs
length(which(is.na(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05)))
length(which(is.na(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05)))


# Correlation between number focal positive MAGs and number of significant nodes.
cor.test(func_rand_summary_POMS_wilcoxon.musicc$num_focal_pos_mags, func_rand_summary_POMS_wilcoxon.musicc$POMS_num_sig_nodes, method = "spearman")
