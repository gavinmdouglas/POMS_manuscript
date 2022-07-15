rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggbeeswarm)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func.based_summary <- readRDS(file = "simulation_summaries/func.based_summary.rds")

ranking_POMS_vs_num.MAG <- ggplot(data = func.based_summary,
                                             aes(y = num_focal_pos_mags,
                                                 x = POMS_rank_0.05)) +
                                                        geom_point(colour = "#1f78b4", size = 2) +
                                                        theme_bw() +
                                                        xlab("") +
                                                        ylab("No. MAGs\nencoding\nfocal gene") +
                                                        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
                                                              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
                                                        ggtitle("POMS") +
                                                        ylim(0, 1600) +
                                                        xlim(0, 5000)

ranking_POMS_vs_num.MAG_marginal <- ggMarginal(ranking_POMS_vs_num.MAG, type = "histogram", size = 10)


ranking_regress_specificity_vs_num.MAG <- ggplot(data = func.based_summary,
                                                  aes(y = num_focal_pos_mags,
                                                      x = regress_specificity_rank_0.05)) +
                        geom_point(colour = "blue", size = 2) +
                        theme_bw() +
                        xlab("") +
                        ylab("") +
                        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
                              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
                        ggtitle("Phylo. regress. (specificity)") +
                        ylim(0, 1600) +
                        xlim(0, 5000)

ranking_regress_specificity_vs_num.MAG_marginal <- ggMarginal(ranking_regress_specificity_vs_num.MAG, type = "histogram", size = 10)


ranking_regress_sig_taxa_vs_num.MAG <- ggplot(data = func.based_summary,
                                                 aes(y = num_focal_pos_mags,
                                                     x = regress_sig_taxa_rank_0.05)) +
        geom_point(colour = "blue4", size = 2) +
        theme_bw() +
        xlab("") +
        ylab("No. MAGs\nencoding\nfocal gene") +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ggtitle("Phylo. regress. (sig. taxa)") +
        ylim(0, 1600) +
        xlim(0, 5000)

ranking_regress_sig_taxa_vs_num.MAG_marginal <- ggMarginal(ranking_regress_sig_taxa_vs_num.MAG, type = "histogram", size = 10)


ranking_aldex2_vs_num.MAG <- ggplot(data = func.based_summary,
                                              aes(y = num_focal_pos_mags,
                                                  x = aldex2_rank_0.05)) +
        geom_point(colour = "firebrick1", size = 2) +
        theme_bw() +
        xlab("") +
        ylab("") +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ggtitle("ALDEx2") +
        ylim(0, 1600) +
        xlim(0, 5000)

ranking_aldex2_vs_num.MAG_marginal <- ggMarginal(ranking_aldex2_vs_num.MAG, type = "histogram", size = 10)

ranking_deseq2_vs_num.MAG <- ggplot(data = func.based_summary,
                                    aes(y = num_focal_pos_mags,
                                        x = deseq2_rank_0.05)) +
        geom_point(colour = "gold2", size = 2) +
        theme_bw() +
        xlab("") +
        ylab("No. MAGs\nencoding\nfocal gene") +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ggtitle("DESeq2") +
        ylim(0, 1600) +
        xlim(0, 5000)

ranking_deseq2_vs_num.MAG_marginal <- ggMarginal(ranking_deseq2_vs_num.MAG, type = "histogram", size = 10)



ranking_limma.voom_vs_num.MAG <- ggplot(data = func.based_summary,
                                    aes(y = num_focal_pos_mags,
                                        x = limma.voom_rank_0.05)) +
        geom_point(colour = "darkslategray3", size = 2) +
        theme_bw() +
        xlab("") +
        ylab("") +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ggtitle("limma-voom") +
        ylim(0, 1600) +
        xlim(0, 5000)

ranking_limma.voom_vs_num.MAG_marginal <- ggMarginal(ranking_limma.voom_vs_num.MAG, type = "histogram", size = 10)


ranking_wilcoxon.musicc_vs_num.MAG <- ggplot(data = func.based_summary,
                                                         aes(y = num_focal_pos_mags,
                                                             x = wilcoxon.musicc_rank_0.05)) +
                                                         geom_point(colour = "mediumorchid", size = 2) +
                                                         theme_bw() +
                                                         xlab("Focal gene ranking") +
                                                         ylab("No. MAGs\nencoding\nfocal gene") +
                                                         theme(plot.title = element_text(hjust = 0.5, vjust = -1),
                                                               axis.title.y = element_text(angle = 0, vjust = 0.5)) +
                                                         ggtitle("Wilcoxon test (corrected)") +
                                                         ylim(0, 1600) +
                                                         xlim(0, 5000)

ranking_wilcoxon.musicc_vs_num.MAG_marginal <- ggMarginal(ranking_wilcoxon.musicc_vs_num.MAG, type = "histogram", size = 10)


ranking_wilcoxon.relab_vs_num.MAG <- ggplot(data = func.based_summary,
                                             aes(y = num_focal_pos_mags,
                                                 x = wilcoxon.relab_rank_0.05)) +
        geom_point(colour = "palegreen", size = 2) +
        theme_bw() +
        xlab("Focal gene ranking") +
        ylab("") +
        theme(plot.title = element_text(hjust = 0.5, vjust = -1),
              axis.title.y = element_text(angle = 0, vjust = 0.5)) +
        ggtitle("Wilcoxon test (relab.)") +
        ylim(0, 1600) +
        xlim(0, 5000)

ranking_wilcoxon.relab_vs_num.MAG_marginal <- ggMarginal(ranking_wilcoxon.relab_vs_num.MAG, type = "histogram", size = 10)





ranking_vs_num_MAG_combined <- plot_grid(ranking_POMS_vs_num.MAG, ranking_regress_specificity_vs_num.MAG,
                                         ranking_regress_sig_taxa_vs_num.MAG, ranking_aldex2_vs_num.MAG, 
                                         ranking_deseq2_vs_num.MAG, ranking_limma.voom_vs_num.MAG,
                                         ranking_wilcoxon.musicc_vs_num.MAG, ranking_wilcoxon.relab_vs_num.MAG,
                                         nrow = 4, ncol = 2, rel_widths = c(1.2, 1), labels = c('', '', '', '', '', '', '', ''))


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking.pdf",
       plot = ranking_vs_num_MAG_combined,
       device = "pdf",
       width = 7,
       height = 8,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_sims_ranking.png",
       plot = ranking_vs_num_MAG_combined,
       device = "png",
       width = 7,
       height = 8,
       dpi = 300)


# Summary statistics
median(func_rand_summary_POMS_wilcoxon.musicc[which(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func_rand_summary_POMS_wilcoxon.musicc[which(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05 > 10), "num_focal_pos_mags"])

median(func_rand_summary_POMS_wilcoxon.musicc[which(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func_rand_summary_POMS_wilcoxon.musicc[which(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05 > 10), "num_focal_pos_mags"])

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
