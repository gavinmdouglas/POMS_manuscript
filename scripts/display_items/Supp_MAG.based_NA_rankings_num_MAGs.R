rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")


MAG.based_POMS_0.05_NA_rank_num_encoding_data <- data.frame(num = func_rand_summary_POMS_wilcoxon.musicc[which(is.na(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05)), "num_focal_pos_mags"])
MAG.based_POMS_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func_rand_summary_POMS_wilcoxon.musicc[which(!is.na(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.05)), "num_focal_pos_mags"])

MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_data <- data.frame(num = func_rand_summary_POMS_wilcoxon.musicc[which(is.na(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05)), "num_focal_pos_mags"])
MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func_rand_summary_POMS_wilcoxon.musicc[which(!is.na(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.05)), "num_focal_pos_mags"])




MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
                                                                        geom_histogram(binwidth = 50, fill = "springgreen4") + 
                                                                        xlab("No. MAGs encoding focal gene") +
                                                                        ylab("No. replicates") +
                                                                        ggtitle("Wilcoxon test - focal gene not significant") +
                                                                        theme(plot.title = element_text(hjust = 0.5)) +
                                                                        coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                        theme_bw() +
                                                                        theme(plot.title = element_text(hjust = 0.5)) +
                                                                        annotate(geom = "text", x = 500, y = 100,
                                                                                 label = paste("n=", as.character(nrow(MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_data)), sep = ""))
                                                                        
MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_hist <- ggplot(MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_data, aes(x = num)) + 
                                                                        geom_histogram(binwidth = 5, fill = "springgreen4") + 
                                                                        xlab("No. MAGs encoding focal gene") +
                                                                        ylab("No. replicates") +
                                                                        ggtitle("Wilcoxon test - focal gene significant") +
                                                                        coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                        theme_bw() +
                                                                        theme(plot.title = element_text(hjust = 0.5)) +
                                                                        annotate(geom = "text", x = 500, y = 100,
                                                                                 label = paste("n=", as.character(nrow(MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_data)), sep = ""))



MAG.based_POMS_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_POMS_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
                                                                geom_histogram(binwidth = 5, fill = "steelblue4") + 
                                                                xlab("No. MAGs encoding focal gene") +
                                                                ylab("No. replicates") +
                                                                ggtitle("POMS - focal gene not significant") +
                                                                coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                theme_bw() +
                                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                                annotate(geom = "text", x = 500, y = 100,
                                                                         label = paste("n=", as.character(nrow(MAG.based_POMS_0.05_NA_rank_num_encoding_data)), sep = ""))


MAG.based_POMS_0.05_nonNA_rank_num_encoding_hist <- ggplot(MAG.based_POMS_0.05_nonNA_rank_num_encoding_data, aes(x = num)) + 
                                                                geom_histogram(binwidth = 5, fill = "steelblue4") + 
                                                                xlab("No. MAGs encoding focal gene") +
                                                                ylab("No. replicates") +
                                                                ggtitle("POMS - focal gene significant") +
                                                                coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                theme_bw() +
                                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                                annotate(geom = "text", x = 500, y = 100,
                                                                         label = paste("n=", as.character(nrow(MAG.based_POMS_0.05_nonNA_rank_num_encoding_data)), sep = ""))



MAG.based_POMS_vs_wilcoxon.musicc_NA_vs_nonNA_ranking_num_MAGs_combined <- plot_grid(MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_hist,
                                                                                     MAG.based_POMS_0.05_nonNA_rank_num_encoding_hist,
                                                                                     MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_hist,
                                                                                     MAG.based_POMS_0.05_NA_rank_num_encoding_hist,
                                                                                     nrow = 2,
                                                                                     ncol = 2,
                                                                                     labels = c('a', 'b', 'c', 'd'))

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_NA_rankings_num_MAGs.pdf",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_NA_vs_nonNA_ranking_num_MAGs_combined,
       device = "pdf",
       width = 8,
       height = 6,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_NA_rankings_num_MAGs.png",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_NA_vs_nonNA_ranking_num_MAGs_combined,
       device = "png",
       width = 8,
       height = 6,
       dpi = 600)
