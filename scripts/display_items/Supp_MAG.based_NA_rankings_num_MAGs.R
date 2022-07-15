rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func.based_summary <- readRDS(file = "simulation_summaries/func.based_summary.rds")


MAG.based_POMS_0.05_NA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(is.na(func.based_summary$POMS_rank_0.05)), "num_focal_pos_mags"])
MAG.based_POMS_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(!is.na(func.based_summary$POMS_rank_0.05)), "num_focal_pos_mags"])

MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(is.na(func.based_summary$wilcoxon.musicc_rank_0.05)), "num_focal_pos_mags"])
MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(!is.na(func.based_summary$wilcoxon.musicc_rank_0.05)), "num_focal_pos_mags"])

MAG.based_regress_specificity_0.05_NA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(is.na(func.based_summary$regress_specificity_rank_0.05)), "num_focal_pos_mags"])
MAG.based_regress_specificity_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(!is.na(func.based_summary$regress_specificity_rank_0.05)), "num_focal_pos_mags"])

MAG.based_regress_sig_taxa_0.05_NA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(is.na(func.based_summary$regress_sig_taxa_rank_0.05)), "num_focal_pos_mags"])
MAG.based_regress_sig_taxa_0.05_nonNA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(!is.na(func.based_summary$regress_sig_taxa_rank_0.05)), "num_focal_pos_mags"])



MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
                                                                        geom_histogram(binwidth = 50, fill = "springgreen4") + 
                                                                        xlab("No. MAGs encoding focal gene") +
                                                                        ylab("No. replicates") +
                                                                        ggtitle("Wilcoxon test - Focal gene not significant") +
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
                                                                        ggtitle("Wilcoxon test - Focal gene significant") +
                                                                        coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                        theme_bw() +
                                                                        theme(plot.title = element_text(hjust = 0.5)) +
                                                                        annotate(geom = "text", x = 500, y = 100,
                                                                                 label = paste("n=", as.character(nrow(MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_data)), sep = ""))



MAG.based_POMS_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_POMS_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
                                                                geom_histogram(binwidth = 5, fill = "steelblue4") + 
                                                                xlab("No. MAGs encoding focal gene") +
                                                                ylab("No. replicates") +
                                                                ggtitle("POMS - Focal gene not significant") +
                                                                coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                theme_bw() +
                                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                                annotate(geom = "text", x = 500, y = 100,
                                                                         label = paste("n=", as.character(nrow(MAG.based_POMS_0.05_NA_rank_num_encoding_data)), sep = ""))


MAG.based_POMS_0.05_nonNA_rank_num_encoding_hist <- ggplot(MAG.based_POMS_0.05_nonNA_rank_num_encoding_data, aes(x = num)) + 
                                                                geom_histogram(binwidth = 5, fill = "steelblue4") + 
                                                                xlab("No. MAGs encoding focal gene") +
                                                                ylab("No. replicates") +
                                                                ggtitle("POMS - Focal gene significant") +
                                                                coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                theme_bw() +
                                                                theme(plot.title = element_text(hjust = 0.5)) +
                                                                annotate(geom = "text", x = 500, y = 100,
                                                                         label = paste("n=", as.character(nrow(MAG.based_POMS_0.05_nonNA_rank_num_encoding_data)), sep = ""))


MAG.based_regress_specificity_0.05_nonNA_rank_num_encoding_hist <- ggplot(MAG.based_regress_specificity_0.05_nonNA_rank_num_encoding_data, aes(x = num)) + 
                                                                           geom_histogram(binwidth = 5, fill = "steelblue4") + 
                                                                           xlab("No. MAGs encoding focal gene") +
                                                                           ylab("No. replicates") +
                                                                           ggtitle("Phylogenetic regression (specificity)\nFocal gene significant") +
                                                                           coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                           theme_bw() +
                                                                           theme(plot.title = element_text(hjust = 0.5)) +
                                                                           annotate(geom = "text", x = 500, y = 100,
                                                                                    label = paste("n=", as.character(nrow(MAG.based_regress_specificity_0.05_nonNA_rank_num_encoding_data)), sep = ""))


MAG.based_regress_specificity_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_regress_specificity_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
                                                                        geom_histogram(binwidth = 50, fill = "steelblue4") + 
                                                                        xlab("No. MAGs encoding focal gene") +
                                                                        ylab("No. replicates") +
                                                                        ggtitle("Phylogenetic regression (specificity)\nFocal gene not significant") +
                                                                        coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
                                                                        theme_bw() +
                                                                        theme(plot.title = element_text(hjust = 0.5)) +
                                                                        annotate(geom = "text", x = 500, y = 100,
                                                                                 label = paste("n=", as.character(nrow(MAG.based_regress_specificity_0.05_NA_rank_num_encoding_data)), sep = ""))


num_encoding_combined <- plot_grid(MAG.based_wilcoxon.musicc_0.05_NA_rank_num_encoding_hist, MAG.based_wilcoxon.musicc_0.05_nonNA_rank_num_encoding_hist,
                                   MAG.based_regress_specificity_0.05_NA_rank_num_encoding_hist, MAG.based_regress_specificity_0.05_nonNA_rank_num_encoding_hist,
                                   MAG.based_POMS_0.05_NA_rank_num_encoding_hist, MAG.based_POMS_0.05_nonNA_rank_num_encoding_hist,
                                   nrow = 3, ncol = 2)

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_NA_rankings_num_MAGs.pdf",
       plot = num_encoding_combined,
       device = "pdf",
       width = 8,
       height = 8,
       dpi = 600)


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Supp_MAG.based_NA_rankings_num_MAGs.png",
       plot = num_encoding_combined,
       device = "png",
       width = 8,
       height = 8,
       dpi = 600)
