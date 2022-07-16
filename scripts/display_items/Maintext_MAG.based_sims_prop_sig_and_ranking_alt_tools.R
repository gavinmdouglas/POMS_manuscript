rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(ggbeeswarm)
library(reshape2)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

taxa.based_summary <- readRDS(file = "simulation_summaries/taxa.based_summary.rds")
clade.based_summary <- readRDS(file = "simulation_summaries/clade.based_summary.rds")
func.based_summary <- readRDS(file = "simulation_summaries/func.based_summary.rds")

qual_col <- c("#1f78b4",
              "blue",
              "blue4",
              "firebrick1",
              "gold2",
              "darkslategray3",
              "mediumorchid",
              "palegreen")

# Now get rankings for each tool (only for Focal gene simulations).

func.based_ranks <- melt(func.based_summary[, grep("_rank_0.05$", colnames(func.based_summary))],
                         variable.name = "Tool",
                         value.name = "Rank")

func.based_ranks$Tool <- as.character(func.based_ranks$Tool)

func.based_ranks <- func.based_ranks[-grep("_rel_rank_0.05$", func.based_ranks$Tool), ]

func.based_ranks[which(func.based_ranks$Tool == "POMS_rank_0.05"), "Tool"] <- "POMS"
func.based_ranks[which(func.based_ranks$Tool == "regress_specificity_rank_0.05"), "Tool"] <- "Phylo. regress.\n(specificity)"
func.based_ranks[which(func.based_ranks$Tool == "regress_sig_taxa_rank_0.05"), "Tool"] <- "Phylo. regress.\n(sig. taxa)"
func.based_ranks[which(func.based_ranks$Tool == "wilcoxon.musicc_rank_0.05"), "Tool"] <- "Wilcoxon test\n(corrected)"
func.based_ranks[which(func.based_ranks$Tool == "wilcoxon.relab_rank_0.05"), "Tool"] <- "Wilcoxon test\n(relab.)"
func.based_ranks[which(func.based_ranks$Tool == "limma.voom_rank_0.05"), "Tool"] <- "limma-voom"
func.based_ranks[which(func.based_ranks$Tool == "aldex2_rank_0.05"), "Tool"] <- "ALDEx2"
func.based_ranks[which(func.based_ranks$Tool == "deseq2_rank_0.05"), "Tool"] <- "DESeq2"

func.based_ranks$Tool <- factor(func.based_ranks$Tool, 
                                  levels = c("POMS",
                                             "Phylo. regress.\n(specificity)",
                                             "Phylo. regress.\n(sig. taxa)",
                                             "ALDEx2",
                                             "DESeq2",
                                             "limma-voom",
                                             "Wilcoxon test\n(corrected)",
                                             "Wilcoxon test\n(relab.)"))

func.based_ranks$Fill <- rep("Differential abundance approach", nrow(func.based_ranks))

func.based_ranks[which(func.based_ranks$Tool %in% c("POMS", "Phylo. regress.\n(specificity)", "Phylo. regress.\n(sig. taxa)")), "Fill"] <- "Phylogenetic-based method"

func.based_ranks$Simulation <- "Focal gene"


func.based_ranks_boxplots <- ggplot(func.based_ranks, aes(x = Tool, y = Rank, fill = Tool)) +
                                        geom_quasirandom(color = "grey") +
                                        geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.5) +
                                        theme_bw() +
                                        facet_grid( ~ Simulation) +
                                        ylab("Focal gene ranking\n(lower is better)") +
                                        xlab("") +
                                        scale_fill_manual(values = qual_col) +
                                        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                                        theme(legend.position = "none")


# Summary of how # MAGs that encode focal gene affects results.
func.based_ranks$num_mags <- func.based_summary$num_focal_pos_mags
func.based_ranks$focal_func <- func.based_summary$focal_func

func.based_ranks$top10 <- NA
func.based_ranks[which(func.based_ranks$Rank <= 10), "top10"] <- "Focal gene ranked <= 10th"
func.based_ranks[which(func.based_ranks$Rank > 10), "top10"] <- "Focal gene ranked > 10th"

func.based_ranks_noNA <- func.based_ranks[-which(is.na(func.based_ranks$Rank)), ]

top10_numMAG_boxplots <- ggplot(data = func.based_ranks_noNA, aes(x = Tool, y = num_mags, fill = Tool)) +
        geom_quasirandom(color = "grey") +
        geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.5) +
        theme_bw() +
        scale_fill_manual(values = qual_col) +
        facet_grid(. ~ top10) +
        ylab("No. MAGs encoding focal gene") +
        xlab("") +
        theme(legend.position = "none") +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))




MAG.based_all_prop_sig <- data.frame(Tool = c(rep("POMS", 2693),

                                              rep("Phylo. regress.\n(specificity)", 2693),

                                              rep("Phylo. regress.\n(sig. taxa)", 2693),

                                              rep("Wilcoxon test\n(corrected)", 2693),

                                              rep("Wilcoxon test\n(relab.)", 2693),

                                              rep("limma-voom", 2693),

                                              rep("ALDEx2", 2693),

                                              rep("DESeq2", 2693)),

                                     Simulation =
                                             c(rep("Random taxa", 1000),
                                               rep("Clade-based", 693),
                                               rep("Focal gene", 1000),
                                               
                                               rep("Random taxa", 1000),
                                               rep("Clade-based", 693),
                                               rep("Focal gene", 1000),
                                               
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Clade-based", 693),
                                              rep("Focal gene", 1000)),
                                              
                                     prop = c(taxa.based_summary$POMS_sig_0.05,
                                              clade.based_summary$POMS_sig_0.05,
                                              func.based_summary$POMS_sig_0.05,
                                              
                                              taxa.based_summary$regress_specificity_func_prop_sig_0.05,
                                              clade.based_summary$regress_specificity_sig_0.05,
                                              func.based_summary$regress_specificity_func_prop_sig_0.05,
                                              
                                              taxa.based_summary$regress_sig_taxa_func_prop_sig_0.05,
                                              clade.based_summary$regress_sig.taxa_sig_0.05,
                                              func.based_summary$regress_sig_taxa_func_prop_sig_0.05,
                                              
                                              taxa.based_summary$wilcoxon.musicc_sig_0.05,
                                              clade.based_summary$wilcoxon.musicc_0.05,
                                              func.based_summary$wilcoxon.musicc_sig_0.05,
                                              
                                              taxa.based_summary$wilcoxon.relab_sig_0.05,
                                              clade.based_summary$wilcoxon.relab_0.05,
                                              func.based_summary$wilcoxon.relab_sig_0.05,
                                              
                                              taxa.based_summary$limma.voom_sig_0.05,
                                              clade.based_summary$limma.voom_0.05,
                                              func.based_summary$limma.voom_sig_0.05,
                                              
                                              taxa.based_summary$aldex2_sig_0.05,
                                              rep(NA, 693), #clade.based_summary$aldex2_sig_0.05,
                                              func.based_summary$aldex2_sig_0.05,
                                              
                                              taxa.based_summary$deseq2_sig_0.05,
                                              rep(NA, 693), #clade.based_summary$deseq2_sig_0.05,
                                              func.based_summary$deseq2_sig_0.05)
                                              
                                    )


MAG.based_all_prop_sig$Simulation <- factor(MAG.based_all_prop_sig$Simulation, 
                                            levels = c("Random taxa", "Clade-based", "Focal gene"))

MAG.based_all_prop_sig$Tool <- factor(MAG.based_all_prop_sig$Tool, 
                                            levels = c("POMS",
                                                       "Phylo. regress.\n(specificity)",
                                                       "Phylo. regress.\n(sig. taxa)",
                                                       "ALDEx2",
                                                       "DESeq2",
                                                       "limma-voom",
                                                       "Wilcoxon test\n(corrected)",
                                                       "Wilcoxon test\n(relab.)"))

MAG.based_all_prop_sig_boxplots <- ggplot(MAG.based_all_prop_sig, aes(x = Tool, y = prop, fill = Tool)) +
                                          geom_quasirandom(color = "grey") +
                                          geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 0.5) +
                                          facet_grid( ~ Simulation) +
                                          theme_bw() +
                                          ylim(0, 1) +
                                          ylab("Proportion sig. gene families\n(lower is better)") +
                                          xlab("") +
                                          scale_fill_manual(values = qual_col) +
                                          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                                          theme(legend.position = "none")


# No. MAGs encoding genes that were not significant.
MAG.based_POMS_0.05_NA_rank_num_encoding_data <- data.frame(num = func.based_summary[which(is.na(func.based_summary$POMS_rank_0.05)), "num_focal_pos_mags"])

MAG.based_POMS_0.05_NA_rank_num_encoding_data$Category <- "POMS - Focal gene not significant"

MAG.based_POMS_0.05_NA_rank_num_encoding_hist <- ggplot(MAG.based_POMS_0.05_NA_rank_num_encoding_data, aes(x = num)) + 
        geom_histogram(binwidth = 25, fill = "steelblue4") + 
        xlab("No. MAGs encoding focal gene") +
        ylab("No. replicates") +
        coord_cartesian(xlim = c(0, 1600), ylim = c(0, 120)) +
        theme_bw() +
        facet_grid( ~ Category) +
        theme(plot.title = element_text(hjust = 0.5)) +
        annotate(geom = "text", x = 500, y = 100,
                 label = paste("n=", as.character(nrow(MAG.based_POMS_0.05_NA_rank_num_encoding_data)), sep = ""))

NA_mag_num_panel <- plot_grid(MAG.based_POMS_0.05_NA_rank_num_encoding_hist, NULL, nrow = 2, rel_heights = c(0.9, 0.1))

focal_gene_row <- plot_grid(func.based_ranks_boxplots,
                            top10_numMAG_boxplots,
                            ncol = 2, rel_widths = c(0.35, 0.65),
                            labels = c('a', 'b'))

combined_figure <- plot_grid(focal_gene_row,
                             MAG.based_all_prop_sig_boxplots,
                             nrow = 2, ncol = 1, labels = c('', 'c'))


ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.pdf",
       plot = combined_figure,
       device = "pdf",
       width = 10,
       height = 7,
       dpi = 600)
       
       
ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.png",
       plot = combined_figure,
       device = "png",
       width = 10,
       height = 7,
       dpi = 300)



# Rank results for main text:
median(func.based_summary$POMS_rank_0.05, na.rm = TRUE)
mean(func.based_summary$POMS_rank_0.05, na.rm = TRUE)
sd(func.based_summary$POMS_rank_0.05, na.rm = TRUE)

median(func.based_summary$deseq2_rank_0.05, na.rm = TRUE)
mean(func.based_summary$deseq2_rank_0.05, na.rm = TRUE)
sd(func.based_summary$deseq2_rank_0.05, na.rm = TRUE)

median(func.based_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
mean(func.based_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
sd(func.based_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)



wilcox.test(func.based_summary$POMS_rank_0.05,
            func.based_summary$deseq2_rank_0.05)


wilcox.test(func.based_summary$POMS_rank_0.05,
            func.based_summary$wilcoxon.musicc_rank_0.05)


# Prop. sig results for main text:
mean(taxa.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)
sd(taxa.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)

mean(taxa.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)
sd(taxa.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)



mean(func.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)
sd(func.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)

mean(func.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)
sd(func.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)


# Percent increases:
100 * ((mean(taxa.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE) / mean(func.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)) / mean(taxa.based_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE))

100 * ((mean(taxa.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE) / mean(func.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)) / mean(taxa.based_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE))

