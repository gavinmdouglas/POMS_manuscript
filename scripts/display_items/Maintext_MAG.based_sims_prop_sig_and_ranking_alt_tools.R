rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(ggbeeswarm)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

taxa_rand_summary <- readRDS(file = "simulation_summaries/taxa_rand_summary.rds")

func_rand_summary <- readRDS(file = "simulation_summaries/func_rand_summary.rds")

qual_col <- c("#1f78b4",
              "firebrick1",
              "gold2",
              "darkslategray3",
              "mediumorchid",
              "palegreen")

# Now get rankings for each tool (only for Focal gene simulations).
MAG.based_all_rank <- data.frame(Tool = c(rep("POMS", 1000),
                                          
                                          rep("Wilcoxon test\n(corrected)", 1000),
                                          
                                          rep("Wilcoxon test\n(relab.)", 1000),
                                          
                                          rep("limma-voom", 1000),
                                          
                                          rep("ALDEx2", 1000),
                                          
                                          rep("DESeq2", 1000)),
                                 
                                 rank = c(func_rand_summary$POMS_rank_0.05,
                                          
                                          func_rand_summary$wilcoxon.musicc_rank_0.05,
                                          
                                          func_rand_summary$wilcoxon.relab_rank_0.05,
                                          
                                          func_rand_summary$limma.voom_rank_0.05,
                                          
                                          func_rand_summary$aldex2_rank_0.05,
                                          
                                          func_rand_summary$deseq2_rank_0.05))


MAG.based_all_rank$Tool <- factor(MAG.based_all_rank$Tool, 
                                  levels = c("POMS", "ALDEx2", "DESeq2", "limma-voom", "Wilcoxon test\n(corrected)", "Wilcoxon test\n(relab.)"))

MAG.based_all_rank$Fill <- rep("Differential abundance approach", nrow(MAG.based_all_rank))
MAG.based_all_rank[which(MAG.based_all_rank$Tool == "POMS"), "Fill"] <- "Balance tree-based approach"




MAG.based_all_rank_boxplots <- ggplot(MAG.based_all_rank, aes(x = Tool, y = rank, fill = Tool)) +
                                        geom_quasirandom(color = "grey") +
                                        geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 1) +
                                        theme_bw() +
                                        ylab("Focal gene ranking\n(lower is better)") +
                                        xlab("") +
                                        scale_fill_manual(values = qual_col) +
                                        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                                        theme(legend.position = "none") 


MAG.based_all_prop_sig <- data.frame(Tool = c(rep("POMS", 2000),
                                              
                                              rep("Wilcoxon test\n(corrected)", 2000),
                                              
                                              rep("Wilcoxon test\n(relab.)", 2000),
                                              
                                              rep("limma-voom", 2000),
                                              
                                              rep("ALDEx2", 2000),
                                              
                                              rep("DESeq2", 2000)),
                                     
                                     Simulation = 
                                              c(rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000)),
                                              
                                     prop = c(taxa_rand_summary$POMS_sig_0.05,
                                              func_rand_summary$POMS_sig_0.05,
                                              
                                              taxa_rand_summary$wilcoxon.musicc_sig_0.05,
                                              func_rand_summary$wilcoxon.musicc_sig_0.05,
                                              
                                              taxa_rand_summary$wilcoxon.relab_sig_0.05,
                                              func_rand_summary$wilcoxon.relab_sig_0.05,
                                              
                                              taxa_rand_summary$limma.voom_sig_0.05,
                                              func_rand_summary$limma.voom_sig_0.05,
                                              
                                              taxa_rand_summary$aldex2_sig_0.05,
                                              func_rand_summary$aldex2_sig_0.05,
                                              
                                              taxa_rand_summary$deseq2_sig_0.05,
                                              func_rand_summary$deseq2_sig_0.05)
                                              
                                    )


MAG.based_all_prop_sig$Simulation <- factor(MAG.based_all_prop_sig$Simulation, 
                                            levels = c("Random taxa", "Focal gene"))

MAG.based_all_prop_sig$Tool <- factor(MAG.based_all_prop_sig$Tool, 
                                            levels = c("POMS", "ALDEx2", "DESeq2", "limma-voom", "Wilcoxon test\n(corrected)", "Wilcoxon test\n(relab.)"))

MAG.based_all_prop_sig_boxplots <- ggplot(MAG.based_all_prop_sig, aes(x = Tool, y = prop, fill = Tool)) +
                                          geom_quasirandom(color = "grey") +
                                          geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 1) +
                                          facet_grid( ~ Simulation) +
                                          theme_bw() +
                                          ylim(0, 1) +
                                          ylab("Proportion sig. gene families\n(lower is better)") +
                                          xlab("") +
                                          scale_fill_manual(values = qual_col) +
                                          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                                          theme(legend.position = "none")

MAG.based_alt_tool_comparisons <- plot_grid(MAG.based_all_rank_boxplots,
                                            MAG.based_all_prop_sig_boxplots,
                                            nrow = 1, ncol = 2, labels = c('a', 'b'), rel_widths = c(1, 1.75))

ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.pdf",
       plot = MAG.based_alt_tool_comparisons,
       device = "pdf",
       width = 10,
       height = 5,
       dpi = 600)
       
       
ggsave(filename = "~/github_repos/POMS_manuscript/display_items/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.png",
       plot = MAG.based_alt_tool_comparisons,
       device = "png",
       width = 10,
       height = 5,
       dpi = 300)


# Rank results for main text:
median(func_rand_summary$POMS_rank_0.05, na.rm = TRUE)
mean(func_rand_summary$POMS_rank_0.05, na.rm = TRUE)
sd(func_rand_summary$POMS_rank_0.05, na.rm = TRUE)

median(func_rand_summary$deseq2_rank_0.05, na.rm = TRUE)
mean(func_rand_summary$deseq2_rank_0.05, na.rm = TRUE)
sd(func_rand_summary$deseq2_rank_0.05, na.rm = TRUE)

median(func_rand_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
mean(func_rand_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)
sd(func_rand_summary$wilcoxon.musicc_rank_0.05, na.rm = TRUE)



wilcox.test(func_rand_summary$POMS_rank_0.05,
            func_rand_summary$deseq2_rank_0.05)


wilcox.test(func_rand_summary$POMS_rank_0.05,
            func_rand_summary$wilcoxon.musicc_rank_0.05)


# Prop. sig results for main text:
mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)
sd(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)

mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)
sd(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)



mean(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)
sd(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)

mean(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)
sd(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)


# Percent increases:
100 * ((mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE) / mean(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE)) / mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.05, na.rm = TRUE))

100 * ((mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE) / mean(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE)) / mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.05, na.rm = TRUE))



