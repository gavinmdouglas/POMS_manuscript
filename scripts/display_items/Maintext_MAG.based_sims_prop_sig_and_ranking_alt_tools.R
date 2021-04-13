rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

taxa_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.musicc.rds")
taxa_rand_summary_POMS_wilcoxon.relab_limma_voom <- readRDS(file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.relab_limma_voom.rds")
taxa_rand_summary_POMS_aldex2_deseq2 <- readRDS(file = "simulation_summaries/taxa_rand_summary_POMS_aldex2_deseq2.rds")

func_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")
func_rand_summary_POMS_wilcoxon.relab_limma_voom <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.relab_limma_voom.rds")
func_rand_summary_POMS_aldex2_deseq2 <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_aldex2_deseq2.rds")

# Create dataframe with combined proportion of significant genes for each tool / simulation type.

MAG.based_all_prop_sig <- data.frame(Tool = c(rep("POMS", 2000),
                                              
                                              rep("Wilcoxon test (MUSiCC)", 2000),
                                              
                                              rep("Wilcoxon test (relab.)", 2000),
                                              
                                              rep("limma-voom", 2000),
                                              
                                              rep("ALDEx2", 20),
                                              
                                              rep("DESeq2", 20)),
                                     
                                     Simulation = 
                                              c(rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 1000),
                                              rep("Focal gene", 1000),
                                              
                                              rep("Random taxa", 10),
                                              rep("Focal gene", 10),
                                              
                                              rep("Random taxa", 10),
                                              rep("Focal gene", 10)),
                                              
                                     prop = c(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25,
                                              func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25,
                                              
                                              taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001,
                                              func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001,
                                              
                                              taxa_rand_summary_POMS_wilcoxon.relab_limma_voom$wilcoxon.relab_sig_0.001,
                                              func_rand_summary_POMS_wilcoxon.relab_limma_voom$wilcoxon.relab_sig_0.001,
                                              
                                              taxa_rand_summary_POMS_wilcoxon.relab_limma_voom$limma.voom_sig_0.001,
                                              func_rand_summary_POMS_wilcoxon.relab_limma_voom$limma.voom_sig_0.001,
                                              
                                              taxa_rand_summary_POMS_aldex2_deseq2$aldex2_sig_0.001,
                                              func_rand_summary_POMS_aldex2_deseq2$aldex2_sig_0.001,
                                              
                                              taxa_rand_summary_POMS_aldex2_deseq2$deseq2_sig_0.001,
                                              func_rand_summary_POMS_aldex2_deseq2$deseq2_sig_0.001)
                                              
                                    )


MAG.based_all_prop_sig$Simulation <- factor(MAG.based_all_prop_sig$Simulation, 
                                            levels = c("Random taxa", "Focal gene"))

MAG.based_all_prop_sig$Tool <- factor(MAG.based_all_prop_sig$Tool, 
                                            levels = c("POMS", "ALDEx2", "DESeq2", "limma-voom", "Wilcoxon test (MUSiCC)", "Wilcoxon test (relab.)"))

MAG.based_all_prop_sig_boxplots <- ggplot(MAG.based_all_prop_sig, aes(x = Tool, y = prop, fill = Simulation)) +
                                          geom_boxplot() +
                                          theme_bw() +
                                          ylim(0, 1) +
                                          ylab("Proportion sig. gene families") +
                                          scale_fill_manual(values = c("skyblue3", "coral2")) +
                                          theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))





# Now get rankings for each tool (only for Focal gene simulations).
MAG.based_all_rank <- data.frame(Tool = c(rep("POMS", 1000),
                                              
                                              rep("Wilcoxon test (MUSiCC)", 1000),
                                              
                                              rep("Wilcoxon test (relab.)", 1000),
                                              
                                              rep("limma-voom", 1000),
                                              
                                              rep("ALDEx2", 10),
                                              
                                              rep("DESeq2", 10)),
                                     
                                     rank = c(func_rand_summary_POMS_wilcoxon.musicc$POMS_rank_0.25,
                                              
                                              func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_rank_0.001,
                                              
                                              func_rand_summary_POMS_wilcoxon.relab_limma_voom$wilcoxon.relab_rank_0.001,
                                              
                                              func_rand_summary_POMS_wilcoxon.relab_limma_voom$limma.voom_rank_0.001,
                                              
                                              func_rand_summary_POMS_aldex2_deseq2$aldex2_rank_0.001,
                                              
                                              func_rand_summary_POMS_aldex2_deseq2$deseq2_rank_0.001))


MAG.based_all_rank$Tool <- factor(MAG.based_all_rank$Tool, 
                                  levels = c("POMS", "ALDEx2", "DESeq2", "limma-voom", "Wilcoxon test (MUSiCC)", "Wilcoxon test (relab.)"))

MAG.based_all_rank$Fill <- rep("Differential abundance approach", nrow(MAG.based_all_rank))
MAG.based_all_rank[which(MAG.based_all_rank$Tool == "POMS"), "Fill"] <- "Balance tree-based approach"




MAG.based_all_rank_boxplots <- ggplot(MAG.based_all_rank, aes(x = Tool, y = rank, fill = Fill)) +
                                      geom_boxplot() +
                                      theme_bw() +
                                      ylab("Focal gene ranking") +
                                      scale_fill_manual(values = c("steelblue4", "springgreen4")) +
                                      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
                                      theme(legend.position = "none")



MAG.based_alt_tool_comparisons <- plot_grid(MAG.based_all_prop_sig_boxplots,
                                            MAG.based_all_rank_boxplots,
                                            nrow = 1, ncol = 2, labels = c('a', 'b'), rel_widths = c(2.5, 1))

ggsave(filename = "~/github_repos/POMS_manuscript/figures/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.pdf",
       plot = MAG.based_alt_tool_comparisons,
       device = "pdf",
       width = 6,
       height = 4,
       dpi = 600)
       
       
ggsave(filename = "~/github_repos/POMS_manuscript/figures/Maintext_MAG.based_sims_prop_sig_and_ranking_alt_tools.png",
       plot = MAG.based_alt_tool_comparisons,
       device = "png",
       width = 6,
       height = 4,
       dpi = 300)
