rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

taxa_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.musicc.rds")

func_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")

MAG.based_sim_taxa_POMS_vs_wilcoxon.musicc_sig_prop <- ggplot(taxa_rand_summary_POMS_wilcoxon.musicc, aes(x = POMS_sig_0.25,
                                                                                                          y = wilcoxon.musicc_sig_0.001)) +
                                                              geom_point(colour = "skyblue3") +
                                                              theme_bw() +
                                                              xlim(0, 1) +
                                                              ylim(0, 1) +
                                                              xlab("POMS proportion sig. gene families") +
                                                              ylab("Wilcoxon test proportion sig. gene families") +
                                                              ggtitle("Random taxa") +
                                                              theme(plot.title = element_text(hjust = 0.5, vjust = -1))

MAG.based_sim_taxa_POMS_vs_wilcoxon.musicc_sig_prop_marginal <- ggMarginal(MAG.based_sim_taxa_POMS_vs_wilcoxon.musicc_sig_prop, type = "histogram", size = 10)

MAG.based_sim_func_POMS_vs_wilcoxon.musicc_sig_prop <- ggplot(func_rand_summary_POMS_wilcoxon.musicc, aes(x = POMS_sig_0.25,
                                                                                                          y = wilcoxon.musicc_sig_0.001)) +
                                                              geom_point(colour = "coral2") +
                                                              theme_bw() +
                                                              xlim(0, 1) +
                                                              ylim(0, 1) +
                                                              xlab("POMS proportion sig. gene families") +
                                                              ylab("Wilcoxon test proportion sig. gene families") +
                                                              ggtitle("Focal gene") +
                                                              theme(plot.title = element_text(hjust = 0.5, vjust = -1))

MAG.based_sim_func_POMS_vs_wilcoxon.musicc_sig_prop_marginal <- ggMarginal(MAG.based_sim_func_POMS_vs_wilcoxon.musicc_sig_prop, type = "histogram", size = 10)


MAG.based_POMS_vs_wilcoxon.musicc_sig_prop_marginal_combined <- plot_grid(MAG.based_sim_taxa_POMS_vs_wilcoxon.musicc_sig_prop_marginal,
                                                                          MAG.based_sim_func_POMS_vs_wilcoxon.musicc_sig_prop_marginal,
                                                                          nrow = 1, ncol = 2, labels = c('a', 'b'))

ggsave(filename = "~/github_repos/POMS_manuscript/figures/Maintext_MAG.based_sims_prop_sig_genes.pdf",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_sig_prop_marginal_combined,
       device = "pdf",
       width = 8,
       height = 4,
       dpi = 600)
       
       
ggsave(filename = "~/github_repos/POMS_manuscript/figures/Maintext_MAG.based_sims_prop_sig_genes.png",
       plot = MAG.based_POMS_vs_wilcoxon.musicc_sig_prop_marginal_combined,
       device = "png",
       width = 8,
       height = 4,
       dpi = 300)


# Calculate summary statistics to report in main-text
mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)
sd(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)

mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)
sd(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)


mean(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)
sd(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)

mean(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)
sd(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)


# Percent increase:
((mean(func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001) - 
  mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001) ) / 
  mean(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)) * 100

((mean(func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25) - 
  mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25) ) / 
  mean(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)) * 100



# Wilcoxon tests:
wilcox.test(taxa_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001,
            func_rand_summary_POMS_wilcoxon.musicc$wilcoxon.musicc_sig_0.001)


wilcox.test(taxa_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25,
            func_rand_summary_POMS_wilcoxon.musicc$POMS_sig_0.25)