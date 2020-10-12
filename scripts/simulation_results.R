rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

# Read in RDS objects that were previously generated.
MAG_sim_taxa.based_summary <- readRDS(file = "MAG_sim_summaries/MAG_sim_taxa.based_summary.rds")
MAG_sim_func.based_summary <- readRDS(file = "MAG_sim_summaries/MAG_sim_func.based_summary.rds")

# Plot scatterplot of prop sig KOs from POMS vs wilcoxon test
MAG_sim_taxa_sig3_other0_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either3_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Taxa-based, cut-off #1") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig3_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig3_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_taxa_sig3_other1_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either3_other1, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Taxa-based, cut-off #2") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig3_other1_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig3_other1_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_taxa_sig2_other0_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either2_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Taxa-based, cut-off #3") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig2_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig2_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig3_other0_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either3_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Gene-based, cut-off #1") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig3_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig3_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig3_other1_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either3_other1, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Gene-based, cut-off #2") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig3_other1_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig3_other1_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig2_other0_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either2_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point() +
  theme_bw() +
  xlim(0, 0.5) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Gene-based, cut-off #3") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig2_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig2_other0_vs_wilcoxon, type = "histogram", size=10)


pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MGS_sim_sig_cutoff2and3_vs_sig_nodes.pdf", width = 7.20472, height=8)
plot_grid(MAG_sim_taxa_sig3_other1_vs_sig_nodes, MAG_sim_func_sig3_other1_vs_sig_nodes,
          MAG_sim_taxa_sig2_other0_vs_sig_nodes, MAG_sim_func_sig2_other0_vs_sig_nodes,
          nrow=2, ncol=2, labels=c('a', 'b', 'c', 'd'))
dev.off()



