rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

# Read in RDS objects that were previously generated.
MAG_sim_unperturbed_summary <- readRDS(file = "MAG_sim_summaries/MAG_sim_unperturbed_summary.rds")

MAG_sim_taxa.based_summary <- readRDS(file = "MAG_sim_summaries/MAG_sim_taxa.based_summary.rds")
MAG_sim_func.based_summary <- readRDS(file = "MAG_sim_summaries/MAG_sim_func.based_summary.rds")

# First plot correlations between num sig genes and num sig balances.
MAG_sim_taxa_sig3_other0_vs_sig_nodes <- ggplot(MAG_sim_taxa.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either3_other0)) +
  geom_point(colour="skyblue3") +
  theme_bw() +
  ylim(0, 0.5) +
  xlim(0, 70) +
  ggtitle("Taxa-based, cut-off #1") +
  xlab("Number of significant nodes") +
  ylab("Proportion of sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_taxa.based_summary[which(MAG_sim_taxa.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)

  
MAG_sim_taxa_sig3_other1_vs_sig_nodes <- ggplot(MAG_sim_taxa.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either3_other1)) +
  geom_point(colour="skyblue3") +
    theme_bw() +
    ylim(0, 0.5) +
    xlim(0, 70) +
    ggtitle("Taxa-based, cut-off #2") +
    xlab("Number of significant nodes") +
    ylab("Proportion of sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_taxa.based_summary[which(MAG_sim_taxa.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)

  
MAG_sim_taxa_sig2_other0_vs_sig_nodes <- ggplot(MAG_sim_taxa.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either2_other0)) +
  geom_point(colour="skyblue3") +
    theme_bw() +
    ylim(0, 0.5) +
    xlim(0, 70) +
    ggtitle("Taxa-based, cut-off #3") +
    xlab("Number of significant nodes") +
    ylab("Proportion of sig. KOs") +
    theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_taxa.based_summary[which(MAG_sim_taxa.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)



MAG_sim_func_sig3_other0_vs_sig_nodes <- ggplot(MAG_sim_func.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either3_other0)) +
  geom_point(colour="coral2") +
  theme_bw() +
  ylim(0, 0.5) +
  xlim(0, 70) +
  ggtitle("Gene-based, cut-off #1") +
  xlab("Number of significant nodes") +
  ylab("Proportion of sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_func.based_summary[which(MAG_sim_func.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)


MAG_sim_func_sig3_other1_vs_sig_nodes <- ggplot(MAG_sim_func.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either3_other1)) +
  geom_point(colour="coral2") +
  theme_bw() +
  ylim(0, 0.5) +
  xlim(0, 70) +
  ggtitle("Gene-based, cut-off #2") +
  xlab("Number of significant nodes") +
  ylab("Proportion of sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_func.based_summary[which(MAG_sim_func.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)


MAG_sim_func_sig2_other0_vs_sig_nodes <- ggplot(MAG_sim_func.based_summary, aes(x=num_sig_balance, y=POMS_total_sig_ko_either2_other0)) +
  geom_point(colour="coral2") +
  theme_bw() +
  ylim(0, 0.5) +
  xlim(0, 70) +
  ggtitle("Gene-based, cut-off #3") +
  xlab("Number of significant nodes") +
  ylab("Proportion of sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_cor(data=MAG_sim_func.based_summary[which(MAG_sim_func.based_summary$num_sig_balance > 5), ],
           method = "spearman", label.x = 30, label.y = 0.45)


pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MAG.based_sig_cutoffs_vs_sig_nodes.pdf", width = 8.5, height=6.5)
plot_grid(MAG_sim_taxa_sig3_other0_vs_sig_nodes, MAG_sim_func_sig3_other0_vs_sig_nodes,
          MAG_sim_taxa_sig3_other1_vs_sig_nodes, MAG_sim_func_sig3_other1_vs_sig_nodes,
          MAG_sim_taxa_sig2_other0_vs_sig_nodes, MAG_sim_func_sig2_other0_vs_sig_nodes,
          nrow=3, labels=c('a', 'b', 'c', 'd', 'e', 'f'))
dev.off()


# Compute additional spearman corr
tmp <- readRDS("MAG_sim_summaries/MAG_sim_func.based_ranking.rds")

cor.test(tmp$all_num_pos_mags, MAG_sim_func.based_summary$POMS_total_sig_ko_either3_other0, method="spearman")


