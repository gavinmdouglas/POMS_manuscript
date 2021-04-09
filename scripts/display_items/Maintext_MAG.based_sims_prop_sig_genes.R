rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

taxa_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/taxa_rand_summary_POMS_wilcoxon.musicc.rds")
func_rand_summary_POMS_wilcoxon.musicc <- readRDS(file = "simulation_summaries/func_rand_summary_POMS_wilcoxon.musicc.rds")

# Plot scatterplot of prop sig KOs from POMS vs wilcoxon test
MAG_sim_taxa_sig3_other0_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either3_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="skyblue3") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Random taxa, cut-off #1") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig3_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig3_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_taxa_sig3_other1_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either3_other1, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="skyblue3") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Random taxa, cut-off #2") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig3_other1_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig3_other1_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_taxa_sig2_other0_vs_wilcoxon <- ggplot(MAG_sim_taxa.based_summary, aes(x=POMS_total_sig_ko_either2_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="skyblue3") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Random taxa, cut-off #3") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_taxa_sig2_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_taxa_sig2_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig3_other0_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either3_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="coral2") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Focal gene, cut-off #1") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig3_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig3_other0_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig3_other1_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either3_other1, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="coral2") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Focal gene, cut-off #2") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig3_other1_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig3_other1_vs_wilcoxon, type = "histogram", size=10)


MAG_sim_func_sig2_other0_vs_wilcoxon <- ggplot(MAG_sim_func.based_summary, aes(x=POMS_total_sig_ko_either2_other0, y=wilcoxon_sig_ko_0.0001)) +
  geom_point(colour="coral2") +
  theme_bw() +
  xlim(0, 0.80) +
  ylim(0, 0.80) +
  xlab("POMS proportion sig. KOs") +
  ylab("Wilcoxon test proportion sig. KOs") +
  ggtitle("Focal gene, cut-off #3") +
  theme(plot.title = element_text(hjust = 0.5))

MAG_sim_func_sig2_other0_vs_wilcoxon_marginal <- ggMarginal(MAG_sim_func_sig2_other0_vs_wilcoxon, type = "histogram", size=10)




pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MAG.based_sim_sig_cutoffs_vs_wilcoxon.pdf",
    width = 6.5, height=8.5)
plot_grid(MAG_sim_taxa_sig3_other0_vs_wilcoxon_marginal, MAG_sim_func_sig3_other0_vs_wilcoxon_marginal,
          MAG_sim_taxa_sig3_other1_vs_wilcoxon_marginal, MAG_sim_func_sig3_other1_vs_wilcoxon_marginal,
          MAG_sim_taxa_sig2_other0_vs_wilcoxon_marginal, MAG_sim_func_sig2_other0_vs_wilcoxon_marginal,
          nrow=3, ncol=2, labels=c('a', 'b', 'c', 'd', 'e', 'f'))
dev.off()



# Calculate summary statistics to report in main-text

mean(MAG_sim_taxa.based_summary$wilcoxon_sig_ko_0.0001)
sd(MAG_sim_taxa.based_summary$wilcoxon_sig_ko_0.0001)

mean(MAG_sim_taxa.based_summary$POMS_total_sig_ko_either3_other0)
sd(MAG_sim_taxa.based_summary$POMS_total_sig_ko_either3_other0)


mean(MAG_sim_func.based_summary$wilcoxon_sig_ko_0.0001)
sd(MAG_sim_func.based_summary$wilcoxon_sig_ko_0.0001)

mean(MAG_sim_func.based_summary$POMS_total_sig_ko_either3_other0)
sd(MAG_sim_func.based_summary$POMS_total_sig_ko_either3_other0)


# Percent increase:
(( mean(MAG_sim_func.based_summary$POMS_total_sig_ko_either3_other0) - 
     mean(MAG_sim_taxa.based_summary$POMS_total_sig_ko_either3_other0) ) / 
    mean(MAG_sim_taxa.based_summary$POMS_total_sig_ko_either3_other0)) * 100



# Wilcoxon tests:
wilcox.test(MAG_sim_func.based_summary$wilcoxon_sig_ko_0.0001,
            MAG_sim_taxa.based_summary$wilcoxon_sig_ko_0.0001)


wilcox.test(MAG_sim_func.based_summary$POMS_total_sig_ko_either3_other0,
            MAG_sim_taxa.based_summary$POMS_total_sig_ko_either3_other0)