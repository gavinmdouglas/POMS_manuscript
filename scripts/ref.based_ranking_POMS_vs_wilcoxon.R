rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

ref_sim_func.based_ranking <- readRDS(file = "ref.genome_sim_summaries/ref.genome_sim_func.based_ranking.rds")

ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.5 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.5, aes(y=all_num_pos_mags, x=wilcoxon_rel_ranks * 100)) +
  geom_point(colour="springgreen4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("Wilcoxon test, Setting 1") +
  ylim(0, 3000) +
  xlim(-2, 100)
         
ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.9 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.9, aes(y=all_num_pos_mags, x=wilcoxon_rel_ranks * 100)) +
  geom_point(colour="springgreen4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("Wilcoxon test, Setting 2") +
  ylim(0, 3000) +
  xlim(-2, 100)
                
ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.99 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.99, aes(y=all_num_pos_mags, x=wilcoxon_rel_ranks * 100)) +
  geom_point(colour="springgreen4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("Wilcoxon test, Setting 3") +
  ylim(0, 3000) +
  xlim(-2, 100)


ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.01_nu_0.99 <- ggplot(ref_sim_func.based_ranking$mu0.01_nu_0.99, aes(y=all_num_pos_mags, x=wilcoxon_rel_ranks * 100)) +
  geom_point(colour="springgreen4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("Wilcoxon test, Setting 4") +
  ylim(0, 3000) +
  xlim(-2, 100)



ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.5 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.5, aes(y=all_num_pos_mags, x=POMS_rel_ranks * 100)) +
  geom_point(colour="steelblue4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("POMS, Setting 1") +
  ylim(0, 3000) +
  xlim(-2, 100)

ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.9 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.9, aes(y=all_num_pos_mags, x=POMS_rel_ranks * 100)) +
  geom_point(colour="steelblue4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("POMS, Setting 2") +
  ylim(0, 3000) +
  xlim(-2, 100)

ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.99 <- ggplot(ref_sim_func.based_ranking$mu0.1_nu_0.99, aes(y=all_num_pos_mags, x=POMS_rel_ranks * 100)) +
  geom_point(colour="steelblue4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("POMS, Setting 3") +
  ylim(0, 3000) +
  xlim(-2, 100)


ref_sim_func.based_POMS_rank_vs_num.genome_mu0.01_nu_0.99 <- ggplot(ref_sim_func.based_ranking$mu0.01_nu_0.99, aes(y=all_num_pos_mags, x=POMS_rel_ranks * 100)) +
  geom_point(colour="steelblue4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5, size=12)) +
  ggtitle("POMS, Setting 4") +
  ylim(0, 3000) +
  xlim(-2, 100)



pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/ref.genome_sim_func.based_ranking_vs_num.genome.pdf", width = 7.20472, height=6)
plot_grid(ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.5, ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.5,
          ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.01_nu_0.99, ref_sim_func.based_POMS_rank_vs_num.genome_mu0.01_nu_0.99,
          nrow=2, ncol=2,
          labels=c('a', 'b', 'c', 'd'))
dev.off()


pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/ref.genome_sim_func.based_ranking_vs_num.genome_additional.pdf", width = 7.20472, height=6)
plot_grid(ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.9, ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.9,
          ref_sim_func.based_wilcoxon_rank_vs_num.genome_mu0.1_nu_0.99, ref_sim_func.based_POMS_rank_vs_num.genome_mu0.1_nu_0.99,
          nrow=2, ncol=2,
          labels=c('a', 'b', 'c', 'd'))
dev.off()



# Summary stats
wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.5`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.5`$POMS_rel_ranks)
V = 3793, p-value = 1.312e-05

wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.9`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.9`$POMS_rel_ranks)
V = 2067, p-value = 0.1156

wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.99`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.99`$POMS_rel_ranks)

wilcox.test(ref_sim_func.based_ranking$`mu=0.01, nu=0.99`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.01, nu=0.99`$POMS_rel_ranks)


min(ref_sim_func.based_ranking$mu0.1_nu_0.9[which(ref_sim_func.based_ranking$mu0.1_nu_0.9$wilcoxon_rel_ranks > 0.1), "all_num_pos_mags"])


mean(ref_sim_func.based_ranking$mu0.1_nu_0.5$wilcoxon_rel_ranks)
sd(ref_sim_func.based_ranking$mu0.1_nu_0.5$wilcoxon_rel_ranks)

mean(ref_sim_func.based_ranking$mu0.1_nu_0.9$wilcoxon_rel_ranks)
sd(ref_sim_func.based_ranking$mu0.1_nu_0.9$wilcoxon_rel_ranks)

wilcox.test(ref_sim_func.based_ranking$mu0.1_nu_0.5$wilcoxon_rel_ranks, ref_sim_func.based_ranking$mu0.1_nu_0.9$wilcoxon_rel_ranks)

cor.test(ref_sim_func.based_ranking$mu0.01_nu_0.99$wilcoxon_rel_ranks, ref_sim_func.based_ranking$mu0.01_nu_0.99$all_num_pos_mags, method="pearson")