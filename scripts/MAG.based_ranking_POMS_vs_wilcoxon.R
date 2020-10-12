rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")


MAG_sim_func.based_ranking <- readRDS(file = "MAG_sim_summaries/MAG_sim_func.based_ranking.rds")

MAG_sim_func.based_ranking_wilcoxon <- ggplot(MAG_sim_func.based_ranking, aes(x=wilcoxon_rel_ranks * 100)) +
                                          geom_histogram(fill="springgreen4") +
                                          xlim(-2, 105) +
                                          xlab("Relative ranking of focal gene (%)") +
                                          ylab("Number of simulation replicates") + 
                                          ggtitle("Wilcoxon test") +
                                          theme_bw() +
                                          theme(plot.title = element_text(hjust = 0.5))
                                          
  

MAG_sim_func.based_ranking_POMS <- ggplot(MAG_sim_func.based_ranking, aes(x=POMS_rel_ranks * 100)) +
                                      geom_histogram(fill="steelblue4") +
                                      xlim(-2, 105) +
                                      xlab("Relative ranking of focal gene (%)") +
                                      ylab("Number of simulation replicates") +
                                      ggtitle("POMS") +
                                      theme_bw() +
                                      theme(plot.title = element_text(hjust = 0.5))
                                      

MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG <- ggplot(MAG_sim_func.based_ranking, aes(y=all_num_pos_mags, x=wilcoxon_rel_ranks * 100)) +
  geom_point(colour="springgreen4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Wilcoxon test") +
  ylim(0, 3000) +
  xlim(-2, 105)

MAG_sim_func.based_ranking_POMS_vs_num.MAG <- ggplot(MAG_sim_func.based_ranking, aes(y=all_num_pos_mags, x=POMS_rel_ranks * 100)) +
  geom_point(colour="steelblue4", size=2) +
  theme_bw() +
  xlab("Relative ranking of focal gene (%)") +
  ylab("No. MAGs encoding focal gene") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("POMS") +
  ylim(0, 3000) +
  xlim(-2, 105)


pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MAG.based_sim_ranking.pdf", width = 7.20472, height=6)
plot_grid(MAG_sim_func.based_ranking_wilcoxon, MAG_sim_func.based_ranking_POMS,
          MAG_sim_func.based_ranking_wilcoxon_vs_num.MAG, MAG_sim_func.based_ranking_POMS_vs_num.MAG,
          nrow=2, ncol=2, labels=c('a', 'b', 'c', 'd'))
dev.off()


# Summary statistics

mean(MAG_sim_func.based_ranking$POMS_rel_ranks, na.rm = TRUE)
sd(MAG_sim_func.based_ranking$POMS_rel_ranks, na.rm = TRUE)

mean(MAG_sim_func.based_ranking$wilcoxon_rel_ranks, na.rm = TRUE)
sd(MAG_sim_func.based_ranking$wilcoxon_rel_ranks, na.rm = TRUE)


wilcox.test(MAG_sim_func.based_ranking$POMS_rel_ranks,
            MAG_sim_func.based_ranking$wilcoxon_rel_ranks)


median(MAG_sim_func.based_ranking[which(MAG_sim_func.based_ranking$wilcoxon_rel_ranks < 0.01), "all_num_pos_mags"])
# Percentages of total MAGs calculated by divided by 2,505