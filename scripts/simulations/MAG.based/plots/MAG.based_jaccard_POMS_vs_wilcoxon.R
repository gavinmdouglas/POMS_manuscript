rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")


MAG_sim_func.based_jaccard <- readRDS(file = "MAG_sim_summaries/MAG_sim_func.based_jaccard.rds")

MAG_sim_func.based_jaccard_cutoff1 <- ggplot(MAG_sim_func.based_jaccard, aes(x=wilcoxon_jaccard_dist_0.0001, y=POMS_jaccard_dist_sig3_other0)) +
                                          geom_point(colour="coral2") +
                                          theme_bw() +
                                          xlim(0.4, 1) +
                                          ylim(0.4, 1) +
                                          ylab("POMS-based Jaccard distance") +
                                          xlab("Wilcoxon test-based Jaccard distance") +
                                          theme(plot.title = element_text(hjust = 0.5)) +
                                          ggtitle("Gene-based, cut-off #1") +
                                          guides(colour=FALSE) +
                                          geom_abline(intercept = 0, slope = 1, color="black", 
                                                      linetype="dashed", size=0.5)

MAG_sim_func.based_jaccard_cutoff1_marginal <- ggMarginal(MAG_sim_func.based_jaccard_cutoff1, type = "histogram", size=10)

MAG_sim_func.based_jaccard_cutoff2 <- ggplot(MAG_sim_func.based_jaccard, aes(x=wilcoxon_jaccard_dist_0.0001, y=POMS_jaccard_dist_sig3_other1)) +
  geom_point(colour="coral2") +
  theme_bw() +
  xlim(0.4, 1) +
  ylim(0.4, 1) +
  ylab("POMS-based Jaccard distance") +
  xlab("Wilcoxon test-based Jaccard distance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Gene-based, cut-off #2") +
  guides(colour=FALSE) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)

MAG_sim_func.based_jaccard_cutoff2_marginal <- ggMarginal(MAG_sim_func.based_jaccard_cutoff2, type = "histogram", size=10)


MAG_sim_func.based_jaccard_cutoff3 <- ggplot(MAG_sim_func.based_jaccard, aes(x=wilcoxon_jaccard_dist_0.0001, y=POMS_jaccard_dist_sig2_other0)) +
  geom_point(colour="coral2") +
  theme_bw() +
  xlim(0.4, 1) +
  ylim(0.4, 1) +
  ylab("POMS-based Jaccard distance") +
  xlab("Wilcoxon test-based Jaccard distance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Gene-based, cut-off #3") +
  guides(colour=FALSE) +
  geom_abline(intercept = 0, slope = 1, color="black", 
              linetype="dashed", size=0.5)

MAG_sim_func.based_jaccard_cutoff3_marginal <- ggMarginal(MAG_sim_func.based_jaccard_cutoff3, type = "histogram", size=10)


pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MGS_sim_sig_cutoff1_jaccard.pdf", width = 7.20472 / 2, height=4)
plot_grid(MAG_sim_func.based_jaccard_cutoff1_marginal, nrow=1, ncol=1)
dev.off()



pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/MGS_sim_sig_cutoff2and3_jaccard.pdf", width = 7.20472, height=4)
plot_grid(MAG_sim_func.based_jaccard_cutoff2_marginal, MAG_sim_func.based_jaccard_cutoff3_marginal,
          nrow=1, ncol=2, labels=c('a', 'b'))
dev.off()


# Summary statistics

mean(MAG_sim_func.based_jaccard$POMS_jaccard_dist_sig3_other0, na.rm = TRUE)
sd(MAG_sim_func.based_jaccard$POMS_jaccard_dist_sig3_other0, na.rm = TRUE)

mean(MAG_sim_func.based_jaccard$wilcoxon_jaccard_dist_0.0001, na.rm = TRUE)
sd(MAG_sim_func.based_jaccard$wilcoxon_jaccard_dist_0.0001, na.rm = TRUE)

wilcox.test(MAG_sim_func.based_jaccard$POMS_jaccard_dist_sig3_other0,
            MAG_sim_func.based_jaccard$wilcoxon_jaccard_dist_0.0001)

ks.test(MAG_sim_func.based_jaccard$POMS_jaccard_dist_sig3_other0,
        MAG_sim_func.based_jaccard$wilcoxon_jaccard_dist_0.0001)