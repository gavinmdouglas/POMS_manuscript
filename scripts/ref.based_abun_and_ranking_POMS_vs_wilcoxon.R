rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

ref_sim_func.based_ranking <- readRDS(file = "ref.genome_sim_summaries/ref.genome_sim_func.based_ranking.rds")

names(ref_sim_func.based_ranking) <- c("Setting 1", "Setting 2", "Setting 3", "Setting 4")

clean_tables <- list()
for(cutoff in names(ref_sim_func.based_ranking)) {
  clean_tables[[cutoff]] <- melt(ref_sim_func.based_ranking[[cutoff]][, -which(colnames(ref_sim_func.based_ranking[[cutoff]]) %in% c("replicate", "func", "all_num_pos_mags"))])
  clean_tables[[cutoff]]$cutoff <- cutoff
}

MAG_sim_func.based_ranking <- readRDS(file = "MAG_sim_summaries/MAG_sim_func.based_ranking.rds")

MAG_sim_func.based_ranking <- MAG_sim_func.based_ranking[, -which(colnames(MAG_sim_func.based_ranking) %in% c("replicate", "func", "all_num_pos_mags"))]
MAG_sim_func.based_ranking$cutoff <- "MAG-based"
clean_tables[["MAG-based"]] <- melt(MAG_sim_func.based_ranking)

ref_ranking <- do.call(rbind, clean_tables)

ref_ranking$variable <- as.character(ref_ranking$variable)

ref_ranking[which(ref_ranking$variable == "wilcoxon_rel_ranks"), "variable"] <- "Wilcoxon test"
ref_ranking[which(ref_ranking$variable == "POMS_rel_ranks"), "variable"] <- "POMS"

ref_ranking$variable <- factor(ref_ranking$variable, levels=c("Wilcoxon test", "POMS"))

ref_ranking$cutoff <- factor(ref_ranking$cutoff, levels=c("Setting 1", "Setting 2", "Setting 3", "Setting 4", "MAG-based"))


ref_sim_func.based_ranking_boxplots <- ggplot(ref_ranking, aes(x=cutoff, y=value * 100, fill=variable)) +
                                              geom_boxplot(outlier.shape = NA) +
                                              ylim(0, 100) +
                                              ylab("Relative ranking of focal gene (%)") +
                                              xlab("Abundance table sim. approach") +
                                              scale_fill_manual(name = "Method", values=c("springgreen4", "steelblue4")) +
                                              theme_bw() +
                                              theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                                    legend.position = c(0.3, 0.7),
                                                    legend.key = element_rect(fill = "white", colour = "white"),
                                                    legend.box.background = element_rect(colour="black"),
                                                    legend.background = element_rect(colour="light grey"))
                                                                                      

genome_percents_boxplots <- readRDS("ref.genome_sim_summaries/sim_abun_prevalence_boxplots.rds")

pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/ref.based_sim_ranking.pdf", width = 7.20472, height=4)
plot_grid(genome_percents_boxplots, ref_sim_func.based_ranking_boxplots, nrow=1, ncol=2, labels=c('a', 'b'))
dev.off()


# Summary stats
###NOTE: THE BELOW NAMES HAVE BEEN CHANGED, BUT THE VALUES HAVE NOT
wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.5`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.5`$POMS_rel_ranks)
V = 3793, p-value = 1.312e-05

wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.9`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.9`$POMS_rel_ranks)
V = 2067, p-value = 0.1156

wilcox.test(ref_sim_func.based_ranking$`mu=0.1, nu=0.99`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.1, nu=0.99`$POMS_rel_ranks)

wilcox.test(ref_sim_func.based_ranking$`mu=0.01, nu=0.99`$wilcoxon_rel_ranks, ref_sim_func.based_ranking$`mu=0.01, nu=0.99`$POMS_rel_ranks)

