rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

func.based_summary <- readRDS(file = "simulation_summaries/func.based_summary.rds")

func.based_summary_rank <- func.based_summary[, grep("rank", colnames(func.based_summary))]
func.based_summary_rank <- func.based_summary_rank[, -grep("rel_rank", colnames(func.based_summary_rank))]

sapply(func.based_summary_rank, mean, na.rm = TRUE)
sort(sapply(func.based_summary_rank, median, na.rm = TRUE))
sapply(func.based_summary_rank, sd, na.rm = TRUE)

wilcox.test(func.based_summary_rank$POMS_rank_0.05, func.based_summary_rank$regress_specificity_rank_0.05)
wilcox.test(func.based_summary_rank$POMS_rank_0.05, func.based_summary_rank$regress_sig_taxa_rank_0.05)
wilcox.test(func.based_summary_rank$regress_sig_taxa_rank_0.05, func.based_summary_rank$regress_specificity_rank_0.05)



# Num MAGs encoding:

median(func.based_summary[which(func.based_summary$POMS_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$regress_specificity_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$regress_sig_taxa_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$aldex2_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$deseq2_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$limma.voom_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$wilcoxon.relab_rank_0.05 <= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$wilcoxon.musicc_rank_0.05 <= 10), "num_focal_pos_mags"])


median(func.based_summary[which(func.based_summary$POMS_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$regress_specificity_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$regress_sig_taxa_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$aldex2_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$deseq2_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$limma.voom_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$wilcoxon.relab_rank_0.05 >= 10), "num_focal_pos_mags"])
median(func.based_summary[which(func.based_summary$wilcoxon.musicc_rank_0.05 >= 10), "num_focal_pos_mags"])



# How often the focal gene was significant.
colSums(! is.na(func.based_summary_rank)) / nrow(func.based_summary_rank)
mean((colSums(! is.na(func.based_summary_rank)) / nrow(func.based_summary_rank))[-1])