rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggExtra)
library(ggplot2)
library(ggpubr)
library(reshape2)

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

MAG.based_abun <- read.table("MAG_sim_starting_files/bwa_depth_min25coverage.tsv.gz", header=TRUE, sep="\t", row.names=1)
MAG.based_samples <- readRDS("MAG_sim_starting_files/almeida_healthy_random_groups.rds")

MAG.based_samples <- c(MAG.based_samples$group1[[1]], MAG.based_samples$group1[[2]], MAG.based_samples$group2[[1]], MAG.based_samples$group2[[2]])
MAG.based_samples <- MAG.based_samples[-which(duplicated(MAG.based_samples))]

MAG.based_abun <- MAG.based_abun[, MAG.based_samples]
MAG.based_abun <- MAG.based_abun[-which(rowSums(MAG.based_abun) == 0), ]


BEZI_tables <- list()

BEZI_tables[["mu0.1_nu_0.5"]] <- read.table(file = "ref.genome_simulated_abundances/BEZI_genome_abun_mu0.1_nu_0.5_sigma1.tsv.gz",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.9"]] <- read.table(file = "ref.genome_simulated_abundances/BEZI_genome_abun_mu0.1_nu_0.9_sigma1.tsv.gz",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.99"]] <- read.table(file = "ref.genome_simulated_abundances/BEZI_genome_abun_mu0.1_nu_0.99_sigma1.tsv.gz",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.01_nu_0.99"]] <- read.table(file = "ref.genome_simulated_abundances/BEZI_genome_abun_mu0.01_nu_0.99_sigma1.tsv.gz",
                                              header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")


genome_percents <- data.frame(mu0.1_nu0.5=c(rowSums(BEZI_tables[["mu0.1_nu_0.5"]] > 0) / ncol(BEZI_tables[["mu0.1_nu_0.5"]]), rep(NA, 3000 - nrow(BEZI_tables[["mu0.1_nu_0.5"]]))),
                              mu0.1_nu0.9=c(rowSums(BEZI_tables[["mu0.1_nu_0.9"]] > 0) / ncol(BEZI_tables[["mu0.1_nu_0.9"]]), rep(NA, 3000 - nrow(BEZI_tables[["mu0.1_nu_0.9"]]))),
                              mu0.1_nu0.99=c(rowSums(BEZI_tables[["mu0.1_nu_0.99"]] > 0) / ncol(BEZI_tables[["mu0.1_nu_0.99"]]), rep(NA, 3000 - nrow(BEZI_tables[["mu0.1_nu_0.99"]]))),
                              mu0.01_nu0.99=c(rowSums(BEZI_tables[["mu0.01_nu_0.99"]] > 0) / ncol(BEZI_tables[["mu0.01_nu_0.99"]]), rep(NA, 3000 - nrow(BEZI_tables[["mu0.01_nu_0.99"]])))) * 100

colnames(genome_percents) <- c("Setting 1", "Setting 2", "Setting 3", "Setting 4")


genome_percents[["MAG-based"]] <- c(rowSums(MAG.based_abun > 0) / ncol(MAG.based_abun), rep(NA, 3000 - nrow(MAG.based_abun))) * 100

genome_percents_melt <- melt(genome_percents)

genome_percents_melt$Dataset <- "Reference genomes"
genome_percents_melt[which(genome_percents_melt$variable == "MAG-based"), "Dataset"] <- "MAGs"

genome_percents_boxplots <- ggplot(genome_percents_melt, aes(x=variable, y=value, fill=Dataset)) +
                                   geom_boxplot() +
                                   scale_fill_manual(values=c("light blue", "grey5")) +
                                   theme_bw() +
                                   ylab("Genome prevalence across samples") +
                                   xlab("Abundance table sim. approach") +
                                   ylim(c(0, 100)) +
                                   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
                                         legend.position = c(0.3, 0.7),
                                         legend.key = element_rect(fill = "white", colour = "white"),
                                         legend.box.background = element_rect(colour="black"),
                                         legend.background = element_rect(colour="light grey"))


# Rather than making plot, save as RDS.
# pdf(file = "/home/gavin/github_repos/POMS_manuscript/figures/sim_abun_tables_prevalence.pdf", width = 7.20472/2, height=4)
# plot_grid(genome_percents_boxplots)
# dev.off()

saveRDS(object = genome_percents_boxplots, file = "ref.genome_sim_summaries/sim_abun_prevalence_boxplots.rds")

# Summary statistics:
mean(rowSums(BEZI_tables[["mu0.01_nu_0.99"]] > 0))

mean(rowSums(BEZI_tables[["mu0.1_nu_0.5"]] > 0))

