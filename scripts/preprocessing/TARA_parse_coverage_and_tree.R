### Parse TARA read mapping files to produce tables of MAG abundances across samples.
### Will output two sets based on calling MAGs as positive based on 25% and 50% coverage.
### Also prep the KO table at the end of this file.

rm(list = ls(all.names = TRUE))

library(ape)

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/TARA/")

MAG_tree <- read.tree("GToTree_output/GToTree_output.tre")

MAG_stats <- read.table("GToTree_output/Genomes_summary_info.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)
MAG_stats <- MAG_stats[MAG_tree$tip.label, ]

MAGs2rm <- rownames(MAG_stats)[which(MAG_stats$perc_comp < 60)]
MAGs2rm <- c(MAGs2rm, rownames(MAG_stats)[which(MAG_stats$perc_redund > 10)])
if (length(which(duplicated(MAGs2rm))) > 0) { MAGs2rm <- MAGs2rm[-which(duplicated(MAGs2rm))] }

TARA_abun <- read.table("NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/abundance.txt", header = TRUE, sep = "\t", row.names = 1)
TARA_coverage <- read.table("NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/mean_coverage.txt", header = TRUE, sep = "\t", row.names = 1)

TARA_abun <- TARA_abun[MAG_tree$tip.label, ]
TARA_coverage <- TARA_coverage[MAG_tree$tip.label, ]

TARA_abun <- TARA_abun[-which(rownames(TARA_abun) %in% MAGs2rm), ]
TARA_coverage <- TARA_coverage[-which(rownames(TARA_coverage) %in% MAGs2rm), ]

TARA_coverage[TARA_coverage < 1] <- 0
TARA_abun[TARA_coverage == 0] <- 0

write.table(x = TARA_abun, file = "NON-REDUNDANT-MAGs-SUMMARY/bins_across_samples/modified/TARA_abundance_min_mean_coverage1.tsv", col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")


MAG_tree <- drop.tip(phy = MAG_tree, tip = MAGs2rm, trim.internal = TRUE)
MAG_tree <- multi2di(MAG_tree)
MAG_tree$node.label <- NULL
write.tree(phy = MAG_tree, file = "GToTree_output/GToTree_output_modified.tre")
