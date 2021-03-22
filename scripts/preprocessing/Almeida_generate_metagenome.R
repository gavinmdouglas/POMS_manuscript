### Multiply MAG and genome annotation tables to get unstratified metagenome abundances per sample.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/Almeida2019_dataset/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

almeida_KO <- data.frame(t(read.table("functional_analyses/modified/kegg_summary.tsv",
                         header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names=FALSE)), check.names=FALSE)

almeida_abun_cov25 <- read.table("mapping_results/modified/bwa_depth_min25coverage.tsv",
                                 header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names=FALSE)

almeida_abun_cov50 <- read.table("mapping_results/modified/bwa_depth_min50coverage.tsv",
                                 header=TRUE, sep="\t", row.names=1, stringsAsFactors = FALSE, check.names=FALSE)

almeida_metagenome_cov25 <- calc_func_abun(in_abun = almeida_abun_cov25, in_func = almeida_KO, ncores = 60)
almeida_metagenome_cov50 <- calc_func_abun(in_abun = almeida_abun_cov50, in_func = almeida_KO, ncores = 60)

almeida_metagenome_cov25 <- almeida_metagenome_cov25[, -which(colSums(almeida_metagenome_cov25) == 0)]
almeida_metagenome_cov50 <- almeida_metagenome_cov50[, -which(colSums(almeida_metagenome_cov50) == 0)]

write.table(x = almeida_metagenome_cov25, file = "functional_analyses/modified/almeida_metagenome_cov25.tsv",
            sep="\t", quote = FALSE, row.names = TRUE, col.names = NA)

write.table(x = almeida_metagenome_cov50, file = "functional_analyses/modified/almeida_metagenome_cov50.tsv",
            sep="\t", quote = FALSE, row.names = TRUE, col.names = NA)

