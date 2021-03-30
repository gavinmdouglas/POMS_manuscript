### Filter to reduced set of high quality reference genomes.
### Simulate metagenomes (i.e. genome abundances across samples) based on the zero-inflated beta distribution for fitting
### Generalized Additive Models for Location Scale and Shape. Then made random sample splits.
### Also get random sample sets and prep tree.

rm(list = ls(all.names = TRUE))

setwd("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/JGI_PICRUSt2_genomes")

library(ape)
library(gamlss.dist)
library(sparseDOSSA)

# First, figure out high quality genomes that are unambiguously in PICRUSt2 database
# (i.e., not clustered in database)

ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)


ref_tree <- read.tree("GToTree_output/GToTree_output.tre")

ref_gtotree_qual <- read.table("GToTree_output/Genomes_summary_info.tsv",
                               header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
ref_gtotree_qual <- ref_gtotree_qual[which(ref_gtotree_qual$perc_comp >= 95), ]
ref_gtotree_qual <- ref_gtotree_qual[which(ref_gtotree_qual$perc_redund <= 5), ]
ref_gtotree_qual <- ref_gtotree_qual[which(rownames(ref_gtotree_qual) %in% rownames(ref_func)), ]

low_qual_tips <- ref_tree$tip.label[which(!ref_tree$tip.label %in% rownames(ref_gtotree_qual))]

ref_tree <- drop.tip(phy = ref_tree, tip = low_qual_tips, trim.internal = TRUE)

# Randomly keep 500 tips.
random_tips_to_drop <- sample(ref_tree$tip.label, size = length(ref_tree$tip.label) - 500)
ref_tree <- drop.tip(phy = ref_tree, tip = random_tips_to_drop, trim.internal = TRUE)

ref_tree <- multi2di(ref_tree)
write.tree(phy = ref_tree, file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/GToTree_output_subset.tre")

num_tips <- length(ref_tree$tip.label)


# First, randomly split the samples into two groups for 1000 replicates.

set.seed(62352)

ref_genomes_random_groups <- list()

BEZI_table <- data.frame(matrix(NA, nrow = num_tips, ncol = 1000))
rownames(BEZI_table) <- ref_tree$tip.label

for (rep in 1:1000) {
  samples <- sample(colnames(BEZI_table))
  ref_genomes_random_groups$group1[[rep]] <- samples[1:500]
  ref_genomes_random_groups$group2[[rep]] <- samples[501:1000]
}

saveRDS(object = ref_genomes_random_groups,
        file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/ref.based_random_groups.rds")


# Simulate tables with contrasting parameters.
BEZI_ensure_nonzero <- function(n_set, mu_set, nu_set, sigma_set, scale_factor=1e6, min_taxa_per_sample=1) {
  BEZI_sim <- c()
  while (length(which(BEZI_sim > 0)) < min_taxa_per_sample) {
    BEZI_sim <- round(rBEZI(n = n_set, mu = mu_set, nu = nu_set, sigma = sigma_set) * scale_factor)
  }
  
  return(BEZI_sim)
}


BEZI_table_mu0.1_nu0.50_sigma1 <- BEZI_table
BEZI_table_mu0.1_nu0.65_sigma1 <- BEZI_table
BEZI_table_mu0.1_nu0.80_sigma1 <- BEZI_table
BEZI_table_mu0.1_nu0.95_sigma1 <- BEZI_table


for (i in 1:1000) {
  BEZI_table_mu0.1_nu0.50_sigma1[, i] <- BEZI_ensure_nonzero(n_set = num_tips, mu_set = 0.1, nu_set = 0.50, sigma_set = 1, min_taxa_per_sample = 5)
  BEZI_table_mu0.1_nu0.65_sigma1[, i] <- BEZI_ensure_nonzero(n_set = num_tips, mu_set = 0.1, nu_set = 0.65, sigma_set = 1, min_taxa_per_sample = 5)
  BEZI_table_mu0.1_nu0.80_sigma1[, i] <- BEZI_ensure_nonzero(n_set = num_tips, mu_set = 0.1, nu_set = 0.80, sigma_set = 1, min_taxa_per_sample = 5)
  BEZI_table_mu0.1_nu0.95_sigma1[, i] <- BEZI_ensure_nonzero(n_set = num_tips, mu_set = 0.1, nu_set = 0.95, sigma_set = 1, min_taxa_per_sample = 5)
}

# Check if any rowSums are 0
which(rowSums(BEZI_table_mu0.1_nu0.50_sigma1) == 0)
which(rowSums(BEZI_table_mu0.1_nu0.65_sigma1) == 0)
which(rowSums(BEZI_table_mu0.1_nu0.80_sigma1) == 0)
which(rowSums(BEZI_table_mu0.1_nu0.95_sigma1) == 0)

# Basic plots contrasting simulations
par(mfrow = c(2, 2))
hist(colSums(BEZI_table_mu0.1_nu0.50_sigma1 > 0))
hist(colSums(BEZI_table_mu0.1_nu0.65_sigma1 > 0))
hist(colSums(BEZI_table_mu0.1_nu0.80_sigma1 > 0))
hist(colSums(BEZI_table_mu0.1_nu0.95_sigma1 > 0))

par(mfrow = c(2, 2))
hist(rowSums(BEZI_table_mu0.1_nu0.50_sigma1 > 0))
hist(rowSums(BEZI_table_mu0.1_nu0.65_sigma1 > 0))
hist(rowSums(BEZI_table_mu0.1_nu0.80_sigma1 > 0))
hist(rowSums(BEZI_table_mu0.1_nu0.95_sigma1 > 0))

write.table(x = BEZI_table_mu0.1_nu0.50_sigma1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/relabun_tables/BEZI_table_mu0.1_nu0.50_sigma1.tsv.gz",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = BEZI_table_mu0.1_nu0.65_sigma1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/relabun_tables/BEZI_table_mu0.1_nu0.65_sigma1.tsv.gz",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = BEZI_table_mu0.1_nu0.80_sigma1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/relabun_tables/BEZI_table_mu0.1_nu0.80_sigma1.tsv.gz",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

write.table(x = BEZI_table_mu0.1_nu0.95_sigma1,
            file = "/home/gavin/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/relabun_tables/BEZI_table_mu0.1_nu0.95_sigma1.tsv.gz",
            col.names = NA, row.names = TRUE, quote = FALSE, sep = "\t")

