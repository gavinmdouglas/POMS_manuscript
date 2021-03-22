### Run POMS and test for significant differences in gene families between 1000 different random sets of healthy samples
### Assumption here is that there is no real differences - which should be true on average, but there may be
### some differences by chance (esp. since these samples are from different studies and populations).

rm(list=ls(all.names=TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(reshape2)
library(ggplot2)
library(parallel)

# Read in pre-determined random sample groupings.
almeida_healthy_random_groups <- readRDS("intermediate_objects/almeida_healthy_random_groups.rds")
ERP003612_healthy_random_groups <- readRDS("intermediate_objects/ERP003612_healthy_random_groups.rds")
SRP095580_healthy_random_groups <- readRDS("intermediate_objects/SRP095580_healthy_random_groups.rds")

almeida_metagenome_cov25 <- read.table("functional_analyses/modified/almeida_metagenome_cov25.tsv",
                                       header=TRUE, sep="\t", row.names=1, check.names = FALSE)

# First determine the FPR for basic Wilcoxon test.
wilcoxon_p_almeida_cov25 <- mclapply(X = 1:1000, FUN = function(x) { wilcoxon_2group_pvalues(almeida_metagenome_cov25, almeida_healthy_random_groups$group1[[x]], almeida_healthy_random_groups$group2[[x]]) }, mc.cores=60)

# Commands for other sample subsets, which weren't used for this particular analysis.
# wilcoxon_p_ERP003612_cov25 <- mclapply(X = 1:1000, FUN = function(x) { wilcoxon_2group_pvalues(almeida_metagenome_cov25, ERP003612_healthy_random_groups$group1[[x]], ERP003612_healthy_random_groups$group2[[x]]) }, mc.cores=60)
# wilcoxon_p_SRP095580_cov25 <- mclapply(X = 1:1000, FUN = function(x) { wilcoxon_2group_pvalues(almeida_metagenome_cov25, SRP095580_healthy_random_groups$group1[[x]], SRP095580_healthy_random_groups$group2[[x]]) }, mc.cores=60)
# wilcoxon_p_ERP003612_cov25_combined <- unlist(wilcoxon_p_ERP003612_cov25)
# wilcoxon_p_SRP095580_cov25_combined <- unlist(wilcoxon_p_SRP095580_cov25)


# Run balance tree approach on this random subsampling.
almeida_func <- read.table(file = "functional_analyses/kegg_summary.csv", sep=",", stringsAsFactors = FALSE, quote="", comment.char = "",
                           header=TRUE, check.names = FALSE, row.names=1)

almeida_func <- data.frame(t(almeida_func), check.names = FALSE)

almeida_tree <- read.tree(file = "phylogenies/raxml_hgr-umgs_phylogeny.nwk")
almeida_abun <- read.table(file = "mapping_results/modified/bwa_depth_min25coverage.tsv", header=TRUE, sep="\t", check.names=FALSE,
                           row.names=1, quote="", comment.char="")

almeida_taxa_umgs <- read.table("taxonomy/taxonomy_umgs.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_umgs) <- almeida_taxa_umgs$MAG_ID
almeida_taxa_umgs <- almeida_taxa_umgs[, -which(colnames(almeida_taxa_umgs) %in% c("UMGS_ID", "MAG_ID"))]

almeida_taxa_hgr <- read.table("taxonomy/taxonomy_hgr.tab", header=TRUE, sep="\t", stringsAsFactors = FALSE)
rownames(almeida_taxa_hgr) <- almeida_taxa_hgr$Genome
almeida_taxa_hgr <- almeida_taxa_hgr[, -which(colnames(almeida_taxa_hgr) == "Genome")]

almeida_taxa <- rbind(almeida_taxa_hgr, almeida_taxa_umgs)

balance_almeida_cov25_reps <- list()

for(rep in 1:1000) {

  print(rep)
  
  group1_samples <- almeida_healthy_random_groups$group1[[rep]]
  group2_samples <- almeida_healthy_random_groups$group2[[rep]]
  
  in_abun <- subset_abun_table(in_abun = almeida_abun, col2keep = c(group1_samples, group2_samples))
  
  in_func <- almeida_func[rownames(in_abun), ]
  
  balance_almeida_cov25_reps[[rep]] <- two_group_balance_tree_pipeline(abun=in_abun,
                                                                            func=in_func,
                                                                            phylogeny=almeida_tree,
                                                                            taxa=almeida_taxa,
                                                                            group1_samples = group1_samples,
                                                                            group2_samples = group2_samples,
                                                                            ncores=60,
                                                                            balance_p_cutoff = 0.05,
                                                                            balance_correction = "none",
                                                                            function_p_cutoff = 0.05,
                                                                            function_correction = "none")
}

# Save objects so that the above commands don't need to be re-run:
# saveRDS(object = balance_almeida_cov25_reps,
#         file = "intermediate_objects/balance_almeida_cov25_reps_random.rds")
# 
# saveRDS(object = list("almeida" = wilcoxon_p_almeida_cov25, "ERP003612" = wilcoxon_p_ERP003612_cov25, "SRP095580" = wilcoxon_p_SRP095580_cov25),
#         file = "intermediate_objects/wilcoxon_cov25_reps_random.rds")

wilcoxon_p_almeida_cov25 <- readRDS(file = "intermediate_objects/wilcoxon_cov25_reps_random.rds")

# Not run to save time:
# aldex2_almeida_cov25 <- mclapply(X = 1:20, FUN = function(x) { run_2group_ALDEx2(almeida_metagenome_cov25, almeida_healthy_random_groups$group1[[x]], almeida_healthy_random_groups$group2[[x]]) }, mc.cores=60)
# corncob_almeida_cov25 <- mclapply(X = 1:20, FUN = function(x) { run_2group_corncob(almeida_metagenome_cov25, almeida_healthy_random_groups$group1[[x]], almeida_healthy_random_groups$group2[[x]]) }, mc.cores=60)


# Plot overall P-value histograms over all replicates.

wilcoxon_p_almeida_cov25_combined <- unlist(wilcoxon_p_almeida_cov25)
wilcoxon_p_almeida_cov25_combined_na_removed <- wilcoxon_p_almeida_cov25_combined[-which(is.na(wilcoxon_p_almeida_cov25_combined))]

# Get raw P values for every gene family in every significant balance.
ALL_balance_almeida_cov25_reps_p <- c()
num_pos_balance_almeida_cov25_reps <- c()
num_neg_balance_almeida_cov25_reps <- c()

for(rep in 1:1000) {
  for(balance in names(balance_almeida_cov25_reps[[rep]]$funcs_per_balance)) {
    ALL_balance_almeida_cov25_reps_p <- c(ALL_balance_almeida_cov25_reps_p, balance_almeida_cov25_reps[[rep]]$funcs_per_balance[[balance]]$P)
  }
  num_pos_balance_almeida_cov25_reps <- c(num_pos_balance_almeida_cov25_reps, balance_almeida_cov25_reps[[rep]]$df$num_sig_balances_pos_enrich)
  num_neg_balance_almeida_cov25_reps <- c(num_neg_balance_almeida_cov25_reps, balance_almeida_cov25_reps[[rep]]$df$num_sig_balances_neg_enrich)
}

# Plotting in base since ggplot2 histograms were being wonky about boundaries (e.g. forcing -0.05 to 0.05 to be a bin so that bin would be really small).
par(mfrow=c(2, 1))

hist(wilcoxon_p_almeida_cov25_combined_na_removed, xlab = "Raw P-value", main="Wilcoxon test", col="grey")
hist(ALL_balance_almeida_cov25_reps_p, xlab = "Raw P-value", main="Post-hoc Fisher's exact test", col="grey", ylim=c(0, 3000000))



nums_balance_almeida_cov25_reps <- data.frame(matrix(0, nrow=21*21, ncol=3))
colnames(nums_balance_almeida_cov25_reps) <- c("pos", "neg", "count")
pos_values <- c()
neg_values <- c()

for(i in 0:20) {
  pos_values <- c(pos_values, rep(i, 21))
  neg_values <- c(neg_values, 0:20)
}

nums_balance_almeida_cov25_reps$pos <- pos_values
nums_balance_almeida_cov25_reps$neg <- neg_values

for(i in 1:length(num_pos_balance_almeida_cov25_reps)) {
  row_i <- which(nums_balance_almeida_cov25_reps$pos == num_pos_balance_almeida_cov25_reps[i] & nums_balance_almeida_cov25_reps$neg == num_neg_balance_almeida_cov25_reps[i])
  nums_balance_almeida_cov25_reps[row_i, "count"] = nums_balance_almeida_cov25_reps[row_i, "count"] + 1
}

nums_balance_almeida_cov25_reps_noZero <- nums_balance_almeida_cov25_reps
nums_balance_almeida_cov25_reps_noZero[1, "count"] <- 0

nums_balance_almeida_cov25_reps_noZero$log_count <- log(nums_balance_almeida_cov25_reps_noZero$count + 1)

nums_balance_almeida_cov25_reps_noZero <- nums_balance_almeida_cov25_reps_noZero[-which(nums_balance_almeida_cov25_reps_noZero$pos > 14), ]
nums_balance_almeida_cov25_reps_noZero <- nums_balance_almeida_cov25_reps_noZero[-which(nums_balance_almeida_cov25_reps_noZero$neg > 14), ]

ggplot(data = nums_balance_almeida_cov25_reps_noZero, aes(x=neg, y=pos, fill=log_count)) +
       geom_tile() +
       xlab("# balances where KO is negatively associated with group 1") +
       ylab("# balances where KO is positively associated with group 1") +
       labs(fill="log(count)")

nums_balance_almeida_cov25_reps_noZero_pos7_neg2 <- nums_balance_almeida_cov25_reps_noZero[which(nums_balance_almeida_cov25_reps_noZero$pos >= 7 & nums_balance_almeida_cov25_reps_noZero$neg <= 2 ), ]
nums_balance_almeida_cov25_reps_noZero_neg7_pos2 <- nums_balance_almeida_cov25_reps_noZero[which(nums_balance_almeida_cov25_reps_noZero$neg >= 7 & nums_balance_almeida_cov25_reps_noZero$pos <= 2 ), ]

sum(nums_balance_almeida_cov25_reps_noZero_pos7_neg2$count)
#[1] 381
sum(nums_balance_almeida_cov25_reps_noZero_neg7_pos2$count)
#[1] 246
sum(nums_balance_almeida_cov25_reps_noZero$count)
#[1] 1566661

#(381/1566661)*100
#0.02431924

#(246/1566661)*100
#0.01570218

### Very low proportion of functiosn categorized as strong signals (<< 1%).

