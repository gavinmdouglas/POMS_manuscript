### Compare results identified when jittering an input table when
### run through basic Wilcoxon test on functions and with POMS.

rm(list=ls(all.names=TRUE))

setwd("/home/gavin/projects/POMS/reference_genome_sim/")
source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

library(cowplot)
library(ggExtra)
library(ggplot2)

# Run wilcoxon tests on all of these simulated tables for comparison.
ref_func <- read.table(file = "/home/gavin/github_repos/picrust_repos/picrust2/picrust2/default_files/prokaryotic/ko.txt.gz",
                       sep="\t", header=TRUE, row.names=1)

func_sim_info_100rep <- readRDS("func_sim_info_100rep.rds")


POMS_out <- readRDS("func_rand_output_POMS_100reps_BEZI_cutoffs.rds")
wilcox_out <- readRDS("func_rand_output_wilcoxon_100reps_BEZI_cutoffs.rds")


BEZI_tables <- list()

BEZI_tables[["mu0.1_nu_0.5"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.5_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.9"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.9_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.99_sigma1.tsv",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.01_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.01_nu_0.99_sigma1.tsv",
                                              header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

cutoff_func <- list()
wilcox_out_BY <- list()
POMS_out_summary <- list()

for(cutoff in names(BEZI_tables)) {
  
  wilcox_out_BY[[cutoff]] <- lapply(wilcox_out[[cutoff]], function(x) { p.adjust(x, "BY") })
  
  cutoff_func[[cutoff]] <- ref_func[rownames(BEZI_tables[[cutoff]]), ]
  
  if(length(which(colSums(cutoff_func[[cutoff]]) == 0)) > 0) {
    cutoff_func[[cutoff]] <- cutoff_func[[cutoff]][, -which(colSums(cutoff_func[[cutoff]]) == 0)]
  }
  
  POMS_out_summary[[cutoff]] <- list()
  
  POMS_out_summary[[cutoff]][["sig_func"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich + x$output$df$num_sig_nodes_neg_enrich) > 0)) })
  POMS_out_summary[[cutoff]][["sig_bal"]] <- sapply(POMS_out[[cutoff]], function(x) { max(x$output$df$num_sig_nodes_pos_enrich + x$output$df$num_sig_nodes_neg_enrich + x$output$df$num_sig_nodes_nonenrich) })
  POMS_out_summary[[cutoff]][["pos3_neg0"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich >= 3 & x$output$df$num_sig_nodes_neg_enrich == 0))) })
  POMS_out_summary[[cutoff]][["neg3_pos0"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich == 0 & x$output$df$num_sig_nodes_neg_enrich >= 3))) })
  POMS_out_summary[[cutoff]][["sig3_other0"]] <- POMS_out_summary[[cutoff]][["pos3_neg0"]] + POMS_out_summary[[cutoff]][["neg3_pos0"]]
  
  
}



# Compute the Jaccard distance in terms of which MAGs significant gene families are found in relative to the perturbed function.
# Also run sanity check that the perturbed function is actually in the significant set.
# The hypothesis is that a higher proportion of gene families in non-perturbed genomes will have been called as significant by Wilcoxon.

POMS_jaccard_summary <- list()
wilcoxon_jaccard_summary <- list()

for(cutoff in names(BEZI_tables)) {
  POMS_jaccard_summary[[cutoff]] <- 1 - POMS_jaccard_wrapper(POMS_output = POMS_out[[cutoff]], func_table = cutoff_func[[cutoff]], num_sig_nodes=3, num_other_nodes=0, num_rep=100)
  wilcoxon_jaccard_summary[[cutoff]] <- 1 - wilcoxon_jaccard_wrapper(POMS_output = POMS_out[[cutoff]], wilcoxon_output=wilcox_out_BY[[cutoff]], func_table = cutoff_func[[cutoff]], p_cutoff=0.0001, num_rep=100)
}


combined_jaccard_summary <- list()
combined_jaccard_summary[["mu0.1_nu_0.5"]] <- data.frame(POMS=POMS_jaccard_summary$mu0.1_nu_0.5,
                                                           wilcoxon=wilcoxon_jaccard_summary$mu0.1_nu_0.5)

combined_jaccard_summary[["mu0.1_nu_0.9"]] <- data.frame(POMS=POMS_jaccard_summary$mu0.1_nu_0.9,
                                                           wilcoxon=wilcoxon_jaccard_summary$mu0.1_nu_0.9)

combined_jaccard_summary[["mu0.1_nu_0.99"]] <- data.frame(POMS=POMS_jaccard_summary$mu0.1_nu_0.99,
                                                            wilcoxon=wilcoxon_jaccard_summary$mu0.1_nu_0.99)

combined_jaccard_summary[["mu0.01_nu_0.99"]] <- data.frame(POMS=POMS_jaccard_summary$mu0.01_nu_0.99,
                                                             wilcoxon=wilcoxon_jaccard_summary$mu0.01_nu_0.99)

### Figure out rank of focal function based on POMS # of enriched nodes and also in terms of Wilcoxon P values.

func_rel_rankings <- list()

for(cutoff in names(BEZI_tables)) {
  func_rel_rankings[[cutoff]] <- POMS_wilcoxon_rank_wrapper(POMS_output = POMS_out[[cutoff]], wilcoxon_output = wilcox_out_BY[[cutoff]],
                                                            func_table = cutoff_func[[cutoff]], num_rep = 100)
}


saveRDS(file = "/home/gavin/github_repos/POMS_manuscript/data/ref.genome_sim_summaries/ref.genome_sim_func.based_ranking.rds", object = func_rel_rankings)
