### First needed to identify focal KOs that would be representative of different POMS signals.
### Did this based on finding focal KOs that consistently showed patterns of low, medium, and high proportions of sig. KOs across all four settings.
### This was looked at based on the ranks of each focal gene for each setting based on the proportion of sig. KOs

rm(list=ls(all.names=TRUE))

devtools::load_all(path = "/home/gavin/github_repos/POMS/")
source("/home/gavin/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")
setwd("/home/gavin/projects/POMS/reference_genome_sim/")

BEZI_tables <- list()

BEZI_tables[["mu0.1_nu_0.5"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.5_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.9"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.9_sigma1.tsv",
                                            header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.1_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.1_nu_0.99_sigma1.tsv",
                                             header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

BEZI_tables[["mu0.01_nu_0.99"]] <- read.table(file = "relabun_tables/BEZI_genome_abun_mu0.01_nu_0.99_sigma1.tsv",
                                              header=TRUE, sep="\t", check.names=FALSE, row.names=1, quote="", comment.char="")

POMS_out <- readRDS(file = "func_rand_output_POMS_100reps_BEZI_cutoffs.rds")



POMS_out_summary <- list()

for(cutoff in names(BEZI_tables)) {
  
  kos_selected <- sapply(POMS_out[[cutoff]], function(x){ x$func })
  
  POMS_out_summary[[cutoff]] <- list()
  POMS_out_summary[[cutoff]][["sig_func"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich + x$output$df$num_sig_nodes_neg_enrich) > 0)) })
  POMS_out_summary[[cutoff]][["sig_bal"]] <- sapply(POMS_out[[cutoff]], function(x) { max(x$output$df$num_sig_nodes_pos_enrich + x$output$df$num_sig_nodes_neg_enrich + x$output$df$num_sig_nodes_nonenrich) })
  POMS_out_summary[[cutoff]][["pos3_neg0"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich >= 3 & x$output$df$num_sig_nodes_neg_enrich == 0))) })
  POMS_out_summary[[cutoff]][["neg3_pos0"]] <- sapply(POMS_out[[cutoff]], function(x) { length(which((x$output$df$num_sig_nodes_pos_enrich == 0 & x$output$df$num_sig_nodes_neg_enrich >= 3))) })
  POMS_out_summary[[cutoff]][["sig3_other0"]] <- POMS_out_summary[[cutoff]][["pos3_neg0"]] + POMS_out_summary[[cutoff]][["neg3_pos0"]]
  
  names(POMS_out_summary[[cutoff]][["sig3_other0"]]) <- kos_selected
  
  POMS_out_summary[[cutoff]][["sig3_other0_rank"]] <- POMS_out_summary[[cutoff]][["sig3_other0"]]
  POMS_out_summary[[cutoff]][["sig3_other0_rank"]] <- rank(POMS_out_summary[[cutoff]][["sig3_other0_rank"]])
}



early_kos <- names(which(POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank < 20))
mid_kos <- names(which(POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank > 35 & POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank < 65))
high_kos <- names(which(POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank > 80))

POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank[early_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0_rank[early_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0_rank[early_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0_rank[early_kos]

POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0[early_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0[early_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0[early_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0[early_kos]


POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank[mid_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0_rank[mid_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0_rank[mid_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0_rank[mid_kos]

POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0[mid_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0[mid_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0[mid_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0[mid_kos]



POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0_rank[high_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0_rank[high_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0_rank[high_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0_rank[high_kos]

POMS_out_summary[["mu0.1_nu_0.5"]]$sig3_other0[high_kos]
POMS_out_summary[["mu0.1_nu_0.9"]]$sig3_other0[high_kos]
POMS_out_summary[["mu0.1_nu_0.99"]]$sig3_other0[high_kos]
POMS_out_summary[["mu0.01_nu_0.99"]]$sig3_other0[high_kos]


# low prop sig KO example: K00480. Ranks: 12, 10.5, 9, 6; Num. sig: 221, 193, 208, 224
# 
# mid prop sig KO example: K01845. Ranks: 57, 52, 51, 52; Num sig: 1407, 1393, 1396, 1366
# 
# high prop sig KO example: K06077. Ranks: 92, 87, 86, 87; Num sig: 3416, 3187, 3132, 3136

