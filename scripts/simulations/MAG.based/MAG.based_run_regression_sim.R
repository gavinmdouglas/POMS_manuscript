rm(list = ls(all.names = TRUE))

# Run regression against functional content based on taxonomic abundances.
# Based on: (1) which taxa are significant based on Wilcoxon relabun and (2) the taxonomic specificities (ala phylogenize).

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")
source("../../../scripts/alt_tool_functions.R")

library(ape)
library(parallel)

# Read in needed prepped files.
almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")

almeida_func_subset <- almeida_func_subset[almeida_tree_subset$tip.label, ]

ptm <- proc.time()

clade.based_sim_regression <-  mclapply(X = 1:693, FUN = function(rep_i) {
 
  clade.based_info <- readRDS(paste("MAG.based_prepped_clade.based_sim_info_sel1.5/", "clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = clade.based_info$taxa_perturb_abun,
                                               group1_samples = clade.based_info$group1,
                                               group2_samples = clade.based_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[almeida_tree_subset$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  almeida_func_subset,
                                                       in_tree = almeida_tree_subset,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(clade.based_info$group1,
                                      clade.based_info$group2),
                             group = c(rep("group1", length(clade.based_info$group1)),
                                       rep("group2", length(clade.based_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = clade.based_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[almeida_tree_subset$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  almeida_func_subset,
                                                          in_tree = almeida_tree_subset,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  return(list(sig_taxa_regress = sig_taxa_regress_out,
              specificity_regress = specificity_regress_out))
   
  }, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = clade.based_sim_regression, file = "regression_out/MAG.based_regression_sim_rand_clade.based_693reps.rds")




ptm <- proc.time()

func_sim_regression <-  mclapply(X = 1:1000, FUN = function(rep_i) {
  
  func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = func_sim_info$taxa_perturb_abun,
                                               group1_samples = func_sim_info$group1,
                                               group2_samples = func_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[almeida_tree_subset$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  almeida_func_subset,
                                                       in_tree = almeida_tree_subset,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(func_sim_info$group1,
                                      func_sim_info$group2),
                             group = c(rep("group1", length(func_sim_info$group1)),
                                       rep("group2", length(func_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = func_sim_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[almeida_tree_subset$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  almeida_func_subset,
                                                          in_tree = almeida_tree_subset,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  return(list(sig_taxa_regress = sig_taxa_regress_out,
              specificity_regress = specificity_regress_out))
  
}, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = func_sim_regression, file = "regression_out/MAG.based_regression_sim_rand_func_1000reps.rds")

rm(func_sim_regression)


ptm <- proc.time()

taxa_sim_regression <-  mclapply(X = 1:1000, FUN = function(rep_i) {
  
  taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = taxa_sim_info$taxa_perturb_abun,
                                               group1_samples = taxa_sim_info$group1,
                                               group2_samples = taxa_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[almeida_tree_subset$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  almeida_func_subset,
                                                       in_tree = almeida_tree_subset,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(taxa_sim_info$group1,
                                      taxa_sim_info$group2),
                             group = c(rep("group1", length(taxa_sim_info$group1)),
                                       rep("group2", length(taxa_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = taxa_sim_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[almeida_tree_subset$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  almeida_func_subset,
                                                          in_tree = almeida_tree_subset,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  return(list(sig_taxa_regress = sig_taxa_regress_out,
              specificity_regress = specificity_regress_out))
  
  }, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = taxa_sim_regression, file = "regression_out/MAG.based_regression_sim_rand_taxa_1000reps.rds")

# Remove regression output object.
rm(taxa_sim_regression)


# Finally, run unperturbed simulations:
almeida_abun_subset <- readRDS("MAG.based_prepped_abun.rds")

ptm <- proc.time()

unperturbed_sim_regression <-  mclapply(X = 1:1000, FUN = function(rep_i) {
  
  taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = almeida_abun_subset,
                                               group1_samples = taxa_sim_info$group1,
                                               group2_samples = taxa_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[almeida_tree_subset$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  almeida_func_subset,
                                                       in_tree = almeida_tree_subset,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(taxa_sim_info$group1,
                                      taxa_sim_info$group2),
                             group = c(rep("group1", length(taxa_sim_info$group1)),
                                       rep("group2", length(taxa_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = almeida_abun_subset,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[almeida_tree_subset$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  almeida_func_subset,
                                                          in_tree = almeida_tree_subset,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  return(list(sig_taxa_regress = sig_taxa_regress_out,
              specificity_regress = specificity_regress_out))
}
, mc.cores = 20)

print(proc.time() - ptm)

saveRDS(object = unperturbed_sim_regression, file = "regression_out/MAG.based_regression_sim_unperturbed_1000reps.rds")
