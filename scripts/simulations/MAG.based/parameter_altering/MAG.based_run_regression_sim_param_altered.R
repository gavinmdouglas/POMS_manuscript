rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

source("../../../../scripts/alt_tool_functions.R")

library(ape)
library(parallel)

parameter_settings <- list()

MAG_nums <- c(1595, 1250, 1000, 750, 500, 250, 100, 50)

pseudocount_settings <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

abun_increase_settings <- c(1.5, 1.3, 1.1, 1.05)

option_count <- 1

for (rep_i in 1:25) {
  for (MAG_num in MAG_nums) {
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    MAG_rep_tree <- readRDS(paste("prepped_trees/subset_",
                                  as.character(MAG_num), "MAGs_tree_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      for (abun_increase_set in abun_increase_settings) {
        
        parameter_settings[[option_count]] <- list()
        parameter_settings[[option_count]][["rep_i"]] <- rep_i
        parameter_settings[[option_count]][["MAG_num"]] <- MAG_num
        parameter_settings[[option_count]][["pseudocount_set"]] <- pseudocount_set
        parameter_settings[[option_count]][["abun_increase_set"]] <- abun_increase_set
        parameter_settings[[option_count]][["MAG_rep_func"]] <- MAG_rep_func
        parameter_settings[[option_count]][["MAG_rep_tree"]] <- MAG_rep_tree
        
        
        option_count <- option_count + 1
      }
    }
  }
}

rm(rep_i)
rm(MAG_num)
rm(pseudocount_set)
rm(abun_increase_set)
rm(MAG_rep_func)
rm(MAG_rep_tree)


func.based_sim_null_out <- mclapply(X = parameter_settings, FUN = function(x) {
  
  rep_i <- x$rep_i
  MAG_num <- x$MAG_num
  pseudocount_set <- x$pseudocount_set
  abun_increase_set <- x$abun_increase_set
  MAG_rep_func <- x$MAG_rep_func
  MAG_rep_tree <- x$MAG_rep_tree
  
  MAG_rep_func <- MAG_rep_func[MAG_rep_tree$tip.label, ]
  
  func.based_sim_info <- readRDS(paste("sim_info/func.based/func_sim_info_",
                                       "rep", as.character(rep_i),
                                       "_MAGs", as.character(MAG_num),
                                       "_pseudo", as.character(pseudocount_set),
                                       "_increase", as.character(abun_increase_set),
                                       ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = func.based_sim_info$taxa_perturb_abun,
                                               group1_samples = func.based_sim_info$group1,
                                               group2_samples = func.based_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[MAG_rep_tree$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  MAG_rep_func,
                                                       in_tree = MAG_rep_tree,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(func.based_sim_info$group1,
                                      func.based_sim_info$group2),
                             group = c(rep("group1", length(func.based_sim_info$group1)),
                                       rep("group2", length(func.based_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = func.based_sim_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[MAG_rep_tree$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  MAG_rep_func,
                                                          in_tree = MAG_rep_tree,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  output <- list(sig_taxa_regress = sig_taxa_regress_out,
                 specificity_regress = specificity_regress_out)
  
  saveRDS(object = output,
          file = paste("regress_out/func.based/regress_func.based_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")
  
}, mc.cores = 20)



taxa.based_sim_null_out <- mclapply(X = parameter_settings, FUN = function(x) {
  
  rep_i <- x$rep_i
  MAG_num <- x$MAG_num
  pseudocount_set <- x$pseudocount_set
  abun_increase_set <- x$abun_increase_set
  MAG_rep_func <- x$MAG_rep_func
  MAG_rep_tree <- x$MAG_rep_tree
  
  MAG_rep_func <- MAG_rep_func[MAG_rep_tree$tip.label, ]
  
  taxa.based_sim_info <- readRDS(paste("sim_info/taxa.based/taxa.based_sim_info_",
                                       "rep", as.character(rep_i),
                                       "_MAGs", as.character(MAG_num),
                                       "_pseudo", as.character(pseudocount_set),
                                       "_increase", as.character(abun_increase_set),
                                       ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = taxa.based_sim_info$taxa_perturb_abun,
                                               group1_samples = taxa.based_sim_info$group1,
                                               group2_samples = taxa.based_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[MAG_rep_tree$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  MAG_rep_func,
                                                       in_tree = MAG_rep_tree,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(taxa.based_sim_info$group1,
                                      taxa.based_sim_info$group2),
                             group = c(rep("group1", length(taxa.based_sim_info$group1)),
                                       rep("group2", length(taxa.based_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = taxa.based_sim_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[MAG_rep_tree$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  MAG_rep_func,
                                                          in_tree = MAG_rep_tree,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  output <- list(sig_taxa_regress = sig_taxa_regress_out,
                 specificity_regress = specificity_regress_out)
  
  saveRDS(object = output,
          file = paste("regress_out/taxa.based/regress_taxa.based_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")
  
}, mc.cores = 20)

clade.based_sim_null_out <- mclapply(X = parameter_settings, FUN = function(x) {
  
  rep_i <- x$rep_i
  MAG_num <- x$MAG_num
  pseudocount_set <- x$pseudocount_set
  abun_increase_set <- x$abun_increase_set
  MAG_rep_func <- x$MAG_rep_func
  
  MAG_rep_func <- MAG_rep_func[clade.based_sim_info$subset_tree$tip.label, ]
  
  clade.based_sim_info <- readRDS(paste("sim_info/clade.based/clade.based_sim_info_",
                                       "rep", as.character(rep_i),
                                       "_MAGs", as.character(MAG_num),
                                       "_pseudo", as.character(pseudocount_set),
                                       "_increase", as.character(abun_increase_set),
                                       ".rds", sep = ""))
  
  rep_wilcox_output <- wilcoxon_2group_pvalues(intable = clade.based_sim_info$taxa_perturb_abun,
                                               group1_samples = clade.based_sim_info$group1,
                                               group2_samples = clade.based_sim_info$group2,
                                               convert_relab = TRUE)
  
  rep_wilcox_output <- rep_wilcox_output[clade.based_sim_info$subset_tree$tip.label, ]
  
  sig_taxa <- rep(0, nrow(rep_wilcox_output))
  sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
  
  sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                       func =  MAG_rep_func,
                                                       in_tree = clade.based_sim_info$subset_tree,
                                                       ncores = 1,
                                                       model_type = "BM")
  
  sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
  
  rep_metadata <- data.frame(samp = c(clade.based_sim_info$group1,
                                      clade.based_sim_info$group2),
                             group = c(rep("group1", length(clade.based_sim_info$group1)),
                                       rep("group2", length(clade.based_sim_info$group2))))
  
  taxa_specificity <- specificity_scores(abun_table = clade.based_sim_info$taxa_perturb_abun,
                                         meta_table = rep_metadata,
                                         focal_var_level = "group1",
                                         var_colname = "group",
                                         sample_colname = "samp",
                                         silence_citation = TRUE)
  
  taxa_specificity <- taxa_specificity$ess[clade.based_sim_info$subset_tree$tip.label]
  
  specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                          func =  MAG_rep_func,
                                                          in_tree = clade.based_sim_info$subset_tree,
                                                          ncores = 1,
                                                          model_type = "BM")
  
  specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
  
  output <- list(sig_taxa_regress = sig_taxa_regress_out,
                 specificity_regress = specificity_regress_out)
  
  saveRDS(object = output,
          file = paste("regress_out/clade.based/regress_clade.based_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")
  
}, mc.cores = 20)
