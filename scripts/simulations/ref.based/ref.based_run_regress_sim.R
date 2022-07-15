rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

source("/home/gavin/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)

cutoffs <- c("nu0.50", "nu0.65", "nu0.80", "nu0.95")

for (cutoff in cutoffs) {
  
  ptm <- proc.time()

  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]
  ref_func_cutoff <- POMS::filter_rare_table_cols(in_tab = ref_func_cutoff, min_nonzero_count = 5, min_nonzero_prop = 0.001)
  all_present_i <- which(colSums(ref_func_cutoff > 0) == nrow(ref_func_cutoff))
  if (length(all_present_i) > 0) {
    ref_func_cutoff <- ref_func_cutoff[, -all_present_i]
  }
  
  tree_cutoff <- readRDS("ref.based_prepped_trees.rds")[[cutoff]]
  tree_cutoff$node.label <- NULL
  tree_cutoff$edge.length <- tree_cutoff$edge.length + 0.001
  tree_cutoff <- POMS::prep_tree(tree_cutoff, tree_cutoff$tip.label)
  
  func_rand_POMS_cutoff <- mclapply(X = 1:1000, FUN = function(rep_i) {

    rep_func_sim_info <- readRDS(paste("func_sim_info_sel1.5/func_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    rep_wilcox_output <- wilcoxon_2group_pvalues(intable = rep_func_sim_info$taxa_perturb_abun,
                                                 group1_samples = rep_func_sim_info$group1,
                                                 group2_samples = rep_func_sim_info$group2,
                                                 convert_relab = TRUE)
    
    rep_wilcox_output <- rep_wilcox_output[tree_cutoff$tip.label, ]
    
    sig_taxa <- rep(0, nrow(rep_wilcox_output))
    sig_taxa[which(rep_wilcox_output$wilcox_p < 0.05)] <- 1
    
    sig_taxa_regress_out <- genome_content_phylo_regress(y = sig_taxa,
                                                         func =  ref_func_cutoff,
                                                         in_tree = tree_cutoff,
                                                         ncores = 1,
                                                         model_type = "BM")
    
    sig_taxa_regress_out$BH <- p.adjust(sig_taxa_regress_out$p, "BH")
    
    rep_metadata <- data.frame(samp = c(rep_func_sim_info$group1,
                                        rep_func_sim_info$group2),
                               group = c(rep("group1", length(rep_func_sim_info$group1)),
                                         rep("group2", length(rep_func_sim_info$group2))))
    
    taxa_specificity <- tryCatch(
      expr = {
        specificity_scores(abun_table = rep_func_sim_info$taxa_perturb_abun,
                           meta_table = rep_metadata,
                           focal_var_level = "group1",
                           var_colname = "group",
                           sample_colname = "samp",
                           silence_citation = TRUE)
      },
      error = function(e) { 
        return("Failed")
      }
    )
    
    if (taxa_specificity[1] != "Failed") {
      
      taxa_specificity <- taxa_specificity$ess[tree_cutoff$tip.label]
      
      specificity_regress_out <- genome_content_phylo_regress(y = taxa_specificity,
                                                              func =  ref_func_cutoff,
                                                              in_tree = tree_cutoff,
                                                              ncores = 1,
                                                              model_type = "BM")
      
      specificity_regress_out$BH <- p.adjust(specificity_regress_out$p, "BH")
    } else {
      specificity_regress_out <- "Failed"
    }
    
    return(list(sig_taxa_regress = sig_taxa_regress_out,
                specificity_regress = specificity_regress_out,
                func = rep_func_sim_info$func))

  }
  , mc.cores = 20)
  
  saveRDS(object = func_rand_POMS_cutoff,
          file = paste("regress_out/ref.based_regress_sim_func.based_1000reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}
