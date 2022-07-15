rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

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
    
    output <- POMS_pipeline(abun = rep_func_sim_info$taxa_perturb_abun,
                            func = ref_func_cutoff,
                            tree = tree_cutoff,
                            group1_samples = rep_func_sim_info$group1,
                            group2_samples = rep_func_sim_info$group2,
                            ncores = 1,
                            BSN_p_cutoff = 0.05,
                            BSN_correction = "none",
                            FSN_p_cutoff = 0.05,
                            FSN_correction = "none",
                            min_func_instances = 5,
                            min_func_prop = 0.001,
                            multinomial_correction = "BH",
                            detailed_output = FALSE,
                            verbose = FALSE)
    
    return(list(func = rep_func_sim_info$func, output = output))
  }
  , mc.cores = 20)
  
  saveRDS(object = func_rand_POMS_cutoff,
          file = paste("POMS_out/ref.based_POMS_sim_func.based_1000reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}
