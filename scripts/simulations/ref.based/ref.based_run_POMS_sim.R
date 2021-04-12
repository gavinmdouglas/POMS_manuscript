rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)

cutoffs <- c("nu0.50", "nu0.65", "nu0.80", "nu0.95")

for (cutoff in cutoffs) {
  
  ptm <- proc.time()

  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]
  
  tree_cutoff <- readRDS("ref.based_prepped_trees.rds")[[cutoff]]
  
  func_rand_POMS_cutoff <- mclapply(X = 1:1000, FUN = function(rep_i) {

    rep_func_sim_info <- readRDS(paste("func_sim_info_sel1.5/func_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    output <- POMS_pipeline(abun = rep_func_sim_info$taxa_perturb_abun,
                            func = ref_func_cutoff,
                            phylogeny = tree_cutoff,
                            group1_samples = rep_func_sim_info$group1,
                            group2_samples = rep_func_sim_info$group2,
                            ncores = 1,
                            balance_p_cutoff = 0.05,
                            balance_correction = "none",
                            function_p_cutoff = 0.05,
                            function_correction = "none",
                            min_func_instances = 0,
                            min_func_prop = 0,
                            run_multinomial_test = TRUE,
                            multinomial_correction = "BH",
                            calc_node_dist = FALSE,
                            detailed_output = FALSE,
                            verbose = FALSE)
    
    return(list(func = rep_func_sim_info$func, output = output))
  }
  , mc.cores = 50)
  
  saveRDS(object = func_rand_POMS_cutoff,
          file = paste("POMS_out/ref.based_POMS_sim_rand_func_1000reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}



for (cutoff in cutoffs) {
  
  ptm <- proc.time()
  
  ref_func_cutoff <- readRDS("ref.based_prepped_func_tables.rds")[[cutoff]]
  
  tree_cutoff <- readRDS("ref.based_prepped_trees.rds")[[cutoff]]
  
  func_rand_POMS_cutoff <- mclapply(X = 1:1000, FUN = function(rep_i) {
    
    rep_taxa_sim_info <- readRDS(paste("taxa_sim_info_sel1.5/taxa_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    output <- POMS_pipeline(abun = rep_taxa_sim_info$taxa_perturb_abun,
                            func = ref_func_cutoff,
                            phylogeny = tree_cutoff,
                            group1_samples = rep_taxa_sim_info$group1,
                            group2_samples = rep_taxa_sim_info$group2,
                            ncores = 1,
                            balance_p_cutoff = 0.05,
                            balance_correction = "none",
                            function_p_cutoff = 0.05,
                            function_correction = "none",
                            min_func_instances = 0,
                            min_func_prop = 0,
                            run_multinomial_test = TRUE,
                            multinomial_correction = "BH",
                            calc_node_dist = FALSE,
                            detailed_output = FALSE,
                            verbose = FALSE)
    
    return(list(func = rep_taxa_sim_info$func, output = output))
  }
  , mc.cores = 30)
  
  saveRDS(object = func_rand_POMS_cutoff,
          file = paste("POMS_out/ref.based_POMS_sim_rand_taxa_1000reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}

