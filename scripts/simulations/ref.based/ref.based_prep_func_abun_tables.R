rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations//")

source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in needed prepped files.
func_subsets <- readRDS("ref.based_prepped_func_tables.rds")

for (cutoff in names(func_subsets)) {
  
  ptm <- proc.time()
  
  func_abun_rand_func_void <- mclapply(X = 1:500, FUN = function(rep_i) {
    
    in_rds <- paste("func_sim_info_sel1.5/", "func_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = "")
    
    rep_sim_info <- readRDS(in_rds)
    
    rep_func_abun <- calc_func_abun(in_abun = rep_sim_info$taxa_perturb_abun, in_func = func_subsets[[cutoff]], ncores = 1)
    
    outfile <- paste("func_abun_tables_func_sim_sel1.5/", "func_abun_tab_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = "")
    
    saveRDS(file = outfile, object = rep_func_abun)
    
  }
  , mc.cores = 30)
  
  print(proc.time() - ptm)
  
}

for (cutoff in names(func_subsets)) {
  
  ptm <- proc.time()
  
  func_abun_rand_taxa_void <- mclapply(X = 1:500, FUN = function(rep_i) {
    
    in_rds <- paste("taxa_sim_info_sel1.5/", "taxa_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = "")
    
    rep_sim_info <- readRDS(in_rds)
    
    rep_func_abun <- calc_func_abun(in_abun = rep_sim_info$taxa_perturb_abun, in_func = func_subsets[[cutoff]], ncores = 1)
    
    outfile <- paste("func_abun_tables_taxa_sim_sel1.5/", "func_abun_tab_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = "")
    
    saveRDS(file = outfile, object = rep_func_abun)
    
  }
  , mc.cores = 30)
  
  print(proc.time() - ptm)
  
}

