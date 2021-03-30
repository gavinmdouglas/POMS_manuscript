rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

# Read in needed prepped files.
almeida_func_subset <- readRDS("MAG.based_prepped_func.rds")
almeida_tree_subset <- readRDS("MAG.based_prepped_tree.rds")


ptm <- proc.time()

mclapply(X = 1:1000, FUN = function(rep_i) {
  
  in_rds <- paste("MAG.based_prepped_func_sim_info_sel1.5/", "func_sim_info_rep", as.character(rep_i), ".rds", sep = "")
  
  rep_func_sim_info <- readRDS(in_rds)
  
  rep_func_abun <- calc_func_abun(in_abun = rep_func_sim_info$taxa_perturb_abun, in_func = almeida_func_subset, ncores = 1)

  outfile <- paste("func_abun_tables_func_sim_sel1.5/", "func_abun_tab_rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(file = outfile, object = rep_func_abun)
  
}
, mc.cores = 30)

print(proc.time() - ptm)


ptm <- proc.time()

mclapply(X = 1:1000, FUN = function(rep_i) {
  
  in_rds <- paste("MAG.based_prepped_taxa_sim_info_sel1.5/", "taxa_sim_info_rep", as.character(rep_i), ".rds", sep = "")
  
  rep_taxa_sim_info <- readRDS(in_rds)
  
  rep_func_abun <- calc_func_abun(in_abun = rep_taxa_sim_info$taxa_perturb_abun, in_func = almeida_func_subset, ncores = 1)
  
  outfile <- paste("func_abun_tables_taxa_sim_sel1.5/", "func_abun_tab_rep", as.character(rep_i), ".rds", sep = "")
  
  saveRDS(file = outfile, object = rep_func_abun)
}
, mc.cores = 30)

print(proc.time() - ptm)


