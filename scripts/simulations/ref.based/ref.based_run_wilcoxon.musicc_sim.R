rm(list = ls(all.names = TRUE))

# Code for running MUSiCC-normalized Wilcoxon tests only.

setwd("~/github_repos/POMS_manuscript/data/intermediates/ref.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

ptm <- proc.time()

cutoffs <- c("nu0.50", "nu0.65", "nu0.80", "nu0.95")

for (cutoff in cutoffs) {
  
  ptm <- proc.time()
  
  func_rand_wilcoxon.musicc_cutoff <- mclapply(X = 1:500, FUN = function(rep_i) {
    
    rep_func_sim_info <- readRDS(paste("func_sim_info_sel1.5/func_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))

    prepped_func_abun <- readRDS(paste("func_abun_tables_func_sim_sel1.5/func_abun_tab_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    wilcoxon.musicc_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                         group1_samples = rep_func_sim_info$group1,
                                         group2_samples = rep_func_sim_info$group2,
                                         USCGs = musicc_uscgs,
                                         tools_to_run = "wilcoxon.musicc")
    
    wilcoxon.musicc_out$func <- rep_func_sim_info$func
    
    return(wilcoxon.musicc_out)
    
  } , mc.cores = 30)
  
  saveRDS(object = func_rand_wilcoxon.musicc_cutoff,
          file = paste("wilcoxon.musicc_out/ref.based_wilcoxon.musicc_sim_rand_func_500reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}



for (cutoff in cutoffs) {
  
  ptm <- proc.time()
  
  taxa_rand_wilcoxon.musicc_cutoff <- mclapply(X = 1:500, FUN = function(rep_i) {
    
    rep_taxa_sim_info <- readRDS(paste("taxa_sim_info_sel1.5/taxa_sim_info_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    prepped_func_abun <- readRDS(paste("func_abun_tables_taxa_sim_sel1.5/func_abun_tab_cutoff_", cutoff, "_rep", as.character(rep_i), ".rds", sep = ""))
    
    wilcoxon.musicc_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                         group1_samples = rep_taxa_sim_info$group1,
                                         group2_samples = rep_taxa_sim_info$group2,
                                         USCGs = musicc_uscgs,
                                         tools_to_run = "wilcoxon.musicc")
    
    wilcoxon.musicc_out$func <- rep_taxa_sim_info$func
    
    return(wilcoxon.musicc_out)
    
    
  } , mc.cores = 30)
  
  saveRDS(object = taxa_rand_wilcoxon.musicc_cutoff,
          file = paste("wilcoxon.musicc_out/ref.based_wilcoxon.musicc_sim_rand_taxa_500reps_cutoff_", as.character(cutoff), ".rds", sep = ""))
  
  print(proc.time() - ptm)
  
}
