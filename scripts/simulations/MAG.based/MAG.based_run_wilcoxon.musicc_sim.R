rm(list = ls(all.names = TRUE))

# Code for running MUSiCC-normalized Wilcoxon tests only.

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

source("~/github_repos/POMS_manuscript/scripts/POMS_manuscript_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

ptm <- proc.time()

func_sim_alt.tools <- mclapply(X = 1:1000, FUN = function(rep_i) {

  prepped_func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))

  prepped_func_abun <- readRDS(paste("func_abun_tables_func_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))

  wilcoxon.musicc_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                       group1_samples = prepped_func_sim_info$group1,
                                       group2_samples = prepped_func_sim_info$group2,
                                       USCGs = musicc_uscgs,
                                       tools_to_run = "wilcoxon.musicc")

  wilcoxon.musicc_out$func <- prepped_func_sim_info$func
  
  return(wilcoxon.musicc_out)
}
, mc.cores = 30)

saveRDS(object = func_sim_alt.tools, file = "MAG.based_wilcoxon.musicc_sim_rand_func_1000reps.rds")

print(proc.time() - ptm)




ptm <- proc.time()

taxa_sim_alt.tools <- mclapply(X = 1:1000, FUN = function(rep_i) {
  
  prepped_taxa_sim_info <- readRDS(paste("MAG.based_prepped_taxa_sim_info_sel1.5/taxa_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  prepped_func_abun <- readRDS(paste("func_abun_tables_taxa_sim_sel1.5/func_abun_tab_rep", as.character(rep_i), ".rds", sep = ""))
  
  wilcoxon.musicc_out <- run_alt.tools(func_abun_table = prepped_func_abun,
                                       group1_samples = prepped_taxa_sim_info$group1,
                                       group2_samples = prepped_taxa_sim_info$group2,
                                       USCGs = musicc_uscgs,
                                       tools_to_run = "wilcoxon.musicc")
  
  wilcoxon.musicc_out$func <- prepped_taxa_sim_info$func
  
  return(wilcoxon.musicc_out)
}
, mc.cores = 30)

saveRDS(object = taxa_sim_alt.tools, file = "MAG.based_wilcoxon.musicc_sim_rand_taxa_1000reps.rds")

print(proc.time() - ptm)
