rm(list = ls(all.names = TRUE))

# Code for running MUSiCC-normalized Wilcoxon tests only.

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

source("~/github_repos/POMS_manuscript/scripts/alt_tool_functions.R")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(ape)
library(parallel)

musicc_uscgs <- read.table("/home/gavin/github_repos/POMS_manuscript/data/key_inputs/MUSiCC_KEGG_single_copy_genes.txt.gz",
                           stringsAsFactors = FALSE)$V1

ptm <- proc.time()

parameter_settings <- list()

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

option_count <- 1

for (rep_i in 1:num_reps) {

  for (MAG_num in MAG_nums) {
    
    MAG_rep_func <- readRDS(paste("prepped_func_tables/subset_",
                                  as.character(MAG_num), "MAGs_func_rep",
                                  as.character(rep_i), ".rds", sep = ""))
    
    for (pseudocount_set in pseudocount_settings) {
      for (abun_increase_set in abun_increase_settings) {
        
        parameter_settings[[option_count]] <- list()
        parameter_settings[[option_count]][["rep_i"]] <- rep_i
        parameter_settings[[option_count]][["MAG_num"]] <- MAG_num
        parameter_settings[[option_count]][["pseudocount_set"]] <- pseudocount_set
        parameter_settings[[option_count]][["abun_increase_set"]] <- abun_increase_set
        parameter_settings[[option_count]][["MAG_rep_func"]] <- MAG_rep_func
        
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


func.based_out <- mclapply(X = parameter_settings, FUN = function(x) {

                     rep_i <- x$rep_i
                     MAG_num <- x$MAG_num
                     pseudocount_set <- x$pseudocount_set
                     abun_increase_set <- x$abun_increase_set
                    
                     func_sim_info <- readRDS(paste("sim_info/func.based/func.based_sim_info_",
                                                    "rep", as.character(rep_i),
                                                    "_MAGs", as.character(MAG_num),
                                                    "_pseudo", as.character(pseudocount_set),
                                                    "_increase", as.character(abun_increase_set),
                                                    ".rds", sep = ""))
                     
                     rep_func_abun <- readRDS(paste("func_abun_tables/func.based/crossprod_abun_",
                                                    "rep", as.character(rep_i),
                                                    "_MAGs", as.character(MAG_num),
                                                    "_pseudo", as.character(pseudocount_set),
                                                    "_increase", as.character(abun_increase_set),
                                                    ".rds", sep = ""))
 
                     wilcoxon.musicc_out <- run_alt.tools(func_abun_table = rep_func_abun,
                                                          group1_samples = func_sim_info$group1,
                                                          group2_samples = func_sim_info$group2,
                                                          USCGs = musicc_uscgs,
                                                          tools_to_run = "wilcoxon.musicc")
                     
                     wilcoxon.musicc_out$func = func_sim_info$func
                     
                     saveRDS(object = wilcoxon.musicc_out,
                             file = paste("wilcoxon_out/func.based/wilcoxon_func.based_out_",
                                          "rep", as.character(rep_i),
                                          "_MAGs", as.character(MAG_num),
                                          "_pseudo", as.character(pseudocount_set),
                                          "_increase", as.character(abun_increase_set),
                                          ".rds", sep = ""))
  
                     return("Success")

  }, mc.cores = 30)





taxa.based_out <- mclapply(X = parameter_settings, FUN = function(x) {
  
  rep_i <- x$rep_i
  MAG_num <- x$MAG_num
  pseudocount_set <- x$pseudocount_set
  abun_increase_set <- x$abun_increase_set
  
  func_sim_info <- readRDS(paste("sim_info/taxa.based/taxa.based_sim_info_",
                                 "rep", as.character(rep_i),
                                 "_MAGs", as.character(MAG_num),
                                 "_pseudo", as.character(pseudocount_set),
                                 "_increase", as.character(abun_increase_set),
                                 ".rds", sep = ""))
  
  rep_func_abun <- readRDS(paste("func_abun_tables/taxa.based/crossprod_abun_",
                                 "rep", as.character(rep_i),
                                 "_MAGs", as.character(MAG_num),
                                 "_pseudo", as.character(pseudocount_set),
                                 "_increase", as.character(abun_increase_set),
                                 ".rds", sep = ""))
  
  wilcoxon.musicc_out <- run_alt.tools(func_abun_table = rep_func_abun,
                                       group1_samples = func_sim_info$group1,
                                       group2_samples = func_sim_info$group2,
                                       USCGs = musicc_uscgs,
                                       tools_to_run = "wilcoxon.musicc")
  
  wilcoxon.musicc_out$func = func_sim_info$func
  
  saveRDS(object = wilcoxon.musicc_out,
          file = paste("wilcoxon_out/taxa.based/wilcoxon_taxa.based_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")

}, mc.cores = 30)






clade.based_out <- mclapply(X = parameter_settings, FUN = function(x) {
  
  rep_i <- x$rep_i
  MAG_num <- x$MAG_num
  pseudocount_set <- x$pseudocount_set
  abun_increase_set <- x$abun_increase_set
  
  func_sim_info <- readRDS(paste("sim_info/clade.based/clade.based_sim_info_",
                                 "rep", as.character(rep_i),
                                 "_MAGs", as.character(MAG_num),
                                 "_pseudo", as.character(pseudocount_set),
                                 "_increase", as.character(abun_increase_set),
                                 ".rds", sep = ""))
  
  rep_func_abun <- readRDS(paste("func_abun_tables/clade.based/crossprod_abun_",
                                 "rep", as.character(rep_i),
                                 "_MAGs", as.character(MAG_num),
                                 "_pseudo", as.character(pseudocount_set),
                                 "_increase", as.character(abun_increase_set),
                                 ".rds", sep = ""))
  
  wilcoxon.musicc_out <- run_alt.tools(func_abun_table = rep_func_abun,
                                       group1_samples = func_sim_info$group1,
                                       group2_samples = func_sim_info$group2,
                                       USCGs = musicc_uscgs,
                                       tools_to_run = "wilcoxon.musicc")
  
  wilcoxon.musicc_out$func = func_sim_info$func
  
  saveRDS(object = wilcoxon.musicc_out,
          file = paste("wilcoxon_out/clade.based/wilcoxon_clade.based_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")

}, mc.cores = 30)


