rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/parameter_altered_files/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

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

                     func.based_sim_info <- readRDS(paste("sim_info/func.based/func_sim_info_",
                                                    "rep", as.character(rep_i),
                                                    "_MAGs", as.character(MAG_num),
                                                    "_pseudo", as.character(pseudocount_set),
                                                    "_increase", as.character(abun_increase_set),
                                                    ".rds", sep = ""))
                     
                     output <- POMS_pipeline(abun = func.based_sim_info$taxa_perturb_abun,
                                             func = MAG_rep_func,
                                             tree = MAG_rep_tree,
                                             group1_samples = func.based_sim_info$group1,
                                             group2_samples = func.based_sim_info$group2,
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
                     
                     output$func = func.based_sim_info$func
                     
                     saveRDS(object = output,
                             file = paste("POMS_out/func.based/POMS_func_out_",
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
          
          taxa.based_sim_info <- readRDS(paste("sim_info/taxa.based/taxa.based_sim_info_",
                                               "rep", as.character(rep_i),
                                               "_MAGs", as.character(MAG_num),
                                               "_pseudo", as.character(pseudocount_set),
                                               "_increase", as.character(abun_increase_set),
                                               ".rds", sep = ""))
          
          output <- POMS_pipeline(abun = taxa.based_sim_info$taxa_perturb_abun,
                                  func = MAG_rep_func,
                                  tree = MAG_rep_tree,
                                  group1_samples = taxa.based_sim_info$group1,
                                  group2_samples = taxa.based_sim_info$group2,
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
          
          output$func = taxa.based_sim_info$func
          
          saveRDS(object = output,
                  file = paste("POMS_out/taxa.based/POMS_func_out_",
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
  
  clade.based_sim_info <- readRDS(paste("sim_info/clade.based/clade.based_sim_info_",
                                       "rep", as.character(rep_i),
                                       "_MAGs", as.character(MAG_num),
                                       "_pseudo", as.character(pseudocount_set),
                                       "_increase", as.character(abun_increase_set),
                                       ".rds", sep = ""))
  
  output <- POMS_pipeline(abun = clade.based_sim_info$taxa_perturb_abun,
                          func = MAG_rep_func,
                          tree = clade.based_sim_info$subset_tree,
                          group1_samples = clade.based_sim_info$group1,
                          group2_samples = clade.based_sim_info$group2,
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
  
  output$func = clade.based_sim_info$func
  
  saveRDS(object = output,
          file = paste("POMS_out/clade.based/POMS_func_out_",
                       "rep", as.character(rep_i),
                       "_MAGs", as.character(MAG_num),
                       "_pseudo", as.character(pseudocount_set),
                       "_increase", as.character(abun_increase_set),
                       ".rds", sep = ""))
  
  return("Success")
  
}, mc.cores = 20)
