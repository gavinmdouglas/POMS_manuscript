# Specifically perform "selection" for "random taxa" replicates.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

library(ape)
devtools::load_all(path = "/home/gavin/github_repos/POMS/")

# Old version of subsetting function, from previous version. Using here to be consistent.
subset_abun_table <- function(in_abun, col2keep) {
  
  in_abun <- in_abun[, which(colnames(in_abun) %in% col2keep)]
  
  missing_rows <- which(rowSums(in_abun) == 0)
  missing_samples <- which(colSums(in_abun) == 0)
  
  if(length(missing_rows) > 0) {
    in_abun <- in_abun[-missing_rows, ]
  }
  
  if(length(missing_samples) > 0) {
    in_abun <- in_abun[, -missing_samples]
  }
  
  return(in_abun)
}

# Read in pre-determined random sample groupings.
random_groups <- readRDS("almeida_healthy_random_groups.rds")

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- subset_abun_table(in_abun = almeida_abun,
                                  col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))
almeida_tree <- read.tree("../../key_inputs/Almeida2019_dataset/phylogenies/raxml_hgr-umgs_phylogeny.nwk")

set.seed(71632)

# Filter rare functions from table.
MAG_nums <- c(1595, 1250, 1000, 750, 500, 250, 100, 50)

pseudocount_settings <- c(0, 0.1, 0.3, 0.5, 0.7, 0.9, 1)

abun_increase_settings <- c(1.5, 1.3, 1.1, 1.05)

for (rep_i in 1:25) {
  
  for (MAG_num in MAG_nums) {
  
    MAG_subset_func <- readRDS(paste("parameter_altered_files/prepped_func_tables/subset_",
                                     as.character(MAG_num),
                                     "MAGs_func_rep",
                                     as.character(rep_i),
                                     ".rds", sep = ""))
    
    representative_sim_info <- readRDS(paste("parameter_altered_files/sim_info/func.based/func_sim_info_rep",
                                             as.character(rep_i),
                                             "_MAGs",
                                             as.character(MAG_num),
                                             "_pseudo0_increase1.05.rds",
                                             sep = ""))
    
    rand_contributors <- sample(x = rownames(MAG_subset_func),
                                size = length(representative_sim_info$orig_contributors))
    
    group1_samples <- representative_sim_info$group1
    group2_samples <- representative_sim_info$group2
  
    # Loop through different selection settings.
    for (pseudocount_set in pseudocount_settings) {
     
       for (abun_increase_set in abun_increase_settings) {
         
         rep_abun_MAG_subset <- almeida_abun[rownames(MAG_subset_func), ]
         
         if (length(which(colSums(rep_abun_MAG_subset) == 0)) > 0) {
           rep_abun_MAG_subset <- rep_abun_MAG_subset[, -which(colSums(rep_abun_MAG_subset) == 0)]
         }

         num_to_bump <- floor(length(rand_contributors) * pseudocount_set)
         if (num_to_bump > 0) {
           for (g1_samp in group1_samples) {
             
             increased_contributors <- sample(rand_contributors, size = num_to_bump, replace = FALSE)
             
             rep_abun_MAG_subset[increased_contributors, g1_samp] <- rep_abun_MAG_subset[increased_contributors, g1_samp] + 1
               
           }
         }

         rep_abun_MAG_subset[rand_contributors, group1_samples] <- rep_abun_MAG_subset[rand_contributors, group1_samples] * abun_increase_set
         
         rep_output <- list(taxa_perturb_abun = rep_abun_MAG_subset,
                            rand_contributors = rand_contributors,
                            group1 = group1_samples,
                            group2 = group2_samples,
                            pseudocount_set = pseudocount_set,
                            abun_increase = abun_increase_set,
                            selected_group = "group1")
         
         outfile <- paste("parameter_altered_files/sim_info/taxa.based/taxa.based_sim_info_",
                          "rep", as.character(rep_i),
                          "_MAGs", as.character(MAG_num),
                          "_pseudo", as.character(pseudocount_set),
                          "_increase", as.character(abun_increase_set),
                          ".rds", sep = "")
         
         saveRDS(file = outfile, object = rep_output)
              
       }
    }
  }
}
