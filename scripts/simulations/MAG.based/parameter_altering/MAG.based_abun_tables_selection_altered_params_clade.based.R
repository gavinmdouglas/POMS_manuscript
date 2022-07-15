# Create simulated tables based on altering the number of MAGs and the "selection pressure".
# This script generates additional simulation replicates where a pseudocount is added to a varying proportion of contributors per sample rather than all or none.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

library(ape)
devtools::load_all(path = "/home/gavin/github_repos/POMS/")

# Read in pre-determined random sample groupings.
random_groups <- readRDS("almeida_healthy_random_groups.rds")

almeida_abun <- read.table(file = "../../key_inputs/Almeida2019_dataset/mapping_results/modified/bwa_depth_min25coverage.tsv.gz",
                           header = TRUE, sep = "\t", check.names = FALSE, row.names = 1, quote = "", comment.char = "")

almeida_abun <- POMS::subset_by_col_and_filt(in_tab = almeida_abun,
                                             col2keep = c(random_groups$group1[[1]], random_groups$group2[[1]]))

set.seed(71632)

MAG_nums <- c(1595, 1000, 500, 250, 100)

pseudocount_settings <- c(0, 0.3, 0.7, 1)

abun_increase_settings <- c(1.5, 1.25, 1.05)

num_reps <- 10

for (rep_i in 1:num_reps) {
  
  for (MAG_num in MAG_nums) {
  
    MAG_subset_func <- readRDS(paste("parameter_altered_files/prepped_func_tables/subset_",
                                     as.character(MAG_num),
                                     "MAGs_func_rep",
                                     as.character(rep_i),
                                     ".rds", sep = ""))

    # Get tips underlying each node.
    tree_subset <- readRDS(paste("parameter_altered_files/prepped_trees/subset_",
                                 as.character(MAG_num),
                                 "MAGs_tree_rep",
                                 as.character(rep_i),
                                 ".rds", sep = ""))
    
    node_features <- parallel::mclapply(tree_subset$node.label,
                                        function(x) { 
                                          tmp <- lhs_rhs_tips(x,
                                                              tree = tree_subset,
                                                              get_node_index = TRUE)
                                          
                                          return(c(tmp$lhs, tmp$rhs))
                                        }, mc.cores=30)
    
    node_feature_lengths <- sapply(node_features, length)
    
    nonnegligible_nodes_i <- which(node_feature_lengths >= 5 & node_feature_lengths < (length(tree_subset$tip.label)))
    
    if (length(nonnegligible_nodes_i) == 0) {
      stop("No nodes remaining after restricting to those with sufficient underlying tips.") 
    }
    
    randomly_selected_node_i <- sample(x = nonnegligible_nodes_i, size = 1)
    randomly_selected_node <- tree_subset$node.label[randomly_selected_node_i]
    underlying_tips <- node_features[[randomly_selected_node_i]]
    
    representative_sim_info <- readRDS(paste("parameter_altered_files/sim_info/func.based/func.based_sim_info_rep",
                                             as.character(rep_i),
                                             "_MAGs",
                                             as.character(MAG_num),
                                             "_pseudo0_increase1.05.rds",
                                             sep = ""))
    
    group1_samples <- representative_sim_info$group1
    group2_samples <- representative_sim_info$group2
  
    # Loop through different selection settings.
    for (pseudocount_set in pseudocount_settings) {
     
       for (abun_increase_set in abun_increase_settings) {
         
         rep_abun_MAG_subset <- almeida_abun[rownames(MAG_subset_func), ]
         rep_abun_MAG_subset <- rep_abun_MAG_subset[, which(colSums(rep_abun_MAG_subset) > 0)]

         num_to_bump <- floor(length(underlying_tips) * pseudocount_set)
         if (num_to_bump > 0) {
           for (g1_samp in group1_samples) {
             
             increased_contributors <- sample(underlying_tips, size = num_to_bump, replace = FALSE)
             
             rep_abun_MAG_subset[increased_contributors, g1_samp] <- rep_abun_MAG_subset[increased_contributors, g1_samp] + 1
               
           }
         }

         rep_abun_MAG_subset[underlying_tips, group1_samples] <- rep_abun_MAG_subset[underlying_tips, group1_samples] * abun_increase_set
         
         rep_output <- list(taxa_perturb_abun = rep_abun_MAG_subset,
                            randomly_selected_node = randomly_selected_node,
                            underlying_tips = underlying_tips,
                            group1 = group1_samples,
                            group2 = group2_samples,
                            pseudocount_set = pseudocount_set,
                            abun_increase = abun_increase_set,
                            selected_group = "group1")
         
         outfile <- paste("parameter_altered_files/sim_info/clade.based/clade.based_sim_info_",
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
