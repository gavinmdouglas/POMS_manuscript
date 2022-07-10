rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

library(parallel)

random_groups <- readRDS("almeida_healthy_random_groups.rds")
almeida_tree <- readRDS(file = "MAG.based_prepped_tree.rds")
almeida_func <- readRDS(file = "MAG.based_prepped_func.rds")
almeida_abun <- readRDS(file = "MAG.based_prepped_abun.rds")

almeida_tree <- prep_tree(almeida_tree, tips2keep = rownames(almeida_func))

node_features <- parallel::mclapply(almeida_tree$node.label,
                                    function(x) { 
                                      tmp <- lhs_rhs_tips(x,
                                                          tree = almeida_tree,
                                                          get_node_index = TRUE)
                                      
                                      return(c(tmp$lhs, tmp$rhs))
                                    }, mc.cores=30)

node_feature_lengths <- sapply(node_features, length)

nonnegligible_nodes_i <- which(node_feature_lengths >= 5 & node_feature_lengths < (length(almeida_tree$tip.label)))

# Note that there are only 693 non-negligible nodes in the tree, so just doing one replicate for each separate node.
clade.based_sim_info <- mclapply(1:length(nonnegligible_nodes_i),
                              
                              function(rep_i) {

                                group1_samples <- random_groups$group1[[rep_i]]
                                group2_samples <- random_groups$group2[[rep_i]]
                                
                                randomly_selected_node_i <- nonnegligible_nodes_i[rep_i]
                                randomly_selected_node <- almeida_tree$node.label[randomly_selected_node_i]
                                underlying_tips <- node_features[[randomly_selected_node_i]]
                                
                                tmp_abun <- almeida_abun
                                
                                # Sanity check that no samples are missing.
                                missing_samples <- which(!c(group1_samples, group2_samples) %in% colnames(tmp_abun))
                                if (length(missing_samples) > 0) { stop("Error - there are samples missing.") }
                                
                                tmp_abun[underlying_tips, group1_samples] <- (tmp_abun[underlying_tips, group1_samples] + 1) * 1.5
                                
                                rep_output <- list(taxa_perturb_abun = tmp_abun,
                                                   tree=almeida_tree,
                                                   randomly_selected_node = randomly_selected_node,
                                                   underlying_tips = underlying_tips,
                                                   group1 = group1_samples,
                                                   group2 = group2_samples,
                                                   selected_group = "group1")
                                
                                outfile <- paste("MAG.based_prepped_clade.based_sim_info_sel1.5/", "clade.based_sim_info_rep", as.character(rep_i), ".rds", sep = "")
                                
                                saveRDS(file = outfile, object = rep_output)
                                
                                return("Success")

                              }, mc.cores = 30)

