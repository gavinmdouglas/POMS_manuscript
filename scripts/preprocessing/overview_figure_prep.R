### Code for figuring out what subset of the data to use as a clear (and small) example for the overview figure.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

library(ggtree)
library(ape)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

test_tree <- readRDS("MAG.based_prepped_tree.rds")

# Try different subsets of 50 tips until you find one with the most non-negligible nodes.

min_num_tips <- 4
range_start <- 1
tip_subset <- test_tree$tip.label[range_start:(range_start + 59)]

nonnegligible_counts <- c()

while (length(which(is.na(tip_subset))) == 0) {
  
  tip_subset <- test_tree$tip.label[range_start:(range_start + 59)]
  test_tree_tmp <- drop.tip(phy = test_tree,
                            tip = test_tree$tip.label[which(!test_tree$tip.label %in% tip_subset)],
                            trim.internal = TRUE)
  test_tree_tmp <- ape::makeNodeLabel(test_tree_tmp, method = "number", prefix = 'n')
  
  node_features <- parallel::mclapply(test_tree_tmp$node.label,
                                      lhs_rhs_asvs,
                                      tree = test_tree_tmp,
                                      get_node_index = TRUE,
                                      mc.cores = 10)
  
  names(node_features) <- test_tree_tmp$node.label
  
  negligible_nodes <- sapply(names(node_features),
                             function(x) {
                               lhs_feat_num <- length(node_features[[x]]$lhs)
                               rhs_feat_num <- length(node_features[[x]]$rhs)
                               if ((lhs_feat_num < min_num_tips) || (rhs_feat_num < min_num_tips)) {
                                 return(TRUE)
                               } else {
                                 return(FALSE) 
                               }
                             })
  
  nonnegligible_counts <- c(nonnegligible_counts, sum(!negligible_nodes))
  
  range_start <- range_start + 1
}


names(node_features) <- nodes2test

which(nonnegligible_counts == max(nonnegligible_counts))

tip_subset <- test_tree$tip.label[661:720]

test_tree <- drop.tip(phy = test_tree,
                      tip = test_tree$tip.label[which(!test_tree$tip.label %in% tip_subset)],
                      trim.internal = TRUE)

test_func <- readRDS("MAG.based_prepped_func.rds")

min_p <- 1
rep_i <- 1

while (min_p > 0.05) {

  print(rep_i)
  
  func_sim_info <- readRDS(paste("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep", as.character(rep_i), ".rds", sep = ""))
  
  # These parameters were set intentionally leniently to help make an example plot.
  output <- POMS_pipeline(abun = func_sim_info$taxa_perturb_abun,
                            func = test_func,
                            phylogeny = test_tree,
                            group1_samples = func_sim_info$group1,
                            group2_samples = func_sim_info$group2,
                            ncores = 10,
                            balance_p_cutoff = 0.05,
                            balance_correction = "none",
                            function_p_cutoff = 0.05,
                            function_correction = "none",
                            min_num_tips = 4,
                            min_func_instances = 0,
                            min_func_prop = 0.2,
                            run_multinomial_test = TRUE,
                            multinomial_correction = "BH",
                            calc_node_dist = FALSE,
                            detailed_output = TRUE,
                            verbose = FALSE)
 
  min_p <- min(output$df$multinomial_p, na.rm = TRUE)
  
  rep_i <- rep_i + 1
  
}

# Determined that rep 76 had at least 1 sig gene with this tree.

func_sim_info <- readRDS("MAG.based_prepped_func_sim_info_sel1.5/func_sim_info_rep76.rds")

# These parameters were set intentionally leniently to help make an example plot.
output <- POMS_pipeline(abun = func_sim_info$taxa_perturb_abun,
                        func = test_func,
                        phylogeny = test_tree,
                        group1_samples = func_sim_info$group1,
                        group2_samples = func_sim_info$group2,
                        ncores = 10,
                        balance_p_cutoff = 0.05,
                        balance_correction = "none",
                        function_p_cutoff = 0.05,
                        function_correction = "none",
                        min_num_tips = 4,
                        min_func_instances = 0,
                        min_func_prop = 0.2,
                        run_multinomial_test = TRUE,
                        multinomial_correction = "BH",
                        calc_node_dist = FALSE,
                        detailed_output = TRUE,
                        verbose = FALSE)


example_enrichment_tree_prepped <- prep_tree_nodes_func(in_list = output, focal_func = "K02036")
example_enrichment_tree_prepped$prepped_tree$edge.length <- NULL

# Add in tips positive for genome as well.
K02036_positive_genomes <- rownames(test_func[which(test_func$K02036 > 0), ])

example_enrichment_tree_prepped$func_positive_tips <- which(example_enrichment_tree_prepped$prepped_tree$tip.label %in% K02036_positive_genomes)

example_enrichment_tree <- ggtree(example_enrichment_tree_prepped$prepped_tree, layout  =  "rectangular") +
                                    geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$tested_nodes)), color = "black", fill = "grey90", alpha = 1, size = 8, pch = 22) +
                                    geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$enriched_nonsig)), color = "black", fill = "white", alpha = 0.75, size = 6, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group1_func_enriched_nodes)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group2_func_enriched_nodes)), color = "black", fill = "dodgerblue3", alpha = 0.75, size = 6, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$sig_balances_nodes)), color = "black", fill = "black", alpha = 1, size = 3, pch = 23) +
                                    geom_tippoint(aes(subset = (node %in% example_enrichment_tree_prepped$func_positive_tips)), pch = 16,size = 2, color = "dodgerblue4") +
                                    theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5)) +
                                    ggtitle("")

relabun <- data.frame(sweep(x = func_sim_info$taxa_perturb_abun, MARGIN = 2,
                            STATS = colSums(func_sim_info$taxa_perturb_abun), FUN = '/')) * 100

mean_relabun <- data.frame("Group 1" = rowMeans(relabun[, func_sim_info$group1]),
                           "Group 2" = rowMeans(relabun[, func_sim_info$group2]),
                           check.names = FALSE)

mean_relabun_log10 <- log10(mean_relabun + 0.01)

example_enrichment_ggtree_heatmap <- gheatmap(example_enrichment_tree, mean_relabun_log10, offset = 0, width = 0.1, font.size = 3, colnames_angle = -45, hjust = 0) +
                                              scale_fill_continuous(name = "log10(% relabun.)", low = "gold", high = "dark green")




example_weak_tree_prepped <- prep_tree_nodes_func(in_list = output, focal_func = "K00174")
example_weak_tree_prepped$prepped_tree$edge.length <- NULL

# Add in tips positive for genome as well.
K00174_positive_genomes <- rownames(test_func[which(test_func$K00174 > 0), ])

example_weak_tree_prepped$func_positive_tips <- which(example_weak_tree_prepped$prepped_tree$tip.label %in% K00174_positive_genomes)

example_weak_tree <- ggtree(example_weak_tree_prepped$prepped_tree, layout  =  "rectangular") +
                            geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$tested_nodes)), color = "black", fill = "grey90", alpha = 1, size = 8, pch = 22) +
                            geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$enriched_nonsig)), color = "black", fill = "white", alpha = 0.75, size = 6, pch = 21) +
                            geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group1_func_enriched_nodes)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
                            geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group2_func_enriched_nodes)), color = "black", fill = "dodgerblue3", alpha = 0.75, size = 6, pch = 21) +
                            geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$sig_balances_nodes)), color = "black", fill = "black", alpha = 1, size = 3, pch = 23) +
                            geom_tippoint(aes(subset = (node %in% example_weak_tree_prepped$func_positive_tips)), pch = 16,size = 2, color = "dodgerblue4") +
                            theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5)) +
                            ggtitle("")


relabun <- data.frame(sweep(x = func_sim_info$taxa_perturb_abun, MARGIN = 2,
                            STATS = colSums(func_sim_info$taxa_perturb_abun), FUN = '/')) * 100

mean_relabun <- data.frame("Group 1" = rowMeans(relabun[, func_sim_info$group1]),
                           "Group 2" = rowMeans(relabun[, func_sim_info$group2]),
                           check.names = FALSE)

mean_relabun_log10 <- log10(mean_relabun + 0.01)

example_weak_ggtree_heatmap <- gheatmap(example_weak_tree, mean_relabun_log10, offset = 0, width = 0.1, font.size = 3, colnames_angle = -45, hjust = 0) +
                                        scale_fill_continuous(name = "log10(% relabun.)", low = "gold", high = "dark green")

ggsave(filename = "../../../display_items/Main_overview_panel_A_RAW.pdf",
       plot = example_enrichment_ggtree_heatmap,
       width = 6, height = 8, dpi = 600, device = "pdf")

ggsave(filename = "../../../display_items/Main_overview_panel_B_RAW.pdf",
       plot = example_weak_ggtree_heatmap,
       width = 6, height = 8, dpi = 600, device = "pdf")
