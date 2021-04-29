### Code for figuring out what subset of the data to use as a clear (and small) example for the overview figure.

rm(list = ls(all.names = TRUE))

setwd("~/github_repos/POMS_manuscript/data/intermediates/MAG.based_simulations/")

library(ggtree)
library(ape)
library(ggplot2)
library(ggbeeswarm)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

test_tree <- readRDS("MAG.based_prepped_tree.rds")

tip_subset <- test_tree$tip.label[661:720]

test_tree <- drop.tip(phy = test_tree,
                      tip = test_tree$tip.label[which(!test_tree$tip.label %in% tip_subset)],
                      trim.internal = TRUE)

test_func <- readRDS("MAG.based_prepped_func.rds")

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


example_enrichment_tree_w_labels <- ggtree(example_enrichment_tree_prepped$prepped_tree, layout  =  "rectangular") +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$tested_nodes)), color = "black", fill = "white", alpha = 1, size = 10, pch = 22) +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group1_balance_dir_i)), color = "black", fill = "red", alpha = 0.3, size = 10, pch = 22) +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group2_balance_dir_i)), color = "black", fill = "blue", alpha = 0.3, size = 10, pch = 22) +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$enriched_nonsig)), color = "black", fill = "grey", alpha = 0.75, size = 6, pch = 21) +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$lhs_balance_function_enrich_dir_i)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
                                          geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$rhs_balance_function_enrich_dir_i)), color = "black", fill = "blue", alpha = 0.75, size = 6, pch = 21) +
                                          geom_tippoint(aes(subset = (node %in% example_enrichment_tree_prepped$func_positive_tips)), pch = 16, size = 6, color = "dodgerblue4") +
                                          theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5)) +
                                          geom_text(aes(label = node), hjust = -.3)


# At this point need to go through each tested node and make sure that the lhs (i.e. top) and rhs (i.e. bottom) side
# of each node are the same in the visualization as they were when used by POMS.
# ***Note that this is an iterative process - when you flip a higher node it can also flip child nodes. So you need to check each time that you make a change.***
# This basic code below was used to do this process manually after each change.
test_node_i <- 77
test_node_label <- example_enrichment_tree_prepped$prepped_tree$node.label[test_node_i - length(example_enrichment_tree_prepped$prepped_tree$tip.label)]
which(example_enrichment_tree_prepped$prepped_tree$tip.label %in% output$balances_info$features[[test_node_label]]$lhs)
which(example_enrichment_tree_prepped$prepped_tree$tip.label %in% output$balances_info$features[[test_node_label]]$rhs)

example_enrichment_tree_w_labels %>% ggtree::rotate(61) %>% ggtree::rotate(77) %>% ggtree::rotate(80) %>% ggtree::rotate(92) %>% ggtree::rotate(93)


# Once all nodes that need to be rotated are finalized, then made actual version:
example_enrichment_tree <- ggtree(example_enrichment_tree_prepped$prepped_tree, layout  =  "rectangular") +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$tested_nodes)), color = "black", fill = "white", alpha = 1, size = 10, pch = 22) +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group1_balance_dir_i)), color = "black", fill = "red", alpha = 0.3, size = 10, pch = 22) +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$group2_balance_dir_i)), color = "black", fill = "blue", alpha = 0.3, size = 10, pch = 22) +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$enriched_nonsig)), color = "black", fill = "grey", alpha = 0.75, size = 6, pch = 21) +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$lhs_balance_function_enrich_dir_i)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
                                  geom_nodepoint(aes(subset = (node %in% example_enrichment_tree_prepped$rhs_balance_function_enrich_dir_i)), color = "black", fill = "blue", alpha = 0.75, size = 6, pch = 21) +
                                  geom_tippoint(aes(subset = (node %in% example_enrichment_tree_prepped$func_positive_tips)), pch = 16, size = 6, color = "dodgerblue4") +
                                  theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5))

example_enrichment_tree <- example_enrichment_tree  %>% ggtree::rotate(61) %>% ggtree::rotate(77) %>% ggtree::rotate(80) %>% ggtree::rotate(92) %>% ggtree::rotate(93)


relabun <- data.frame(sweep(x = func_sim_info$taxa_perturb_abun, MARGIN = 2,
                            STATS = colSums(func_sim_info$taxa_perturb_abun), FUN = '/')) * 100

mean_relabun <- data.frame("Group 1" = rowMeans(relabun[, func_sim_info$group1]),
                           "Group 2" = rowMeans(relabun[, func_sim_info$group2]),
                           check.names = FALSE)

mean_relabun_log10 <- log10(mean_relabun + 0.01)

example_enrichment_ggtree_heatmap <- gheatmap(example_enrichment_tree, mean_relabun_log10, offset = 0, width = 0.1, font.size = 3, colnames_angle = -45, hjust = 0) +
                                              scale_fill_continuous(name = "log10(% relabun.)", low = "gold", high = "dark green")




example_weak_tree_prepped <- prep_tree_nodes_func(in_list = output, focal_func = "K07106")
example_weak_tree_prepped$prepped_tree$edge.length <- NULL

# Add in tips positive for genome as well.
K07106_positive_genomes <- rownames(test_func[which(test_func$K07106 > 0), ])

example_weak_tree_prepped$func_positive_tips <- which(example_weak_tree_prepped$prepped_tree$tip.label %in% K07106_positive_genomes)

example_weak_tree <- ggtree(example_weak_tree_prepped$prepped_tree, layout  =  "rectangular") +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$tested_nodes)), color = "black", fill = "white", alpha = 1, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group1_balance_dir_i)), color = "black", fill = "red", alpha = 0.3, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group2_balance_dir_i)), color = "black", fill = "blue", alpha = 0.3, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$enriched_nonsig)), color = "black", fill = "grey", alpha = 0.75, size = 6, pch = 21) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$lhs_balance_function_enrich_dir_i)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$rhs_balance_function_enrich_dir_i)), color = "black", fill = "blue", alpha = 0.75, size = 6, pch = 21) +
  geom_tippoint(aes(subset = (node %in% example_weak_tree_prepped$func_positive_tips)), pch = 16, size = 6, color = "dodgerblue4") +
  theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5))

# NOTE THAT HAD TO MANUALLY SPECIFY ONE NODE TO BE ENRICHED ON TOP SIDE BASED ON VISUALIZATION (RATHER THAN KEEPING GREY CIRCLE CATEGORY)

example_weak_tree_prepped$enriched_nonsig <- c()
example_weak_tree_prepped$lhs_balance_function_enrich_dir_i <- c(example_weak_tree_prepped$lhs_balance_function_enrich_dir_i, 93)

example_weak_tree <- ggtree(example_weak_tree_prepped$prepped_tree, layout  =  "rectangular") +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$tested_nodes)), color = "black", fill = "white", alpha = 1, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group1_balance_dir_i)), color = "black", fill = "red", alpha = 0.3, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$group2_balance_dir_i)), color = "black", fill = "blue", alpha = 0.3, size = 10, pch = 22) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$enriched_nonsig)), color = "black", fill = "grey", alpha = 0.75, size = 6, pch = 21) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$lhs_balance_function_enrich_dir_i)), color = "black", fill = "red", alpha = 0.75, size = 6, pch = 21) +
  geom_nodepoint(aes(subset = (node %in% example_weak_tree_prepped$rhs_balance_function_enrich_dir_i)), color = "black", fill = "blue", alpha = 0.75, size = 6, pch = 21) +
  geom_tippoint(aes(subset = (node %in% example_weak_tree_prepped$func_positive_tips)), pch = 16, size = 6, color = "dodgerblue4") +
  theme(plot.title  =  element_text(hjust  =  0.5, vjust = -27.5))

example_weak_tree <- example_weak_tree  %>% ggtree::rotate(61) %>% ggtree::rotate(77) %>% ggtree::rotate(80) %>% ggtree::rotate(92) %>% ggtree::rotate(93)


relabun <- data.frame(sweep(x = func_sim_info$taxa_perturb_abun, MARGIN = 2,
                            STATS = colSums(func_sim_info$taxa_perturb_abun), FUN = '/')) * 100

mean_relabun <- data.frame("Group 1" = rowMeans(relabun[, func_sim_info$group1]),
                           "Group 2" = rowMeans(relabun[, func_sim_info$group2]),
                           check.names = FALSE)

mean_relabun_log10 <- log10(mean_relabun + 0.01)

example_weak_ggtree_heatmap <- gheatmap(example_weak_tree, mean_relabun_log10, offset = 0, width = 0.1, font.size = 3, colnames_angle = -45, hjust = 0) +
                                        scale_fill_continuous(name = "log10(% relabun.)", low = "gold", high = "dark green")


# Print example sample balances at one node to illustrate
# Because this is just a plot to introduce the tool rather than a result,
# I simulated normally distributed data to plot to make it prettier then the simulated data the plot is based on.

ex_balances <- data.frame(balances = c(rnorm(n = 25, mean = mean(output$balances_info$balances$n17[output$input_param$group1_samples]), sd = 0.25),
                                       rnorm(n = 25, mean = mean(output$balances_info$balances$n17[output$input_param$group2_samples]), sd = 0.25)),
                          group = c(rep(x = "Group 1", times = 25), rep(x = "Group 2", times = 25)))

example_balances_boxplots <- ggplot(ex_balances, aes(x = group, y = balances)) +
                                          geom_quasirandom(color = "black") +
                                          geom_boxplot(outlier.shape = NA, alpha = 0.5, lwd = 1, fill = "grey") +
                                          theme_bw() +
                                          ylab("ILR(top taxa rel. abun. / bottom taxa rel. abun.)") +
                                          xlab("") +
                                          theme(legend.position = "none")


ggsave(filename = "../../../display_items/Maintext_overview_panel_A_RAW.pdf",
       plot = example_enrichment_ggtree_heatmap,
       width = 6, height = 8, dpi = 600, device = "pdf")

ggsave(filename = "../../../display_items/Maintext_overview_panel_B_RAW.pdf",
       plot = example_weak_ggtree_heatmap,
       width = 6, height = 8, dpi = 600, device = "pdf")

ggsave(filename = "../../../display_items/Maintext_overview_panel_C_RAW.pdf",
       plot = example_balances_boxplots,
       width = 4, height = 3, dpi = 600, device = "pdf")
