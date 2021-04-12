rm(list = ls(all.names = TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

almeida_taxa <- readRDS(file = "key_inputs/Almeida2019_dataset/taxonomy/taxa_table.rds")

ERP003612_out <- readRDS(file = "results/Almeida_2019_POMS_output/ERP003612_POMS_out.rds")

ERP003612_out_K01993_tree_prep <- prep_tree_nodes_func(in_list = ERP003612_out, focal_func = "K01993", taxa_table = almeida_taxa)

ERP003612_out_K01993_tree <- ggtree(ERP003612_out_K01993_tree_prep$prepped_tree, layout  =  "circular") +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_K01993_tree_prep$nonenriched_sig)), color = "black", fill = "black", alpha = 0.75, size = 4, pch = 23) +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_K01993_tree_prep$nonenriched_nonsig)), color = "black", fill = "white", , alpha = 0.75, size = 4, pch = 23) +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_K01993_tree_prep$enriched_nonsig)), color = "black", fill = "grey", alpha = 0.75, size = 4, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_K01993_tree_prep$enriched)), color = "black", fill = "red", alpha = 0.75, size = 4, pch = 21) +
                                    geom_nodepoint(aes(subset = (node %in% ERP003612_out_K01993_tree_prep$depleted)), color = "black", fill = "blue", alpha = 0.75, size = 4, pch = 21) +
                                    theme(plot.margin = unit(c(-35,-10,-25,-25), "mm"),
                                          plot.title  =  element_text(hjust  =  0.5, vjust = -27.5)) +
                                    ggtitle("K01993 enrichment pattern")


ERP003612_out_K01993_multinomial <- data.frame(Type = c(rep("Observed", 3), rep("Expected", 3)),
                                               
                                               Count = c(as.numeric(ERP003612_out$df["K01993", c("num_sig_nodes_group1_enrich", "num_sig_nodes_group2_enrich", "num_nonsig_nodes_enrich")]),
                                                          ERP003612_out$multinomial_exp_prop[1] * ERP003612_out$df["K01993", "num_nodes_enriched"],
                                                          ERP003612_out$multinomial_exp_prop[2] * ERP003612_out$df["K01993", "num_nodes_enriched"],
                                                          ERP003612_out$multinomial_exp_prop[3] * ERP003612_out$df["K01993", "num_nodes_enriched"]),
                                               
                                               Direction = c("Obese", "Control", "Non-sig.", "Obese", "Control", "Non-sig."))

ERP003612_out_K01993_multinomial$Direction <- factor(ERP003612_out_K01993_multinomial$Direction, levels = c("Control", "Obese", "Non-sig."))

ERP003612_out_K01993_multinomial_plot <- ggplot(ERP003612_out_K01993_multinomial, aes(x = Direction, y = Count, fill = Type)) +
                                                geom_bar(stat = "identity", position = position_dodge(), colour = "black") +
                                                scale_fill_manual(name = "Data type", values = c("white", "gold")) +
                                                theme_bw() +
                                                theme(plot.title = element_text(hjust = 0.5, vjust = -0.5)) +
                                                xlab("Gene family enrichment direction") +
                                                ylab("Number of enriched nodes") +
                                                ggtitle("K01993 node enrichment summary")


ERP003612_K01993_example_combined <- plot_grid(ERP003612_out_K01993_tree,
                                               ERP003612_out_K01993_multinomial_plot,
                                               nrow = 1, labels = c('a', 'b'))

ggsave(filename = "../display_items/Maintext_obesity1_K01993_example_RAW.pdf",
       plot = ERP003612_K01993_example_combined,
       width = 10, height = 5, dpi = 600, device = "pdf")
