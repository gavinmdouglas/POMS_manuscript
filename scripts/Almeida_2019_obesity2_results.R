rm(list=ls(all.names=TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

# Note that plot.margin gives parameters from the top clockwise.

almeida_taxa <- readRDS(file = "Almeida_2019_POMS_output/taxa_table.rds")

ERP003612_out <- readRDS(file = "Almeida_2019_POMS_output/ERP003612_POMS_out.rds")

ERP003612_pseudo_null_sig <- readRDS(file = "Almeida_2019_POMS_output/ERP003612_pseudo.null_sig.rds")

ERP003612_K13038_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP003612_K13038_pseudo.null.rds")

ERP003612_out_tree_prep <- prep_tree_sig_nodes(in_list = ERP003612_out, taxa_table = almeida_taxa)
ERP003612_sig_nodes_tree <- ggtree(ERP003612_out_tree_prep$prepped_tree, layout = "circular") +
                                    geom_nodepoint(aes(subset=(node %in% ERP003612_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=4) +
                                    theme(plot.margin=unit(c(-35,-10,-25, -25), "mm")) +
                                    ggtitle("All significant nodes") +
                                    theme(plot.title = element_text(hjust = 0.5, vjust=-27.5))

ERP003612_out_K13038_tree_prep <- prep_tree_nodes_func(in_list = ERP003612_out, focal_func = "K13038", taxa_table = almeida_taxa)
ERP003612_out_K13038_tree <- ggtree(ERP003612_out_K13038_tree_prep$prepped_tree, layout = "circular") +
                                      geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$nonenriched_sig)), color="black", fill="white", alpha=0.75, size=4, pch=21) +
                                      geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$nonenriched_nonsig)), color="black", alpha=0.75, size=4, pch=20) +
                                      geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$enriched_nonsig)), color="grey", alpha=0.75, size=4, pch=18) +
                                      geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$enriched)), color="black", fill="red", alpha=0.75, size=4, pch=21) +
                                      geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$depleted)), color="black", fill="blue", alpha=0.75, size=4, pch=21) +
                                      theme(plot.margin=unit(c(-35,-10,-25,-25), "mm"),
                                            plot.title = element_text(hjust = 0.5, vjust=-27.5)) +
                                      ggtitle("K13038 enrichment pattern")


ERP003612_out_num_enriched <- ERP003612_out$df[, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]

ERP003612_out_num_enriched <- ddply(ERP003612_out_num_enriched, .(num_sig_nodes_pos_enrich, num_sig_nodes_neg_enrich), nrow)

colnames(ERP003612_out_num_enriched)[3] <- "Count"

ERP003612_out_num_enriched_heatmap <- ggplot(ERP003612_out_num_enriched, aes(x=num_sig_nodes_pos_enrich, y=num_sig_nodes_neg_enrich, fill=log10(Count))) +
                                              geom_tile() +
                                              xlab("No. positively enriched nodes") +
                                              ylab("No. negatively enriched nodes") +
                                              labs(fill=expression("log" [10]*"(Count)")) +
                                              scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 15)) +
                                              scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 12)) +
                                              theme_bw() +
                                              ggtitle("All KO enrichments") +
                                              theme(legend.position = c(0.865, 0.7),
                                                    legend.box.background = element_rect(colour="black"),
                                                    legend.title = element_text(size = 7), 
                                                    legend.text = element_text(size = 7),
                                                    legend.key.size = unit(0.9, "lines"),
                                                    plot.title = element_text(hjust = 0.5)) +
                                              guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                     color = guide_legend(override.aes = list(size = 0.5)))
                                            

K13038_pseudo_null <- ddply(ERP003612_K13038_pseudo.null, .(pos_nodes, neg_nodes), nrow)
colnames(K13038_pseudo_null)[3] <- "Count"

K13038_pseudo_null_heatmap <- ggplot(K13038_pseudo_null, aes(x=pos_nodes, y=neg_nodes, fill=Count, width=1, height=1)) +
                                      geom_tile() +
                                      xlab("No. positively enriched nodes") +
                                      ylab("No. negatively enriched nodes") +
                                      labs(fill="Count") +
                                      scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 15)) +
                                      scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 12)) +
                                      theme_bw() +
                                      ggtitle("K13038 enrich. pseudo-null distribution") +
                                      theme(legend.position = c(0.85, 0.7),
                                            legend.box.background = element_rect(colour="black"),
                                            legend.title = element_text(size = 7),
                                            legend.text = element_text(size = 7),
                                            legend.key.size = unit(0.9, "lines"),
                                            plot.title = element_text(hjust = 0.5)) +
                                      guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                             color = guide_legend(override.aes = list(size = 0.5))) +
                                      geom_point(aes(x=13.5, y=1.5), colour="gold3", pch=4, size=6)


obesity2_sig_nodes_tree_and_heatmap <- plot_grid(ERP003612_sig_nodes_tree,
                                                 ERP003612_out_num_enriched_heatmap,
                                                  nrow=1, labels=c('a', 'b'))

obesity2_K13038_tree_and_pseudo.null_heatmap <- plot_grid(ERP003612_out_K13038_tree,
                                                              K13038_pseudo_null_heatmap,
                                                              nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/obesity2_K13038_tree_and_pseudo.null_heatmap.pdf", plot = obesity2_K13038_tree_and_pseudo.null_heatmap,
       width = 7.20472, height=3.5)
