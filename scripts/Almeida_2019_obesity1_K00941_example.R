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

ERP002061_out <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_POMS_out.rds")

ERP002061_K00941_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_K00941_pseudo.null.rds")

ERP002061_out_K00941_tree_prep <- prep_tree_nodes_func(in_list = ERP002061_out, focal_func = "K00941", taxa_table = almeida_taxa)

ERP002061_out_K00941_tree <- ggtree(ERP002061_out_K00941_tree_prep$prepped_tree, layout = "circular") +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$nonenriched_sig)), color="black", fill="white", alpha=0.75, size=4, pch=21) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$nonenriched_nonsig)), color="black", alpha=0.75, size=4, pch=20) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$enriched_nonsig)), color="grey", alpha=0.75, size=4, pch=18) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$enriched)), color="black", fill="red", alpha=0.75, size=4, pch=21) +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$depleted)), color="black", fill="blue", alpha=0.75, size=4, pch=21) +
                                    theme(plot.margin=unit(c(-35,-10,-25,-25), "mm"),
                                          plot.title = element_text(hjust = 0.5, vjust=-27.5)) +
                                    ggtitle("K00941 enrichment pattern")


K00941_pseudo_null <- ddply(ERP002061_K00941_pseudo.null, .(pos_nodes, neg_nodes), nrow)
colnames(K00941_pseudo_null)[3] <- "Count"

K00941_pseudo_null_heatmap <- ggplot(K00941_pseudo_null, aes(x=pos_nodes, y=neg_nodes, fill=Count, width=1, height=1)) +
                                        geom_tile() +
                                        xlab("No. positively enriched nodes") +
                                        ylab("No. negatively enriched nodes") +
                                        labs(fill="Count") +
                                        scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 15)) +
                                        scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 12)) +
                                        theme_bw() +
                                        ggtitle("K00941 enrich. pseudo-null distribution") +
                                        theme(legend.position = c(0.2, 0.7),
                                              legend.box.background = element_rect(colour="black"),
                                              legend.title = element_text(size = 7),
                                              legend.text = element_text(size = 7),
                                              plot.title = element_text(hjust = 0.5)) +
                                        guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                               color = guide_legend(override.aes = list(size = 0.5))) +
                                        geom_point(aes(x=14.25, y=1), colour="gold3", pch=4, size=6)
                                            

ERP002061_K00941_tree_and_pseudo.null_heatmap <- plot_grid(ERP002061_out_K00941_tree, K00941_pseudo_null_heatmap, nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/ERP002061_K00941_tree_and_pseudo.null_heatmap.pdf", plot = ERP002061_K00941_tree_and_pseudo.null_heatmap,
       width = 7.20472, height=3.5)



