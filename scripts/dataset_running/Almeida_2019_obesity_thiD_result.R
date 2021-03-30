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

ERP002061_out_tree_prep <- prep_tree_sig_nodes(in_list = ERP002061_out, taxa_table = almeida_taxa)


ERP002061_out_K00941_tree_prep <- prep_tree_nodes_func(in_list = ERP002061_out, focal_func = "K00941", taxa_table = almeida_taxa)

ERP002061_out_K00941_tree <- ggtree(ERP002061_out_K00941_tree_prep$prepped_tree, layout = "circular") +
                                     geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$enriched)), color="red", alpha=0.75, size=4) +
                                     geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$depleted)), color="blue", alpha=0.75, size=4) +
                                     geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$nonsig)), color="grey", alpha=0.75, size=4)+
                                     theme(plot.margin=unit(c(-35,-10,-25,-25), "mm"))


ERP002061_out_K00941_pseudo_null <- readRDS("Almeida_2019_POMS_output/ERP002061_K00941_pseudo.null.rds")

K00941_pseudo_null <- ddply(ERP002061_out_K00941_pseudo_null, .(pos_nodes, neg_nodes), nrow)

colnames(K00941_pseudo_null)[3] <- "Count"

K00941_pseudo_null_heatmap <- ggplot(K00941_pseudo_null, aes(x=pos_nodes, y=neg_nodes, fill=Count, width=1, height=1)) +
                                              geom_tile() +
                                              xlab("Nom. positively enriched nodes") +
                                              ylab("Nom. negatively enriched nodes") +
                                              labs(fill="Count") +
                                              scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 15)) +
                                              scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 12)) +
                                              theme_bw() +
                                              theme(legend.position = c(0.2, 0.725),
                                                    legend.box.background = element_rect(colour="black"),
                                                    legend.title = element_text(size = 7), 
                                                    legend.text = element_text(size = 7)) +
                                              guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                     color = guide_legend(override.aes = list(size = 0.5))) +
                                              geom_point(aes(x=14.25, y=1), colour="gold3", pch=4, size=6)
                                            

K00941_plot <- plot_grid(ERP002061_out_K00941_tree, K00941_pseudo_null_heatmap, nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/ERP002061_K00941_tree_and_pseudo.null_heatmap.pdf", plot = K00941_plot,
       width = 7.20472, height=3.5)




