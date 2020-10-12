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

ERP012177_out <- readRDS(file = "Almeida_2019_POMS_output/ERP012177_POMS_out.rds")

ERP012177_pseudo_null_sig <- readRDS(file = "Almeida_2019_POMS_output/ERP012177_pseudo.null_sig.rds")


ERP012177_K00702_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP012177_K00702_pseudo.null.rds")

ERP012177_out_K00702_tree_prep <- prep_tree_nodes_func(in_list = ERP012177_out, focal_func = "K00702", taxa_table = almeida_taxa)
ERP012177_out_K00702_tree <- ggtree(ERP012177_out_K00702_tree_prep$prepped_tree, layout = "circular") +
                                      geom_nodepoint(aes(subset=(node %in% ERP012177_out_K00702_tree_prep$nonenriched_sig)), color="black", fill="white", alpha=0.75, size=4, pch=21) +
                                      geom_nodepoint(aes(subset=(node %in% ERP012177_out_K00702_tree_prep$nonenriched_nonsig)), color="black", alpha=0.75, size=4, pch=20) +
                                      geom_nodepoint(aes(subset=(node %in% ERP012177_out_K00702_tree_prep$enriched_nonsig)), color="grey", alpha=0.75, size=4, pch=18) +
                                      geom_nodepoint(aes(subset=(node %in% ERP012177_out_K00702_tree_prep$enriched)), color="black", fill="red", alpha=0.75, size=4, pch=21) +
                                      geom_nodepoint(aes(subset=(node %in% ERP012177_out_K00702_tree_prep$depleted)), color="black", fill="blue", alpha=0.75, size=4, pch=21) +
                                      theme(plot.margin=unit(c(-35,-10,-25,-25), "mm"),
                                            plot.title = element_text(hjust = 0.5, vjust=-27.5)) +
                                      ggtitle("K00702 enrichment pattern")
                                                                                

K00702_pseudo_null <- ddply(ERP012177_K00702_pseudo.null, .(pos_nodes, neg_nodes), nrow)
colnames(K00702_pseudo_null)[3] <- "Count"

K00702_pseudo_null_heatmap <- ggplot(K00702_pseudo_null, aes(x=pos_nodes, y=neg_nodes, fill=Count, width=1, height=1)) +
                                      geom_tile() +
                                      xlab("No. positively enriched nodes") +
                                      ylab("No. negatively enriched nodes") +
                                      labs(fill="Count") +
                                      scale_x_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 5)) +
                                      scale_y_continuous(expand = expansion(mult = c(0, 0)), limits=c(0, 5)) +
                                      theme_bw() +
                                      ggtitle("K00702 enrich. pseudo-null distribution") +
                                      theme(legend.position = c(0.85, 0.7),
                                            legend.box.background = element_rect(colour="black"),
                                            legend.title = element_text(size = 7),
                                            legend.text = element_text(size = 7),
                                            legend.key.size = unit(0.9, "lines"),
                                            plot.title = element_text(hjust = 0.5)) +
                                      guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                             color = guide_legend(override.aes = list(size = 0.5))) +
                                      geom_point(aes(x=4.5, y=1), colour="gold3", pch=4, size=6)


crc_K00702_tree_and_pseudo.null_heatmap <- plot_grid(ERP012177_out_K00702_tree,
                                                     K00702_pseudo_null_heatmap,
                                                     nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/crc_K00702_tree_and_pseudo.null_heatmap.pdf", plot = crc_K00702_tree_and_pseudo.null_heatmap,
       width = 7.20472, height=3.5)






### K09691 negatively enriched nodes;
"n74"   "n921"  "n1028" "n1030"

node_taxa(lhs_features = ERP012177_out$balances_info$features[["n74"]]$lhs,
          rhs_features = ERP012177_out$balances_info$features[["n74"]]$rhs,
          taxa = almeida_taxa, threshold = 0.5)

node_taxa(lhs_features = ERP012177_out$balances_info$features[["n921"]]$lhs,
          rhs_features = ERP012177_out$balances_info$features[["n921"]]$rhs,
          taxa = almeida_taxa)

node_taxa(lhs_features = ERP012177_out$balances_info$features[["n1028"]]$lhs,
          rhs_features = ERP012177_out$balances_info$features[["n1028"]]$rhs,
          taxa = almeida_taxa)

node_taxa(lhs_features = ERP012177_out$balances_info$features[["n1030"]]$lhs,
          rhs_features = ERP012177_out$balances_info$features[["n1030"]]$rhs,
          taxa = almeida_taxa)

