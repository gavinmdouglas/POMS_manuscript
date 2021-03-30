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


ERP002061_sig_nodes_tree <- ggtree(ERP002061_out_tree_prep$prepped_tree, layout = "circular") +
                                    geom_nodepoint(aes(subset=(node %in% ERP002061_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=4) +
                                    theme(plot.margin=unit(c(-35,-10,-25,-25), "mm"))


ERP002061_out_num_enriched <- ERP002061_out$df[, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]

ERP002061_out_num_enriched <- ddply(ERP002061_out_num_enriched, .(num_sig_nodes_pos_enrich, num_sig_nodes_neg_enrich), nrow)

colnames(ERP002061_out_num_enriched)[3] <- "Count"

ERP002061_out_num_enriched_heatmap <- ggplot(ERP002061_out_num_enriched, aes(x=num_sig_nodes_pos_enrich, y=num_sig_nodes_neg_enrich, fill=log10(Count))) +
                                              geom_tile() +
                                              xlab("Nom. positively enriched balances") +
                                              ylab("Nom. negatively enriched balances") +
                                              labs(fill=expression("log" [10]*"(Count)")) +
                                              theme_bw() +
                                              theme(legend.position = c(0.865, 0.725),
                                                    legend.box.background = element_rect(colour="black"),
                                                    legend.title = element_text(size = 7), 
                                                    legend.text = element_text(size = 7)) +
                                              guides(shape = guide_legend(override.aes = list(size = 0.5)),
                                                     color = guide_legend(override.aes = list(size = 0.5)))
                                            

ERP002061_sig_nodes_tree_and_heatmap <- plot_grid(ERP002061_sig_nodes_tree, ERP002061_out_num_enriched_heatmap, nrow=1, labels=c('a', 'b'))

ggsave(filename = "../figures/ERP002061_sig_nodes_tree_and_heatmap.pdf", plot = ERP002061_sig_nodes_tree_and_heatmap,
       width = 7.20472, height=3.5)

ERP003612_out <- readRDS("balance_tree_out/ERP003612_almeida_out.rds")
ERP003612_out_tree_prep <- prep_tree_sig_balances(ERP003612_out, taxa_table = almeida_taxa)






ggtree(ERP003612_out_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP003612_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=5) +
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=2)



ERP012177_out <- readRDS("balance_tree_out/ERP012177_almeida_out.rds")
ERP012177_out_tree_prep <- prep_tree_sig_balances(ERP012177_out, taxa_table = almeida_taxa)

ggtree(ERP012177_out_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP012177_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=5) +
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=2)


SRP091570_out <- readRDS("balance_tree_out/SRP091570_almeida_out.rds")
SRP091570_out_tree_prep <- prep_tree_sig_balances(SRP091570_out, taxa_table = almeida_taxa)

ggtree(SRP091570_out_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% SRP091570_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=5) +
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=2)



### K09765 - epoxyqueuosine reductase, which is last step in biosynthesis of queuosine, that has been linked to obesity.

ERP002061_out_K09765_tree_prep <- prep_tree_balances_func(in_list = ERP002061_out, focal_func = "K09765", taxa_table = almeida_taxa)

ggtree(ERP002061_out_K09765_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K09765_tree_prep$enriched)), color="red", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K09765_tree_prep$depleted)), color="blue", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K09765_tree_prep$nonsig)), color="grey", alpha=0.75, size=5)


### K00091 - E1.1.1.219; dihydroflavonol-4-reductase

ERP002061_out_K00091_tree_prep <- prep_tree_balances_func(in_list = ERP002061_out, focal_func = "K00091", taxa_table = almeida_taxa)

ggtree(ERP002061_out_K00091_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00091_tree_prep$enriched)), color="red", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00091_tree_prep$depleted)), color="blue", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00091_tree_prep$nonsig)), color="grey", alpha=0.75, size=5)


### K09779 - uncharacterized protein

### Enriched in genomes found within obese samples
### K00941 - hydroxymethylpyrimidine/phosphomethylpyrimidine kinase
### Involved in thiamine (vitamin B1) metabolism.

ERP002061_out_K00941_tree_prep <- prep_tree_balances_func(in_list = ERP002061_out, focal_func = "K00941", taxa_table = almeida_taxa)

ggtree(ERP002061_out_K00941_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$enriched)), color="red", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$depleted)), color="blue", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00941_tree_prep$nonsig)), color="grey", alpha=0.75, size=5)


### Mixed signal of enriched and depleted in genomes found within obese samples
### K00128 - aldehyde dehydrogenase (NAD+)

ERP002061_out_K00128_tree_prep <- prep_tree_balances_func(in_list = ERP002061_out, focal_func = "K00128", taxa_table = almeida_taxa)

ggtree(ERP002061_out_K00128_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00128_tree_prep$enriched)), color="red", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00128_tree_prep$depleted)), color="blue", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP002061_out_K00128_tree_prep$nonsig)), color="grey", alpha=0.75, size=5)





### Make heatmap of # positively and negatively enriched nodes and point out outliers

ERP002061_out_num_enriched <- ERP002061_out$df[, c("num_sig_balances_pos_enrich", "num_sig_balances_neg_enrich")]

ERP002061_out_num_enriched <- ddply(ERP002061_out_num_enriched, .(num_sig_balances_pos_enrich, num_sig_balances_neg_enrich), nrow)

colnames(ERP002061_out_num_enriched)[3] <- "Count"

ggplot(ERP002061_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=Count)) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")

ggplot(ERP002061_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=log10(Count))) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")







### Obesity dataset - ERP003612
ERP003612_out_num_enriched <- ERP003612_out$df[, c("num_sig_balances_pos_enrich", "num_sig_balances_neg_enrich")]

ERP003612_out_num_enriched <- ddply(ERP003612_out_num_enriched, .(num_sig_balances_pos_enrich, num_sig_balances_neg_enrich), nrow)

colnames(ERP003612_out_num_enriched)[3] <- "Count"

ggplot(ERP003612_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=Count)) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")

ggplot(ERP003612_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=log10(Count))) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")


### K13038 - involved in coenzyme A production

ERP003612_out_K13038_tree_prep <- prep_tree_balances_func(in_list = ERP003612_out, focal_func = "K13038", taxa_table = almeida_taxa)

ggtree(ERP003612_out_K13038_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$enriched)), color="red", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$depleted)), color="blue", alpha=0.75, size=5) +
  geom_nodepoint(aes(subset=(node %in% ERP003612_out_K13038_tree_prep$nonsig)), color="grey", alpha=0.75, size=5)






### CRC dataset - ERP012177

ggtree(ERP012177_out_tree_prep$prepped_tree, layout = "circular") +
  geom_nodepoint(aes(subset=(node %in% ERP012177_out_tree_prep$nodes2plot)), color="#b5e521", alpha=0.75, size=5) +
  geom_text2(aes(label=label, subset=!isTip), hjust=-.2, size=2)



ERP012177_out_num_enriched <- ERP012177_out$df[, c("num_sig_balances_pos_enrich", "num_sig_balances_neg_enrich")]

ERP012177_out_num_enriched <- ddply(ERP012177_out_num_enriched, .(num_sig_balances_pos_enrich, num_sig_balances_neg_enrich), nrow)

colnames(ERP012177_out_num_enriched)[3] <- "Count"

ggplot(ERP012177_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=Count)) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")

ggplot(ERP012177_out_num_enriched, aes(x=num_sig_balances_pos_enrich, y=num_sig_balances_neg_enrich, fill=log10(Count))) +
  geom_tile() +
  xlab("# Positively enriched balances") +
  ylab("# Negatively enriched balances")


ERP012177_out$df[which(ERP012177_out$df$num_sig_balances_neg_enrich == 5 & ERP012177_out$df$num_sig_balances_pos_enrich == 0), ]
ERP012177_out$df[which(ERP012177_out$df$num_sig_balances_neg_enrich == 0 & ERP012177_out$df$num_sig_balances_pos_enrich == 4), ]

