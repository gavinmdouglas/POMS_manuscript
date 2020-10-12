rm(list=ls(all.names=TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/")

ERP003612_out <- readRDS(file = "Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
ERP003612_out_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP003612_POMS_pseudo_null.rds")
ERP003612_pseudo_null_sig_KOs <- names(which(p.adjust(ERP003612_out_pseudo.null$pseudo_null_p, "BY") < 0.05))
ERP003612_pseudo_null_sig_KO_enrich <- ERP003612_out$df[ERP003612_pseudo_null_sig_KOs, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]
ERP003612_pseudo_null_sig_KO_enrich$enrich_diff <- ERP003612_pseudo_null_sig_KO_enrich$num_sig_nodes_pos_enrich - ERP003612_pseudo_null_sig_KO_enrich$num_sig_nodes_neg_enrich
ERP003612_pseudo_null_sig_KOs_up <- ERP003612_pseudo_null_sig_KOs[which(ERP003612_pseudo_null_sig_KO_enrich$enrich_diff > 0)]
ERP003612_pseudo_null_sig_KOs_down <- ERP003612_pseudo_null_sig_KOs[which(ERP003612_pseudo_null_sig_KO_enrich$enrich_diff < 0)]
ERP003612_pseudo_null_enriched_pathways <- gene_set_enriched_pathways(up_genes = ERP003612_pseudo_null_sig_KOs_up,
                                                                           down_genes = ERP003612_pseudo_null_sig_KOs_down,
                                                                           gene_background = rownames(ERP003612_out$df))

ERP003612_pseudo_null_sig <- list(KO_up=ERP003612_pseudo_null_sig_KOs_up, KO_down=ERP003612_pseudo_null_sig_KOs_down, pathways=ERP003612_pseudo_null_enriched_pathways)

ERP003612_K13038_pseudo.null <- gene_pseudo_null_dist(ERP003612_out_pseudo.null$pseudo_null_dist, "K13038")



ERP012177_out <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
ERP012177_out_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP012177_POMS_pseudo_null.rds")
ERP012177_pseudo_null_sig_KOs <- names(which(p.adjust(ERP012177_out_pseudo.null$pseudo_null_p, "BY") < 0.05))
ERP012177_pseudo_null_sig_KO_enrich <- ERP012177_out$df[ERP012177_pseudo_null_sig_KOs, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]
ERP012177_pseudo_null_sig_KO_enrich$enrich_diff <- ERP012177_pseudo_null_sig_KO_enrich$num_sig_nodes_pos_enrich - ERP012177_pseudo_null_sig_KO_enrich$num_sig_nodes_neg_enrich
ERP012177_pseudo_null_sig_KOs_up <- ERP012177_pseudo_null_sig_KOs[which(ERP012177_pseudo_null_sig_KO_enrich$enrich_diff > 0)]
ERP012177_pseudo_null_sig_KOs_down <- ERP012177_pseudo_null_sig_KOs[which(ERP012177_pseudo_null_sig_KO_enrich$enrich_diff < 0)]
ERP012177_pseudo_null_enriched_pathways <- gene_set_enriched_pathways(up_genes = ERP012177_pseudo_null_sig_KOs_up,
                                                                      down_genes = ERP012177_pseudo_null_sig_KOs_down,
                                                                      gene_background = rownames(ERP012177_out$df))

ERP012177_pseudo_null_sig <- list(KO_up=ERP012177_pseudo_null_sig_KOs_up, KO_down=ERP012177_pseudo_null_sig_KOs_down, pathways=ERP012177_pseudo_null_enriched_pathways)

ERP012177_K00702_pseudo.null <- gene_pseudo_null_dist(ERP012177_out_pseudo.null$pseudo_null_dist, "K00702")
ERP012177_K09691_pseudo.null <- gene_pseudo_null_dist(ERP012177_out_pseudo.null$pseudo_null_dist, "K09691")


saveRDS(object = ERP012177_pseudo_null_sig, file = "Almeida_2019_POMS_output/ERP012177_pseudo.null_sig.rds")
saveRDS(object = ERP003612_pseudo_null_sig, file = "Almeida_2019_POMS_output/ERP003612_pseudo.null_sig.rds")

saveRDS(object = ERP012177_K00702_pseudo.null, file = "Almeida_2019_POMS_output/ERP012177_K00702_pseudo.null.rds")
saveRDS(object = ERP012177_K09691_pseudo.null, file = "Almeida_2019_POMS_output/ERP012177_K09691_pseudo.null.rds")
saveRDS(object = ERP003612_K13038_pseudo.null, file = "Almeida_2019_POMS_output/ERP003612_K13038_pseudo.null.rds")

