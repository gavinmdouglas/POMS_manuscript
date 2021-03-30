rm(list=ls(all.names=TRUE))

library(ggtree)
library(ggplot2)
library(ggplotify)
library(reshape2)
library(cowplot)
library(plyr)

devtools::load_all(path = "/home/gavin/github_repos/POMS/")

setwd("/home/gavin/github_repos/POMS_manuscript/data/")


ERP002061_out <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
ERP002061_out_pseudo.null <- readRDS(file = "Almeida_2019_POMS_output/ERP002061_POMS_pseudo_null.rds")
ERP002061_pseudo_null_sig_KOs <- names(which(p.adjust(ERP002061_out_pseudo.null$pseudo_null_p, "BY") < 0.05))


ERP002061_K00941_pseudo.null <- gene_pseudo_null_dist(ERP002061_out_pseudo.null$pseudo_null_dist, "K00941")

ERP002061_pseudo_null_sig_KOs <- names(which(p.adjust(ERP002061_out_pseudo.null$pseudo_null_p, "BY") < 0.05))
ERP002061_pseudo_null_sig_KO_enrich <- ERP002061_out$df[ERP002061_pseudo_null_sig_KOs, c("num_sig_nodes_pos_enrich", "num_sig_nodes_neg_enrich")]
ERP002061_pseudo_null_sig_KO_enrich$enrich_diff <- ERP002061_pseudo_null_sig_KO_enrich$num_sig_nodes_pos_enrich - ERP002061_pseudo_null_sig_KO_enrich$num_sig_nodes_neg_enrich

ERP002061_pseudo_null_sig_KOs_up <- ERP002061_pseudo_null_sig_KOs[which(ERP002061_pseudo_null_sig_KO_enrich$enrich_diff > 0)]
ERP002061_pseudo_null_sig_KOs_down <- ERP002061_pseudo_null_sig_KOs[which(ERP002061_pseudo_null_sig_KO_enrich$enrich_diff < 0)]

ERP002061_pseudo_null_enriched_pathways <- gene_set_enriched_pathways(up_genes = ERP002061_pseudo_null_sig_KOs_up,
                                                                           down_genes = ERP002061_pseudo_null_sig_KOs_down,
                                                                           gene_background = rownames(ERP002061_out$df))

ERP002061_pseudo_null_sig <- list(KO_up=ERP002061_pseudo_null_sig_KOs_up, KO_down=ERP002061_pseudo_null_sig_KOs_down, pathways=ERP002061_pseudo_null_enriched_pathways)

saveRDS(object = ERP002061_K00941_pseudo.null, file = "Almeida_2019_POMS_output/ERP002061_K00941_pseudo.null.rds")

saveRDS(object = ERP002061_pseudo_null_sig, file = "Almeida_2019_POMS_output/ERP002061_pseudo.null_sig.rds")


# Write out table of significant hits.
ERP002061_pseudo_null_sig_KO_enrich$description <- ERP002061_out$df[rownames(ERP002061_pseudo_null_sig_KO_enrich), "description"]
ERP002061_pseudo_null_sig_KO_enrich$abs_enrich_diff <- abs(ERP002061_pseudo_null_sig_KO_enrich$enrich_diff)
ERP002061_pseudo_null_sig_KO_enrich$gene <- rownames(ERP002061_pseudo_null_sig_KO_enrich)
write.table(x = ERP002061_pseudo_null_sig_KO_enrich, file = "/home/gavin/github_repos/POMS_manuscript/data/Almeida_2019_tables/ERP002061_obesity1_pseudo_null_sig_KO_enrich.tsv",
            sep="\t", quote=FALSE, col.names = TRUE, row.names = FALSE)

