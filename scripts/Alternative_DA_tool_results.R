rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)
library(ggVennDiagram)

setwd("/home/gavin/github_repos/POMS_manuscript/data")
source("../scripts/POMS_manuscript_functions.R")


DA_tool_BY_sig_features <- function(DA_tool_out, correction="BY") {
                                DA_tool_BY <- list()
                                DA_tool_BY[["limma.voom"]] <- rownames(DA_tool_out$limma.voom)[which(p.adjust(DA_tool_out$limma.voom$wi.ep, correction) < 0.05)]
                                DA_tool_BY[["DESeq2"]] <- rownames(DA_tool_out$deseq2)[which(p.adjust(DA_tool_out$deseq2$pvalue, correction) < 0.05)]
                                DA_tool_BY[["Limma-Voom"]] <- rownames(DA_tool_out$limma.voom)[which(p.adjust(DA_tool_out$limma.voom$P.Value, correction) < 0.05)]
                                DA_tool_BY[["Wilcoxon (relab.)"]] <- rownames(DA_tool_out$wilcoxon.relab)[which(p.adjust(DA_tool_out$wilcoxon.relab$wilcox_p, correction) < 0.05)]
                                DA_tool_BY[["Wilcoxon (USCG-norm.)"]] <- rownames(DA_tool_out$wilcoxon.musicc)[which(p.adjust(DA_tool_out$wilcoxon.musicc$wilcox_p, correction) < 0.05)]
                                return(DA_tool_BY)
}

# DA_tool_BY_sig <- function(DA_tool_out) {
#                                   DA_tool_BY <- list()
#                                   DA_tool_BY[["limma.voom"]] <- p.adjust(DA_tool_out$limma.voom$wi.ep, "BY")
#                                   DA_tool_BY[["DESeq2"]] <- p.adjust(DA_tool_out$deseq2$pvalue, "BY")
#                                   DA_tool_BY[["Limma-Voom"]] <- p.adjust(DA_tool_out$limma.voom$P.Value, "BY")
#                                   DA_tool_BY[["Wilcoxon (relab.)"]] <- p.adjust(DA_tool_out$wilcoxon.relab$wilcox_p, "BY")
#                                   DA_tool_BY[["Wilcoxon (USCG-norm.)"]] <- p.adjust(DA_tool_out$wilcoxon.musicc$wilcox_p, "BY")
#                                   return(DA_tool_BY)
# }
# 

almeida_DA_out <- readRDS(file = "Almeida_2019_DA_tool_output.rds")
almeida_DA_out_features_BY_0.05 <- lapply(almeida_DA_out, DA_tool_BY_sig_features)



almeida_POMS_out_ERP002061 <- readRDS("Almeida_2019_POMS_output/ERP002061_pseudo.null_sig.rds")
almeida_DA_out_features_BY_0.05$ERP002061$POMS <- c(almeida_POMS_out_ERP002061$KO_up, almeida_POMS_out_ERP002061$KO_down)
almeida_DA_out_features_BY_0.05$ERP002061[["Wilcoxon (USCG-norm.)"]] <- NULL
almeida_DA_out_features_BY_0.05$ERP002061[["DESeq2"]] <- NULL
ERP002061_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP002061) +
                                scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
                                ggtitle("Primary obesity dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5),
                                      plot.margin=unit(c(0, 0, 0, 0), "mm"))


almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP003612_pseudo.null_sig.rds")
almeida_DA_out_features_BY_0.05$ERP003612$POMS <- c(almeida_POMS_out_ERP003612$KO_up, almeida_POMS_out_ERP003612$KO_down)
almeida_DA_out_features_BY_0.05$ERP003612[["Wilcoxon (USCG-norm.)"]] <- NULL
almeida_DA_out_features_BY_0.05$ERP003612[["DESeq2"]] <- NULL
ERP003612_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP003612) +
  scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
  ggtitle("Secondary obesity dataset sig. KOs") +
  theme(plot.title = element_text(hjust = 0.5))


almeida_POMS_out_ERP012177 <- readRDS("Almeida_2019_POMS_output/ERP012177_pseudo.null_sig.rds")
almeida_DA_out_features_BY_0.05$ERP012177$POMS <- c(almeida_POMS_out_ERP012177$KO_up, almeida_POMS_out_ERP012177$KO_down)
almeida_DA_out_features_BY_0.05$ERP012177[["Wilcoxon (USCG-norm.)"]] <- NULL
almeida_DA_out_features_BY_0.05$ERP012177[["DESeq2"]] <- NULL
ERP012177_venn <- ggVennDiagram(almeida_DA_out_features_BY_0.05$ERP012177) +
                                scale_fill_gradient(name="Count", low="white",high = "firebrick3") +
                                ggtitle("Colorectal cancer dataset sig. KOs") +
                                theme(plot.title = element_text(hjust = 0.5))



bottom_row <- plot_grid(ggplot() + theme_void(), ERP012177_venn, ggplot() + theme_void(),
                        rel_widths = c(0.25, 5, 0.25), nrow=1, ncol=3, labels=c('', 'c', ''), label_x = 0.225, label_y=1.05)
top_row <- plot_grid(ERP002061_venn, ERP003612_venn, labels=c('a', 'b'))

almeida_venn_plot <- plot_grid(top_row, bottom_row,nrow=2, ncol=1)

ggsave(filename = "../figures/almeida_DA_tool_venn_plot.pdf", plot = almeida_venn_plot,
       width = 10, height=8)


### Test for enriched pathways

devtools::load_all(path = "/home/gavin/github_repos/POMS/")


almeida_2019_ERP002061_background <- rownames(almeida_DA_out$ERP002061$wilcoxon.relab)


almeida_DA_out$ERP002061$aldex2_sig <- almeida_DA_out$ERP002061$aldex2[almeida_DA_out_features_BY_0.05$ERP002061$ALDEx2, ]
almeida_2019_ERP002061_aldex2_up <- rownames(almeida_DA_out$ERP002061$aldex2_sig)[which(almeida_DA_out$ERP002061$aldex2_sig$diff.btw < 0)]
almeida_2019_ERP002061_aldex2_down <- rownames(almeida_DA_out$ERP002061$aldex2_sig)[which(almeida_DA_out$ERP002061$aldex2_sig$diff.btw > 0)]


almeida_2019_ERP002061_aldex2_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP002061_aldex2_up,
                                                                     down_genes=almeida_2019_ERP002061_aldex2_down,
                                                                     gene_background=almeida_2019_ERP002061_background,
                                                                     p_corr_method="BY",
                                                                     corr_P_cutoff=0.05,
                                                                     pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                     pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")
                                            
  


almeida_DA_out$ERP002061$limma.voom_sig <- almeida_DA_out$ERP002061$limma.voom[almeida_DA_out_features_BY_0.05$ERP002061$`Limma-Voom`, ]
almeida_2019_ERP002061_limma.voom_up <- rownames(almeida_DA_out$ERP002061$limma.voom_sig)[which(almeida_DA_out$ERP002061$limma.voom_sig$logFC < 0)]
almeida_2019_ERP002061_limma.voom_down <- rownames(almeida_DA_out$ERP002061$limma.voom_sig)[which(almeida_DA_out$ERP002061$limma.voom_sig$logFC > 0)]

almeida_2019_ERP002061_limma.voom_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP002061_limma.voom_up,
                                                                     down_genes=almeida_2019_ERP002061_limma.voom_down,
                                                                     gene_background=almeida_2019_ERP002061_background,
                                                                     p_corr_method="BY",
                                                                     corr_P_cutoff=0.05,
                                                                     pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                     pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")


almeida_DA_out$ERP002061$wilcoxon.relab_sig <- almeida_DA_out$ERP002061$wilcoxon.relab[almeida_DA_out_features_BY_0.05$ERP002061$`Wilcoxon (relab.)`, ]
almeida_2019_ERP002061_wilcoxon.relab_up <- rownames(almeida_DA_out$ERP002061$wilcoxon.relab_sig)[which(almeida_DA_out$ERP002061$wilcoxon.relab_sig$mean_group1 > almeida_DA_out$ERP002061$wilcoxon.relab_sig$mean_group2)]
almeida_2019_ERP002061_wilcoxon.relab_down <- rownames(almeida_DA_out$ERP002061$wilcoxon.relab_sig)[which(almeida_DA_out$ERP002061$wilcoxon.relab_sig$mean_group1 < almeida_DA_out$ERP002061$wilcoxon.relab_sig$mean_group2)]

almeida_2019_ERP002061_wilcoxon.relab_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP002061_wilcoxon.relab_up,
                                                                             down_genes=almeida_2019_ERP002061_wilcoxon.relab_down,
                                                                              gene_background=almeida_2019_ERP002061_background,
                                                                              p_corr_method="BY",
                                                                              corr_P_cutoff=0.05,
                                                                              pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                              pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")
                                                                  





almeida_2019_ERP003612_background <- rownames(almeida_DA_out$ERP003612$wilcoxon.relab)


almeida_DA_out$ERP003612$aldex2_sig <- almeida_DA_out$ERP003612$aldex2[almeida_DA_out_features_BY_0.05$ERP003612$ALDEx2, ]
almeida_2019_ERP003612_aldex2_up <- rownames(almeida_DA_out$ERP003612$aldex2_sig)[which(almeida_DA_out$ERP003612$aldex2_sig$diff.btw < 0)]
almeida_2019_ERP003612_aldex2_down <- rownames(almeida_DA_out$ERP003612$aldex2_sig)[which(almeida_DA_out$ERP003612$aldex2_sig$diff.btw > 0)]


almeida_2019_ERP003612_aldex2_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP003612_aldex2_up,
                                                                     down_genes=almeida_2019_ERP003612_aldex2_down,
                                                                     gene_background=almeida_2019_ERP003612_background,
                                                                     p_corr_method="BY",
                                                                     corr_P_cutoff=0.05,
                                                                     pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                     pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")




almeida_DA_out$ERP003612$limma.voom_sig <- almeida_DA_out$ERP003612$limma.voom[almeida_DA_out_features_BY_0.05$ERP003612$`Limma-Voom`, ]
almeida_2019_ERP003612_limma.voom_up <- rownames(almeida_DA_out$ERP003612$limma.voom_sig)[which(almeida_DA_out$ERP003612$limma.voom_sig$logFC < 0)]
almeida_2019_ERP003612_limma.voom_down <- rownames(almeida_DA_out$ERP003612$limma.voom_sig)[which(almeida_DA_out$ERP003612$limma.voom_sig$logFC > 0)]

almeida_2019_ERP003612_limma.voom_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP003612_limma.voom_up,
                                                                         down_genes=almeida_2019_ERP003612_limma.voom_down,
                                                                         gene_background=almeida_2019_ERP003612_background,
                                                                         p_corr_method="BY",
                                                                         corr_P_cutoff=0.05,
                                                                         pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                         pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")


almeida_DA_out$ERP003612$wilcoxon.relab_sig <- almeida_DA_out$ERP003612$wilcoxon.relab[almeida_DA_out_features_BY_0.05$ERP003612$`Wilcoxon (relab.)`, ]
almeida_2019_ERP003612_wilcoxon.relab_up <- rownames(almeida_DA_out$ERP003612$wilcoxon.relab_sig)[which(almeida_DA_out$ERP003612$wilcoxon.relab_sig$mean_group1 > almeida_DA_out$ERP003612$wilcoxon.relab_sig$mean_group2)]
almeida_2019_ERP003612_wilcoxon.relab_down <- rownames(almeida_DA_out$ERP003612$wilcoxon.relab_sig)[which(almeida_DA_out$ERP003612$wilcoxon.relab_sig$mean_group1 < almeida_DA_out$ERP003612$wilcoxon.relab_sig$mean_group2)]

almeida_2019_ERP003612_wilcoxon.relab_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP003612_wilcoxon.relab_up,
                                                                             down_genes=almeida_2019_ERP003612_wilcoxon.relab_down,
                                                                             gene_background=almeida_2019_ERP003612_background,
                                                                             p_corr_method="BY",
                                                                             corr_P_cutoff=0.05,
                                                                             pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                             pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")







almeida_2019_ERP012177_background <- rownames(almeida_DA_out$ERP012177$wilcoxon.relab)


almeida_DA_out$ERP012177$aldex2_sig <- almeida_DA_out$ERP012177$aldex2[almeida_DA_out_features_BY_0.05$ERP012177$ALDEx2, ]
almeida_2019_ERP012177_aldex2_up <- rownames(almeida_DA_out$ERP012177$aldex2_sig)[which(almeida_DA_out$ERP012177$aldex2_sig$diff.btw < 0)]
almeida_2019_ERP012177_aldex2_down <- rownames(almeida_DA_out$ERP012177$aldex2_sig)[which(almeida_DA_out$ERP012177$aldex2_sig$diff.btw > 0)]


almeida_2019_ERP012177_aldex2_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP012177_aldex2_up,
                                                                     down_genes=almeida_2019_ERP012177_aldex2_down,
                                                                     gene_background=almeida_2019_ERP012177_background,
                                                                     p_corr_method="BY",
                                                                     corr_P_cutoff=0.05,
                                                                     pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                     pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")




almeida_DA_out$ERP012177$limma.voom_sig <- almeida_DA_out$ERP012177$limma.voom[almeida_DA_out_features_BY_0.05$ERP012177$`Limma-Voom`, ]
almeida_2019_ERP012177_limma.voom_up <- rownames(almeida_DA_out$ERP012177$limma.voom_sig)[which(almeida_DA_out$ERP012177$limma.voom_sig$logFC < 0)]
almeida_2019_ERP012177_limma.voom_down <- rownames(almeida_DA_out$ERP012177$limma.voom_sig)[which(almeida_DA_out$ERP012177$limma.voom_sig$logFC > 0)]

almeida_2019_ERP012177_limma.voom_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP012177_limma.voom_up,
                                                                         down_genes=almeida_2019_ERP012177_limma.voom_down,
                                                                         gene_background=almeida_2019_ERP012177_background,
                                                                         p_corr_method="BY",
                                                                         corr_P_cutoff=0.05,
                                                                         pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                         pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")


almeida_DA_out$ERP012177$wilcoxon.relab_sig <- almeida_DA_out$ERP012177$wilcoxon.relab[almeida_DA_out_features_BY_0.05$ERP012177$`Wilcoxon (relab.)`, ]
almeida_2019_ERP012177_wilcoxon.relab_up <- rownames(almeida_DA_out$ERP012177$wilcoxon.relab_sig)[which(almeida_DA_out$ERP012177$wilcoxon.relab_sig$mean_group1 > almeida_DA_out$ERP012177$wilcoxon.relab_sig$mean_group2)]
almeida_2019_ERP012177_wilcoxon.relab_down <- rownames(almeida_DA_out$ERP012177$wilcoxon.relab_sig)[which(almeida_DA_out$ERP012177$wilcoxon.relab_sig$mean_group1 < almeida_DA_out$ERP012177$wilcoxon.relab_sig$mean_group2)]

almeida_2019_ERP012177_wilcoxon.relab_pathways <- gene_set_enriched_pathways(up_genes=almeida_2019_ERP012177_wilcoxon.relab_up,
                                                                             down_genes=almeida_2019_ERP012177_wilcoxon.relab_down,
                                                                             gene_background=almeida_2019_ERP012177_background,
                                                                             p_corr_method="BY",
                                                                             corr_P_cutoff=0.05,
                                                                             pathway2func_map = "/home/gavin/projects/POMS/KEGG_mappings/prepped/KO_pathways_22Aug2019.tsv",
                                                                             pathway_descrip_infile = "/home/gavin/projects/POMS/KEGG_mappings/prepped/path_descrip_22Aug2019.tsv")




