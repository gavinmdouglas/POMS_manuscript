rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)

setwd("/home/gavin/github_repos/POMS_manuscript/data/results")

POMS_out <- list()

POMS_out[["Obesity 1"]] <- readRDS("Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
POMS_out[["Obesity 2"]] <- readRDS("Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
POMS_out[["Colorectal cancer"]] <- readRDS("Almeida_2019_POMS_output/ERP012177_POMS_out.rds")
POMS_out[["Mean salinity"]] <- readRDS("TARA_POMS_out.rds")[["Mean_Salinity"]]
POMS_out[["PO4"]] <- readRDS("TARA_POMS_out.rds")[["PO4"]]

tmp <- readRDS("TARA_POMS_out.rds")

for (dataset in names(POMS_out)) {
  
  
  
}

almeida_POMS_out_ERP002061 <- readRDS("Almeida_2019_POMS_output/ERP002061_POMS_out.rds")
almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP003612_POMS_out.rds")
almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP012177_POMS_out.rds")




almeida_DA_out <- readRDS(file = "Almeida_2019_DA_tool_output.rds")

ERP002061_sig_KOs <- rownames(almeida_POMS_out_ERP002061$df)[which(almeida_POMS_out_ERP002061$df$multinomial_corr < 0.25)]
ERP003612_sig_KOs <- rownames(almeida_POMS_out_ERP003612$df)[which(almeida_POMS_out_ERP003612$df$multinomial_corr < 0.25)]

ERP002061_sig_enrich_vs_abun <- data.frame(matrix(NA, nrow = length(ERP002061_sig_KOs), ncol = 8))
colnames(ERP002061_sig_enrich_vs_abun) <- c("dataset", "gene", "group1_enrich", "group2_enrich", "group1_relabun", "group2_relabun", "group1_clr", "group2_clr")
ERP002061_sig_enrich_vs_abun$dataset <- "Obesity 1"
ERP002061_sig_enrich_vs_abun$gene <- ERP002061_sig_KOs
ERP002061_sig_enrich_vs_abun$group1_enrich <- almeida_POMS_out_ERP002061$df[ERP002061_sig_enrich_vs_abun$gene, "num_sig_nodes_group1_enrich"]
ERP002061_sig_enrich_vs_abun$group2_enrich <- almeida_POMS_out_ERP002061$df[ERP002061_sig_enrich_vs_abun$gene, "num_sig_nodes_group2_enrich"]
ERP002061_sig_enrich_vs_abun$group1_relabun <- almeida_DA_out$ERP002061$wilcoxon.relab[ERP002061_sig_enrich_vs_abun$gene, "mean_group1"]
ERP002061_sig_enrich_vs_abun$group2_relabun <- almeida_DA_out$ERP002061$wilcoxon.relab[ERP002061_sig_enrich_vs_abun$gene, "mean_group2"]
ERP002061_sig_enrich_vs_abun$group1_clr <- almeida_DA_out$ERP002061$aldex2[ERP002061_sig_enrich_vs_abun$gene, "rab.win.group1"]
ERP002061_sig_enrich_vs_abun$group2_clr <- almeida_DA_out$ERP002061$aldex2[ERP002061_sig_enrich_vs_abun$gene, "rab.win.group2"]

ERP003612_sig_enrich_vs_abun <- data.frame(matrix(NA, nrow = length(ERP003612_sig_KOs), ncol = 8))
colnames(ERP003612_sig_enrich_vs_abun) <- c("dataset", "gene", "group1_enrich", "group2_enrich", "group1_relabun", "group2_relabun", "group1_clr", "group2_clr")
ERP003612_sig_enrich_vs_abun$dataset <- "Obesity 2"
ERP003612_sig_enrich_vs_abun$gene <- ERP003612_sig_KOs
ERP003612_sig_enrich_vs_abun$group1_enrich <- almeida_POMS_out_ERP003612$df[ERP003612_sig_enrich_vs_abun$gene, "num_sig_nodes_group1_enrich"]
ERP003612_sig_enrich_vs_abun$group2_enrich <- almeida_POMS_out_ERP003612$df[ERP003612_sig_enrich_vs_abun$gene, "num_sig_nodes_group2_enrich"]
ERP003612_sig_enrich_vs_abun$group1_relabun <- almeida_DA_out$ERP003612$wilcoxon.relab[ERP003612_sig_enrich_vs_abun$gene, "mean_group1"]
ERP003612_sig_enrich_vs_abun$group2_relabun <- almeida_DA_out$ERP003612$wilcoxon.relab[ERP003612_sig_enrich_vs_abun$gene, "mean_group2"]
ERP003612_sig_enrich_vs_abun$group1_clr <- almeida_DA_out$ERP003612$aldex2[ERP003612_sig_enrich_vs_abun$gene, "rab.win.group1"]
ERP003612_sig_enrich_vs_abun$group2_clr <- almeida_DA_out$ERP003612$aldex2[ERP003612_sig_enrich_vs_abun$gene, "rab.win.group2"]


sig_enrich_vs_abun <- rbind(ERP002061_sig_enrich_vs_abun, ERP003612_sig_enrich_vs_abun)

sig_enrich_vs_abun$diff_enrich <- sig_enrich_vs_abun$group1_enrich - sig_enrich_vs_abun$group2_enrich
sig_enrich_vs_abun$log2_relabun <- log2(sig_enrich_vs_abun$group1_relabun / sig_enrich_vs_abun$group2_relabun)
sig_enrich_vs_abun$diff_relabun <- sig_enrich_vs_abun$group1_relabun - sig_enrich_vs_abun$group2_relabun
sig_enrich_vs_abun$diff_clr <- sig_enrich_vs_abun$group1_clr - sig_enrich_vs_abun$group2_clr

sig_relabun_vs_enrich_w_legend <- ggplot(sig_enrich_vs_abun, aes(x = log2_relabun, y = diff_enrich, fill = dataset)) +
                                          geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
                                          geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 1) +
                                          geom_point(size = 4, shape = 21) +
                                          scale_fill_manual(name = "Dataset", values = c("gray45", "coral1", "royalblue1")) +
                                          ylab("POMS enrichment diff.") +
                                          xlab(expression(paste("log"[2]*"-fold difference in rel. abun."))) +
                                          theme_bw()
                                

sig_relabun_vs_enrich <- sig_relabun_vs_enrich_w_legend + theme(legend.position = "none")

sig_clr_vs_enrich <- ggplot(sig_enrich_vs_abun, aes(x = diff_clr, y = diff_enrich, fill = dataset)) +
                             geom_hline(yintercept = 0, linetype = "dashed", color = "grey", size = 1) +
                             geom_vline(xintercept = 0, linetype = "dashed", color = "grey", size = 1) +
                             geom_point(size = 4, shape = 21) +
                             scale_fill_manual(values = c("gray45", "coral1", "royalblue1")) +
                             ylab("POMS enrichment diff.") +
                             xlab("Diff. in CLR-transformed rel. abun.") +
                             theme_bw() +
                             theme(legend.position = "none")



legend <- get_legend(sig_relabun_vs_enrich_w_legend)

sig_abun_vs_enrich_plot <- plot_grid(sig_relabun_vs_enrich,
                                     sig_clr_vs_enrich,
                                     legend,
                                     labels = c('a', 'b', ''),
                                     nrow = 1,
                                     rel_widths = c(0.42, 0.42, 0.16))

ggsave(filename = "../../display_items/Supp_obesity_compare_enrichment_vs_relabun.pdf",
       plot = sig_abun_vs_enrich_plot,
       width = 8.5, height = 3.5, device = "pdf", dpi = 600)


ggsave(filename = "../../display_items/Supp_obesity_compare_enrichment_vs_relabun.png",
       plot = sig_abun_vs_enrich_plot,
       width = 8.5, height = 3.5, device = "png", dpi = 300)

