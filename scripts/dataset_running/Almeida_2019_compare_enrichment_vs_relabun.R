rm(list=ls(all.names=TRUE))

library(cowplot)
library(ggplot2)

setwd("/home/gavin/github_repos/POMS_manuscript/data")
source("../scripts/POMS_manuscript_functions.R")


almeida_DA_out <- readRDS(file = "Almeida_2019_DA_tool_output.rds")

almeida_POMS_out_ERP002061_sig <- readRDS("Almeida_2019_POMS_output/ERP002061_pseudo.null_sig.rds")
almeida_POMS_out_ERP002061 <- readRDS("Almeida_2019_POMS_output/ERP002061_POMS_out.rds")

almeida_POMS_out_ERP003612_sig <- readRDS("Almeida_2019_POMS_output/ERP003612_pseudo.null_sig.rds")
almeida_POMS_out_ERP003612 <- readRDS("Almeida_2019_POMS_output/ERP003612_POMS_out.rds")

almeida_POMS_out_ERP012177_sig <- readRDS("Almeida_2019_POMS_output/ERP012177_pseudo.null_sig.rds")
almeida_POMS_out_ERP012177 <- readRDS("Almeida_2019_POMS_output/ERP012177_POMS_out.rds")


ERP002061_all_sig <- c(almeida_POMS_out_ERP002061_sig$KO_up, almeida_POMS_out_ERP002061_sig$KO_down)
ERP002061_sig_enrich_vs_abun <- data.frame(matrix(NA, nrow=length(ERP002061_all_sig), ncol=8))
colnames(ERP002061_sig_enrich_vs_abun) <- c("dataset", "gene", "pos_enrich", "neg_enrich", "group1_relabun", "group2_relabun", "group1_clr", "group2_clr")
ERP002061_sig_enrich_vs_abun$dataset <- "Primary obesity"
ERP002061_sig_enrich_vs_abun$gene <- c(almeida_POMS_out_ERP002061_sig$KO_up, almeida_POMS_out_ERP002061_sig$KO_down)
ERP002061_sig_enrich_vs_abun$pos_enrich <- almeida_POMS_out_ERP002061$df[ERP002061_sig_enrich_vs_abun$gene, "num_sig_nodes_pos_enrich"]
ERP002061_sig_enrich_vs_abun$neg_enrich <- almeida_POMS_out_ERP002061$df[ERP002061_sig_enrich_vs_abun$gene, "num_sig_nodes_neg_enrich"]
ERP002061_sig_enrich_vs_abun$group1_relabun <- almeida_DA_out$ERP002061$wilcoxon.relab[ERP002061_sig_enrich_vs_abun$gene, "mean_group1"]
ERP002061_sig_enrich_vs_abun$group2_relabun <- almeida_DA_out$ERP002061$wilcoxon.relab[ERP002061_sig_enrich_vs_abun$gene, "mean_group2"]
ERP002061_sig_enrich_vs_abun$group1_clr <- almeida_DA_out$ERP002061$aldex2[ERP002061_sig_enrich_vs_abun$gene, "rab.win.group1"]
ERP002061_sig_enrich_vs_abun$group2_clr <- almeida_DA_out$ERP002061$aldex2[ERP002061_sig_enrich_vs_abun$gene, "rab.win.group2"]



ERP003612_all_sig <- c(almeida_POMS_out_ERP003612_sig$KO_up, almeida_POMS_out_ERP003612_sig$KO_down)
ERP003612_sig_enrich_vs_abun <- data.frame(matrix(NA, nrow=length(ERP003612_all_sig), ncol=8))
colnames(ERP003612_sig_enrich_vs_abun) <- c("dataset", "gene", "pos_enrich", "neg_enrich", "group1_relabun", "group2_relabun", "group1_clr", "group2_clr")
ERP003612_sig_enrich_vs_abun$dataset <- "Secondary obesity"
ERP003612_sig_enrich_vs_abun$gene <- c(almeida_POMS_out_ERP003612_sig$KO_up, almeida_POMS_out_ERP003612_sig$KO_down)
ERP003612_sig_enrich_vs_abun$pos_enrich <- almeida_POMS_out_ERP003612$df[ERP003612_sig_enrich_vs_abun$gene, "num_sig_nodes_pos_enrich"]
ERP003612_sig_enrich_vs_abun$neg_enrich <- almeida_POMS_out_ERP003612$df[ERP003612_sig_enrich_vs_abun$gene, "num_sig_nodes_neg_enrich"]
ERP003612_sig_enrich_vs_abun$group1_relabun <- almeida_DA_out$ERP003612$wilcoxon.relab[ERP003612_sig_enrich_vs_abun$gene, "mean_group1"]
ERP003612_sig_enrich_vs_abun$group2_relabun <- almeida_DA_out$ERP003612$wilcoxon.relab[ERP003612_sig_enrich_vs_abun$gene, "mean_group2"]
ERP003612_sig_enrich_vs_abun$group1_clr <- almeida_DA_out$ERP003612$aldex2[ERP003612_sig_enrich_vs_abun$gene, "rab.win.group1"]
ERP003612_sig_enrich_vs_abun$group2_clr <- almeida_DA_out$ERP003612$aldex2[ERP003612_sig_enrich_vs_abun$gene, "rab.win.group2"]


ERP012177_all_sig <- c(almeida_POMS_out_ERP012177_sig$KO_up, almeida_POMS_out_ERP012177_sig$KO_down)
ERP012177_sig_enrich_vs_abun <- data.frame(matrix(NA, nrow=length(ERP012177_all_sig), ncol=8))
colnames(ERP012177_sig_enrich_vs_abun) <- c("dataset", "gene", "pos_enrich", "neg_enrich", "group1_relabun", "group2_relabun", "group1_clr", "group2_clr")
ERP012177_sig_enrich_vs_abun$dataset <- "Colorectal cancer"
ERP012177_sig_enrich_vs_abun$gene <- c(almeida_POMS_out_ERP012177_sig$KO_up, almeida_POMS_out_ERP012177_sig$KO_down)
ERP012177_sig_enrich_vs_abun$pos_enrich <- almeida_POMS_out_ERP012177$df[ERP012177_sig_enrich_vs_abun$gene, "num_sig_nodes_pos_enrich"]
ERP012177_sig_enrich_vs_abun$neg_enrich <- almeida_POMS_out_ERP012177$df[ERP012177_sig_enrich_vs_abun$gene, "num_sig_nodes_neg_enrich"]
ERP012177_sig_enrich_vs_abun$group1_relabun <- almeida_DA_out$ERP012177$wilcoxon.relab[ERP012177_sig_enrich_vs_abun$gene, "mean_group1"]
ERP012177_sig_enrich_vs_abun$group2_relabun <- almeida_DA_out$ERP012177$wilcoxon.relab[ERP012177_sig_enrich_vs_abun$gene, "mean_group2"]
ERP012177_sig_enrich_vs_abun$group1_clr <- almeida_DA_out$ERP012177$aldex2[ERP012177_sig_enrich_vs_abun$gene, "rab.win.group1"]
ERP012177_sig_enrich_vs_abun$group2_clr <- almeida_DA_out$ERP012177$aldex2[ERP012177_sig_enrich_vs_abun$gene, "rab.win.group2"]


sig_enrich_vs_abun <- rbind(ERP002061_sig_enrich_vs_abun, ERP003612_sig_enrich_vs_abun, ERP012177_sig_enrich_vs_abun)

sig_enrich_vs_abun$diff_enrich <- sig_enrich_vs_abun$pos_enrich - sig_enrich_vs_abun$neg_enrich
sig_enrich_vs_abun$log2_relabun <- log2(sig_enrich_vs_abun$group1_relabun / sig_enrich_vs_abun$group2_relabun)
sig_enrich_vs_abun$diff_relabun <- sig_enrich_vs_abun$group1_relabun - sig_enrich_vs_abun$group2_relabun
sig_enrich_vs_abun$diff_clr <- sig_enrich_vs_abun$group1_clr - sig_enrich_vs_abun$group2_clr

sig_relabun_vs_enrich_w_legend <- ggplot(sig_enrich_vs_abun, aes(x=log2_relabun, y=diff_enrich, fill=dataset)) +
                                          geom_hline(yintercept=0, linetype="dashed", color = "grey", size=1) +
                                          geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
                                          geom_point(size=4, shape=21) +
                                          scale_fill_manual(name="Dataset", values=c("gray45", "coral1", "royalblue1")) +
                                          ylab("POMS enrichment diff.") +
                                          xlab(expression(paste("log"[2]*"-fold difference in rel. abun."))) +
                                          theme_bw()
                                

sig_relabun_vs_enrich <- sig_relabun_vs_enrich_w_legend + theme(legend.position = "none")

sig_clr_vs_enrich <- ggplot(sig_enrich_vs_abun, aes(x=diff_clr, y=diff_enrich, fill=dataset)) +
                             geom_hline(yintercept=0, linetype="dashed", color = "grey", size=1) +
                             geom_vline(xintercept=0, linetype="dashed", color = "grey", size=1) +
                             geom_point(size=4, shape=21) +
                             scale_fill_manual(values=c("gray45", "coral1", "royalblue1")) +
                             ylab("POMS enrichment diff.") +
                             xlab("Diff. in CLR-transformed rel. abun.") +
                             theme_bw() +
                             theme(legend.position = "none")



legend <- get_legend(sig_relabun_vs_enrich_w_legend)

sig_abun_vs_enrich_plot <- plot_grid(sig_relabun_vs_enrich,
                                     sig_clr_vs_enrich,
                                     legend,
                                     labels=c('a', 'b', ''),
                                     nrow=1,
                                     rel_widths = c(0.42, 0.42, 0.16))

ggsave(filename = "../figures/Almeida_sig_KOs_abun_vs_enrich_plot.pdf", plot = sig_abun_vs_enrich_plot,
       width = 8.5, height=3.5)
